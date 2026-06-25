%% build_transducer_library  Build or refresh the full transducer calibration library
%
% For each equipment combo in config/equipment/ this script either:
%   - Runs an empirical geometric-baseline build (geometric phases + velocity fit)
%     when axial characterisation data (char_data_path) is available, or
%   - Computes a pure geometry-only entry from ring dimensions alone.
%
% Both paths produce a global_model section in config/calibration/{combo}.yaml
% that load_transducer_from_library can interpolate at any focal depth.
%
% QUICK START
%   1. Set CHAR_DATA_ROOT to the folder that contains your axial profile CSVs
%      (or leave '' if char_data_path entries in the equipment YAML are absolute).
%   2. Optionally restrict processing with COMBOS_TO_RUN.
%   3. Run the script. Results are written to config/calibration/ inside the repo.
%
% ADDING A NEW TRANSDUCER
%   Add a YAML to config/equipment/, set char_data_path if you have hydrophone
%   data, then re-run this script. The new combo is processed automatically.
%
% PER-DEPTH REFINEMENTS
%   Use examples/calibration_standalone.m to layer single-depth refinements
%   (calibrate_transducer Mode='single_ref') on top of the library baseline.
%
% See also: LOAD_TRANSDUCER_FROM_LIBRARY, UPDATE_TRANSDUCER_LIBRARY,
%           VALIDATE_TRANSDUCER_LIBRARY, CALIBRATION_STANDALONE

close all;

%% ── Paths ────────────────────────────────────────────────────────────────────
repo_root    = fullfile(fileparts(mfilename('fullpath')), '..');
equip_path   = fullfile(repo_root, 'config', 'equipment');
library_path = fullfile(repo_root, 'config', 'calibration');

addpath(genpath(fullfile(repo_root, 'functions')));
addpath(genpath(fullfile(repo_root, 'external')));

%% ── CONFIGURATION ────────────────────────────────────────────────────────────
% Variables are only set here when not already defined in the caller's workspace.
% Site-specific wrapper scripts (e.g. build_transducer_library_donders.m) can
% pre-set any of these before calling run() on this file to override defaults.

% Which combos to process.
%   'all'  — every combo found in config/equipment/
%   cell   — subset, e.g. {'IS_PCD15473_01001_IGT_32_ch_comb_10_ch'}
if ~exist('COMBOS_TO_RUN',      'var'); COMBOS_TO_RUN      = 'all'; end

% Root folder for axial-profile CSVs referenced by char_data_path in the
% equipment YAML. Set to '' if those paths are already absolute.
if ~exist('CHAR_DATA_ROOT',     'var'); CHAR_DATA_ROOT     = ''; end

% Root folder for manufacturer phase table files referenced by phase_table in
% the equipment YAML. Set to '' if those paths are already absolute.
if ~exist('PHASE_TABLE_ROOT',   'var'); PHASE_TABLE_ROOT   = ''; end

% Intensities [W/cm²] to calibrate amplitude for (empirical combos only).
if ~exist('DESIRED_INTENSITIES','var'); DESIRED_INTENSITIES = [30]; end

% Scratch folder for intermediate per-depth YAMLs.
% Defaults to a temp subfolder; set to a persistent path to keep them.
if ~exist('OUTPUT_YAML_FOLDER', 'var')
    OUTPUT_YAML_FOLDER = fullfile(tempdir, 'prestus_library_build');
end

% Amplitude fit settings
if ~exist('OPT_UPPER_VELOCITY', 'var'); OPT_UPPER_VELOCITY = 0.2;  end
if ~exist('OPT_WEIGHTS',        'var'); OPT_WEIGHTS        = 0;    end
if ~exist('SKIP_FRONT_PEAK_MM', 'var'); SKIP_FRONT_PEAK_MM = 10;   end
if ~exist('INITIAL_VELOCITY',   'var'); INITIAL_VELOCITY   = 0.05; end

% Water sound speed [m/s] used for geometric phase computation.
if ~exist('SOUND_SPEED',        'var'); SOUND_SPEED        = 1500; end

% ── END CONFIGURATION ─────────────────────────────────────────────────────────

if ~exist(OUTPUT_YAML_FOLDER, 'dir'); mkdir(OUTPUT_YAML_FOLDER); end

%% ── Load equipment config ────────────────────────────────────────────────────
equip_param = load_equipment_config(equip_path);
all_serials = fieldnames(equip_param.trans);

% Build default parameters struct from repo config
config_default_path = fullfile(repo_root, 'config', 'config_default.yaml');
if isfile(config_default_path)
    parameters = yaml.loadFile(config_default_path, 'ConvertToArray', true);
else
    parameters = struct();
    warning('build_transducer_library: config_default.yaml not found at %s — using empty base parameters.', ...
        config_default_path);
end

%% ── Build combo list ─────────────────────────────────────────────────────────
combo_list = {};   % {combo_name, tran_serial, ds_serial, has_char}

for si = 1:numel(all_serials)
    tran_serial = all_serials{si};
    tran        = equip_param.trans.(tran_serial);

    if ~isfield(tran, 'combos') || isempty(tran.combos)
        combo_name = tran_serial;
        if ~strcmp(COMBOS_TO_RUN, 'all') && ~any(strcmp(COMBOS_TO_RUN, combo_name))
            continue;
        end
        combo_list(end+1,:) = {combo_name, tran_serial, '', false}; %#ok<AGROW>
        continue;
    end

    ds_keys = fieldnames(tran.combos);
    for di = 1:numel(ds_keys)
        ds_serial  = ds_keys{di};
        combo_name = [tran_serial '_' ds_serial];

        if ~strcmp(COMBOS_TO_RUN, 'all') && ~any(strcmp(COMBOS_TO_RUN, combo_name))
            continue;
        end

        has_char = isfield(tran.combos.(ds_serial), 'char_data_path') && ...
                   ~isempty(tran.combos.(ds_serial).char_data_path);
        combo_list(end+1,:) = {combo_name, tran_serial, ds_serial, has_char}; %#ok<AGROW>
    end
end

if isempty(combo_list)
    fprintf('No combos matched COMBOS_TO_RUN — nothing to do.\n');
    return;
end

n_total     = size(combo_list, 1);
n_empirical = sum([combo_list{:,4}]);
n_geometric = n_total - n_empirical;

fprintf('\n=== build_transducer_library: %d combo(s) — %d empirical, %d geometric ===\n\n', ...
    n_total, n_empirical, n_geometric);

%% ── Process each combo ───────────────────────────────────────────────────────
results = struct('combo', {}, 'level', {}, 'status', {}, 'msg', {});

for ci = 1:n_total
    combo_name  = combo_list{ci,1};
    tran_serial = combo_list{ci,2};
    ds_serial   = combo_list{ci,3};
    has_char    = combo_list{ci,4};
    tran        = equip_param.trans.(tran_serial);

    level = 'geometric';
    if has_char; level = 'empirical'; end

    fprintf('[%d/%d] %s  (%s)\n', ci, n_total, combo_name, level);
    rec = struct('combo', combo_name, 'level', level, 'status', '', 'msg', '');

    try
        if strcmp(level, 'geometric')
            build_geometric_entry(tran_serial, ds_serial, tran, library_path, SOUND_SPEED);
        else
            combo_entry = tran.combos.(ds_serial);
            build_empirical_entry(combo_name, tran_serial, ds_serial, tran, combo_entry, ...
                equip_param, parameters, ...
                CHAR_DATA_ROOT, PHASE_TABLE_ROOT, ...
                DESIRED_INTENSITIES, OPT_UPPER_VELOCITY, OPT_WEIGHTS, ...
                SKIP_FRONT_PEAK_MM, INITIAL_VELOCITY, ...
                OUTPUT_YAML_FOLDER, library_path);
        end
        rec.status = 'ok';
        fprintf('  [OK]\n');
    catch ME
        rec.status = 'error';
        rec.msg    = ME.message;
        fprintf('  [ERROR] %s  (%s line %d)\n', ME.message, ME.stack(1).file, ME.stack(1).line);
    end

    results(end+1) = rec; %#ok<AGROW>
end

%% ── Post-build validation ────────────────────────────────────────────────────
fprintf('\n--- Validating library ---\n');
try
    validate_transducer_library(library_path, equip_path);
catch ME
    fprintf('[WARN] Validation error: %s\n', ME.message);
end

%% ── Update library_version.yaml ─────────────────────────────────────────────
manifest_path = fullfile(library_path, 'library_version.yaml');
if isfile(manifest_path)
    fprintf('\n--- Updating library_version.yaml ---\n');
    manifest = yaml.loadFile(manifest_path, 'ConvertToArray', true);
    n_added  = 0;
    for ci = 1:numel(results)
        r = results(ci);
        if ~strcmp(r.status, 'ok') || isfield(manifest.combos, r.combo); continue; end
        manifest.combos.(r.combo) = struct( ...
            'calibration_type', r.level, ...
            'added_date',       datestr(now, 'yyyy-mm-dd')); %#ok<TNOW1,DATST>
        n_added = n_added + 1;
    end
    if n_added > 0
        yaml.dumpFile(manifest_path, manifest, 'block');
        fprintf('[OK]  %d new combo(s) added. Bump patch version before committing.\n', n_added);
    else
        fprintf('[OK]  library_version.yaml already up to date.\n');
    end
end

%% ── Summary ──────────────────────────────────────────────────────────────────
n_ok  = sum(strcmp({results.status}, 'ok'));
n_err = sum(strcmp({results.status}, 'error'));

fprintf('\n=== Done: %d ok, %d error(s). Library at: %s ===\n', n_ok, n_err, library_path);
if n_err > 0
    fprintf('\nFailed combos:\n');
    failed = results(strcmp({results.status}, 'error'));
    for fi = 1:numel(failed)
        fprintf('  %-50s %s\n', failed(fi).combo, failed(fi).msg);
    end
end

%% ═══════════════════════ Local helpers ════════════════════════════════════════

function build_geometric_entry(tran_serial, ds_serial, tran, library_path, sound_speed)
% Compute geometric phases from ring dimensions and write library YAML.

    if isempty(ds_serial)
        combo_key = tran_serial;
    else
        combo_key = [tran_serial '_' ds_serial];
    end

    % Auto depth grid: 10 steps from ~15% to ~95% of curvature radius
    r_c = tran.transducer.annular.curv_radius_mm;
    lo  = r_c * 0.15;
    hi  = r_c * 0.95;
    if isfield(tran, 'min_foc') && ~isempty(tran.min_foc); lo = tran.min_foc; end
    if isfield(tran, 'max_foc') && ~isempty(tran.max_foc); hi = tran.max_foc; end
    depths_mm = round(linspace(lo, hi, 10), 1);
    fprintf('  Depth grid (auto): %s mm\n', mat2str(depths_mm));

    % Minimal parameters struct for set_real_phases
    params = struct();
    params.medium_properties.water.sound_speed = sound_speed;
    params.medium_properties.water.density     = 994;
    params.transducer                           = tran.transducer;
    params.transducer.n_elem                   = tran.transducer.annular.elem_n;
    params = load_parameters(params);

    n_elem   = tran.transducer.annular.elem_n;
    n_depths = numel(depths_mm);
    phase_mat = zeros(n_depths, n_elem);
    dist_bowl_offset = tran.transducer.annular.curv_radius_mm - ...
                       tran.transducer.annular.dist_geom_ep_mm;

    for j = 1:n_depths
        d = depths_mm(j);
        params.transducer.focal_distance_ep   = d;
        params.transducer.focal_distance_bowl = d + dist_bowl_offset;
        phase_mat(j, :) = set_real_phases(tran, d, params);
    end

    out_file = fullfile(library_path, [combo_key '.yaml']);
    if isfile(out_file)
        lib = yaml.loadFile(out_file, 'ConvertToArray', true);
    else
        lib = struct();
    end

    lib.meta.tran_serial      = tran_serial;
    lib.meta.created_at       = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    lib.meta.calibration_type = 'geometric';
    try; lib.meta.prestus_version = prestus_version(); catch; lib.meta.prestus_version = 'unknown'; end

    lib.global_model.depths_ep_mm   = depths_mm;
    lib.global_model.elem_phase_deg = num2cell(phase_mat, 2)';

    yaml.dumpFile(out_file, lib);
    fprintf('  Written: %s\n', out_file);
end


function build_empirical_entry(combo_name, tran_serial, ds_serial, tran, combo_entry, ...
        equip_param, parameters, ...
        char_data_root, phase_table_root, ...
        desired_intensities, opt_upper_velocity, opt_weights, ...
        skip_front_peak_mm, initial_velocity, ...
        output_yaml_folder, library_path)
% Geometric phases + per-depth velocity fit from axial hydrophone profiles.

    parameters.transducer = tran.transducer;
    parameters.transducer.annular.elem_phase_deg = zeros(1, tran.transducer.annular.elem_n);
    parameters.transducer.annular.elem_p_amp       = repmat(200000, 1, tran.transducer.annular.elem_n);
    parameters.medium_properties.water.sound_speed = 1500;
    parameters.medium_properties.water.density     = 1000;
    if ~isstruct(parameters.calibration); parameters.calibration = struct(); end
    parameters.calibration.opt_upper_velocity  = opt_upper_velocity;
    parameters.calibration.opt_weights         = opt_weights;
    parameters.calibration.skip_front_peak_mm  = skip_front_peak_mm;
    parameters.calibration.forward_model       = 'rayleigh';

    tmp_out = fullfile(tempdir, ['prestus_calib_' combo_name]);
    if ~exist(tmp_out, 'dir'); mkdir(tmp_out); end
    parameters.io.outputs_folder = tmp_out;

    % Resolve char data path
    char_rel = combo_entry.char_data_path;
    if isempty(char_data_root)
        char_path = char_rel;
    else
        char_path = fullfile(char_data_root, char_rel);
    end
    if ~isfile(char_path)
        warning('build_transducer_library:missingCharData', ...
            'Characterisation data not found: %s\nFalling back to geometry-only entry.\n(Set CHAR_DATA_ROOT if paths are relative.)', char_path);
        build_geometric_entry(tran_serial, ds_serial, tran, library_path, 1500);
        return;
    end

    charac_data         = readmatrix(char_path, 'NumHeaderLines', 1);
    available_depths_ep = round(charac_data(1, 2:end), 2);
    dist_ep             = charac_data(2:end, 1);
    intens_data         = charac_data(2:end, 2:end);
    fprintf('  Depths: %s mm\n', mat2str(available_depths_ep));

    dist_bowl_offset = tran.transducer.annular.curv_radius_mm - ...
                       tran.transducer.annular.dist_geom_ep_mm;
    dist_bowl        = dist_ep + dist_bowl_offset;

    % Optionally load manufacturer phase table (ignored for geometric baseline)
    phase_table = [];
    if isfield(combo_entry, 'phase_table') && ~isempty(combo_entry.phase_table)
        pt_rel = combo_entry.phase_table;
        pt_path = pt_rel;
        if ~isempty(phase_table_root)
            pt_path = fullfile(phase_table_root, pt_rel);
        end
        if ~isfile(pt_path)
            fprintf('  [INFO] phase_table not found (%s) — using geometry-only phases.\n', pt_path);
        elseif strcmpi(tran.manufact, 'Sonic Concepts')
            phase_table = readtable(pt_path);
        elseif strcmpi(tran.manufact, 'Imasonic')
            phase_table = read_ini_file(pt_path);
        else
            fprintf('  [INFO] Unknown manufacturer ''%s'' — ignoring phase_table.\n', tran.manufact);
        end
    end

    n_depths = numel(available_depths_ep);
    fprintf('  Computing geometric phases and fitting velocity at %d depths...\n', n_depths);

    for j = 1:n_depths
        depth_j  = available_depths_ep(j);
        params_j = parameters;
        params_j.transducer.focal_distance_ep          = depth_j;
        params_j.transducer.focal_distance_bowl        = depth_j + dist_bowl_offset;
        params_j.calibration.desired_focal_distance_ep = depth_j;
        params_j = load_parameters(params_j);

        geo_deg_j  = set_real_phases(tran, depth_j, params_j);
        phases_rad = geo_deg_j * pi / 180;
        params_j.transducer.annular.elem_phase_deg = geo_deg_j;
        params_j.transducer.annular.elem_phase_rad = phases_rad;

        emp_j = intens_data(:, j);
        emp_j(emp_j < 0) = 0;

        for ii = 1:numel(desired_intensities)
            isppa  = desired_intensities(ii);
            emp_sc = emp_j * (isppa / max(emp_j + eps));

            profile_struct.axial_intensity     = emp_sc;
            profile_struct.axial_distance_bowl = dist_bowl;
            params_j.calibration.desired_intensity = isppa;

            [vel_fit, ~, ~] = fit_velocity_to_intensity( ...
                params_j, profile_struct, phases_rad, initial_velocity, isppa);

            water_Z      = params_j.medium_properties.water.density * ...
                           params_j.medium_properties.water.sound_speed;
            elem_p_amp_val = round(vel_fit * water_Z);

            foc_str   = strrep(num2str(depth_j), '.', 'p');
            int_str   = strrep(num2str(isppa),   '.', 'p');
            yaml_name = sprintf('%s-F%smm-I%swpercm2.yaml', combo_name, foc_str, int_str);
            yaml_path = fullfile(output_yaml_folder, yaml_name);

            data = struct();
            data.meta.tran_serial = tran_serial;
            data.meta.created_at  = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
            try; data.meta.prestus_version = prestus_version(); catch; data.meta.prestus_version = 'unknown'; end
            data.transducer.focal_distance_ep       = depth_j;
            data.transducer.set_intensity_w_per_cm2 = isppa;
            data.transducer.annular.elem_phase_deg  = geo_deg_j;
            data.transducer.annular.elem_phase_rad  = phases_rad;
            data.transducer.annular.elem_p_amp        = repmat(elem_p_amp_val, 1, tran.transducer.annular.elem_n);

            yaml.dumpFile(yaml_path, data);
            fprintf('  Saved: %s (amp=%d)\n', yaml_name, elem_p_amp_val);
        end
    end

    update_transducer_library(combo_name, output_yaml_folder, library_path);
end
