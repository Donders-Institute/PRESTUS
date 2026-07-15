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
repo_root    = fullfile(fileparts(mfilename('fullpath')), '..', '..');
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

function build_geometric_entry(tran_serial, ds_serial, tran, library_path, sound_speed) %#ok<INUSD>
% Write a geometry-only library entry (no empirical phase correction).
%
% Stores an empty global_model.elem_phase_correction_deg so that
% load_transducer_from_library reconstructs pure geometric steering phases
% from ring geometry at the requested depth. No amplitude calibration is
% available without characterisation data.

    if isempty(ds_serial)
        combo_key = tran_serial;
    else
        combo_key = [tran_serial '_' ds_serial];
    end

    out_file = fullfile(library_path, [combo_key '.yaml']);
    if isfile(out_file)
        lib = yaml.loadFile(out_file, 'ConvertToArray', true);
    else
        lib = struct();
    end

    lib.meta.tran_serial = tran_serial;
    if ~isempty(ds_serial); lib.meta.driving_system = ds_serial; end
    if ~isfield(lib, 'meta') || ~isfield(lib.meta, 'created_at') || isempty(lib.meta.created_at)
        lib.meta.created_at = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    end
    lib.meta.last_updated     = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    lib.meta.calibration_type = 'geometric';
    try; lib.meta.prestus_version = prestus_version(); catch; lib.meta.prestus_version = 'unknown'; end

    % Drop any legacy per-depth phase tables; phases are computed at load time.
    if isfield(lib, 'global_model')
        for legacy = {'depths_ep_mm', 'elem_phase_deg'}
            if isfield(lib.global_model, legacy{1})
                lib.global_model = rmfield(lib.global_model, legacy{1});
            end
        end
    end
    lib.global_model.elem_phase_correction_deg = [];

    yaml.dumpFile(out_file, lib);
    fprintf('  Written (geometric): %s\n', out_file);
end


function build_empirical_entry(combo_name, tran_serial, ds_serial, tran, combo_entry, ...
        equip_param, parameters, ...
        char_data_root, phase_table_root, ...
        desired_intensities, opt_upper_velocity, opt_weights, ...
        skip_front_peak_mm, initial_velocity, ...
        output_yaml_folder, library_path) %#ok<INUSD>
% Joint multi-depth phase correction + amplitude scaling from axial profiles.
%
% Fits a single depth-independent per-element phase correction
% (elem_phase_correction_deg) jointly across all characterisation depths
% via perform_joint_depth_fit, then derives a source amplitude per requested
% intensity. Only the correction and amplitude scaling are written to the
% library; geometric phases are reconstructed at load time by
% load_transducer_from_library.

    parameters.transducer = tran.transducer;
    parameters.transducer.annular.elem_phase_deg = zeros(1, tran.transducer.annular.elem_n);
    parameters.transducer.annular.elem_p_amp       = repmat(200000, 1, tran.transducer.annular.elem_n);
    parameters.medium_properties.water.sound_speed = 1500;
    parameters.medium_properties.water.density     = 1000;
    if ~isstruct(parameters.calibration); parameters.calibration = struct(); end
    parameters.calibration.opt_upper_velocity  = opt_upper_velocity;
    parameters.calibration.opt_weights         = opt_weights;
    parameters.calibration.skip_front_peak_mm  = skip_front_peak_mm;
    parameters.calibration.initial_velocity    = initial_velocity;
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

    % --- Joint multi-depth phase correction (depth-independent) ---
    ref_intensity = desired_intensities(1);
    fit_opts = struct('ForwardModel', parameters.calibration.forward_model);
    fprintf('  Fitting joint phase correction across %d depths...\n', numel(available_depths_ep));
    [opt_delta_rad, opt_velocity_ref, ~] = perform_joint_depth_fit( ...
        parameters, intens_data, available_depths_ep, dist_bowl, dist_bowl_offset, ...
        ref_intensity, tran, fit_opts);
    elem_phase_correction_deg = mod(opt_delta_rad(:)' * 180/pi, 360);

    % --- Amplitude scaling per intensity ---
    % Source amplitude = particle velocity x water impedance; velocity scales
    % as sqrt(I) since I ∝ v². Anchored on the joint-fit reference velocity.
    water_Z    = parameters.medium_properties.water.density * ...
                 parameters.medium_properties.water.sound_speed;
    amp_struct = struct();
    for ii = 1:numel(desired_intensities)
        isppa = desired_intensities(ii);
        vel_i = opt_velocity_ref * sqrt(isppa / ref_intensity);
        amp_i = round(vel_i * water_Z);
        key   = sprintf('i%s', strrep(num2str(isppa), '.', 'p'));
        amp_struct.(key) = amp_i;
        fprintf('  I=%g W/cm² -> amp=%d\n', isppa, amp_i);
    end

    % --- Write library YAML (correction-only global model) ---
    out_file = fullfile(library_path, [combo_name '.yaml']);
    if isfile(out_file)
        lib = yaml.loadFile(out_file, 'ConvertToArray', true);
    else
        lib = struct();
    end

    lib.meta.tran_serial = tran_serial;
    if ~isempty(ds_serial); lib.meta.driving_system = ds_serial; end
    if ~isfield(lib, 'meta') || ~isfield(lib.meta, 'created_at') || isempty(lib.meta.created_at)
        lib.meta.created_at = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    end
    lib.meta.last_updated     = datestr(now, 'yyyy-mm-dd'); %#ok<TNOW1,DATST>
    lib.meta.calibration_type = 'empirical';
    try; lib.meta.prestus_version = prestus_version(); catch; lib.meta.prestus_version = 'unknown'; end

    % Drop any legacy per-depth phase tables so the entry is correction-only.
    if isfield(lib, 'global_model')
        for legacy = {'depths_ep_mm', 'elem_phase_deg'}
            if isfield(lib.global_model, legacy{1})
                lib.global_model = rmfield(lib.global_model, legacy{1});
            end
        end
    end
    lib.global_model.elem_phase_correction_deg = elem_phase_correction_deg;
    lib.global_model.amplitude_scaling         = amp_struct;

    yaml.dumpFile(out_file, lib);
    fprintf('  Written (correction-only): %s\n', out_file);
end
