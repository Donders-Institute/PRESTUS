function sensor_data = babelbrain_wrapper( ...
        parameters, medium_masks, kwave_medium, focus_pos, trans_pos, planimg, equip_param)
% BABELBRAIN_WRAPPER  Run BabelBrain acoustic simulation as PRESTUS backend
%
% Drop-in replacement for ACOUSTIC_WRAPPER that uses BabelBrain's
% BabelViscoFDTD solver instead of k-Wave. Accepts the same pre-processed
% PRESTUS state (medium masks, transducer parameters, focus position) and
% returns sensor_data in the same format as acoustic_wrapper so that all
% downstream stages (acoustic_analysis, thermal_simulation) remain unchanged.
%
% Pipeline:
%   1. Map PRESTUS transducer → BabelBrain preset (babelbrain_transducer_map)
%   2. Write BabelViscoInput NIfTI mask (write_babelbrain_input)
%   3. Write OptimizedWeightsFile HDF5 (export_babelbrain_weights)
%   4. Invoke BabelBrain Python bridge via system() call
%   5. Read *DataForSim.h5 output (read_babelbrain_output)
%
% Required config fields (parameters.babelbrain):
%   .tx_system          - BabelBrain preset name, e.g. 'CTX_500'.
%                         If empty, babelbrain_transducer_map is called to
%                         auto-detect. Set explicitly when using a custom
%                         RUN_SIM subclass.
%   .computing_backend  - 1=CUDA, 2=OpenCL, 3=Metal (Apple Silicon), 4=CPU
%   .device_name        - device name string, e.g. 'Apple' or 'RTX3090'
%   .ppw                - points-per-wavelength (default 6)
%   .do_refocusing      - logical; use BabelBrain skull refocusing (default false)
%   .z_beyond_focus_mm  - axial truncation past focus (default 40, mm)
%   .python_env         - path to Python venv/conda env with BabelBrain installed;
%                         empty string = use system Python
%   .tx_mech_adjust_z_m - transducer mechanical Z offset into skin [m] (default 0: placement
%                         geometry is already encoded in BabelViscoInput.nii.gz)
%   .z_steering_m       - ZSteering [m]: axial offset of electronic focus from preset FocalLength.
%                         Empty (default): auto-derived as focal_distance_mm - FocalLength when
%                         no calibrated weights file is present; 0 when weights are provided.
%
% Use as:
%   sensor_data = babelbrain_wrapper(parameters, medium_masks, kwave_medium, ...
%                     focus_pos, trans_pos, planimg)
%   sensor_data = babelbrain_wrapper(..., equip_param)
%
% Input:
%   parameters   - PRESTUS config struct with .babelbrain sub-struct
%   medium_masks - [Nx x Ny x Nz] medium label volume (from preproc_medium_mask)
%   kwave_medium - k-Wave medium struct (not used for simulation; kept for
%                  interface compatibility with acoustic_wrapper)
%   focus_pos    - [1x3] focus voxel in grid coordinates
%   trans_pos    - [1x3] transducer voxel in grid coordinates (for provenance)
%   planimg      - planimg struct from grid_tissue_setup (affine info)
%   equip_param  - (optional) equipment config struct for transducer lookup
%
% Output:
%   sensor_data  - struct with p_max_all [Nx x Ny x Nz] [Pa], matching the
%                  format returned by acoustic_wrapper.
%
% See also: ACOUSTIC_WRAPPER, WRITE_BABELBRAIN_INPUT, BABELBRAIN_TRANSDUCER_MAP,
%           READ_BABELBRAIN_OUTPUT, EXPORT_BABELBRAIN_WEIGHTS

arguments
    parameters   (1,1) struct
    medium_masks (:,:,:) {mustBeNumeric}
    kwave_medium (1,1) struct  %#ok<INUSA>  kept for interface symmetry
    focus_pos    (1,3) {mustBeNumeric}
    trans_pos    (1,3) {mustBeNumeric}  %#ok<INUSA>
    planimg      (1,1) struct
    equip_param        struct  = struct()
end

%% Resolve BabelBrain config with defaults
bb = struct();
if isfield(parameters, 'babelbrain')
    bb = parameters.babelbrain;
end

bb = apply_default(bb, 'tx_system',         '');
bb = apply_default(bb, 'computing_backend',  3);    % Metal
bb = apply_default(bb, 'device_name',        'Apple');
bb = apply_default(bb, 'ppw',                6);
bb = apply_default(bb, 'do_refocusing',      false);
bb = apply_default(bb, 'z_beyond_focus_mm',  40);
bb = apply_default(bb, 'python_env',         '');
bb = apply_default(bb, 'tx_mech_adjust_z_m', 0.0);  % 0: geometry already in NIfTI
bb = apply_default(bb, 'z_steering_m',       []);   % empty = auto-derived below
bb = apply_default(bb, 'root',               '');

%% Locate BabelBrain Python package and runner script
% The runner script lives in functions/babelbrain/ alongside this file.
% BabelBrain source (external/BabelBrain submodule) is located separately.
prestus_root = get_prestus_path();
this_dir     = fileparts(mfilename('fullpath'));
bb_runner    = fullfile(this_dir, 'run_babelbrain_case.py');

if ~isfile(bb_runner)
    error(['babelbrain_wrapper: run_babelbrain_case.py not found at:\n  %s\n' ...
        'This file should be in functions/babelbrain/ alongside babelbrain_wrapper.m.'], ...
        bb_runner);
end

% BabelBrain source root (contains TranscranialModeling/)
bb_search = {
    bb.root, ...
    fullfile(prestus_root, 'external', 'BabelBrain'), ...
};

bb_root = '';
for k = 1:numel(bb_search)
    cand = bb_search{k};
    if isempty(cand); continue; end
    if exist(fullfile(cand, 'TranscranialModeling'), 'dir')
        bb_root = cand;
        break
    end
end

if isempty(bb_root)
    error(['babelbrain_wrapper: BabelBrain source not found.\n' ...
        'Initialise the submodule:\n' ...
        '  git submodule update --init external/BabelBrain\n' ...
        'Or set parameters.babelbrain.root explicitly.']);
end


%% Select Python executable
if ~isempty(bb.python_env)
    python_bin = fullfile(bb.python_env, 'bin', 'python');
    if ~isfile(python_bin)
        python_bin = fullfile(bb.python_env, 'python');  % Windows / conda root
    end
else
    python_bin = 'python3';
end

%% Query BabelBrain version (before any simulation work)
bb_version = query_babelbrain_version(python_bin, bb_runner, bb_root);
fprintf('[babelbrain_wrapper] BabelBrain %s  (%s)\n', ...
    bb_version.version_txt, bb_version.git_describe);

%% I/O directories
out_dir   = parameters.io.dir_cache;
subj_id   = sprintf('sub-%03d', parameters.subject_id);
target    = strrep(parameters.io.output_affix, '_', '');
if isempty(target)
    target = 'target';
end

tr = parameters.transducer(1);

%% Step 1 — Transducer preset
if isempty(bb.tx_system)
    [bb.tx_system, csv_path] = babelbrain_transducer_map(parameters, equip_param, out_dir);
    if isempty(bb.tx_system)
        error(['babelbrain_wrapper: No BabelBrain preset matched and auto-detection ' ...
            'returned empty preset.\n' ...
            'Set parameters.babelbrain.tx_system explicitly or register the ' ...
            'transducer as described in babelbrain_transducer_map.']);
    end
    if ~isempty(csv_path)
        fprintf(['[babelbrain_wrapper] Custom element CSV generated.\n' ...
            '  You must register the preset in BabelBrain before proceeding.\n' ...
            '  CSV: %s\n'], csv_path);
        error('babelbrain_wrapper: Custom preset registration required. See above.');
    end
end
fprintf('[babelbrain_wrapper] Using preset: %s\n', bb.tx_system);

%% Derive ZSteering from focal distance vs preset FocalLength
% ZSteering shifts BabelBrain's Rayleigh pre-step target axially (metres).
% When calibrated phases are provided, the phases already encode any steering
% offset, so ZSteering should stay 0. When no weights file is available,
% derive ZSteering from the difference between the PRESTUS focal distance and
% the preset's natural focal length so that geometric focusing is correct.
if isempty(bb.z_steering_m)
    has_calib_weights = isfield(tr, 'annular') && isfield(tr.annular, 'elem_phase_rad') ...
        && ~isempty(tr.annular.elem_phase_rad);
    if ~has_calib_weights
        % Attempt to read FocalLength from the preset default.yaml
        preset_yaml = fullfile(bb_root, sprintf('Babel_%s', bb.tx_system), 'default.yaml');
        preset_focal_length_m = [];
        if isfile(preset_yaml)
            try
                preset_cfg = read_yaml(preset_yaml);
                if isfield(preset_cfg, 'FocalLength')
                    preset_focal_length_m = preset_cfg.FocalLength;
                end
            catch
            end
        end
        if ~isempty(preset_focal_length_m) && isfield(tr, 'focal_distance_mm')
            focal_dist_m = tr.focal_distance_mm / 1000;
            z_steering_m = focal_dist_m - preset_focal_length_m;
            fprintf('[babelbrain_wrapper] ZSteering=%.4f m (focal_dist=%.3f m, preset FL=%.3f m)\n', ...
                z_steering_m, focal_dist_m, preset_focal_length_m);
        else
            z_steering_m = 0.0;
        end
    else
        z_steering_m = 0.0;  % calibrated phases absorb the steering delta
    end
else
    z_steering_m = bb.z_steering_m;
end

%% Step 2 — Write BabelViscoInput NIfTI
fprintf('[babelbrain_wrapper] Writing input mask...\n');
nii_path = write_babelbrain_input(medium_masks, focus_pos, parameters, planimg, out_dir, target);

%% Step 3 — Write OptimizedWeightsFile (calibrated phases)
weights_path = fullfile(out_dir, subj_id, [target '_CalibWeights.h5']);

if isfield(tr, 'annular') && isfield(tr.annular, 'elem_phase_rad') && isfield(tr.annular, 'elem_amp')
    phases_rad = tr.annular.elem_phase_rad;
    elem_amp   = tr.annular.elem_amp(1);

    % Geometric phases at the calibration focal depth for delta computation.
    % If available from the calibration pipeline, pass them through.
    geo_phases_rad = [];
    if isfield(tr.annular, 'geo_phase_rad')
        geo_phases_rad = tr.annular.geo_phase_rad;
    end

    export_babelbrain_weights(phases_rad, elem_amp, parameters, weights_path, geo_phases_rad);
else
    weights_path = '';
    fprintf('[babelbrain_wrapper] No element phases found — BabelBrain will use geometric focusing.\n');
end

%% Step 4 — Invoke Python bridge
fprintf('[babelbrain_wrapper] Launching BabelBrain simulation...\n');

freq_hz  = tr.freq_hz;
res_mm   = parameters.grid.resolution_mm;
c_min    = 1500;
if isfield(parameters, 'medium_properties') && isfield(parameters.medium_properties, 'brain')
    c_min = parameters.medium_properties.brain.sound_speed;
end
ppw_auto = round(c_min / (freq_hz * res_mm * 1e-3));
ppw_use  = max(bb.ppw, ppw_auto);  % never go below geometry-required PPW

dry_run = isfield(bb, 'dry_run') && bb.dry_run;

cmd_args = sprintf([ ...
    '"%s" "%s" '  ...
    '--basedir "%s" '  ...
    '--id "%s" '  ...
    '--target "%s" '  ...
    '--tx_system "%s" '  ...
    '--freq %.0f '  ...
    '--ppw %d '  ...
    '--computing_backend %d '  ...
    '--device_name "%s" '  ...
    '--tx_mech_adjust_z %.6f '  ...
    '--z_steering %.6f '  ...
    '--z_beyond_focus_mm %.1f '  ...
    '%s '  ...  % optional babelbrain_root flag
    '%s '  ...  % optional weights flag
    '%s '  ...  % optional refocusing flag
    '%s '], ...  % optional dry_run flag
    python_bin, bb_runner, ...
    out_dir, subj_id, target, bb.tx_system, ...
    freq_hz, ppw_use, bb.computing_backend, bb.device_name, ...
    bb.tx_mech_adjust_z_m, z_steering_m, bb.z_beyond_focus_mm, ...
    bb_root_flag(bb_root), ...
    bb_weights_flag(weights_path), ...
    bb_refocus_flag(bb.do_refocusing), ...
    bb_dryrun_flag(dry_run));

fprintf('[babelbrain_wrapper] Command:\n  %s\n', cmd_args);

t0     = tic;
status = system(cmd_args);
elapsed = toc(t0);

if status ~= 0
    error('babelbrain_wrapper: BabelBrain Python process exited with status %d.', status);
end
fprintf('[babelbrain_wrapper] BabelBrain finished in %.1f s.\n', elapsed);

%% Step 5 — Locate and read output HDF5
if dry_run
    fprintf('[babelbrain_wrapper] Dry run — skipping HDF5 output read.\n');
    sensor_data = struct('p_max_all', [], 'babelbrain_source', 'dry_run', ...
                         'babelbrain_version', bb_version);
    return
end

% BabelBrain names output: <basedir>/<ID>/<target>_<freq>kHz_<PPW>PPW_DataForSim.h5
freq_kHz = round(freq_hz / 1e3);
h5_name  = sprintf('%s_%dkHz_%dPPW_DataForSim.h5', target, freq_kHz, ppw_use);
h5_path  = fullfile(out_dir, subj_id, h5_name);

if ~isfile(h5_path)
    % Try glob in case BabelBrain version uses a slightly different naming
    candidates = dir(fullfile(out_dir, subj_id, '*DataForSim.h5'));
    if isempty(candidates)
        error('babelbrain_wrapper: Output HDF5 not found.\n  Expected: %s', h5_path);
    end
    h5_path = fullfile(candidates(end).folder, candidates(end).name);
    fprintf('[babelbrain_wrapper] Using output: %s\n', h5_path);
end

sensor_data = read_babelbrain_output(h5_path, bb.do_refocusing);
sensor_data.babelbrain_version = bb_version;

end

%% Local helpers

function s = apply_default(s, field, val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = val;
    end
end

function flag = bb_weights_flag(weights_path)
    if ~isempty(weights_path) && isfile(weights_path)
        flag = sprintf('--weights_file "%s"', weights_path);
    else
        flag = '';
    end
end

function flag = bb_root_flag(bb_root)
    if ~isempty(bb_root)
        flag = sprintf('--babelbrain_root "%s"', bb_root);
    else
        flag = '';
    end
end

function flag = bb_refocus_flag(do_refocus)
    if do_refocus
        flag = '--do_refocusing';
    else
        flag = '';
    end
end

function flag = bb_dryrun_flag(dry_run)
    if dry_run
        flag = '--dry_run';
    else
        flag = '';
    end
end

function ver = query_babelbrain_version(python_bin, bb_runner, bb_root)
% Call run_babelbrain_case.py --get_version; parse KEY=VALUE stdout lines.
    ver = struct( ...
        'version_txt',  '', ...
        'git_describe', '', ...
        'git_commit',   '', ...
        'bb_root',      bb_root);

    root_flag = '';
    if ~isempty(bb_root)
        root_flag = sprintf('--babelbrain_root "%s"', bb_root);
    end
    cmd = sprintf('"%s" "%s" --get_version %s 2>/dev/null', python_bin, bb_runner, root_flag);
    [status, out] = system(cmd);
    if status ~= 0
        fprintf('[babelbrain_wrapper] Warning: could not query BabelBrain version.\n');
        return
    end

    lines = strsplit(strtrim(out), newline);
    for k = 1:numel(lines)
        parts = regexp(lines{k}, '^(\w+)=(.*)$', 'tokens', 'once');
        if numel(parts) == 2
            switch parts{1}
                case 'BABELBRAIN_VERSION_TXT';  ver.version_txt  = strtrim(parts{2});
                case 'BABELBRAIN_GIT_DESCRIBE';  ver.git_describe = strtrim(parts{2});
                case 'BABELBRAIN_GIT_COMMIT';    ver.git_commit   = strtrim(parts{2});
                case 'BABELBRAIN_ROOT';          ver.bb_root      = strtrim(parts{2});
            end
        end
    end
end
