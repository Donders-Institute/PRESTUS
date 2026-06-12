function sensor_data = babelbrain_wrapper( ...
        parameters, medium_masks, kwave_medium, focus_pos, trans_pos, planimg, equip_param, pseudoCT)
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
%                         'auto' (default): detect Metal on macOS, CUDA on Linux
%                         with nvidia-smi, else CPU.
%   .device_name        - device name string, e.g. 'Apple' or 'RTX3090'.
%                         'auto' (default): derived from computing_backend.
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
%   .ct_mapping         - BabelBrain HU→property mapping method. Non-empty activates the CT
%                         path (bUseCT=True) and passes the pseudoCT HU volume to BabelBrain.
%                         Valid values: 'Webb-Marsac' (default when non-empty), 'Aubry',
%                         'Pichardo', 'McDannold', 'Marsac-Aubry', 'Pichardo-Marsac',
%                         'McDannold-Marsac'. Empty (default): use label-based MatFreq path.
%                         Note: only set this when passing raw HU values from a UTE/PETRA scan
%                         for BabelBrain to map itself. If PRESTUS's pCT pipeline has already
%                         converted HU to acoustic properties, use the default label path.
%   .ct_use_petra       - logical; set true when the HU values in pseudoCT were derived from a
%                         PETRA sequence (uses Björn's density formula in BabelBrain).
%                         Default false (appropriate for ZTE/UTE-derived pCT).
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
%   pseudoCT     - (optional) [Nx x Ny x Nz] Hounsfield-unit pseudo-CT in
%                  grid space (from grid_tissue_setup). Required when
%                  parameters.babelbrain.ct_mapping is non-empty.
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
    pseudoCT           {mustBeNumeric} = []
end

%% Resolve BabelBrain config with defaults
bb = struct();
if isfield(parameters, 'babelbrain')
    bb = parameters.babelbrain;
end

bb = apply_default(bb, 'tx_system',         '');
bb = apply_default(bb, 'computing_backend',  'auto');
bb = apply_default(bb, 'device_name',        'auto');
bb = apply_default(bb, 'ppw',                6);
bb = apply_default(bb, 'do_refocusing',      false);
bb = apply_default(bb, 'z_beyond_focus_mm',  40);
bb = apply_default(bb, 'python_env',         '');
bb = apply_default(bb, 'tx_mech_adjust_z_m', 0.0);  % 0: geometry already in NIfTI
bb = apply_default(bb, 'z_steering_m',       []);   % empty = auto-derived below
bb = apply_default(bb, 'root',               '');
bb = apply_default(bb, 'ct_mapping',         '');   % empty = label-based path
bb = apply_default(bb, 'ct_use_petra',       false);

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
python_bin = resolve_python_env(bb.python_env, 'python');

%% Ensure HDF5 plugin path is valid (pass python_bin to locate blosc filter)
ensure_hdf5_plugin_path(python_bin);

%% Query BabelBrain version (before any simulation work)
bb_version = query_babelbrain_version(python_bin, bb_runner, bb_root);
fprintf('[babelbrain_wrapper] BabelBrain %s  (%s)\n', ...
    bb_version.version_txt, bb_version.git_describe);

%% I/O directories
out_dir   = parameters.io.dir_cache;
if out_dir(end) ~= filesep, out_dir(end+1) = filesep; end
subj_id   = sprintf('sub-%03d', parameters.subject_id);
target    = strrep(parameters.io.output_affix, '_', '');
if isempty(target)
    target = 'target';
end

%% Multi-transducer coherent combination
% When more than one transducer is present and async mode is not active,
% simulate each transducer independently and sum the complex Fourier pressure
% fields (coherent superposition). Each recursive call reduces to the
% single-transducer path below.
n_tx = numel(parameters.transducer);
is_async = isfield(parameters, 'simulation') && ...
           isfield(parameters.simulation, 'transducer_coupling') && ...
           strcmp(parameters.simulation.transducer_coupling, 'async');

if n_tx > 1 && ~is_async
    fprintf('[babelbrain_wrapper] %d transducers — running coherent combination.\n', n_tx);
    p_combined  = [];
    last_sd     = struct();
    for ti = 1:n_tx
        fprintf('[babelbrain_wrapper] === Transducer %d / %d ===\n', ti, n_tx);
        p_ti           = parameters;
        p_ti.transducer = parameters.transducer(ti);
        % Prevent re-entry on recursive call.
        p_ti.simulation.transducer_coupling = 'coherent';
        % Give each sub-run a unique affix so output files don't collide.
        base_affix = parameters.io.output_affix;
        p_ti.io.output_affix = [base_affix, trans_suffix(ti, n_tx)];
        focus_ti = parameters.transducer(ti).focus_pos;
        trans_ti = parameters.transducer(ti).trans_pos;
        sd_ti = babelbrain_wrapper(p_ti, medium_masks, kwave_medium, ...
            focus_ti, trans_ti, planimg, equip_param, pseudoCT);
        if ~isfield(sd_ti, 'p_complex') || isempty(sd_ti.p_complex)
            error(['babelbrain_wrapper: p_complex missing from transducer %d output.\n' ...
                'Ensure the BabelBrain DataForSim.h5 contains /p_complex.'], ti);
        end
        % Shape consistency across transducers
        if ~isempty(p_combined) && ~isequal(size(sd_ti.p_complex), size(p_combined))
            error(['babelbrain_wrapper: p_complex size mismatch at transducer %d.\n' ...
                'Got [%d %d %d], expected [%d %d %d].\n' ...
                'All transducers must use the same simulation grid.'], ...
                ti, size(sd_ti.p_complex,1), size(sd_ti.p_complex,2), size(sd_ti.p_complex,3), ...
                size(p_combined,1), size(p_combined,2), size(p_combined,3));
        end
        % Finite-value check
        n_bad = nnz(~isfinite(sd_ti.p_complex(:)));
        if n_bad > 0
            error('babelbrain_wrapper: p_complex from transducer %d contains %d NaN/Inf values.', ti, n_bad);
        end
        if isempty(p_combined)
            p_combined = sd_ti.p_complex;
        else
            p_combined = p_combined + sd_ti.p_complex;
        end
        last_sd = sd_ti;
    end
    % Final combined field check
    n_bad = nnz(~isfinite(p_combined(:)));
    if n_bad > 0
        error('babelbrain_wrapper: combined p_complex contains %d NaN/Inf values after summing %d transducers.', n_bad, n_tx);
    end
    sensor_data                    = last_sd;
    sensor_data.p_complex          = p_combined;
    sensor_data.p_max_all          = abs(p_combined);
    sensor_data.coherent_n_tx      = n_tx;
    return
end

tr = parameters.transducer(1);

%% Resolve frequency and PPW (needed for NIfTI filenames in Steps 2–3b)
freq_hz  = tr.freq_hz;
res_mm   = parameters.grid.resolution_mm;
c_min    = 1500;
if isfield(parameters, 'medium_properties') && isfield(parameters.medium_properties, 'brain')
    c_min = parameters.medium_properties.brain.sound_speed;
end
ppw_auto = round(c_min / (freq_hz * res_mm * 1e-3));
ppw_use  = max(bb.ppw, ppw_auto);

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
nii_path = write_babelbrain_input(medium_masks, focus_pos, parameters, planimg, out_dir, target, ppw_use);

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

%% Step 3b — Write pseudo-CT NIfTI for BabelBrain CT path (optional)
ct_nii_path = '';
if ~isempty(bb.ct_mapping)
    if isempty(pseudoCT)
        error(['babelbrain_wrapper: ct_mapping is set but no pseudoCT was supplied.\n' ...
            'Pass the pseudoCT output of grid_tissue_setup as the 8th argument.']);
    end
    % BabelBrain expects CT.nii.gz alongside BabelViscoInput.nii.gz with the
    % same frequency/PPW naming prefix.
    ct_fname = sprintf('%s_%dkHz_%dPPW_CT.nii.gz', target, round(freq_hz/1e3), ppw_use);
    ct_nii_path = fullfile(out_dir, subj_id, ct_fname);
    if ~isfolder(fileparts(ct_nii_path))
        mkdir(fileparts(ct_nii_path));
    end
    niftiwrite(single(pseudoCT), strrep(ct_nii_path, '.nii.gz', ''), 'Compressed', true);
    fprintf('[babelbrain_wrapper] CT NIfTI written: %s\n', ct_nii_path);
end

%% Step 4 — Invoke Python bridge
[bb.computing_backend, bb.device_name] = resolve_bb_backend(bb.computing_backend, bb.device_name);

% Checkpoint: verify input files are present and log their sizes
d_nii = dir(nii_path);
if isempty(d_nii)
    error('babelbrain_wrapper: NIfTI input not found before launching simulation:\n  %s', nii_path);
end
fprintf('[babelbrain_wrapper] Input NIfTI:   %s  (%.1f MB)\n', nii_path, d_nii.bytes/1e6);
if ~isempty(weights_path) && isfile(weights_path)
    d_wt = dir(weights_path);
    fprintf('[babelbrain_wrapper] Weights HDF5:  %s  (%.1f MB)\n', weights_path, d_wt.bytes/1e6);
end

% Per-run log file (Python stdout+stderr appended here for post-hoc debugging)
log_dir = fullfile(out_dir, subj_id, 'log');
if ~isfolder(log_dir), mkdir(log_dir); end
log_file = fullfile(log_dir, sprintf('%s_%s_%dkHz_babelbrain_acoustic.log', subj_id, target, round(freq_hz/1e3)));
fprintf('[babelbrain_wrapper] Log file: %s\n', log_file);

fprintf('[babelbrain_wrapper] Launching BabelBrain simulation...\n');

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
    '--log_file "%s" '  ...
    '%s '  ...  % optional babelbrain_root flag
    '%s '  ...  % optional weights flag
    '%s '  ...  % optional refocusing flag
    '%s '  ...  % optional dry_run flag
    '%s '  ...  % optional use_ct flag
    '%s '  ...  % optional ct_mapping flag
    '%s '], ...  % optional use_petra flag
    python_bin, bb_runner, ...
    out_dir, subj_id, target, bb.tx_system, ...
    freq_hz, ppw_use, bb.computing_backend, bb.device_name, ...
    bb.tx_mech_adjust_z_m, z_steering_m, bb.z_beyond_focus_mm, ...
    log_file, ...
    bb_root_flag(bb_root), ...
    bb_weights_flag(weights_path), ...
    bb_refocus_flag(bb.do_refocusing), ...
    bb_dryrun_flag(dry_run), ...
    bb_usect_flag(~isempty(bb.ct_mapping)), ...
    bb_ctmapping_flag(bb.ct_mapping), ...
    bb_petra_flag(bb.ct_use_petra));

fprintf('[babelbrain_wrapper] Command:\n  %s\n', cmd_args);

t0     = tic;
status = system(cmd_args);
elapsed = toc(t0);

if status ~= 0
    error('babelbrain_wrapper: BabelBrain Python process exited with status %d.\nCommand was:\n  %s', ...
        status, cmd_args);
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

rewrite_hdf5_no_blosc(h5_path, python_bin, this_dir);

sensor_data = read_babelbrain_output(h5_path, bb.do_refocusing, parameters, planimg);
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

function [backend, device] = resolve_bb_backend(backend, device)
% Resolve 'auto' computing_backend / device_name to concrete values.
% Priority: Metal (macOS) > CUDA (nvidia-smi) > CPU.
    if ~(ischar(backend) && strcmp(backend, 'auto'))
        % Already a numeric backend — only resolve device if still 'auto'
        if ischar(device) && strcmp(device, 'auto')
            device = default_device_for_backend(backend);
        end
        return
    end

    if ismac
        backend = 3;   % Metal
        device  = 'Apple';
        fprintf('[babelbrain_wrapper] Auto-detected backend: Metal (macOS)\n');
        return
    end

    % Linux/Windows: probe for CUDA
    [status, out] = system('nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null');
    if status == 0
        lines = strtrim(strsplit(strtrim(out), newline));
        lines = lines(~cellfun(@isempty, lines));
        if ~isempty(lines)
            backend = 1;   % CUDA
            device  = lines{1};
            fprintf('[babelbrain_wrapper] Auto-detected backend: CUDA, device="%s"\n', device);
            return
        end
    end

    % Fallback: CPU — warn so users on GPU machines notice driver/detection issues
    backend = 4;
    device  = 'CPU';
    warning('babelbrain_wrapper:cpuFallback', ...
        '[babelbrain_wrapper] No GPU detected (nvidia-smi failed or returned no devices); falling back to CPU. Simulation may be slow.');
end

function device = default_device_for_backend(backend)
% Return a safe default device string for a known numeric backend.
    switch backend
        case 3;  device = 'Apple';
        otherwise; device = 'CPU';
    end
end

function flag = bb_dryrun_flag(dry_run)
    if dry_run
        flag = '--dry_run';
    else
        flag = '';
    end
end

function flag = bb_usect_flag(use_ct)
    if use_ct
        flag = '--use_ct';
    else
        flag = '';
    end
end

function flag = bb_ctmapping_flag(ct_mapping)
    if ~isempty(ct_mapping)
        flag = sprintf('--ct_mapping "%s"', ct_mapping);
    else
        flag = '';
    end
end

function flag = bb_petra_flag(use_petra)
    if use_petra
        flag = '--use_petra';
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

function rewrite_hdf5_no_blosc(h5_path, python_bin, this_dir)
% Rewrite h5_path in-place replacing blosc compression with gzip so MATLAB
% can read it without a blosc filter plugin.
rewriter = fullfile(this_dir, 'rewrite_hdf5_no_blosc.py');
if ~isfile(rewriter)
    warning('babelbrain_wrapper: rewrite_hdf5_no_blosc.py not found at %s — blosc read may fail.', rewriter);
    return
end
cmd = sprintf('"%s" "%s" "%s"', python_bin, rewriter, h5_path);
[status, out] = system(cmd);
if status ~= 0
    warning('babelbrain_wrapper: HDF5 rewrite failed (exit %d):\n%s\n  h5read may fail if blosc plugin is missing.', status, out);
end
end
