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
%   4. Invoke BabelBrain Python bridge via run_command()
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
%   .ppw                - BabelBrain points-per-wavelength (default 12, matching
%                         the BabelViscoFDTD benchmark/recommended fidelity;
%                         BabelBrain's native step = SmallestSOS/(freq*ppw)).
%                         Lower values (e.g. 6) run faster but add numerical
%                         dispersion. The mask is upsampled to this native step
%                         regardless, so ppw need not track grid.resolution_mm.
%   .do_refocusing      - logical; use BabelBrain skull refocusing (default false)
%   (domain size)       - controlled by parameters.grid.acoustic_fov_length_mm (not babelbrain.*):
%                         when set, BabelBrain simulates the same FOV-cropped box as k-Wave/jwave.
%   .tight_narrow_beam  - logical; use BabelBrain's narrow-beam domain (default false). Saves
%                         compute on full-brain runs, but the pre-focal field is not returned.
%                         Ignored when a PRESTUS FOV crop is active.
%   .python_env         - path to Python venv/conda env with BabelBrain installed;
%                         empty string = use system Python
%   .tx_mech_adjust_z_m - additional transducer mechanical Z offset into skin [m] (default 0).
%                         A placement correction of (preset FocalLength - |focus_pos - trans_pos|)
%                         is ALWAYS added on top of this value: BabelBrain builds its domain
%                         with the apex at preset-FL distance from the target voxel, whereas
%                         PRESTUS places it at trans_pos — without the correction the
%                         resampled field is rigidly shifted by the difference. Falls back to
%                         focal_distance_bowl when trans_pos is unavailable.
%   .z_steering_m       - ZSteering [m]: axial offset of the electronic aim from the preset
%                         FocalLength (i.e. aim_from_apex = preset_FL + ZSteering).
%                         Empty (default): auto-derived — 0 when explicit element phases exist
%                         (the weights file carries the full drive); otherwise
%                         aim_bowl - preset_FL (aim = calibrated steer D* or focal_distance_bowl).
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
%   .water_amp_calibration - logical (default true). Scale the output pressure by
%                         k = sqrt(target_isppa_wcm2 / achieved free-water ISPPA), so the
%                         free-water ISPPA matches the transducer's intended value. Applied only
%                         when a BabelBrain free-water run exists; k is stored in
%                         sensor_data.amp_calibration.
%   .use_rayleigh_for_water - logical (default false). true: compute water-only runs with
%                         BabelBrain's Rayleigh integral (faster); false: full-wave FDTD, which
%                         is comparable to k-Wave/jwave. Tissue runs always use FDTD.
%   .focal_calibration  - logical (default true). Without explicit element phases, choose the
%                         geometric steer so the achieved focus lands at focal_distance_ep
%                         (CALIBRATE_FOCAL_STEERING); with explicit phases, only report it.
%   .focal_calibration_target - 'flhm_center' (default) | 'peak'; quantity placed at
%                         focal_distance_ep by focal_calibration.
%   .profile_calibration - logical (default false). Re-fit annular drive phases to the whole
%                         axial profile (CALIBRATE_BABELBRAIN_PHASE_PROFILE).
%   .root               - BabelBrain source root; empty (default) = auto-detect.
%   .cleanup_cache      - logical (default true); delete BabelBrain's cache folder after reading.
%   .source_amp_pa      - BabelBrain Rayleigh surface pressure amplitude [Pa].
%                         Empty (default): auto-populated from elem_p_amp when available,
%                         because elem_p_amp = rho_water * c_water * v_opt is the physical
%                         face pressure derived by the O'Neill/Rayleigh calibration — the
%                         same quantity BabelBrain expects here. This gives consistent
%                         free-water focal Isppa between kWave and BabelBrain at the
%                         calibrated depth without any additional cross-calibration step.
%                         Set explicitly to override the auto-populated value.
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
%   trans_pos    - [1x3] transducer (bowl rear-centre / apex) voxel in grid coordinates;
%                  sets BabelBrain's Tx placement and the beam direction
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
    trans_pos    (1,3) {mustBeNumeric}
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
bb = apply_default(bb, 'ppw',                12);
bb = apply_default(bb, 'do_refocusing',      false);
bb = apply_default(bb, 'tight_narrow_beam',  []);   % [] = auto (off unless no FOV crop AND explicitly requested); see below
% Focal calibration: choose the geometric steer so the achieved FWHM-centre
% (or peak) lands at focal_distance_ep, rather than steering geometrically to
% focal_distance_ep directly (which places the peak proximal via the diffraction
% focal shift).  See calibrate_focal_steering.m.  Always on by default; the flag
% is an escape hatch.
bb = apply_default(bb, 'focal_calibration',        true);
bb = apply_default(bb, 'focal_calibration_target', 'flhm_center');  % 'flhm_center' | 'peak'
% Whole-profile phase calibration (annular explicit-phase drives only): re-fit
% the BabelBrain drive phases so the discrete-ring Rayleigh field matches the
% continuous-bowl O'Neil ideal across the whole axial profile, rather than only
% placing the FLHM centre. Corrects the ring-gap discretisation bias in
% BabelBrain's Rayleigh water field. See calibrate_babelbrain_phase_profile.m.
bb = apply_default(bb, 'profile_calibration',      false);

% Domain size: derived from parameters.grid.acoustic_fov_length_mm (not user-facing
% on the babelbrain sub-struct).
%
% BabelBrain's narrow-beam mode (bTightNarrowBeamDomain) is off by default: it
% solves only the region around and beyond the focus and does not emit the
% pre-focal field, so the returned axial profile would be zero up to ~the
% focus. Enable it with babelbrain.tight_narrow_beam = true to save compute on
% full-brain runs. It stays off whenever a PRESTUS FOV crop is active, since the
% mask is then already the beam-centred box used for k-Wave/jwave and a second
% crop can cut away skull near the transducer. See doc_babelbrain.md.
fov_len_mm = [];
if isfield(parameters, 'grid') && isfield(parameters.grid, 'acoustic_fov_length_mm') && ...
        isnumeric(parameters.grid.acoustic_fov_length_mm)
    fov_len_mm = parameters.grid.acoustic_fov_length_mm;
end

if ~isempty(bb.tight_narrow_beam)
    tight_narrow_beam = logical(bb.tight_narrow_beam) && isempty(fov_len_mm);
    if logical(bb.tight_narrow_beam) && ~isempty(fov_len_mm)
        fprintf(['[babelbrain_wrapper] babelbrain.tight_narrow_beam requested but a PRESTUS FOV crop ' ...
            'is active (fov=%.1f mm) -- ignoring (narrow-beam would double-crop).\n'], fov_len_mm);
    end
else
    tight_narrow_beam = false;   % default: full (non-narrow) domain, so the pre-focal field is emitted
end

z_beyond_focus_mm = 40;
if tight_narrow_beam
    fprintf(['[babelbrain_wrapper] BabelBrain narrow-beam domain ENABLED (babelbrain.tight_narrow_beam=true): ' ...
        'z_beyond_focus=%.1f mm. NOTE: the pre-focal field is not emitted -- the returned axial profile will be ' ...
        'zero from the transducer to ~the focus.\n'], z_beyond_focus_mm);
elseif isempty(fov_len_mm)
    fprintf(['[babelbrain_wrapper] Full (non-narrow) BabelBrain domain, no PRESTUS FOV crop -- BabelBrain sizes the ' ...
        'domain from the full input mask (may be large). Set a grid.acoustic_fov_length_mm crop to bound it, or ' ...
        'babelbrain.tight_narrow_beam=true to use BabelBrain''s narrow domain (drops the pre-focal field).\n']);
else
    fprintf(['[babelbrain_wrapper] PRESTUS-side FOV crop active (fov=%.1f mm) -- using BabelBrain''s regular ' ...
        '(non-narrow) domain sizing so it matches the already-cropped mask, with no further BabelBrain-side shrink.\n'], fov_len_mm);
end
bb = apply_default(bb, 'python_env',         '');
bb = apply_default(bb, 'tx_mech_adjust_z_m', 0.0);  % 0: geometry already in NIfTI
bb = apply_default(bb, 'z_steering_m',       []);   % empty = auto-derived below
bb = apply_default(bb, 'root',               '');
bb = apply_default(bb, 'ct_mapping',         '');   % empty = label-based path
bb = apply_default(bb, 'ct_use_petra',       false);
bb = apply_default(bb, 'cleanup_cache',      true);   % delete BabelBrain cache dir after reading output
bb = apply_default(bb, 'source_amp_pa',     []);    % [] = use BabelBrain built-in default (~50 kPa)
bb = apply_default(bb, 'water_amp_calibration', true);  % scale output to the intended free-water ISPPA
% Water-field solver. false (default) => full-wave FDTD for the water-only field,
% matching the other PRESTUS backends (k-Wave/jwave) so BabelBrain is comparable
% on a like-for-like basis. true => BabelBrain's native Rayleigh–Sommerfeld
% integral (faster, exact linear free-field, but focuses slightly tighter/deeper
% than a gridded full-wave solve → reads as a rigid axial offset + higher peak
% vs the other backends). Only affects water-only runs; the transcranial field
% always uses FDTD.
bb = apply_default(bb, 'use_rayleigh_for_water', false);

% Water-medium mode: run BabelBrain's --water_only path as the primary simulation.
% medium_masks will be all-water (zeros in BabelBrain label space); the NIfTI is
% still written for transducer geometry, but tissue acoustics are bypassed.
is_water_medium = isfield(parameters, 'simulation') && ...
    isfield(parameters.simulation, 'medium') && ...
    strcmpi(parameters.simulation.medium, 'water');

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

% BabelBrain's native step is SmallestSOS/(freq*PPW), where SmallestSOS is
% normally its Pichardo shear-speed floor (1102.515 m/s; same formula as in
% write_babelbrain_input.m). ppw_auto makes that step match grid.resolution_mm;
% use it or babelbrain.ppw, whichever is finer. write_babelbrain_input resamples
% the mask to the resulting step.
pichardo_shear_floor_ms = 1000.0 * 0.422 + 680.515;
smallest_sos_ms = pichardo_shear_floor_ms;
if isfield(parameters, 'medium_properties')
    mp_names_ = fieldnames(parameters.medium_properties);
    for mi_ = 1:numel(mp_names_)
        mp_ = parameters.medium_properties.(mp_names_{mi_});
        if isfield(mp_, 'sound_speed') && mp_.sound_speed > 0
            smallest_sos_ms = min(smallest_sos_ms, mp_.sound_speed);
        end
    end
end
ppw_auto = round(smallest_sos_ms / (freq_hz * res_mm * 1e-3));
ppw_use  = max(bb.ppw, ppw_auto);

dx_bb_est_mm = (smallest_sos_ms / freq_hz) / ppw_use * 1e3;
if abs(res_mm - dx_bb_est_mm) > 0.01
    fprintf(['[babelbrain_wrapper] Resolution note:\n' ...
        '  PRESTUS grid:   %.4f mm (parameters.grid.resolution_mm)\n' ...
        '  BabelBrain native: %.4f mm (%d PPW at %.0f Hz, SmallestSOS=%.1f m/s)\n' ...
        '  These differ — write_babelbrain_input.m resamples the mask to %.4f mm\n' ...
        '  before writing it, so BabelBrain''s voxel-count-based domain sizing still\n' ...
        '  covers the intended physical extent; no action needed. BabelBrain runs at\n' ...
        '  its own native step (PPW=%d, the BabelViscoFDTD benchmark fidelity) and the\n' ...
        '  result is resampled back onto the PRESTUS grid; the gap is expected.\n'], ...
        res_mm, dx_bb_est_mm, ppw_use, freq_hz, smallest_sos_ms, dx_bb_est_mm, ppw_use);
end

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

%% Focal calibration — locate or re-steer the geometric aim
% PRESTUS's focal_distance_ep is the intended FWHM-centre of the achieved field,
% but steering is geometric (peak lands proximal via the diffraction focal shift).
% Drive modes:
%   explicit phases (tr.annular.elem_phase_rad) — the per-element drive is fixed
%     and BabelBrain must reproduce it exactly (Step 3a exports it via the
%     weights file; ZSteering stays 0). Calibration runs locate-only, purely to
%     REPORT the predicted focus of that drive — it never replaces it.
%   no explicit phases — calibrate_focal_steering root-finds the geometric steer
%     D* whose achieved FLHM-centre (or peak) equals focal_distance_ep, using
%     the in-repo Rayleigh forward model; D* feeds the ZSteering block below.
% Skipped when the user pins z_steering_m.
calib = struct('active', false, 'mode', 'off');
has_explicit_phases = isfield(tr, 'annular') && isfield(tr.annular, 'elem_phase_rad') && ...
    ~isempty(tr.annular.elem_phase_rad);
is_annular_calib = isfield(tr, 'annular') && ...
    all(isfield(tr.annular, {'elem_n','elem_id_mm','elem_od_mm','curv_radius_mm', ...
                             'dist_geom_ep_mm','elem_p_amp'})) && ...
    isfield(tr, 'focal_distance_ep') && isfield(tr, 'freq_hz');
if isfield(tr, 'focal_distance_offset') && ~isempty(tr.focal_distance_offset)
    bowl_off_mm = tr.focal_distance_offset;
elseif isfield(tr, 'annular') && all(isfield(tr.annular, {'curv_radius_mm','dist_geom_ep_mm'}))
    bowl_off_mm = tr.annular.curv_radius_mm - tr.annular.dist_geom_ep_mm;
else
    bowl_off_mm = 0;
end
if bb.focal_calibration && is_annular_calib && isempty(bb.z_steering_m)
    copts = struct('target', bb.focal_calibration_target);
    if has_explicit_phases
        copts.phases_rad = tr.annular.elem_phase_rad;   % locate-only: drive is fixed
    end
    try
        [steer_ep, cdiag] = calibrate_focal_steering(tr, tr.focal_distance_ep, parameters, copts);
        calib.active        = true;
        calib.mode          = cdiag.mode;
        calib.steer_ep_mm   = steer_ep;
        calib.steer_bowl_mm = steer_ep + bowl_off_mm;
        calib.diag          = cdiag;
        fprintf(['[babelbrain_wrapper] Focal calibration (%s, target=%s): ' ...
            'steer %.2f mm EP (naive %.2f) -> achieved FLHM %.2f mm, peak %.2f mm, span %.2f mm\n'], ...
            calib.mode, bb.focal_calibration_target, calib.steer_ep_mm, tr.focal_distance_ep, ...
            cdiag.achieved_flhm_ep_mm, cdiag.achieved_peak_ep_mm, cdiag.flhm_span_mm);
    catch ME
        warning('babelbrain_wrapper:focalCalib', ...
            'Focal calibration failed (%s); falling back to naive geometric steer.', ME.message);
        calib.active = false;
    end
elseif bb.focal_calibration && ~is_annular_calib
    fprintf('[babelbrain_wrapper] Focal calibration skipped (non-annular or missing geometry).\n');
end

%% Transducer placement (TxMechanicalAdjustmentZ) and electronic aim (ZSteering)
% BabelBrain places the transducer apex at its natural focal distance (preset
% FocalLength) from the target, whereas PRESTUS places it at trans_pos, i.e.
% |focus_pos - trans_pos| away (trans_pos is the bowl rear centre = apex).
% That distance equals focal_distance_bowl only when the target sits at the
% nominal focal depth; with planned/pFOV placements it generally does not (e.g.
% 71 mm vs 92 mm), and using focal_distance_bowl then shifts the whole
% BabelBrain field axially relative to the head. Electronic steering cannot
% make up for that offset, so correct it mechanically:
%   TxMechanicalAdjustmentZ = preset_FL - |focus_pos - trans_pos| * res_mm
% (fallback without a valid trans_pos: preset_FL - focal_distance_bowl).
% This moves the transducer and BabelBrain's phase-reference aim together, so
% with ZSteering = 0 the aim stays at the natural focus. Then:
%   explicit per-element phases -> ZSteering = 0; the weights file carries the
%     full drive, so BabelBrain reproduces the phases the other solvers use.
%   no explicit phases          -> geometric steering from the apex:
%     ZSteering = aim_bowl - preset_FL, with aim_bowl = D* (calibrated) or
%     focal_distance_bowl (uncalibrated).
% See doc_babelbrain.md (Transducer placement, ZSteering, and the drive phases).

% Read FocalLength from the preset default.yaml. BabelBrain's folder name can
% differ from the preset string (e.g. 'CTX_500' lives in Babel_CTX500/), so try
% the canonical folder name first, then the preset string as-is.
preset_focal_length_m = [];   % also used by focal_length_mm_for_resample
preset_folder_name = bb_preset_folder(bb.tx_system);
preset_yaml_candidates = {
    fullfile(bb_root, 'BabelBrain', sprintf('Babel_%s', preset_folder_name), 'default.yaml'), ...
    fullfile(bb_root, 'BabelBrain', sprintf('Babel_%s', bb.tx_system),        'default.yaml'), ...
};
for k = 1:numel(preset_yaml_candidates)
    preset_yaml = preset_yaml_candidates{k};
    if isfile(preset_yaml)
        try
            preset_cfg = yaml.loadFile(preset_yaml, 'ConvertToArray', true);
            if isfield(preset_cfg, 'FocalLength')
                preset_focal_length_m = preset_cfg.FocalLength;
                fprintf('[babelbrain_wrapper] Preset yaml: %s\n', preset_yaml);
            end
        catch ME
            warning('babelbrain_wrapper: failed to read preset yaml %s: %s', preset_yaml, ME.message);
        end
        break
    end
end

% Focal distance from the bowl apex [mm]
focal_distance_bowl_mm = [];
if isfield(tr, 'focal_distance_bowl') && ~isempty(tr.focal_distance_bowl)
    focal_distance_bowl_mm = tr.focal_distance_bowl;
elseif isfield(tr, 'focal_distance_ep') && ~isempty(tr.focal_distance_ep)
    focal_distance_bowl_mm = tr.focal_distance_ep + bowl_off_mm;
end

% Mechanical placement correction (see block comment above), composed with any
% user-requested standoff.
[geom_corr_m, apex_to_target_mm, placement_src] = babelbrain_placement_correction( ...
    tr, focus_pos, trans_pos, res_mm, preset_focal_length_m);
switch placement_src
    case 'trans_pos'
        fprintf(['[babelbrain_wrapper] Tx placement correction: %+.1f mm ' ...
            '(preset FL %.1f mm - apex-to-target |focus_pos - trans_pos| %.2f mm)\n'], ...
            geom_corr_m*1e3, preset_focal_length_m*1e3, apex_to_target_mm);
        if ~isempty(focal_distance_bowl_mm) && abs(apex_to_target_mm - focal_distance_bowl_mm) > 1
            fprintf(['[babelbrain_wrapper]   note: target is %+.1f mm from the nominal focal depth ' ...
                '(focal_distance_bowl %.1f mm); BabelBrain apex placed at trans_pos.\n'], ...
                apex_to_target_mm - focal_distance_bowl_mm, focal_distance_bowl_mm);
        end
    case 'focal_distance_bowl'
        fprintf(['[babelbrain_wrapper] Tx placement correction: %+.1f mm ' ...
            '(preset FL %.1f mm - focal_distance_bowl %.1f mm; trans_pos unavailable, ' ...
            'assuming the target sits at the nominal focal depth)\n'], ...
            geom_corr_m*1e3, preset_focal_length_m*1e3, apex_to_target_mm);
end
tx_mech_adjust_total_m = bb.tx_mech_adjust_z_m + geom_corr_m;

if isempty(bb.z_steering_m)
    if has_explicit_phases
        % Drive is carried in full by the weights file (delta from the
        % natural-focus reference) — no electronic re-aim.
        z_steering_m = 0.0;
        fprintf('[babelbrain_wrapper] ZSteering=0 (explicit element phases carried by weights file)\n');
    elseif ~isempty(preset_focal_length_m) && ~isempty(focal_distance_bowl_mm)
        % ZSteering = aim_bowl - preset_FL (aim measured from the apex).
        if calib.active
            aim_bowl_mm   = calib.steer_bowl_mm;
            focal_ref_tag = sprintf('focal_calib_%s', calib.mode);
        else
            aim_bowl_mm   = focal_distance_bowl_mm;
            focal_ref_tag = 'focal_distance_bowl (uncalibrated)';
        end
        z_steering_m = aim_bowl_mm / 1000 - preset_focal_length_m;
        fprintf('[babelbrain_wrapper] ZSteering=%.4f m (%s: aim %.2f mm bowl, preset FL=%.3f m)\n', ...
            z_steering_m, focal_ref_tag, aim_bowl_mm, preset_focal_length_m);
    else
        z_steering_m = 0.0;
    end
else
    z_steering_m = bb.z_steering_m;
end

%% Pre-simulation parameter report and sanity checks
% Flags values that are unusual or that have caused bugs in the past so
% problems are visible before a multi-hour simulation runs to completion.

% Resolve computing backend early so it appears in the summary below.
[bb.computing_backend, bb.device_name] = resolve_bb_backend(bb.computing_backend, bb.device_name);

backend_names = {'CUDA', 'OpenCL', 'Metal', 'CPU'};
if isnumeric(bb.computing_backend) && bb.computing_backend >= 1 && bb.computing_backend <= 4
    backend_str = sprintf('%s (%s)', backend_names{bb.computing_backend}, bb.device_name);
else
    backend_str = sprintf('%s (%s)', num2str(bb.computing_backend), bb.device_name);
end

fprintf('\n');
fprintf('[babelbrain_wrapper] ---- BabelBrain simulation summary ----\n');

% --- Transducer & frequency ---
fprintf('[babelbrain_wrapper]  backend   : %s\n', backend_str);
fprintf('[babelbrain_wrapper]  transducers: %d\n', n_tx);
fprintf('[babelbrain_wrapper]  Preset:           %s\n', bb.tx_system);
fprintf('[babelbrain_wrapper]  Frequency:        %.0f Hz (%.0f kHz)\n', freq_hz, freq_hz/1e3);
fprintf('[babelbrain_wrapper]  PRESTUS grid:     %.4f mm\n', res_mm);
fprintf('[babelbrain_wrapper]  BabelBrain est.:  %.4f mm  (%d PPW)\n', dx_bb_est_mm, ppw_use);

% --- Focal geometry ---
if isfield(tr, 'focal_distance_ep')
    fprintf('[babelbrain_wrapper]  focal_distance_ep:  %.1f mm\n', tr.focal_distance_ep);
end
if ~isempty(preset_focal_length_m)
    fprintf('[babelbrain_wrapper]  Preset FocalLength: %.1f mm (from yaml)\n', preset_focal_length_m*1e3);
else
    fprintf('[babelbrain_wrapper]  Preset FocalLength: NOT FOUND  ← fallback path active\n');
end
fprintf('[babelbrain_wrapper]  ZSteering:        %.4f m (%.1f mm)\n', z_steering_m, z_steering_m*1e3);
fprintf('[babelbrain_wrapper]  TxMechAdjustZ:    %.4f m (%.1f mm = %.1f mm placement corr + %.1f mm user)\n', ...
    tx_mech_adjust_total_m, tx_mech_adjust_total_m*1e3, geom_corr_m*1e3, bb.tx_mech_adjust_z_m*1e3);
if ~isempty(preset_focal_length_m)
    fprintf('[babelbrain_wrapper]  Resample anchor:  /TargetLocation tissue marker (fallback: preset FL %.1f mm)\n', ...
        preset_focal_length_m*1e3);
elseif isfield(tr, 'focal_distance_ep') && ~isempty(tr.focal_distance_ep)
    fprintf('[babelbrain_wrapper]  Resample anchor:  /TargetLocation tissue marker (fallback: focal_distance_bowl %.1f mm — preset FL unknown)\n', ...
        tr.focal_distance_ep + bowl_off_mm);
else
    fprintf('[babelbrain_wrapper]  Resample anchor:  /TargetLocation tissue marker (no fallback available)\n');
end

% --- Domain ---
if tight_narrow_beam
    fprintf('[babelbrain_wrapper]  Domain:           narrow-beam  (z_beyond_focus=%.1f mm)\n', z_beyond_focus_mm);
else
    fprintf('[babelbrain_wrapper]  Domain:           full BabelBrain default\n');
end

% --- Medium & backend ---
if is_water_medium
    fprintf('[babelbrain_wrapper]  Medium:           water-only  (--water_only)\n');
else
    fprintf('[babelbrain_wrapper]  Medium:           layered/tissue\n');
end
if ~isempty(bb.ct_mapping)
    fprintf('[babelbrain_wrapper]  CT mapping:       %s  (pseudoCT path active)\n', bb.ct_mapping);
end

% --- Focus position in grid ---
gs = size(medium_masks);
fprintf('[babelbrain_wrapper]  focus_pos (grid): [%d %d %d]  grid size: [%d %d %d]\n', ...
    round(focus_pos), gs(1), gs(2), gs(3));

fprintf('[babelbrain_wrapper] ---- Checks ----\n');
any_warn = false;

% CHECK: preset yaml found?
if isempty(preset_focal_length_m)
    fprintf('[babelbrain_wrapper]  ⚠  WARN: Preset default.yaml not found — Tx placement\n');
    fprintf('[babelbrain_wrapper]           correction, ZSteering, and resample anchor are ALL\n');
    fprintf('[babelbrain_wrapper]           falling back. Results will be misregistered if\n');
    fprintf('[babelbrain_wrapper]           preset FL ≠ focal_distance_bowl.\n');
    any_warn = true;
end

% CHECK: large ZSteering
if abs(z_steering_m) > 0.020
    fprintf('[babelbrain_wrapper]  ⚠  WARN: |ZSteering| = %.1f mm > 20 mm — large axial steering.\n', ...
        abs(z_steering_m)*1e3);
    fprintf('[babelbrain_wrapper]           Verify focal_distance_ep (%.1f mm) and preset FL (%.1f mm).\n', ...
        tr.focal_distance_ep, preset_focal_length_m*1e3);
    any_warn = true;
end

% CHECK: resampling needed?
if abs(res_mm - dx_bb_est_mm) > 0.01
    fprintf('[babelbrain_wrapper]  ⚠  NOTE: Grid mismatch — resampling from %.4f → %.4f mm.\n', ...
        dx_bb_est_mm, res_mm);
    any_warn = true;
end

% CHECK: focus outside medium mask
fp = round(focus_pos);
if any(fp < 1) || any(fp > gs)
    fprintf('[babelbrain_wrapper]  ✗  ERROR: focus_pos [%d %d %d] is outside grid [%d %d %d]!\n', ...
        fp, gs(1), gs(2), gs(3));
    any_warn = true;
else
    focus_label = medium_masks(fp(1), fp(2), fp(3));
    medium_labels_local = fieldnames(parameters.medium_properties);
    if focus_label >= 1 && focus_label <= numel(medium_labels_local)
        focus_tissue = medium_labels_local{focus_label};
    else
        focus_tissue = sprintf('label=%d (background/water)', focus_label);
    end
    if contains(focus_tissue, 'skull') || contains(focus_tissue, 'bone')
        fprintf('[babelbrain_wrapper]  ⚠  WARN: focus_pos voxel is inside "%s" — not in brain/water.\n', focus_tissue);
        any_warn = true;
    else
        fprintf('[babelbrain_wrapper]  ✓  focus_pos tissue: %s\n', focus_tissue);
    end
end

% CHECK: focal distance beyond domain
if ~isempty(fov_len_mm) && isfield(tr, 'focal_distance_ep')
    if tr.focal_distance_ep > fov_len_mm - 5
        fprintf('[babelbrain_wrapper]  ⚠  WARN: focal_distance_ep (%.1f mm) is within 5 mm of FOV end (%.1f mm).\n', ...
            tr.focal_distance_ep, fov_len_mm);
        any_warn = true;
    end
end

% CHECK: ZSteering consistency.  With the mechanical placement correction the
% geometry inside BabelBrain is physical (apex at trans_pos, |focus_pos -
% trans_pos| from the target), so the auto-derived invariant is (aim measured from the apex):
%   explicit phases -> ZSteering = 0 (drive carried by the weights file)
%   calibrated      -> ZSteering = calib.steer_bowl_mm/1e3 - preset_FL
%   uncalibrated    -> ZSteering = focal_distance_bowl/1e3 - preset_FL
% A mismatch is expected only when z_steering_m was pinned manually.
if isempty(bb.z_steering_m) && ~isempty(preset_focal_length_m) && ~isempty(focal_distance_bowl_mm)
    if has_explicit_phases
        expected_z_steering = 0.0;
        ref_dist_tag = '0 (explicit phases in weights file)';
    elseif calib.active
        expected_z_steering = calib.steer_bowl_mm / 1e3 - preset_focal_length_m;
        ref_dist_tag = 'calib.steer_bowl - preset FL';
    else
        expected_z_steering = focal_distance_bowl_mm / 1e3 - preset_focal_length_m;
        ref_dist_tag = 'focal_distance_bowl - preset FL';
    end
    if abs(z_steering_m - expected_z_steering) > 1e-4
        fprintf('[babelbrain_wrapper]  ⚠  WARN: ZSteering (%.4f m) ≠ %s (%.4f m).\n', ...
            z_steering_m, ref_dist_tag, expected_z_steering);
        fprintf('[babelbrain_wrapper]           This is only valid if z_steering_m was set manually.\n');
        any_warn = true;
    end
end

if ~any_warn
    fprintf('[babelbrain_wrapper]  ✓  All checks passed.\n');
end
fprintf('[babelbrain_wrapper] ------------------------------------\n\n');

%% Step 2 — Write BabelViscoInput NIfTI (and companion CT NIfTI, if configured)
fprintf('[babelbrain_wrapper] Writing input mask...\n');
if ~isempty(bb.ct_mapping) && isempty(pseudoCT)
    error(['babelbrain_wrapper: ct_mapping is set but no pseudoCT was supplied.\n' ...
        'Pass the pseudoCT output of grid_tissue_setup as the 8th argument.']);
end
pseudoCT_for_write = [];
if ~isempty(bb.ct_mapping)
    pseudoCT_for_write = pseudoCT;
end
[nii_path, ct_nii_path] = write_babelbrain_input(medium_masks, focus_pos, parameters, planimg, out_dir, target, ppw_use, pseudoCT_for_write, trans_pos);

%% Step 2b — Export PRESTUS medium properties for BabelBrain material injection
% BabelBrain's non-CT path uses its own built-in material database (MatFreq).
% Dump PRESTUS's configured properties to JSON so run_babelbrain_case.py can
% overwrite MatFreq with them — giving k-Wave / jwave / BabelBrain identical
% tissue acoustics. Failure here is non-fatal: BabelBrain falls back to its
% built-in database.
medium_json_path = fullfile(out_dir, subj_id, sprintf('%s_medium_properties.json', target));
try
    fid = fopen(medium_json_path, 'w');
    if fid == -1
        error('cannot open %s for writing', medium_json_path);
    end
    fwrite(fid, jsonencode(parameters.medium_properties), 'char');
    fclose(fid);
    fprintf('[babelbrain_wrapper] Wrote medium properties JSON: %s\n', medium_json_path);
catch ME
    if exist('fid', 'var') && fid ~= -1
        fclose(fid);
    end
    warning('babelbrain_wrapper:mediumJson', ...
        'Failed to write medium properties JSON (%s). BabelBrain will use its built-in material database.', ...
        ME.message);
    medium_json_path = '';
end

%% Step 3 — Write OptimizedWeightsFile (calibrated phases)
weights_path = fullfile(out_dir, subj_id, [target '_CalibWeights.h5']);

% CHECK: source amplitude available when weights exist (must follow weights_path definition)
if ~isempty(weights_path) && isfile(weights_path) && ~exist('elem_p_amp', 'var')
    fprintf('[babelbrain_wrapper]  ⚠  WARN: Weights found but elem_p_amp not set — source amplitude flag will be omitted.\n');
    any_warn = true;
end

if isfield(tr, 'annular') && isfield(tr.annular, 'elem_phase_rad') && isfield(tr.annular, 'elem_p_amp')
    phases_rad = tr.annular.elem_phase_rad;
    elem_p_amp   = tr.annular.elem_p_amp(1);

    % Whole-profile Rayleigh calibration of the BabelBrain drive phases.
    % Replaces single-point FLHM-centre matching: instead of reproducing the
    % (O'Neil-calibrated) library phases verbatim -- which BabelBrain's
    % discrete-ring Rayleigh field renders with a residual proximal bias and
    % axial broadening -- the phases are re-fit so that the discrete-ring
    % Rayleigh profile matches the continuous-bowl O'Neil ideal across the
    % whole focal region. BabelBrain-specific: only the weights delta below
    % uses these; the shared library phases (k-Wave / jwave) are untouched.
    if bb.profile_calibration && is_annular_calib
        try
            [prof_phases, ~, pdiag] = calibrate_babelbrain_phase_profile(tr, parameters, ...
                struct('outputs_folder', fullfile(out_dir, subj_id)));
            fprintf(['[babelbrain_wrapper] Whole-profile phase calibration: ' ...
                'library peak %.1f mm EP -> fit peak %.1f mm EP (O''Neil target %.1f mm, ' ...
                'fit FLHM %.1f mm; %d pts, err %.3g).\n'], ...
                pdiag.library_peak_ep_mm, pdiag.fit_peak_ep_mm, pdiag.target_peak_ep_mm, ...
                pdiag.fit_flhm_ep_mm, pdiag.n_pts, pdiag.min_err);
            phases_rad = prof_phases;
        catch ME
            warning('babelbrain_wrapper:profileCalib', ...
                ['Whole-profile phase calibration failed (%s); ' ...
                 'falling back to library drive phases.'], ME.message);
        end
    end

    % Delta reference: BabelBrain multiplies the weights onto ring phases it
    % derives internally (AllPhi) by Rayleigh back-propagation from its steering
    % target, and the delta must be referenced to the geometric phases at THAT
    % aim so that AllPhi + delta reproduces phases_rad. The exact aim BabelBrain
    % uses is entangled with its voxelised grid setup (see the analysis in
    % doc_babelbrain.md / CHANGELOG) and is NOT cleanly reproduced by an ideal-
    % geometry compute_phases call: sweeping the aim, the geometric focus
    % (preset_FL - bowl_off = dist_geom_ep) is the closest single reference
    % (water focal centre lands at ~84.5 mm vs the k-Wave/jwave 85 mm), but a
    % residual ~2-3 mm proximal peak bias persists. The durable fix is to
    % reference the delta to BabelBrain's OWN emitted AllPhi rather than an
    % analytic reconstruction; until that is wired in, use the geometric-focus
    % aim as the best available approximation.
    if ~isempty(preset_focal_length_m)
        aim_ep_ref_mm = (preset_focal_length_m + z_steering_m) * 1e3 - bowl_off_mm;
    else
        aim_ep_ref_mm = tr.annular.dist_geom_ep_mm + z_steering_m * 1e3;
    end
    try
        tran_geo = struct('n_elem', tr.annular.elem_n, 'manufact', '');
        tran_geo.transducer.freq_hz = tr.freq_hz;
        tran_geo.transducer.annular = tr.annular;
        geo_phases_deg = compute_phases( ...
            parameters.medium_properties.water.sound_speed, ...
            tran_geo, aim_ep_ref_mm);
        geo_phases_rad = geo_phases_deg * pi / 180;
        fprintf(['[babelbrain_wrapper] Weights delta referenced to BabelBrain aim ' ...
            '%.1f mm EP (geometric focus; preset FL %+.1f mm ZSteering).\n'], ...
            aim_ep_ref_mm, z_steering_m*1e3);
    catch ME
        error('babelbrain_wrapper:geoPhases', ...
            ['Could not compute the delta-reference geometric phases at %.1f mm EP: %s\n' ...
             'Ensure tr.annular has elem_id_mm, elem_od_mm, curv_radius_mm, dist_geom_ep_mm,\n' ...
             'tr.freq_hz, and parameters.medium_properties.water.sound_speed.'], ...
            aim_ep_ref_mm, ME.message);
    end

    elem_amp_for_weights = [];
    if isfield(tr.annular, 'elem_amp') && ~isempty(tr.annular.elem_amp)
        elem_amp_for_weights = tr.annular.elem_amp;
    end
    export_babelbrain_weights(phases_rad, elem_p_amp, parameters, weights_path, geo_phases_rad, elem_amp_for_weights);
else
    weights_path = '';
    fprintf('[babelbrain_wrapper] No element phases found — BabelBrain will use geometric focusing.\n');
end

%% Step 4 — Invoke Python bridge
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

% BabelBrain source amplitude:
% elem_p_amp = rho_water * c_water * v_opt is the physical face pressure derived
% by the O'Neill/Rayleigh calibration fit.  BabelBrain's --source_amp_pa is the
% same quantity (Rayleigh surface pressure amplitude), so elem_p_amp transfers
% directly — no cross-calibration is needed for free-water simulations.
% An explicit bb.source_amp_pa overrides the auto-populated value.
resolved_source_amp_pa = bb.source_amp_pa;
if isempty(resolved_source_amp_pa) || resolved_source_amp_pa <= 0
    if exist('elem_p_amp', 'var') && elem_p_amp > 0
        resolved_source_amp_pa = elem_p_amp;
        fprintf('[babelbrain_wrapper] source_amp_pa auto-set from elem_p_amp (O''Neill face pressure): %.1f Pa\n', resolved_source_amp_pa);
    end
end

source_amp_pa_flag = {};
if ~isempty(resolved_source_amp_pa) && resolved_source_amp_pa > 0
    source_amp_pa_flag = {'--source_amp_pa', sprintf('%.2f', resolved_source_amp_pa)};
    if ~isempty(bb.source_amp_pa) && bb.source_amp_pa > 0
        fprintf('[babelbrain_wrapper] BabelBrain source amplitude (explicit override): %.1f Pa\n', resolved_source_amp_pa);
    end
else
    fprintf('[babelbrain_wrapper]  ⚠  WARN: source_amp_pa not set and elem_p_amp unavailable — BabelBrain will use built-in default (~50 kPa).\n');
    fprintf('[babelbrain_wrapper]           Isppa output will NOT be comparable to kWave.\n');
end

% run_command invokes python directly via argv (no shell), so every path,
% subject and target string reaches the process as a literal argument.
babelbrain_argv = [{python_bin, bb_runner, ...
    '--basedir', out_dir, ...
    '--id', subj_id, ...
    '--target', target, ...
    '--tx_system', bb.tx_system, ...
    '--freq', sprintf('%.0f', freq_hz), ...
    '--ppw', num2str(ppw_use), ...
    '--computing_backend', num2str(bb.computing_backend), ...
    '--device_name', bb.device_name, ...
    '--tx_mech_adjust_z', sprintf('%.6f', tx_mech_adjust_total_m), ...
    '--z_steering', sprintf('%.6f', z_steering_m), ...
    '--z_beyond_focus_mm', sprintf('%.1f', z_beyond_focus_mm), ...
    '--log_file', log_file}, ...
    bb_root_flag(bb_root), ...
    bb_weights_flag(weights_path), ...
    source_amp_pa_flag, ...
    bb_refocus_flag(bb.do_refocusing), ...
    bb_dryrun_flag(dry_run), ...
    bb_no_tight_flag(tight_narrow_beam), ...
    bb_usect_flag(~isempty(bb.ct_mapping)), ...
    bb_ctmapping_flag(bb.ct_mapping), ...
    bb_petra_flag(bb.ct_use_petra), ...
    bb_medium_json_flag(medium_json_path), ...
    bb_wateronly_flag(is_water_medium), ...
    bb_rayleigh_water_flag(bb.use_rayleigh_for_water), ...
    bb_debug_flag(parameters)];

if is_water_medium
    fprintf('[babelbrain_wrapper] Water medium: running in --water_only mode.\n');
end
fprintf('[babelbrain_wrapper] Command:\n  %s\n', strjoin(babelbrain_argv, ' '));

% Pre-run cache invalidation: BabelBrain's RunCases() skips simulation if the
% DataForSim.h5 already exists. Delete it when overwrite_files='always' so that
% changed parameters (source_amp_pa, z_steering, etc.) actually take effect.
if isfield(parameters, 'io') && isfield(parameters.io, 'overwrite_files') && ...
        any(strcmp(parameters.io.overwrite_files, {'always', 'force'}))
    freq_kHz_pre = round(freq_hz / 1e3);
    if is_water_medium
        h5_pre = sprintf('%s_%dkHz_%dPPW_Water_DataForSim.h5', target, freq_kHz_pre, ppw_use);
    else
        h5_pre = sprintf('%s_%dkHz_%dPPW_DataForSim.h5', target, freq_kHz_pre, ppw_use);
    end
    h5_pre_path = fullfile(out_dir, subj_id, h5_pre);
    if isfile(h5_pre_path)
        delete(h5_pre_path);
        fprintf('[babelbrain_wrapper] Deleted stale DataForSim.h5 cache (overwrite_files=always):\n  %s\n', h5_pre_path);
    end
end

% MATLAB and the BabelBrain child share one SLURM cgroup memory limit; log
% MATLAB's share so an OOM kill can be attributed (the child logs its own).
log_memory_usage('before BabelBrain launch');

t0     = tic;
% Echo=true streams output live to the console, matching the original
% single-output `status = system(cmd_args)` (which prints stdout/stderr to
% the Command Window in real time when its output isn't captured).
[status, cmd_out] = run_command(babelbrain_argv, 'Echo', true);
elapsed = toc(t0);
mem_after = log_memory_usage('after BabelBrain exit');
log_memory_summary('BabelBrain', cmd_out, mem_after);

if status ~= 0
    error('babelbrain_wrapper: BabelBrain Python process exited with status %d.\nCommand was:\n  %s\nOutput:\n%s', ...
        status, strjoin(babelbrain_argv, ' '), cmd_out);
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
% In water medium mode, --water_only is passed so BabelBrain writes the water variant.
freq_kHz = round(freq_hz / 1e3);
if is_water_medium
    h5_name = sprintf('%s_%dkHz_%dPPW_Water_DataForSim.h5', target, freq_kHz, ppw_use);
else
    h5_name = sprintf('%s_%dkHz_%dPPW_DataForSim.h5', target, freq_kHz, ppw_use);
end
h5_path  = fullfile(out_dir, subj_id, h5_name);

if ~isfile(h5_path)
    % Fallback: glob for DataForSim.h5 files (BabelBrain version may name differently)
    candidates = dir(fullfile(out_dir, subj_id, '*DataForSim.h5'));
    if is_water_medium
        % In water mode the primary output IS the water H5
        candidates = candidates(contains({candidates.name}, 'Water'));
    else
        % Exclude water baseline files from the acoustic output search
        candidates = candidates(~contains({candidates.name}, 'Water'));
    end
    if isempty(candidates)
        error('babelbrain_wrapper: Output HDF5 not found.\n  Expected: %s\n  No *DataForSim.h5 files in: %s', ...
            h5_path, fullfile(out_dir, subj_id));
    end
    [~, newest] = max([candidates.datenum]);
    h5_path = fullfile(candidates(newest).folder, candidates(newest).name);
    warning('babelbrain_wrapper:outputNameMismatch', ...
        'Expected output not found; using most-recently-modified DataForSim.h5:\n  %s', h5_path);
end

rewrite_hdf5_no_blosc(h5_path, python_bin, this_dir);

% focal_length_mm: FALLBACK z anchor for the output resample.
% read_babelbrain_output anchors on /TargetLocation (the label-5 marker
% write_babelbrain_input placed at focus_pos), which maps every returned volume
% back onto the PRESTUS grid. This value is used only when /TargetLocation is
% missing from the HDF5: BabelBrain places the marker at z = preset FocalLength,
% so that is the fallback; focal_distance_bowl is a last resort (the reader
% warns). See doc_babelbrain.md (Output alignment).
focal_length_mm_for_resample = [];
if ~isempty(preset_focal_length_m)
    focal_length_mm_for_resample = preset_focal_length_m * 1e3;   % marker z by construction
    fprintf('[babelbrain_wrapper] Fallback resample anchor: preset FocalLength = %.1f mm (marker z).\n', ...
        focal_length_mm_for_resample);
elseif isfield(tr, 'focal_distance_ep') && ~isempty(tr.focal_distance_ep)
    focal_length_mm_for_resample = tr.focal_distance_ep + bowl_off_mm;   % = focal_distance_bowl (last resort)
    fprintf(['[babelbrain_wrapper] Preset FocalLength unavailable; fallback resample anchor set to ' ...
        'focal_distance_bowl = %.1f mm. Only exact if it equals the marker z (preset FL); ' ...
        'normally /TargetLocation in the HDF5 makes this fallback unused.\n'], focal_length_mm_for_resample);
end
if isempty(focal_length_mm_for_resample)
    error('babelbrain_wrapper: cannot determine focal_length_mm for resampling. Provide focal_distance_ep or preset FocalLength.');
end
sensor_data = read_babelbrain_output(h5_path, bb.do_refocusing, parameters, planimg, focus_pos, focal_length_mm_for_resample, trans_pos);
sensor_data.babelbrain_version   = bb_version;
sensor_data.focus_pos            = focus_pos;
sensor_data.focal_length_mm      = focal_length_mm_for_resample;
sensor_data.trans_pos            = trans_pos;   % beam-axis direction for the thermal resample

% Derive free-water ISPPA.
% In water medium mode, the primary simulation IS the water simulation, so
% freefield_isppa_wcm2 is derived directly from sensor_data without a separate run.
% In layered/phantom mode, a dedicated water-only BabelBrain run provides the baseline.
p_bb_water_peak_pa = [];   % BabelBrain free-water focal peak, for amplitude calibration below
if is_water_medium
    p_max_pa = max(sensor_data.p_max_all(:));
    p_bb_water_peak_pa = p_max_pa;
    sensor_data.freefield_isppa_wcm2 = p_max_pa^2 / (2 * 1000 * 1500) / 1e4;
    fprintf('[babelbrain_wrapper] Water medium — freefield ISPPA: %.4f W/cm²\n', ...
        sensor_data.freefield_isppa_wcm2);
else
    % The water baseline output is named <prefix>_Water_DataForSim.h5 in the same
    % directory. If absent (e.g. run_water_baseline=0 or external workflow), the
    % water-only simulation is triggered here before reading. The ISPPA is needed
    % downstream for babelbrain_thermal_wrapper.
    water_h5 = strrep(h5_path, '_DataForSim.h5', '_Water_DataForSim.h5');
    % Same cache invalidation as the main h5 above (line ~847): without this,
    % a stale Water_DataForSim.h5 from an earlier/buggy run is silently reused
    % forever, since the block below only triggers a re-run when the file is
    % entirely absent. A stale water baseline contaminates the free-water ISPPA
    % normalization used for the skull-attenuation ratio and downstream thermal
    % calibration.
    if isfile(water_h5) && isfield(parameters, 'io') && isfield(parameters.io, 'overwrite_files') && ...
            any(strcmp(parameters.io.overwrite_files, {'always', 'force'}))
        delete(water_h5);
        fprintf('[babelbrain_wrapper] Deleted stale Water_DataForSim.h5 cache (overwrite_files=always):\n  %s\n', water_h5);
    end
    if ~isfile(water_h5)
        % Water file missing — run the water-only simulation now.
        % This can happen when the acoustic run was done without the water baseline
        % (e.g. run_water_baseline=0 in config, or via an external workflow).
        fprintf('[babelbrain_wrapper] Water_DataForSim.h5 not found — running water-only simulation...\n');
        water_log = strrep(log_file, '_babelbrain_acoustic.log', '_babelbrain_water.log');
        % run_command invokes python directly via argv, no shell at all.
        water_argv = [{python_bin, bb_runner, ...
            '--basedir', out_dir, ...
            '--id', subj_id, ...
            '--target', target, ...
            '--tx_system', bb.tx_system, ...
            '--freq', sprintf('%.0f', freq_hz), ...
            '--ppw', num2str(ppw_use), ...
            '--computing_backend', num2str(bb.computing_backend), ...
            '--device_name', bb.device_name, ...
            '--tx_mech_adjust_z', sprintf('%.6f', tx_mech_adjust_total_m), ...
            '--z_steering', sprintf('%.6f', z_steering_m), ...
            '--z_beyond_focus_mm', sprintf('%.1f', z_beyond_focus_mm), ...
            '--log_file', water_log, ...
            '--water_only'}, ...
            bb_root_flag(bb_root), ...
            bb_weights_flag(weights_path), ...  % needed for z-steering geometry match
            source_amp_pa_flag, ...             % same drive as the transcranial run, else BabelBrain's ~50 kPa default
            bb_refocus_flag(bb.do_refocusing), ...
            bb_no_tight_flag(tight_narrow_beam), ...
            bb_usect_flag(~isempty(bb.ct_mapping)), ...
            bb_ctmapping_flag(bb.ct_mapping), ...
            bb_petra_flag(bb.ct_use_petra), ...
            bb_rayleigh_water_flag(bb.use_rayleigh_for_water), ...
            bb_debug_flag(parameters)];
        fprintf('[babelbrain_wrapper] Water-only command:\n  %s\n', strjoin(water_argv, ' '));
        t0_w = tic;
        [status_w, water_out] = run_command(water_argv, 'Echo', true);
        elapsed_w = toc(t0_w);
        log_memory_summary('BabelBrain water-only', water_out, ...
            log_memory_usage('after BabelBrain water-only exit'));
        if status_w ~= 0
            warning('babelbrain_wrapper:waterSimFailed', ...
                'Water-only simulation exited with status %d (%.1f s). ISPPA normalization may be inaccurate.\nOutput:\n%s', ...
                status_w, elapsed_w, water_out);
        else
            fprintf('[babelbrain_wrapper] Water-only simulation finished in %.1f s.\n', elapsed_w);
        end
    end
    if ~isfile(water_h5)
        candidates_w = dir(fullfile(fileparts(h5_path), '*Water_DataForSim.h5'));
        if ~isempty(candidates_w)
            [~, newest_w] = max([candidates_w.datenum]);
            water_h5 = fullfile(candidates_w(newest_w).folder, candidates_w(newest_w).name);
            warning('babelbrain_wrapper:waterNameMismatch', ...
                'Expected water output not found; using most-recently-modified Water_DataForSim.h5:\n  %s', water_h5);
        end
    end
    if isfile(water_h5)
        rewrite_hdf5_no_blosc(water_h5, python_bin, this_dir);
        try
            p_water = double(h5read(water_h5, '/p_amp'));
            p_bb_water_peak_pa = max(p_water(:));
            % ISPPA = p_peak^2 / (2 * rho_water * c_water), converted Pa→W/cm²
            sensor_data.freefield_isppa_wcm2 = p_bb_water_peak_pa^2 / (2 * 1000 * 1500) / 1e4;
            fprintf('[babelbrain_wrapper] Free-water ISPPA: %.4f W/cm²\n', ...
                sensor_data.freefield_isppa_wcm2);
        catch ME
            warning('babelbrain_wrapper:waterIsppa', ...
                'Could not read water ISPPA from %s: %s', water_h5, ME.message);
        end
    else
        warning('babelbrain_wrapper:noWaterH5', ...
            'Water_DataForSim.h5 not found — babelbrain thermal will need isppa_wcm2 set manually.');
    end
end  % if is_water_medium / else

%% Step 5b — Free-water intensity normalization (intended vs achieved ISPPA)
% Scale the returned pressure so the free-water ISPPA matches the transducer's
% intended value, as assemble_limited_fov_fields does:
%   k = sqrt(target_isppa_wcm2 / achieved_freewater_isppa)
% Runs only when a BabelBrain free-water baseline was simulated and
% target_isppa_wcm2 is set; otherwise the raw solver field is kept. Only
% sensor_data is scaled: DataForSim.h5 is left untouched for the thermal
% backend, whose BHTE renormalizes to its own requested Isppa. The factor is
% stored in sensor_data.amp_calibration. See doc_babelbrain.md (Free-water
% amplitude calibration).
if bb.water_amp_calibration
    target_isppa = NaN;
    if isfield(tr, 'target_isppa_wcm2') && isscalar(tr.target_isppa_wcm2) ...
            && isfinite(tr.target_isppa_wcm2) && tr.target_isppa_wcm2 > 0
        target_isppa = tr.target_isppa_wcm2;
    end
    achieved_isppa = NaN;
    if isfield(sensor_data, 'freefield_isppa_wcm2') && isscalar(sensor_data.freefield_isppa_wcm2) ...
            && isfinite(sensor_data.freefield_isppa_wcm2) && sensor_data.freefield_isppa_wcm2 > 0
        achieved_isppa = sensor_data.freefield_isppa_wcm2;
    end

    if isempty(p_bb_water_peak_pa) || ~isfinite(achieved_isppa)
        fprintf(['[babelbrain_wrapper] Free-water intensity normalization skipped: ' ...
            'no BabelBrain water baseline / achieved free-water ISPPA unavailable. ' ...
            'Raw solver field kept.\n']);
    elseif ~isfinite(target_isppa)
        fprintf(['[babelbrain_wrapper] Free-water intensity normalization skipped: ' ...
            'transducer.target_isppa_wcm2 not set. Raw solver field kept ' ...
            '(achieved free-water ISPPA = %.2f W/cm²).\n'], achieved_isppa);
    else
        k_amp = sqrt(target_isppa / achieved_isppa);
        fprintf(['[babelbrain_wrapper] Free-water intensity normalization: ' ...
            'target %.2f / achieved %.2f W/cm² -> k = %.4f (pressure)\n'], ...
            target_isppa, achieved_isppa, k_amp);
        if k_amp > 4 || k_amp < 0.25
            warning('babelbrain_wrapper:ampCalibLarge', ...
                ['Free-water normalization factor k = %.3f (>4× or <0.25× pressure). ' ...
                 'BabelBrain''s free-water run achieved only %.2f W/cm² against an intended ' ...
                 '%.2f W/cm² — a gap this large points to a setup problem (e.g. the free-water ' ...
                 'focus is mis-placed / defocused rather than landing at the intended focal ' ...
                 'depth), not ordinary drive scaling. Inspect the free-water focal profile.'], ...
                k_amp, achieved_isppa, target_isppa);
        end
        sensor_data.p_max_all = sensor_data.p_max_all * k_amp;
        if isfield(sensor_data, 'p_complex')
            sensor_data.p_complex = sensor_data.p_complex * k_amp;
        end
        sensor_data.freefield_isppa_wcm2 = achieved_isppa * k_amp^2;   % == target_isppa
        sensor_data.amp_calibration = struct( ...
            'k_pressure',                    k_amp, ...
            'target_isppa_wcm2',             target_isppa, ...
            'achieved_freewater_isppa_wcm2', achieved_isppa, ...
            'p_bb_water_peak_pa',            p_bb_water_peak_pa, ...
            'reference', 'target_isppa_wcm2 / achieved free-water ISPPA (BabelBrain water run)');
    end
end

%% Step 6 — Optional cache cleanup
if bb.cleanup_cache
    bb_cache_dir = fullfile(out_dir, subj_id);
    if isfolder(bb_cache_dir)
        % Move DataForSim.h5 and Water_DataForSim.h5 out of the cache before
        % deleting it so that babelbrain_thermal_wrapper can still find them.
        persistent_h5 = fullfile(out_dir, h5_name);
        try
            movefile(h5_path, persistent_h5);
            sensor_data.babelbrain_source = persistent_h5;
            fprintf('[babelbrain_wrapper] Moved acoustic HDF5 to: %s\n', persistent_h5);
        catch ME
            warning('babelbrain_wrapper:h5MoveFailed', ...
                'Could not move DataForSim.h5 before cache cleanup — skipping cleanup: %s', ME.message);
            bb.cleanup_cache = false;
        end
        if bb.cleanup_cache && exist('water_h5', 'var') && isfile(water_h5)
            water_h5_name    = sprintf('%s_%dkHz_%dPPW_Water_DataForSim.h5', target, freq_kHz, ppw_use);
            persistent_water = fullfile(out_dir, water_h5_name);
            try
                movefile(water_h5, persistent_water);
                water_h5 = persistent_water;
                fprintf('[babelbrain_wrapper] Moved water HDF5 to: %s\n', persistent_water);
            catch ME
                warning('babelbrain_wrapper:waterH5MoveFailed', ...
                    'Could not move Water_DataForSim.h5 before cache cleanup: %s', ME.message);
            end
        end
        if bb.cleanup_cache
            try
                rmdir(bb_cache_dir, 's');
                fprintf('[babelbrain_wrapper] Cache cleaned up: %s\n', bb_cache_dir);
            catch ME
                warning('babelbrain_wrapper:rmdirFailed', ...
                    'Cache cleanup failed (NFS lock or permissions): %s\n  %s', ...
                    bb_cache_dir, ME.message);
            end
        end
    end
end

end

%% Local helpers

function sfx = trans_suffix(ti, n_tx)
% Return a unique affix for transducer ti of n_tx (e.g. '_tx1of3').
    sfx = sprintf('_tx%dof%d', ti, n_tx);
end

function s = apply_default(s, field, val)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = val;
    end
end

function flag = bb_weights_flag(weights_path)
    if ~isempty(weights_path) && isfile(weights_path)
        flag = {'--weights_file', weights_path};
    else
        flag = {};
    end
end

function flag = bb_root_flag(bb_root)
    if ~isempty(bb_root)
        flag = {'--babelbrain_root', bb_root};
    else
        flag = {};
    end
end

function flag = bb_refocus_flag(do_refocus)
    if do_refocus
        flag = {'--do_refocusing'};
    else
        flag = {};
    end
end

function flag = bb_dryrun_flag(dry_run)
    if dry_run
        flag = {'--dry_run'};
    else
        flag = {};
    end
end

function flag = bb_usect_flag(use_ct)
    if use_ct
        flag = {'--use_ct'};
    else
        flag = {};
    end
end

function flag = bb_ctmapping_flag(ct_mapping)
    if ~isempty(ct_mapping)
        flag = {'--ct_mapping', ct_mapping};
    else
        flag = {};
    end
end

function flag = bb_petra_flag(use_petra)
    if use_petra
        flag = {'--use_petra'};
    else
        flag = {};
    end
end

function flag = bb_no_tight_flag(tight_narrow_beam)
    if ~tight_narrow_beam
        flag = {'--no_tight_narrow_beam'};
    else
        flag = {};
    end
end

function flag = bb_wateronly_flag(is_water)
    if is_water
        flag = {'--water_only'};
    else
        flag = {};
    end
end

function flag = bb_rayleigh_water_flag(use_rayleigh)
    % Pass --no_rayleigh_for_water (FDTD water) unless Rayleigh is requested.
    if ~use_rayleigh
        flag = {'--no_rayleigh_for_water'};
    else
        flag = {};
    end
end

function flag = bb_medium_json_flag(json_path)
    if ~isempty(json_path)
        flag = {'--medium_properties_json', json_path};
    else
        flag = {};
    end
end

function flag = bb_debug_flag(parameters)
% Keep BabelBrain debug-only outputs (currently: *_Sub_NORM.nii.gz, a
% brain-masked/normalized FullElasticSolution useful for visualizing the
% brain/skull mask boundary) when running in PRESTUS debug mode.
    if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'debug') ...
            && parameters.simulation.debug == 1
        flag = {'--debug'};
    else
        flag = {};
    end
end

function ver = query_babelbrain_version(python_bin, bb_runner, bb_root)
% Call run_babelbrain_case.py --get_version; parse KEY=VALUE stdout lines.
    ver = struct( ...
        'version_txt',  '', ...
        'git_describe', '', ...
        'git_commit',   '', ...
        'bb_root',      bb_root);

    root_flag = {};
    if ~isempty(bb_root)
        root_flag = {'--babelbrain_root', bb_root};
    end
    % run_command invokes python directly via argv, no shell at all.
    [status, out] = run_command([{python_bin, bb_runner, '--get_version'}, root_flag]);
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

function folder = bb_preset_folder(preset_name)
% Map a BabelBrain preset string to the on-disk Babel_<folder> directory name.
% Mirrors PRESET_MODULE_MAP in run_babelbrain_case.py.
    folder_map = struct( ...
        'CTX_500',    'CTX500', ...
        'CTX_250',    'CTX250', ...
        'CTX_250_2ch','CTX250_2ch', ...
        'DPX_500',    'DPX500', ...
        'DPXPC_300',  'DPXPC300', ...
        'Single',     'SingleTx', ...
        'BSonix',     'SingleTx' ...
    );
    if isfield(folder_map, preset_name)
        folder = folder_map.(preset_name);
    else
        folder = preset_name;   % most presets use the preset name directly
    end
end

function log_memory_summary(label, child_out, mem_matlab)
% One [memory_summary] line per BabelBrain run, so the pipeline log states the
% split without reading the per-sample [memory_bb]/[memory] lines. Prints
% nothing where neither side could measure (non-Linux).
bb = parse_memory_tags(child_out);
cg_peak = max([bb.cgroup_peak_gb, mem_matlab.cg_peak_gb]);   % max ignores NaN
if all(isnan([bb.bb_peak_gb, mem_matlab.hwm_gb, cg_peak]))
    return
end
fprintf('[memory_summary] %s: BabelBrain peak=%.1f GB | MATLAB peak=%.1f GB | cgroup peak=%.1f GB, limit=%.1f GB\n', ...
    label, bb.bb_peak_gb, mem_matlab.hwm_gb, cg_peak, mem_matlab.cg_limit_gb);
end

function rewrite_hdf5_no_blosc(h5_path, python_bin, this_dir)
% Rewrite h5_path in-place replacing blosc compression with gzip so MATLAB
% can read it without a blosc filter plugin.
rewriter = fullfile(this_dir, 'rewrite_hdf5_no_blosc.py');
if ~isfile(rewriter)
    warning('babelbrain_wrapper: rewrite_hdf5_no_blosc.py not found at %s — blosc read may fail.', rewriter);
    return
end
% run_command invokes python directly via argv, no shell at all.
[status, out] = run_command({python_bin, rewriter, h5_path});
if status ~= 0
    warning('babelbrain_wrapper: HDF5 rewrite failed (exit %d):\n%s\n  h5read may fail if blosc plugin is missing.', status, out);
end
end
