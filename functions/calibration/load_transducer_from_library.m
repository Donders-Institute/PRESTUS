function parameters = load_transducer_from_library(combo_name, focal_distance_ep, desired_intensity, equip_param, library_path, condition)
% LOAD_TRANSDUCER_FROM_LIBRARY  Load calibrated transducer parameters from library
%
% Reads the calibration library YAML for the given equipment combination and
% the geometry from equip_param (load_equipment_config), and returns element
% phases and amplitude for the requested focal depth and intensity.
%
% Phases are always geometric steering plus a per-element offset:
%
%   phase(d) = set_real_phases(d) + offset(d)
%
%   offset(d) = static offset  (global_model.elem_phase_correction_deg;
%                               depth-independent, shipped default)
%               blended towards
%               refinement offsets (refinement.elem_phase_correction_deg at
%                               refinement.depths_ep_mm; optional,
%                               depth-specific fits to hydrophone data)
%
% Within the span of refined depths the refinement offsets are interpolated
% (per element, on the unit circle). Outside it the nearest refinement offset
% fades back to the static offset over refinement.blend_mm (default 5 mm),
% so phases are continuous in depth.
%
% Legacy per-depth forms are read as refinements: global_model.depths_ep_mm +
% elem_phase_deg (absolute phases) and calibration.focal_depths.f* entries
% (linear/monotonic parametric phases). Their absolute phases are converted
% to offsets relative to geometric steering at each depth.
%
% Use as:
%   parameters = load_transducer_from_library(combo_name, focal_distance_ep, ...
%                    desired_intensity, equip_param)
%   parameters = load_transducer_from_library(combo_name, focal_distance_ep, ...
%                    desired_intensity, equip_param, library_path)
%   parameters = load_transducer_from_library(..., library_path, condition)
%
% With condition (EXPERIMENTAL, not validated) the phases, gains and drive
% come from the dephased condition lib.conditions.<condition>
% (LIBRARY_CONDITION) instead of a focal depth; focal_distance_ep may then be
% empty and defaults to the condition's nominal setpoint.
%
% Input:
%   combo_name         - equipment combination key, e.g.
%                        'IS_PCD15473_01001_IGT_32_ch_comb_10_ch'
%   focal_distance_ep  - desired focal distance from exit plane [mm]
%   desired_intensity  - target ISPPA [W/cm²]
%   equip_param        - struct from load_equipment_config()
%   library_path       - (optional) path to transducer library folder;
%                        defaults to config/calibration/ under PRESTUS root
%
% Output:
%   parameters - struct with:
%     .transducer              — geometry + calibrated phases and amplitude
%     .transducer.tx_calibration.isppa_wcm2        — ISPPA the elem_p_amp produces [W/cm²]
%     .transducer.tx_calibration.focal_depth_ep_mm — requested depth [mm]
%     .transducer.tx_calibration.phase_source      — 'library'
%     .transducer.tx_calibration.phase_model       — 'geometric' | 'static' | 'refined' | 'blend' | 'condition'
%     .transducer.tx_calibration.condition         — dephased condition label ('' = focal depth)
%     .transducer.tx_calibration.refinement_weight — 0 (static) … 1 (refined)
%     .transducer.tx_calibration.isppa_source      — 'library'
%     .transducer.tx_calibration.amp_kwave_calibrated — true iff elem_p_amp at this depth
%                                                     was calibrated by KWAVE_CALIBRATE_DRIVE
%                                                     to a free-water k-Wave simulation,
%                                                     rather than being analytical-only
%     .transducer.tx_calibration.amp_kwave_nearest_depth_mm/.amp_kwave_r2/.amp_kwave_focal_err_mm
%                                                   — nearest k-Wave drive calibration, for context even
%                                                     when amp_kwave_calibrated is false
%     (.amp_kwave_verified / .amp_kwave_ref_depth_mm: deprecated aliases of
%      .amp_kwave_calibrated / .amp_kwave_nearest_depth_mm, to be removed)
%     .transducer.tx_calibration.library_key       — combo name used for lookup
%     .transducer.tx_calibration.library_meta      — provenance from the library YAML
%
% A warning (load_transducer_from_library:uncalibratedAmplitude) is issued when an
% amplitude is resolved but was never k-Wave-calibrated — see KWAVE_CALIBRATE_DRIVE.
%
% Static-model fields written by FIT_TRANSDUCER_MODEL (all optional):
%   global_model.steering_sound_speed — sound speed of the steering law [m/s]
%                                       (geometric delays; default 1500)
%   global_model.elem_amp_weight      — per-element gain (mean 1) applied to elem_p_amp
%   refinement.elem_p_amp, .isppa_wcm2 — per-depth drive of the refinements
%                                       (LIBRARY_ADD_REFINEMENTS); interpolated in
%                                       depth and blended with the global amplitude
%                                       by the phase weight (tx_calibration.amplitude_source)
%   global_model.amplitude_by_depth   — {depths_ep_mm, elem_p_amp, isppa_wcm2}:
%                                       drive that gave isppa_wcm2 in water at each
%                                       measured depth; interpolated in depth and
%                                       scaled as sqrt(I). Without it, the
%                                       depth-independent amplitude_scaling is used.
%
% The depth logic itself is LIBRARY_MODEL_AT_DEPTH (in-memory entry, no file
% lookup); this function adds the file lookup, validation and bookkeeping.
%
% See also: LIBRARY_MODEL_AT_DEPTH, SET_REAL_PHASES, UPDATE_TRANSDUCER_LIBRARY,
%           REFINE_TRANSDUCER_LIBRARY, BUILD_TRANSDUCER_LIBRARY

    if nargin < 5 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'calibration');
    end
    if nargin < 6; condition = ''; end
    condition = char(condition);

    %% Load calibration library entry
    yaml_path = fullfile(library_path, [combo_name, '.yaml']);
    if ~isfile(yaml_path)
        error('load_transducer_from_library: library file not found:\n  %s', yaml_path);
    end
    lib = yaml.loadFile(yaml_path, 'ConvertToArray', true);
    if ~isstruct(lib) || ~isfield(lib, 'meta')
        error('load_transducer_from_library:missingMeta', ...
            'Library file is missing a ''meta'' section:\n  %s', yaml_path);
    end
    if ~isfield(lib.meta, 'tran_serial')
        error('load_transducer_from_library:missingSerial', ...
            'Library file meta section is missing ''tran_serial'':\n  %s', yaml_path);
    end

    %% Load geometry from equipment config
    tran_serial = lib.meta.tran_serial;
    if ~isfield(equip_param.trans, tran_serial)
        error('load_transducer_from_library: transducer ''%s'' not found in equip_param.', tran_serial);
    end
    tran        = equip_param.trans.(tran_serial);
    tran        = library_transducer(tran, lib);   % sub-rings if the entry has a ring profile
    n_elem      = tran.transducer.annular.elem_n;
    tran.n_elem = n_elem;
    parameters.transducer = tran.transducer;

    %% Phases, gains and drive at the requested depth
    if isempty(condition)
        m = library_model_at_depth(lib, tran, focal_distance_ep, desired_intensity, yaml_path);
    else
        assert(~isfield(tran.transducer.annular, 'elem_weight'), 'load_transducer_from_library:conditionProfile', ...
            'Dephased conditions are not supported for entries with a ring profile (%s).', yaml_path);
        m = library_condition(lib, tran, condition, desired_intensity, yaml_path);
        if isempty(focal_distance_ep) || isnan(focal_distance_ep)
            focal_distance_ep = m.setpoint_ep_mm;
        end
    end
    parameters.transducer.annular.elem_phase_rad = m.elem_phase_rad;
    parameters.transducer.annular.elem_phase_deg = m.elem_phase_deg;
    if ~isempty(m.elem_p_amp)
        parameters.transducer.annular.elem_p_amp = m.elem_p_amp;
    end
    % Ring profile: each sub-ring is driven with its channel's phase and drive
    % times its fixed weight, so the simulated source vibrates as the profile
    % model does (the analytical models apply annular.elem_weight instead)
    ring_profile = [];
    if isfield(parameters.transducer.annular, 'ring_profile')
        rp = parameters.transducer.annular.ring_profile;
        ring_profile = struct('n_sub', rp.n_sub, 'smooth', field_or(rp, 'smooth', NaN), ...
            'r2_profile_median', field_or(rp, 'r2_profile_median', NaN));
    end
    parameters.transducer.annular = ring_profile_drive(parameters.transducer.annular, ...
        m.elem_phase_deg, m.elem_p_amp);
    phase_model = m.phase_model;
    w           = m.refinement_weight;
    amp_source  = m.amplitude_source;
    elem_p_amp  = m.drive_pa;
    amp_kwave_calibrated = m.amp_kwave_calibrated;

    %% Focal distance bookkeeping
    parameters.transducer.focal_distance_ep   = focal_distance_ep;
    parameters.transducer.focal_distance_bowl = focal_distance_ep + ...
        annular_bowl_depth_mm(tran.transducer.annular);
    parameters.transducer.tx_calibration.isppa_wcm2        = desired_intensity;
    parameters.transducer.tx_calibration.focal_depth_ep_mm = focal_distance_ep;
    parameters.transducer.tx_calibration.condition         = condition;
    parameters.transducer.tx_calibration.phase_source      = 'library';
    parameters.transducer.tx_calibration.phase_model       = phase_model;
    parameters.transducer.tx_calibration.ring_profile      = ring_profile;   % [] = uniform rings
    parameters.transducer.tx_calibration.refinement_weight = w;
    parameters.transducer.tx_calibration.isppa_source      = 'library';
    parameters.transducer.tx_calibration.amplitude_source  = amp_source;   % 'global' | 'refinement' | 'blend'
    parameters.transducer.tx_calibration.amp_kwave_calibrated       = amp_kwave_calibrated;
    parameters.transducer.tx_calibration.amp_kwave_nearest_depth_mm = m.amp_kwave_nearest_depth_mm;
    parameters.transducer.tx_calibration.amp_kwave_r2           = m.amp_kwave_r2;
    parameters.transducer.tx_calibration.amp_kwave_focal_err_mm = m.amp_kwave_focal_err_mm;
    parameters.transducer.tx_calibration.library_key       = combo_name;
    parameters.transducer.tx_calibration.library_meta      = lib.meta;
    % Deprecated aliases (former names), kept for one release for scripts
    % that read them; prefer amp_kwave_calibrated / amp_kwave_nearest_depth_mm
    parameters.transducer.tx_calibration.amp_kwave_verified     = amp_kwave_calibrated;
    parameters.transducer.tx_calibration.amp_kwave_ref_depth_mm = m.amp_kwave_nearest_depth_mm;

    if ~isempty(elem_p_amp) && ~amp_kwave_calibrated && isempty(condition)
        if isfinite(m.amp_kwave_nearest_depth_mm)
            nearest_str = sprintf(['nearest k-Wave drive calibration is at %.1f mm from the exit plane ' ...
                '(R²=%.3f, focal error %.2f mm) but does not cover the requested depth'], ...
                m.amp_kwave_nearest_depth_mm, m.amp_kwave_r2, m.amp_kwave_focal_err_mm);
        else
            nearest_str = 'no depth in this entry has a k-Wave-calibrated drive';
        end
        warning('load_transducer_from_library:uncalibratedAmplitude', ...
            ['%s @ %.1f mm: elem_p_amp is analytical only — it was fit to the O''Neil/Rayleigh ' ...
             'model and measured hydrophone profile but never calibrated against a free-water ' ...
             'k-Wave simulation, so its absolute intensity may be biased by the analytical-vs-' ...
             'simulation offset (%s). Run KWAVE_CALIBRATE_DRIVE on this entry to calibrate it.'], ...
            combo_name, focal_distance_ep, nearest_str);
    end

    cal_type = '';
    if isfield(lib.meta, 'calibration_type'), cal_type = char(lib.meta.calibration_type); end
    if isempty(elem_p_amp)
        amp_str = 'no amplitude calibration';
    else
        amp_str = sprintf('amp %d @ %.1f W/cm²', round(elem_p_amp), desired_intensity);
    end
    fprintf('Library loaded (%s, %s): %s | focal %.1f mm | %s (w=%.2f) | %s\n', ...
        phase_model, cal_type, combo_name, focal_distance_ep, phase_model, w, amp_str);
end

function v = field_or(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
