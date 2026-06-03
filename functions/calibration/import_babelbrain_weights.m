function correction_deg = import_babelbrain_weights(...
    h5_path, equipment_yaml_path, tran, ref_depth_ep, parameters)
% IMPORT_BABELBRAIN_WEIGHTS  Import a BabelBrain calibration into PRESTUS
%
% Reads per-element phase weights from a BabelBrain OptimizedWeightsFile
% (.h5), derives the per-element hardware correction relative to the
% geometric phases at a reference depth, and writes the correction into
% the transducer equipment YAML via SAVE_ELEM_CORRECTION so it is picked
% up automatically by load_equipment_config / load_transducer_from_library.
%
% BabelBrain stores calibrated element weights as complex64 phasors where
% angle(w_i) encodes a DELTA correction relative to geometric phases, not
% an absolute phase. BabelBrain's TxCalibration builds a transfer matrix A
% with geometric Rayleigh phases already applied as source excitation, and
% optimises correction weights b starting from zero phases. Consequently:
%   angle(w_i) = opt_phase_i - geo_phase_i    [hardware delta, rad]
%   abs(w_i)   = amplitude scaling (ignored — PRESTUS uses elem_amp)
%
% To derive the PRESTUS correction (also a delta from geo), the absolute
% phases are first reconstructed:
%   abs_opt_deg = mod(geo_deg + angle(w)*180/pi, 360)
% then passed to SAVE_ELEM_CORRECTION which computes:
%   correction_deg = mod(unwrap(abs_opt) - unwrap(geo), 360) = angle(w)*180/pi
%
% The two-step round-trip is equivalent to using angle(w) directly but
% keeps the code path consistent with PRESTUS global-search calibrations.
%
% Use as:
%   correction_deg = import_babelbrain_weights( ...
%       h5_path, equipment_yaml_path, tran, ref_depth_ep, parameters)
%
% Input:
%   h5_path             - path to BabelBrain OptimizedWeightsFile (.h5)
%                         Must contain /CALIBRATION/real and /CALIBRATION/imag
%                         datasets, each [N_elem x 1] float32.
%   equipment_yaml_path - path to the transducer equipment YAML to update
%                         (e.g. {equip_param.path}/{tran_serial}.yaml)
%   tran                - transducer struct from load_equipment_config;
%                         used to compute geometric phases via
%                         generate_tran_ini_from_geometry / set_real_phases.
%   ref_depth_ep        - focal depth at which the BabelBrain calibration
%                         was performed, measured from the exit plane [mm].
%                         This must match the depth used in BabelBrain's
%                         TxCalibration run.
%   parameters          - PRESTUS config struct; uses
%                         medium_properties.water.sound_speed.
%
% Output:
%   correction_deg - [1 x N_elem] per-element phase correction [°];
%                    also written to equipment_yaml_path.
%
% See also: EXPORT_BABELBRAIN_WEIGHTS, SAVE_ELEM_CORRECTION,
%           CALIBRATION_EXTRACT_CORRECTION, SET_REAL_PHASES

arguments
    h5_path             (1,:) char
    equipment_yaml_path (1,:) char
    tran                (1,1) struct
    ref_depth_ep        (1,1) {mustBeNumeric, mustBePositive}
    parameters          (1,1) struct
end

    if ~isfile(h5_path)
        error('import_babelbrain_weights: file not found: %s', h5_path);
    end
    if ~isfile(equipment_yaml_path)
        error('import_babelbrain_weights: equipment YAML not found: %s', equipment_yaml_path);
    end

    % ── Read phasors from BabelBrain HDF5 ────────────────────────────────
    try
        w_real = h5read(h5_path, '/CALIBRATION/real');   % [N_elem x 1] single
        w_imag = h5read(h5_path, '/CALIBRATION/imag');   % [N_elem x 1] single
    catch ME
        error('import_babelbrain_weights: could not read /CALIBRATION/real|imag from %s.\n%s', ...
            h5_path, ME.message);
    end

    n_elem_h5 = numel(w_real);
    n_elem_tran = tran.n_elem;
    if n_elem_h5 ~= n_elem_tran
        error(['import_babelbrain_weights: HDF5 has %d elements but transducer ' ...
            'expects %d. Check that h5_path matches tran.'], n_elem_h5, n_elem_tran);
    end

    phasors    = double(w_real(:)') + 1i * double(w_imag(:)');
    opt_phases_rad = angle(phasors);              % [1 x N_elem], range (-pi, pi]
    opt_phases_deg = opt_phases_rad * 180 / pi;

    fprintf('BabelBrain weights loaded: %d elements from %s\n', n_elem_h5, h5_path);
    fprintf('  Phase range: [%.1f, %.1f] deg\n', min(opt_phases_deg), max(opt_phases_deg));
    amplitude = abs(phasors);
    if any(abs(amplitude - 1) > 0.05)
        fprintf(['  Note: phasors have non-unit amplitude (range [%.3f, %.3f]).\n' ...
                 '  Amplitude scaling is ignored — PRESTUS uses elem_amp for drive level.\n'], ...
            min(amplitude), max(amplitude));
    end

    % ── Compute geometric phases at the reference depth ───────────────────
    tran_ini_data  = generate_tran_ini_from_geometry(tran);
    geo_phases_deg = set_real_phases(tran_ini_data, tran, ref_depth_ep, parameters);

    if numel(geo_phases_deg) ~= n_elem_tran
        error(['import_babelbrain_weights: set_real_phases returned %d phases ' ...
            'but transducer has %d elements.'], numel(geo_phases_deg), n_elem_tran);
    end

    % ── Derive and save correction ────────────────────────────────────────
    % BabelBrain stores weights as corrections relative to geometric phases:
    % the transfer matrix A is built with geometric phases already applied,
    % so angle(b) is a delta from geometric, not an absolute phase.
    % Reconstruct absolute phases before calling save_elem_correction, which
    % internally subtracts geo — otherwise geometric phases would be subtracted twice.
    opt_phases_abs_deg = mod(geo_phases_deg + opt_phases_deg, 360);
    save_elem_correction(equipment_yaml_path, opt_phases_abs_deg, geo_phases_deg, ref_depth_ep);

    geo_unwrapped  = unwrap(geo_phases_deg       * pi/180) * 180/pi;
    opt_unwrapped  = unwrap(opt_phases_abs_deg   * pi/180) * 180/pi;
    correction_deg = mod(opt_unwrapped - geo_unwrapped, 360);

end
