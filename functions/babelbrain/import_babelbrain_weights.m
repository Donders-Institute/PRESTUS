function correction_deg = import_babelbrain_weights(...
    h5_path, equipment_yaml_path, tran, ref_depth_ep, parameters, elem_counts)
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
%   abs(w_i)   = amplitude scaling (ignored — PRESTUS uses elem_p_amp)
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
%   correction_deg = import_babelbrain_weights(..., elem_counts)
%
% Input:
%   h5_path             - path to BabelBrain OptimizedWeightsFile (.h5).
%                         Accepted /CALIBRATION layouts:
%                           - BabelBrain native (SaveToH5py, TxCalibration.py):
%                             1-D compound {r, i} dataset (h5py complex64 or
%                             complex128), or a real dataset (FitType 'Amp',
%                             read as zero phase);
%                           - legacy PRESTUS export: group with float32
%                             datasets /CALIBRATION/real and /CALIBRATION/imag.
%                         BabelBrain compresses arrays above 1 KiB with Blosc,
%                         which MATLAB cannot read: convert such files first
%                         with functions/babelbrain/rewrite_hdf5_no_blosc.py.
%   equipment_yaml_path - path to the transducer equipment YAML to update
%                         (e.g. {equip_param.path}/{tran_serial}.yaml)
%   tran                - transducer struct from load_equipment_config;
%                         used to compute geometric phases via set_real_phases.
%   ref_depth_ep        - focal depth at which the BabelBrain calibration
%                         was performed, measured from the exit plane [mm].
%                         This must match the depth used in BabelBrain's
%                         TxCalibration run.
%   parameters          - PRESTUS config struct; uses
%                         medium_properties.water.sound_speed.
%   elem_counts         - (optional) [1 x N_ring] sub-elements per ring, i.e.
%                         BabelBrain's Tx['elemdims'][:,0]. TxCalibration fits
%                         one weight per SUB-element; when the file holds
%                         sum(elem_counts) weights they are reduced to one per
%                         ring by the phase of their complex sum. Not needed
%                         when the file holds one weight per ring.
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
    elem_counts         (1,:) {mustBeNumeric, mustBeInteger, mustBePositive} = []
end

    if ~isfile(h5_path)
        error('import_babelbrain_weights: file not found: %s', h5_path);
    end
    if ~isfile(equipment_yaml_path)
        error('import_babelbrain_weights: equipment YAML not found: %s', equipment_yaml_path);
    end

    % ── Read phasors from BabelBrain HDF5 ────────────────────────────────
    phasors = read_calibration_phasors(h5_path);     % [1 x N] double complex

    n_elem_h5 = numel(phasors);
    n_elem_tran = tran.n_elem;
    if n_elem_h5 ~= n_elem_tran && ~isempty(elem_counts) && n_elem_h5 == sum(elem_counts)
        % Per-sub-element weights (TxCalibration) → one phasor per ring.
        if numel(elem_counts) ~= n_elem_tran
            error(['import_babelbrain_weights: elem_counts has %d rings but transducer ' ...
                'has %d elements.'], numel(elem_counts), n_elem_tran);
        end
        ring_of_sub = repelem(1:n_elem_tran, elem_counts);
        phasors = accumarray(ring_of_sub(:), phasors(:)).' ./ elem_counts;
        fprintf('  Reduced %d sub-element weights to %d rings.\n', n_elem_h5, n_elem_tran);
        n_elem_h5 = n_elem_tran;
    end
    if n_elem_h5 ~= n_elem_tran
        error(['import_babelbrain_weights: HDF5 has %d elements but transducer ' ...
            'expects %d. Check that h5_path matches tran; for per-sub-element ' ...
            'BabelBrain weights pass elem_counts.'], n_elem_h5, n_elem_tran);
    end
    opt_phases_rad = angle(phasors);              % [1 x N_elem], range (-pi, pi]
    opt_phases_deg = opt_phases_rad * 180 / pi;

    fprintf('BabelBrain weights loaded: %d elements from %s\n', n_elem_h5, h5_path);
    fprintf('  Phase range: [%.1f, %.1f] deg\n', min(opt_phases_deg), max(opt_phases_deg));
    amplitude = abs(phasors);
    if any(abs(amplitude - 1) > 0.05)
        fprintf(['  Note: phasors have non-unit amplitude (range [%.3f, %.3f]).\n' ...
                 '  Amplitude scaling is ignored — PRESTUS uses elem_p_amp for drive level.\n'], ...
            min(amplitude), max(amplitude));
    end

    % ── Compute geometric phases at the reference depth ───────────────────
    geo_phases_deg = set_real_phases(tran, ref_depth_ep, parameters);

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


function phasors = read_calibration_phasors(h5_path)
% Read /CALIBRATION as a [1 x N] complex double row, accepting the BabelBrain
% native layout (compound r/i or real dataset) and the legacy PRESTUS group.
    try
        info = h5info(h5_path, '/');
    catch ME
        error('import_babelbrain_weights: could not open %s.\n%s', h5_path, ME.message);
    end
    is_group = ~isempty(info.Groups)   && any(strcmp({info.Groups.Name}, '/CALIBRATION'));
    is_dset  = ~isempty(info.Datasets) && any(strcmp({info.Datasets.Name}, 'CALIBRATION'));
    try
        if is_group
            w_real = h5read(h5_path, '/CALIBRATION/real');
            w_imag = h5read(h5_path, '/CALIBRATION/imag');
        elseif is_dset
            raw = h5read(h5_path, '/CALIBRATION');
            if isstruct(raw)
                % h5py complex: compound with fields 'r' and 'i'
                if ~all(isfield(raw, {'r', 'i'}))
                    error('compound /CALIBRATION has fields %s, expected r and i.', ...
                        strjoin(fieldnames(raw), ', '));
                end
                w_real = raw.r;
                w_imag = raw.i;
            else
                w_real = raw;
                w_imag = zeros(size(raw));
            end
        else
            error('no /CALIBRATION dataset or group.');
        end
    catch ME
        if contains(ME.message, {'filter', 'Blosc', 'blosc'}, 'IgnoreCase', true)
            error(['import_babelbrain_weights: %s is Blosc-compressed, which MATLAB ' ...
                'cannot read. Convert it first: python rewrite_hdf5_no_blosc.py <file>\n%s'], ...
                h5_path, ME.message);
        end
        error('import_babelbrain_weights: could not read /CALIBRATION from %s.\n%s', ...
            h5_path, ME.message);
    end
    phasors = double(w_real(:)') + 1i * double(w_imag(:)');
end
