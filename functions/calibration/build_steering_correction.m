function [sc, T] = build_steering_correction(parameters, tr, meas_file, opts)
% BUILD_STEERING_CORRECTION  Steering correction of a matrix array from hydrophone measurements
%
% Reads measured focal peaks at a set of steering targets
% (READ_STEERING_MEASUREMENTS), models the focal peak at the same targets
% and element amplitude (MATRIX_FOCAL_PEAK; or the file's
% model_pressure_mpa column), fits measured / modelled
% (FIT_STEERING_CORRECTION) and optionally writes a library entry that
% RESOLVE_TRANSDUCER_FROM_SERIAL attaches to the transducer, so that
% CREATE_MATRIX_KARRAY scales the drive by k(target).
%
% The element amplitude used for the model (tr.matrix.elem_p_amp(1)) and
% the measurement drive voltage are stored as elem_p_amp_ref and
% drive_ref_v: elem_p_amp = elem_p_amp_ref * V / drive_ref_v is the nominal
% amplitude for drive voltage V. Rows measured at different voltages are
% scaled to the median voltage first (pressure proportional to voltage).
%
% Use as:
%   [sc, T] = build_steering_correction(parameters, tr, meas_file)
%   [sc, T] = build_steering_correction(parameters, tr, meas_file, opts)
%
% Input:
%   parameters - PRESTUS config (medium_properties.water)
%   tr         - matrix transducer struct (geometry; freq_hz; elem_p_amp)
%   meas_file  - measurement CSV (READ_STEERING_MEASUREMENTS)
%   opts       - optional: .fit (options for FIT_STEERING_CORRECTION),
%                .peak (options for MATRIX_FOCAL_PEAK), .out_file (library
%                YAML to write, e.g. config/calibration/<serial>[_<ds>].yaml),
%                .serial, .driving_system (written to the entry's meta)
%
% Output:
%   sc - steering correction (FIT_STEERING_CORRECTION) plus elem_p_amp_ref,
%        drive_ref_v, freq_hz, model_source, measurements, created_at
%   T  - measurement table with model_pressure_pa and ratio (measured /
%        modelled, after scaling to drive_ref_v)
%
% See also: FIT_STEERING_CORRECTION, MATRIX_FOCAL_PEAK, RESOLVE_TRANSDUCER_FROM_SERIAL

arguments
    parameters
    tr (1,1) struct
    meas_file {mustBeTextScalar}
    opts (1,1) struct = struct()
end

    if nargin < 4; opts = struct(); end
    T = read_steering_measurements(meas_file);
    X = [T.target_x_mm, T.target_y_mm, T.target_z_mm];
    p_meas = T.peak_pressure_mpa * 1e6;

    drive_ref_v = NaN;
    if ismember('drive_v', T.Properties.VariableNames) && all(isfinite(T.drive_v))
        drive_ref_v = median(T.drive_v);
        p_meas = p_meas .* drive_ref_v ./ T.drive_v;
    end

    if ismember('model_pressure_mpa', T.Properties.VariableNames)
        p_model = T.model_pressure_mpa * 1e6;
        model_source = 'file';
    else
        peak_opts = struct(); if isfield(opts, 'peak'); peak_opts = opts.peak; end
        M = matrix_focal_peak(parameters, tr, X, peak_opts);
        p_model = M.peak_pressure_pa;
        model_source = 'rayleigh';
    end
    T.model_pressure_pa = p_model;
    T.ratio = p_meas ./ p_model;

    fit_opts = struct(); if isfield(opts, 'fit'); fit_opts = opts.fit; end
    sc = fit_steering_correction(X, p_model, p_meas, fit_opts);
    sc.elem_p_amp_ref = tr.matrix.elem_p_amp(1);
    sc.drive_ref_v = drive_ref_v;
    sc.freq_hz = tr.freq_hz;
    sc.model_source = model_source;
    [~, fn, ext] = fileparts(meas_file);
    sc.measurements = [fn ext];
    sc.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd'));

    fprintf(['Steering correction (%d targets, measured/modelled %.3f-%.3f): leave-one-out median error ' ...
        '%.1f%% with one scale factor, %.1f%% with a quadratic map -> %s.\n'], sc.n_targets, sc.ratio_range, ...
        sc.cv.scale_med_pct, sc.cv.quadratic_med_pct, sc.model);

    if isfield(opts, 'out_file') && ~isempty(opts.out_file)
        lib = struct();
        lib.meta = struct('tran_serial', get_opt(opts, 'serial', ''), ...
            'driving_system', get_opt(opts, 'driving_system', ''), ...
            'calibration_type', 'steering_map', 'created_at', sc.created_at);
        lib.steering_correction = sc;
        yaml.dumpFile(opts.out_file, lib);
        fprintf('Steering correction written to %s\n', opts.out_file);
    end
end

function v = get_opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
