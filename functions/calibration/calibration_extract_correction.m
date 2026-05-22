function correction_deg = calibration_extract_correction(...
    calibrated_yaml_path, parameters, tran, phase_table, ref_depth_ep)
% CALIBRATION_EXTRACT_CORRECTION  Retroactively extract per-element correction from a prior calibration
%
% For users who have already completed a full calibration at a reference
% depth: loads the optimised phases from the saved YAML, recomputes the
% geometric phases at that depth, and saves the per-element correction to
% a YAML file suitable for use with the geometric steering calibration mode.
%
% For new calibrations, set parameters.calibration.save_elem_correction = true
% in calibration_standalone.m — the correction is then extracted automatically
% during the reference-depth run without calling this function.
%
% Use as:
%   correction_deg = calibration_extract_correction( ...
%       calibrated_yaml_path, parameters, tran, phase_table, ref_depth_ep)
%
% Input:
%   calibrated_yaml_path - path to a YAML written by save_optimized_values
%                          (contains transducer.annular.elem_phase_deg)
%   parameters           - PRESTUS config; uses medium_properties.water.sound_speed
%                          and calibration.path_output_profiles / equipment_name
%   tran                 - transducer struct (from load_transducer_from_library)
%   phase_table          - manufacturer phase data (table or ini struct)
%   ref_depth_ep         - reference focal depth w.r.t. exit plane [mm]
%
% Output:
%   correction_deg - per-element correction [°, 1xN]; also written to
%                    {path_output_profiles}/{equipment_name}_elem_correction.yaml
%
% See also: SAVE_ELEM_CORRECTION, CALIBRATION_TRANSDUCER, COMPUTE_PHASES, SET_REAL_PHASES

arguments
    calibrated_yaml_path (1,:) char
    parameters           (1,1) struct
    tran                 (1,1) struct
    phase_table
    ref_depth_ep         (1,1) {mustBeNumeric}
end

    if ~isfile(calibrated_yaml_path)
        error('calibration_extract_correction: file not found: %s', calibrated_yaml_path);
    end

    % Load optimised phases from the prior calibration YAML
    cal_data       = yaml.loadFile(calibrated_yaml_path, 'ConvertToArray', true);
    opt_phases_deg = cal_data.transducer.annular.elem_phase_deg;

    % Recompute geometric phases at the same reference depth
    geo_phases_deg = set_real_phases(phase_table, tran, ref_depth_ep, parameters);

    % Delegate to save_elem_correction for consistent computation and output
    save_elem_correction(parameters, opt_phases_deg, geo_phases_deg, ref_depth_ep);

    geo_unwrapped  = unwrap(geo_phases_deg * pi/180) * 180/pi;
    opt_unwrapped  = unwrap(opt_phases_deg * pi/180) * 180/pi;
    correction_deg = mod(opt_unwrapped - geo_unwrapped, 360);

end
