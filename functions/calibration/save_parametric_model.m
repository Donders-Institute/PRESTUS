function save_parametric_model(parameters, opt_params_raw, precession_mode)
% SAVE_PARAMETRIC_MODEL  Save compact precession-model calibration output
%
% Called by calibration_transducer when calibration_mode='parametric'.
% Writes a YAML sidecar containing the fitted precession parameters and
% amplitude, alongside the standard per-element CSV entry so the output
% remains compatible with downstream PRESTUS loaders.
%
% Use as:
%   save_parametric_model(parameters, opt_params_raw, precession_mode)
%
% Input:
%   parameters      - PRESTUS config (same as passed to save_optimized_values)
%   opt_params_raw  - raw optimisation vector from perform_global_search:
%                     'linear':   [phase_start_rad, phase_step_rad, velocity_m_s]
%                     'monotonic':[phase_start_rad, delta_1..delta_{N-1}, velocity_m_s]
%   precession_mode - 'linear' or 'monotonic'
%
% See also: SAVE_OPTIMIZED_VALUES, PERFORM_GLOBAL_SEARCH, CALIBRATION_TRANSDUCER

arguments
    parameters      (1,1) struct
    opt_params_raw  (1,:) {mustBeNumeric}
    precession_mode (1,:) char
end

    %% Build compact model struct
    model.precession_mode          = precession_mode;
    model.focal_distance_ep_mm     = parameters.calibration.desired_focal_distance_ep;
    model.desired_intensity_w_cm2  = parameters.calibration.desired_intensity;
    model.elem_amp                 = double(parameters.transducer.annular.elem_amp(1));
    model.calibration_date         = datestr(now, 'yyyy-mm-dd');

    switch precession_mode
        case 'linear'
            model.phase_start_deg = opt_params_raw(1) / pi * 180;
            model.phase_step_deg  = opt_params_raw(2) / pi * 180;
            model.velocity_m_s    = opt_params_raw(end);
        case 'monotonic'
            n_elem = parameters.transducer.annular.elem_n;
            model.phase_start_deg    = opt_params_raw(1) / pi * 180;
            model.phase_increments_deg = opt_params_raw(2:n_elem) / pi * 180;
            model.velocity_m_s       = opt_params_raw(end);
    end

    %% Write YAML sidecar
    if mod(parameters.calibration.desired_focal_distance_ep, 1) == 0
        yaml_file = sprintf('%s-F%.0fmm-I%.0fwpercm2_parametric.yaml', ...
            parameters.calibration.equipment_name, ...
            parameters.calibration.desired_focal_distance_ep, ...
            parameters.calibration.desired_intensity);
    else
        yaml_file = sprintf('%s-F%.1fmm-I%.0fwpercm2_parametric.yaml', ...
            parameters.calibration.equipment_name, ...
            parameters.calibration.desired_focal_distance_ep, ...
            parameters.calibration.desired_intensity);
    end
    yaml_path = fullfile(parameters.calibration.path_output_profiles, yaml_file);

    data = struct('parametric_model', model, 'transducer', parameters.transducer);
    yaml.dumpFile(yaml_path, data);
    fprintf('Parametric model saved: %s\n', yaml_path);

    %% Also write per-element CSV entry for compatibility with set_real_phases
    save_optimized_values(parameters);

end
