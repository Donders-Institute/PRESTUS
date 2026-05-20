function [corrected_velocity, I_peak_before, I_peak_after] = fit_velocity_to_intensity(...
    parameters, profile_analytical, opt_phases, opt_velocity, desired_intensity)
% FIT_VELOCITY_TO_INTENSITY  Analytically correct particle velocity to match desired peak intensity
%
% After global search optimises profile shape (phases and velocity), the
% peak intensity may not perfectly match the desired_intensity.
%
% If only amplitude needs to be calibrated, this also does not require a
% recalibration of phases.
%
% Since I ∝ v², velocity for a target amplitude can be computed analytically:
% v_new = v_old * sqrt(I_desired / I_peak).
%
% The forward model is selected by parameters.calibration.forward_model:
%   'oneil'    (default) — O'Neil closed-form via focusedAnnulusONeil
%   'rayleigh'           — Rayleigh–Sommerfeld via rayleigh_axial_intensity
%
% Use as:
%   [corrected_velocity, I_peak_before, I_peak_after] = ...
%       fit_velocity_to_intensity(parameters, profile_analytical, opt_phases, ...
%                                 opt_velocity, desired_intensity)
%
% Input:
%   parameters          - PRESTUS config with transducer.annular geometry,
%                         medium_properties.water, calibration.skip_front_peak_mm [mm],
%                         calibration.forward_model ('oneil'|'rayleigh', optional)
%   profile_analytical  - struct with axial_distance_bowl [mm]
%   opt_phases          - optimised phases per element [rad]
%   opt_velocity        - optimised particle velocity from global search [m/s]
%   desired_intensity   - target peak intensity [W/cm²]
%
% Output:
%   corrected_velocity - velocity adjusted to yield desired peak intensity [m/s]
%   I_peak_before      - peak intensity before correction [W/cm²]
%   I_peak_after       - peak intensity after correction [W/cm²]
%
% See also: PERFORM_GLOBAL_SEARCH, COMPUTE_ANALYTICAL_SOLUTION, CALIBRATION_TRANSDUCER,
%           RAYLEIGH_AXIAL_INTENSITY

arguments
    parameters          (1,1) struct
    profile_analytical  (1,1) struct
    opt_phases          (1,:) {mustBeNumeric}
    opt_velocity        (1,1) {mustBeNumeric}
    desired_intensity   (1,1) {mustBeNumeric}
end

    if isfield(parameters.calibration, 'forward_model')
        forward_model = parameters.calibration.forward_model;
    else
        forward_model = 'oneil';
    end

    axial_position = profile_analytical.axial_distance_bowl;

    % Evaluate forward model with current optimised parameters
    if strcmp(forward_model, 'rayleigh')
        I_axial = rayleigh_axial_intensity(opt_phases, opt_velocity, parameters, axial_position);
    elseif strcmp(forward_model, 'oneil')
        p_axial = focusedAnnulusONeil( ...
            parameters.transducer.annular.curv_radius_mm / 1e3, ...
            [parameters.transducer.annular.elem_id_mm; parameters.transducer.annular.elem_od_mm] / 1e3, ...
            repmat(opt_velocity, 1, parameters.transducer.annular.elem_n), ...
            [opt_phases], ...
            parameters.transducer.freq_hz, ...
            parameters.medium_properties.water.sound_speed, ...
            parameters.medium_properties.water.density, ...
            (axial_position - 0.5) * 1e-3);
        I_axial = p_axial .^ 2 / (2 * parameters.medium_properties.water.sound_speed * parameters.medium_properties.water.density) * 1e-4;
    else
        error('fit_velocity_to_intensity: unknown forward_model ''%s''; expected ''oneil'' or ''rayleigh''.', forward_model);
    end

    % Exclude near-field if requested
    I_filtered = I_axial;
    axial_pos_filtered = axial_position;
    if parameters.calibration.skip_front_peak_mm ~= 0
        i_remove = axial_pos_filtered <= parameters.calibration.skip_front_peak_mm;
        axial_pos_filtered(i_remove) = [];
        I_filtered(i_remove) = [];
    end

    % Find peak intensity
    I_peak_before = max(I_filtered);

    % Guard against zero/near-zero peak
    if I_peak_before < eps
        warn('fit_velocity_to_intensity: Peak intensity is near zero. Returning original velocity.');
        corrected_velocity = opt_velocity;
        I_peak_after = I_peak_before;
        return;
    end

    % Analytical target
    analytical_target = desired_intensity;

    % Correct velocity analytically: I ∝ v² => v_new = v_old * sqrt(I_target / I_peak)
    correction_factor = sqrt(analytical_target / I_peak_before);
    corrected_velocity = opt_velocity * correction_factor;
    I_peak_after = analytical_target;

    % Warn if corrected velocity exceeds upper bound
    if isfield(parameters.calibration, 'opt_upper_velocity') && ...
            corrected_velocity > parameters.calibration.opt_upper_velocity
        warn('fit_velocity_to_intensity: Corrected velocity (%.4f m/s) exceeds opt_upper_velocity (%.4f m/s).', ...
            corrected_velocity, parameters.calibration.opt_upper_velocity);
    end

    fprintf('Velocity correction: %.4f -> %.4f m/s (factor: %.4f)\n', ...
        opt_velocity, corrected_velocity, correction_factor);
    fprintf('Analytical target: %.2f W/cm^2 (desired_intensity = %.2f)\n', ...
        analytical_target, desired_intensity);
    fprintf('Expected Isppa in simulation: %.2f W/cm^2\n', desired_intensity);

end
