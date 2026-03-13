function [corrected_velocity, I_peak_before, I_peak_after] = fit_velocity_to_intensity(...
    parameters, profile_oneil, opt_phases, opt_velocity, desired_intensity, simulated_analytical_scaling)
% fit_velocity_to_intensity - Correct velocity so peak intensity matches desired intensity.
%
% After global search optimizes profile shape (phases + velocity), the peak
% intensity may not exactly equal desired_intensity. Since I ∝ v², we can
% analytically correct velocity: v_new = v_old * sqrt(I_desired / I_peak).
%
% Because opt_source_amp divides by simulated_analytical_scaling, the
% analytical target must be desired_intensity * scaling so that the final
% simulation's max_Isppa equals desired_intensity.
%
% Arguments:
%   parameters                   - Structure with transducer and medium parameters.
%   profile_oneil                - Structure with .axial_distance_bowl (mm from bowl).
%   opt_phases                   - Optimized phases for each element [rad].
%   opt_velocity                 - Optimized particle velocity from global search [m/s].
%   desired_intensity            - Target peak intensity [W/cm^2].
%   simulated_analytical_scaling - Ratio of simulated to analytical peak intensity.
%
% Returns:
%   corrected_velocity - Velocity adjusted to yield desired peak intensity [m/s].
%   I_peak_before      - Peak intensity before correction [W/cm^2].
%   I_peak_after       - Peak intensity after correction [W/cm^2].

    axial_position = profile_oneil.axial_distance_bowl;

    % Compute analytical pressure profile with current optimized parameters
    p_axial = focusedAnnulusONeil(...
        parameters.transducer.curv_radius_mm / 1e3, ...
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm] / 1e3, ...
        repmat(opt_velocity, 1, parameters.transducer.n_elements), ...
        [opt_phases], ...
        parameters.transducer.source_freq_hz, ...
        parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, ...
        (axial_position - 0.5) * 1e-3);

    % Convert pressure to intensity [W/cm^2]
    I_axial = p_axial .^ 2 / (2 * parameters.medium.water.sound_speed * parameters.medium.water.density) * 1e-4;

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
        warning('fit_velocity_to_intensity: Peak intensity is near zero. Returning original velocity.');
        corrected_velocity = opt_velocity;
        I_peak_after = I_peak_before;
        return;
    end

    % Analytical target accounts for scaling: simulation divides velocity by
    % scaling, so analytical peak must overshoot by that factor.
    analytical_target = desired_intensity * simulated_analytical_scaling;

    % Correct velocity analytically: I ∝ v² => v_new = v_old * sqrt(I_target / I_peak)
    correction_factor = sqrt(analytical_target / I_peak_before);
    corrected_velocity = opt_velocity * correction_factor;
    I_peak_after = analytical_target;

    % Warn if corrected velocity exceeds upper bound
    if isfield(parameters.calibration, 'opt_upper_velocity') && ...
            corrected_velocity > parameters.calibration.opt_upper_velocity
        warning('fit_velocity_to_intensity: Corrected velocity (%.4f m/s) exceeds opt_upper_velocity (%.4f m/s).', ...
            corrected_velocity, parameters.calibration.opt_upper_velocity);
    end

    fprintf('Velocity correction: %.4f -> %.4f m/s (factor: %.4f)\n', ...
        opt_velocity, corrected_velocity, correction_factor);
    fprintf('Analytical target: %.2f W/cm^2 (desired_intensity * scaling = %.2f * %.4f)\n', ...
        analytical_target, desired_intensity, simulated_analytical_scaling);
    fprintf('Expected max_Isppa in simulation: %.2f W/cm^2\n', desired_intensity);

end
