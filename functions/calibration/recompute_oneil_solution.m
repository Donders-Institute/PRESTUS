function profile_oneil_opt = recompute_oneil_solution(parameters, profile_oneil, profile_target, opt_phases, opt_velocity)
    % Recalculate analytical solution based on optimized phases and velocity.
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % -     .calibration.desired_focal_distance_ep: Focal depth with respect to the transducer exit plane [mm].
    % -     .calibration.desired_intensity: Target intensity for optimization [W/cm^2].
    % -     .calibration.equipment_name: Name of the equipment used.
    % - profile_oneil.axial_distance_bowl: Axial position vector [mm from bowl].
    % - profile_oneil.axial_intensity: Initial O'Neil solution for intensity [W/cm2].
    % - profile_target.axial_intensity: Adjusted desired intensity profile [W/cm^2].
    % - opt_phases: Optimized phases for each transducer element [rad].
    % - opt_velocity: Optimized particle velocity [m/s].
    %
    % Returns:
    % - profile_oneil_opt.axial_intensity: Optimized O'Neil solution for intensity along the beam axis [W/cm2].
    % - profile_oneil_opt.axial_distance_bowl: Axial position vector [mm from bowl].

    % Get axial position
    axial_position = profile_oneil.axial_distance_bowl;

    % Compute optimized analytical pressure profile
    p_axial_oneil_opt = focusedAnnulusONeil(...
        parameters.transducer.curv_radius_mm / 1e3, ...
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm] / 1e3, ...
        repmat(opt_velocity, 1, parameters.transducer.n_elements), ...
        [opt_phases], ...
        parameters.transducer.source_freq_hz, ...
        parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, ...
        (axial_position - 0.5) * 1e-3);
    
    % Convert pressure to intensity
    i_axial_oneil = profile_oneil.axial_intensity;
    i_axial_oneil_opt = p_axial_oneil_opt .^ 2 / (2 * parameters.medium.water.sound_speed * parameters.medium.water.density) * 1e-4;

    % Plot comparison of profiles
    figure('Position', [10, 10, 900, 500]);
    plot(axial_position, i_axial_oneil, ...
        'LineWidth', 1, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Original Profile (Analytical)');
    hold on;
    plot(axial_position, i_axial_oneil_opt, ...
        'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'Optimized Profile (Analytical)');
    plot(axial_position, profile_target.axial_intensity, ...
        'LineWidth', 2, 'Color', [1 0 0], 'DisplayName', 'Desired Profile');
    hold off;
    xline(parameters.expected_focal_distance_bowl, '--', 'DisplayName', 'Focal Point wrt Mid-Bowl');
    yline(parameters.calibration.desired_intensity, '--', 'DisplayName', 'Desired Intensity');

    xlabel('Distance wrt Mid-Bowl of Transducer [mm]');
    ylabel('Intensity [W/cm^2]');
    legend('Location', 'NorthEast');
    legend('boxoff')
    title('Pressure Along the Beam Axis');
    ylim([0 inf]);
    xlim([-5 inf]);
    
    % Save the profile comparison figure
    fig_path = fullfile(parameters.outputs_folder, ...
        strcat('Recalculated_oneil_at_F_', num2str(parameters.calibration.desired_focal_distance_ep), ...
        '_at_I_', num2str(parameters.calibration.desired_intensity), ...
        '_', parameters.calibration.equipment_name, '.png'));
    saveas(gcf, fig_path);
    close(gcf);

    % exclude the near-field if requested
    i_axial_oneil_opt_summary = i_axial_oneil_opt;
    if parameters.calibration.skip_front_peak_mm ~=0
        i_remove = axial_position <= parameters.calibration.skip_front_peak_mm;
        axial_position(i_remove) = [];
        i_axial_oneil_opt_summary(i_remove) = [];
    end

    fprintf('Estimated distance to the point of maximum intensity: %.2f mm\n', ...
        axial_position(i_axial_oneil_opt_summary == max(i_axial_oneil_opt_summary)))
    try
        fprintf('Estimated distance to the center of half-maximum range: %.2f mm\n', ...
            get_flhm_center_position(axial_position, i_axial_oneil_opt_summary))
    catch
        warning('Could not estimated distance to the center of half-maximum range');
    end

    %% Collect in structure

    profile_oneil_opt.axial_intensity = i_axial_oneil_opt;
    profile_oneil_opt.axial_distance_bowl = axial_position;

    
end