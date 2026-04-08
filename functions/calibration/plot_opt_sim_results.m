function plot_opt_sim_results(parameters, profile_target, profile_oneil, profile_oneil_opt, profile_sim, profile_sim_opt, min_err)
    % Plot optimized simulation results and compare with desired profiles
    %
    % Arguments:
    % - parameters: Structure containing optimized parameters.
    %       .calibration.path_output: Directory for saving output.
    %       .calibration.save_in_calibration_folder: Option to save data in the general PRESTUS output folder.
    %       .calibration.desired_focal_distance_ep: Focal depth with respect to the transducer exit plane [mm].
    %       .calibration.desired_intensity: Target intensity for optimization [W/cm^2].
    %       .calibration.equipment_name: Name of the equipment used.
    %       .outputs_folder: Directory containing output simulation results.
    % - profile_target
    %       .axial_distance_bowl: Axial position vector [mm from bowl].
    %       .axial_intensity: Adjusted desired intensity profile [W/cm^2].
    % - profile_oneil
    %       .axial_intensity: Original O'Neil solution for intensity [W/cm^2].
    % - profile_oneil_opt
    %       .axial_intensity: Optimized O'Neil solution for intensity [W/cm^2].
    % - profile_sim
    %       .axial_intensity: Original simulation intensity profile [W/cm^2].
    % - profile_sim_opt
    %       .axial_intensity: Optimized simulation intensity profile [W/cm^2].
    % - min_err: Minimum optimization error.
    
    %% Plot comparison of profiles

    figure('Position', [10, 10, 900, 500]);
    hold on;
    plot(profile_target.axial_distance_bowl, profile_target.axial_intensity, ...
        'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Target Profile');
    plot(profile_oneil.axial_distance_bowl, profile_oneil.axial_intensity, ...
        'LineWidth', 1, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Original (Analytical))');
    plot(profile_oneil_opt.axial_distance_bowl, profile_oneil_opt.axial_intensity, ...
        'LineWidth', 1, 'Color', [0 0 0], 'LineStyle', ':', 'DisplayName', 'Optimized (Analytical)');
    plot(profile_sim.axial_distance_bowl, profile_sim.axial_intensity, ...
        'LineWidth', 1, 'Color', [0.75 0.75 0.75], 'DisplayName', 'Original (Simulated))');
    plot(profile_sim_opt.axial_distance_bowl, profile_sim_opt.axial_intensity, ...
        'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'Optimized (Simulated)');
    hold off;

    % Add focus and intensity reference lines
    xline(parameters.transducer(1).focal_distance_bowl, '--', 'DisplayName', 'Focal Point wrt mid-bowl');
    yline(parameters.calibration.desired_intensity, '--', 'DisplayName', 'Desired Intensity');

    xlabel('Distance wrt Mid-Bowl of Transducer [mm]');
    ylabel('Intensity [W/cm^2]');
    legend('Location', 'best');
    title(sprintf('Desired vs Optimized Profiles - Optimization Error: %.4f', min_err));
    ylim([0 inf]);
    xlim([0 inf]);
    
    % Save the profile comparison figure
    fig_path = fullfile(parameters.io.outputs_folder, ...
        strcat('Opt_simulation_at_F_', num2str(parameters.calibration.desired_focal_distance_ep), ...
            '_at_I_', num2str(parameters.calibration.desired_intensity), ...
            '_', parameters.calibration.equipment_name, '.png'));
    saveas(gcf, fig_path);
    close(gcf);

    %% Plot a simplified version and save in the virtual transducer directory

    figure('Position', [10, 10, 900, 500]);
    hold on;
    plot(profile_sim_opt.axial_distance_bowl, profile_sim_opt.axial_intensity, ...
        'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'Optimized (Simulated)');
    plot(profile_target.axial_distance_bowl, profile_target.axial_intensity, ...
        'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Desired Profile');
    hold off;
    % Add focus and intensity reference lines
    xline(parameters.transducer(1).focal_distance_bowl, '--', 'DisplayName', 'Focal Point wrt mid-bowl');
    yline(parameters.calibration.desired_intensity, '--', 'DisplayName', 'Desired Intensity');
    % add expected distance
    if isfield(parameters.transducer(1), 'focal_distance_bowl')
        xline(parameters.transducer(1).focal_distance_bowl, '--', ...
            'LineWidth', 1.2, 'DisplayName', 'Expected Focal Distance (mm from bowl)', 'Color', [1 0 0]);
    end
    % add exit plane
    if isfield(parameters.transducer(1), 'focal_distance_offset')
        xline(parameters.transducer(1).focal_distance_offset, '--', ...
            'LineWidth', 1.2, 'DisplayName', 'Exit Plane');
    end
    xlabel('Distance wrt Mid-Bowl of Transducer [mm]');
    ylabel('Intensity [W/cm^2]');
    legend('Location', 'best');
    ylim([0 inf]);
    xlim([0 inf]);

    % Save the profile comparison figure
    fig_path = fullfile(parameters.calibration.path_output_profiles, ...
        strcat('OptProfile_at_F_', num2str(parameters.calibration.desired_focal_distance_ep), ...
            '_at_I_', num2str(parameters.calibration.desired_intensity), ...
            '_', parameters.calibration.equipment_name, '.png'));
    saveas(gcf, fig_path);
    close(gcf);

    %% Display summary

    parameters = focal_distance_calculation(parameters);

    max_intensity_index = find(profile_sim_opt.axial_intensity == max(profile_sim_opt.axial_intensity), 1);
    fprintf('Expected distance from bowl to the point of maximum intensity: %.2f mm\n', ...
        parameters.transducer(1).focal_distance_bowl);
    fprintf('Estimated distance from bowl to the point of maximum intensity: %.2f mm\n', ...
        profile_sim_opt.axial_distance_bowl(max_intensity_index));
    
    fprintf('Expected distance from exit plane to the point of maximum intensity: %.2f mm\n', ...
        parameters.transducer(1).focal_distance_ep);
    fprintf('Estimated distance from exit plane to the point of maximum intensity: %.2f mm\n', ...
        profile_sim_opt.axial_distance_bowl(max_intensity_index)-parameters.transducer(1).focal_distance_offset);

    
end