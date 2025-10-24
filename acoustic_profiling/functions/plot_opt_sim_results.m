function plot_opt_sim_results(opt_param, outputs_folder, sim_id, axial_position, dist_exit_plane, adjusted_profile_focus, p_axial_oneil_opt, p_axial_oneil, focus_wrt_exit_plane, desired_intensity, equipment_name, min_err)
    % Plot optimized simulation results and compare with desired profiles
    %
    % Arguments:
    % - opt_param: Structure containing optimized parameters.
    %   opt_param.calibration.path_output: Directory for saving output.
    %   opt_param.calibration.save_in_calibration_folder: Option to save data in the general PRESTUS output folder.
    % - outputs_folder: Directory containing output simulation results.
    % - sim_id: Simulation ID for loading specific results.
    % - axial_position: Axial position vector [mm].
    % - dist_exit_plane: Distance vector for desired profile focus [mm].
    % - adjusted_profile_focus: Adjusted desired intensity profile [W/cm^2].
    % - p_axial_oneil_opt: Optimized O'Neil solution for pressure [Pa].
    % - p_axial_oneil: Original O'Neil solution for pressure [Pa].
    % - focus_wrt_exit_plane: Focal depth with respect to the transducer exit plane [mm].
    % - desired_intensity: Target intensity for optimization [W/cm^2].
    % - equipment_name: Name of the equipment used.
    % - min_err: Minimum optimization error.
    
    % Load optimized simulation results    
    opt_res = load(sprintf('%s/sub-%03d_water_results%s.mat', outputs_folder, sim_id, opt_param.results_filename_affix),'sensor_data','parameters');

    % Extract maximum pressure data
    p_max = gather(opt_res.sensor_data.p_max_all);

    % Plot 2D intensity map for the focal plane
    figure;
    imagesc((1:size(p_max, 1)) * opt_res.parameters.grid_step_mm, ...
           (1:size(p_max, 3)) * opt_res.parameters.grid_step_mm, ...
           squeeze(p_max(:, opt_res.parameters.transducer.pos_grid(2), :))');
    axis image;
    colormap(getColorMap());
    xlabel('Lateral Position [mm]');
    ylabel('Axial Position [mm]');
    colorbar;
    title('Pressure for the Focal Plane')
    
    % Save the 2D intensity map figure
    if opt_param.calibration.save_in_calibration_folder
        fig_path = fullfile(opt_param.calibration.path_output, strcat('Opt_intensity_map_2D_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    else
        fig_path = fullfile(opt_param.output_location, strcat('Opt_intensity_map_2D_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    end
    saveas(gcf, fig_path);

    % Simulated pressure along the focal axis
    pred_axial_pressure_opt = squeeze(p_max(opt_res.parameters.transducer.pos_grid(1), opt_res.parameters.transducer.pos_grid(2), :));
  
    % Compute focal position relative to the mid-bowl of the transducer
    focus_wrt_mid_bowl = focus_wrt_exit_plane + (opt_res.parameters.transducer.curv_radius_mm - opt_res.parameters.transducer.dist_to_plane_mm);

    % Plot comparison of profiles
    figure('Position', [10, 10, 900, 500]);
    hold on;
    plot(axial_position, p_axial_oneil .^ 2 / (2 * opt_param.medium.water.sound_speed * opt_param.medium.water.density) * 1e-4, 'DisplayName', 'Original Simulation');
    plot(axial_position, p_axial_oneil_opt .^ 2 / (2 * opt_param.medium.water.sound_speed * opt_param.medium.water.density) * 1e-4, 'DisplayName', 'Optimized (Analytical)');
    
    % Adjust axial position for simulated results
    sim_res_axial_position = axial_position - (opt_res.parameters.transducer.pos_grid(3) - 1) * 0.5;
    plot(sim_res_axial_position, pred_axial_pressure_opt .^ 2 / (2 * opt_param.medium.water.sound_speed * opt_param.medium.water.density) * 1e-4, 'DisplayName', 'Optimized (Simulated)');
    plot(dist_exit_plane, adjusted_profile_focus, 'DisplayName', 'Desired Profile');
    hold off;

    % Add focus and intensity reference lines
    xline(focus_wrt_mid_bowl, '--', 'DisplayName', 'Focal Point wrt mid-bowl');
    yline(desired_intensity, '--', 'DisplayName', 'Desired Intensity');

    xlabel('Distance wrt Mid-Bowl of Transducer [mm]');
    ylabel('Intensity [W/cm^2]');
    legend('Location', 'best');
    title(sprintf('Desired vs Optimized Profiles - Optimization Error: %.4f', min_err));
    
    % Save the profile comparison figure
    if opt_param.calibration.save_in_calibration_folder
        fig_path = fullfile(opt_param.calibration.path_output, strcat('Opt_simulation_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    else
        fig_path = fullfile(opt_param.output_location, strcat('Opt_simulation_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    end

    saveas(gcf, fig_path);

    max_pressure_index = find(pred_axial_pressure_opt == max(pred_axial_pressure_opt), 1);
    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n', sim_res_axial_position(max_pressure_index));
    
end