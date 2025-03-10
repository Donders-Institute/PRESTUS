function p_axial_oneil_opt = recalculate_analytical_sol(parameters, p_axial_oneil, opt_phases, opt_velocity, dist_from_tran, adjusted_profile_focus, axial_position, focus_wrt_exit_plane, desired_intensity, prestus_dir, equipment_name, save_in_general_folder)
    % Recalculate analytical solution based on optimized phases and velocity.
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - p_axial_oneil: Initial O'Neil solution for pressure [Pa].
    % - opt_phases: Optimized phases for each transducer element [rad].
    % - opt_velocity: Optimized particle velocity [m/s].
    % - dist_from_tran: Distance vector for desired profile focus [mm].
    % - adjusted_profile_focus: Adjusted desired intensity profile [W/cm^2].
    % - axial_position: Axial position vector [mm].
    % - focus_wrt_exit_plane: Focal depth with respect to the transducer exit plane [mm].
    % - desired_intensity: Target intensity for optimization [W/cm^2].
    % - prestus_dir: Directory for saving output.
    % - equipment_name: Name of the equipment used.
    % - save_in_general_folder: Option to save data in the general PRESTUS output folder.
    %
    % Returns:
    % - p_axial_oneil_opt: Optimized O'Neil solution for pressure along the beam axis [Pa].

    % Compute optimized analytical pressure profile
    p_axial_oneil_opt = focusedAnnulusONeil(parameters.transducer.curv_radius_mm / 1e3, ...
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm] / 1e3, ...
        repmat(opt_velocity, 1, parameters.transducer.n_elements), ...
        [0 opt_phases], parameters.transducer.source_freq_hz, parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, (axial_position - 0.5) * 1e-3);
    
    % Convert pressure to intensity
    i_axial_oneil = p_axial_oneil .^ 2 / (2 * parameters.medium.water.sound_speed * parameters.medium.water.density) * 1e-4;
    i_axial_oneil_opt = p_axial_oneil_opt .^ 2 / (2 * parameters.medium.water.sound_speed * parameters.medium.water.density) * 1e-4;
    
    % Compute focal position relative to the mid-bowl of the transducer
    focus_wrt_mid_bowl = focus_wrt_exit_plane + (parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);

    % Plot comparison of profiles
    figure('Position', [10, 10, 900, 500]);
    plot(axial_position, i_axial_oneil, 'DisplayName', 'Original Analytical Profile');
    hold on;
    plot(axial_position, i_axial_oneil_opt, 'DisplayName', 'Optimized Analytical Profile');
    plot(dist_from_tran, adjusted_profile_focus, 'DisplayName', 'Desired Profile');
    hold off;
    xline(focus_wrt_mid_bowl, '--', 'DisplayName', 'Focal Point wrt Mid-Bowl');
    yline(desired_intensity, '--', 'DisplayName', 'Desired Intensity');

    xlabel('Distance wrt Mid-Bowl of Transducer [mm]');
    ylabel('Intensity [W/cm^2]');
    legend('Location', 'best');
    title('Pressure Along the Beam Axis');
    
    % Save the profile comparison figure
    if save_in_general_folder
        fig_path = fullfile(prestus_dir, strcat('Recalculated_oneil_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    else
        fig_path = fullfile(parameters.output_location, strcat('Recalculated_oneil_at_F_', num2str(focus_wrt_exit_plane), '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    end

    saveas(gcf, fig_path);

    fprintf('Estimated distance to the point of maximum pressure: %.2f mm\n', axial_position(p_axial_oneil_opt == max(p_axial_oneil_opt)))
    fprintf('Estimated distance to the center of half-maximum range: %.2f mm\n', get_flhm_center_position(axial_position, p_axial_oneil_opt))
    
end