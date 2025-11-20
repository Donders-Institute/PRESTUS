function [p_axial_oneil, simulated_grid_adj_factor, velocity, axial_position] = ...
    compute_oneil_solution(parameters, pred_axial_pressure, dist_exit_plane, ...
    adjusted_profile_focus, focus_wrt_exit_plane, desired_intensity, equipment_name)
    % Compute O'Neil solution and plot it along with comparisons
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    %   parameters.calibration.path_output: Directory for saving results and figures.
    % - pred_axial_pressure: Predicted pressure along the beam axis [Pa].
    % - dist_exit_plane: Axial distance from the transducer exit plane [mm].
    % - adjusted_profile_focus: Adjusted intensity profile for the focus.
    % - focus_wrt_exit_plane: Focal distance relative to the exit plane [mm].
    % - desired_intensity: Desired intensity at the focal point [W/cm^2].
    % - equipment_name: Name of the equipment for labeling results.
    %
    % Returns:
    % - p_axial_oneil: Computed O'Neil solution for pressure along the beam axis [Pa].
    % - simulated_grid_adj_factor: Adjustment factor to align simulated pressure with analytical solution.
    % - velocity: Particle velocity [m/s].
    % - axial_position: Axial position vector [mm].

    % Compute particle velocity [m/s]
    velocity = parameters.transducer.source_amp(1) / ...
               (parameters.medium.water.density * parameters.medium.water.sound_speed);
    
    % TODO: should the resolution be retrieved from somewhere else?
    % Define the axial position vector [mm]
    % Assumes grid spacing of 0.5 mm (default resolution).
    axial_position = (1:parameters.default_grid_dims(end)) * 0.5;
    
    % Compute O'Neil analytical solution for pressure along the beam axis [Pa]
    p_axial_oneil = focusedAnnulusONeil(...
        parameters.transducer.curv_radius_mm / 1e3, ... % Convert radius to meters
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm] / 1e3, ... % Element dimensions in meters
        repmat(velocity, 1, parameters.transducer.n_elements), ... % Velocity array
        parameters.transducer.source_phase_rad, ... % Source phases [radians]
        parameters.transducer.source_freq_hz, ... % Source frequency [Hz]
        parameters.medium.water.sound_speed, ... % Sound speed in water [m/s]
        parameters.medium.water.density, ... % Water density [kg/m^3]
        (axial_position - 0.5) * 1e-3); % Axial positions (adjusted, in meters)
    
    % Convert pressures to intensities [W/cm^2]
    i_axial_oneil = p_axial_oneil .^ 2 / (2 * parameters.medium.water.sound_speed * parameters.medium.water.density) * 1e-4;
    pred_axial_intensity = pred_axial_pressure .^ 2 / (2 * parameters.medium.water.sound_speed * parameters.medium.water.density) * 1e-4;

    % Plot intensity along the beam axis
    figure('Position', [10, 10, 900, 500]);
    plot(axial_position, i_axial_oneil, ...
        'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'O''Neil Analytical Solution');
    hold on;
    plot(axial_position - (parameters.transducer.pos_grid(end) - 1) * 0.5, pred_axial_intensity, ...
        '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Inital Simulated Intensity');
    plot(dist_exit_plane, adjusted_profile_focus, ...
        'LineWidth', 2, 'Color', [1 0 0], 'DisplayName', 'Desired Profile');
    if isfield(parameters, 'expected_focal_distance_EP_mm')
        xline(parameters.expected_focal_distance_EP_mm, '--', ...
            'LineWidth', 1.2, 'DisplayName', 'Expected Focal Distance (mm from EP)');
    end
    hold off;

    xlabel('Distance w.r.t. Exit Plane [mm]');
    ylabel('Intensity [W/cm^2]');
    legend('show');
    title('Intensity Along the Beam Axis');
    grid on;
    ylim([0 inf]);
    xlim([-5 inf]);

    % Save the figure
    fig_path = fullfile(parameters.outputs_folder, sprintf('Initial_Simulation_F_%.2f_I_%.2f_%s.png', ...
        focus_wrt_exit_plane, desired_intensity, equipment_name));
    saveas(gcf, fig_path);
    close(gcf);

    % Report the estimated distance to the point of maximum pressure
    [~, max_idx] = max(p_axial_oneil);
    fprintf('Estimated distance to maximum pressure: %.2f mm\n', axial_position(max_idx));

    % Compute adjustment factor to align simulated and analytical pressures
    simulated_grid_adj_factor = max(pred_axial_pressure(:)) / max(p_axial_oneil(:));
    
end