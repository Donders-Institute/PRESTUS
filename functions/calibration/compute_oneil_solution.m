function [profile_oneil, simulated_oneil_scaling] = compute_oneil_solution(parameters, profile_sim, profile_target)

    % Compute O'Neil solution and plot it along with comparisons
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    %       .calibration.path_output: Directory for saving results and figures.
    %       .calibration.desired_focal_distance_ep: Focal distance relative to the exit plane [mm]. Will be used only for labelling.
    %       .calibration.desired_intensity: Desired intensity at the focal point [W/cm^2].
    %       .calibration.equipment_name: Name of the equipment for labeling results.
    % - profile_sim.axial_intensity: Simulated axial intensity [W/cm2].
    % - profile_sim.axial_distance_bowl: Axial distance of simulated pressure [mm from bowl]
    % - profile_sim.velocity: Particle profile_sim.velocity [m/s].
    % - profile_target.axial_intensity: Adjusted intensity profile for the focus.
    %
    % Returns:
    % - profile_oneil.axial_intensity: Computed O'Neil solution for pressure along the beam axis [Pa].
    % - profile_oneil.axial_distance_bowl: Axial position vector [mm from bowl].
    % - simulated_oneil_scaling: Adjustment factor to align simulated intensity with analytical solution.
    
    % Define the axial position vector [mm]
    axial_position = profile_sim.axial_distance_bowl;

    % Flip the desired profile
    profile_target.axial_intensity = profile_target.axial_intensity';
    
    % Compute O'Neil analytical solution for pressure along the beam axis [Pa]
    p_axial_oneil = focusedAnnulusONeil(...
        parameters.transducer.curv_radius_mm / 1e3, ... % Convert radius to meters
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm] / 1e3, ... % Element dimensions in meters
        repmat(profile_sim.velocity, 1, parameters.transducer.n_elements), ... % Velocity array
        parameters.transducer.source_phase_rad, ... % Source phases [radians]
        parameters.transducer.source_freq_hz, ... % Source frequency [Hz]
        parameters.medium.water.sound_speed, ... % Sound speed in water [m/s]
        parameters.medium.water.density, ... % Water density [kg/m^3]
        (axial_position - 0.5) * 1e-3); % Axial positions (adjusted, in meters)
    
    % Convert pressure to intensities [W/cm^2]
    i_axial_oneil = p_axial_oneil .^ 2 / (2 * parameters.medium.water.sound_speed * parameters.medium.water.density) * 1e-4;

    % Plot intensity along the beam axis
    figure('Position', [10, 10, 900, 500]);
    plot(axial_position, i_axial_oneil, ...
        'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'O''Neil Analytical Solution');
    hold on;
    plot(axial_position, profile_sim.axial_intensity, ...
        '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Inital Simulated Intensity');
    plot(profile_target.axial_distance_bowl, profile_target.axial_intensity, ...
        'LineWidth', 2, 'Color', [1 0 0], 'DisplayName', 'Desired Profile');
    if isfield(parameters, 'expected_focal_distance_bowl')
        xline(parameters.expected_focal_distance_bowl, '--', ...
            'LineWidth', 1.2, 'DisplayName', 'Expected Focal Distance (mm from bowl)', 'Color', [1 0 0]);
    end
    if isfield(parameters.transducer, 'focal_distance_offset')
        xline(parameters.transducer.focal_distance_offset, '--', ...
            'LineWidth', 1.2, 'DisplayName', 'Exit Plane');
    end
    hold off;

    xlabel('Distance w.r.t. Transducer Bowl [mm]');
    ylabel('Intensity [W/cm^2]');
    legend('show');
    title('Intensity Along the Beam Axis');
    grid on;
    ylim([0 inf]);
    xlim([0 inf]);

    % Save the figure
    fig_path = fullfile(parameters.outputs_folder, sprintf('Initial_Simulation_F_%.2f_I_%.2f_%s.png', ...
        parameters.calibration.desired_focal_distance_ep, parameters.calibration.desired_intensity, parameters.calibration.equipment_name));
    saveas(gcf, fig_path);
    close(gcf);

    % Report the estimated distance to the point of maximum intensity
    [~, max_idx] = max(i_axial_oneil);
    fprintf('Estimated distance to maximum intensity: %.2f mm\n', axial_position(max_idx));

    % Compute adjustment factor to align simulated and analytical intensities
    simulated_oneil_scaling = max(profile_sim.axial_intensity(:)) / max(i_axial_oneil(:));

    %% Collect outputs
    profile_oneil.axial_intensity = i_axial_oneil;
    profile_oneil.axial_distance_bowl = axial_position;
    
end