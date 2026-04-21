function [profile_oneil, simulated_oneil_scaling] = compute_oneil_solution(parameters, profile_sim, profile_target)
% COMPUTE_ONEIL_SOLUTION  Compute the O'Neil analytical solution and compare with simulation
%
% Evaluates the focusedAnnulusONeil axial pressure model, converts to
% intensity, and plots the result alongside the simulated and target profiles.
% Returns the O'Neil profile and a scaling factor for aligning simulation to
% the analytical solution.
%
% Use as:
%   [profile_oneil, simulated_oneil_scaling] = compute_oneil_solution(parameters, profile_sim, profile_target)
%
% Input:
%   parameters   - PRESTUS config; uses transducer.annular geometry,
%                  medium_properties.water, calibration.desired_focal_distance_ep [mm],
%                  calibration.desired_intensity [W/cm²], calibration.equipment_name
%   profile_sim  - struct with axial_intensity [W/cm²], axial_distance_bowl [mm], velocity [m/s]
%   profile_target - struct with axial_intensity and axial_distance_bowl [mm]
%
% Output:
%   profile_oneil          - struct with axial_intensity [W/cm²] and axial_distance_bowl [mm]
%   simulated_oneil_scaling - scaling factor to align simulated intensity with analytical solution
%
% See also: CALIBRATION_TRANSDUCER, FIT_VELOCITY_TO_INTENSITY, RECOMPUTE_ONEIL_SOLUTION

arguments
    parameters   (1,1) struct
    profile_sim  (1,1) struct
    profile_target (1,1) struct
end
    
    % Define the axial position vector [mm]
    axial_position = profile_sim.axial_distance_bowl;

    % Flip the desired profile
    profile_target.axial_intensity = profile_target.axial_intensity';
    
    % Compute O'Neil analytical solution for pressure along the beam axis [Pa]
    p_axial_oneil = focusedAnnulusONeil(...
        parameters.transducer.annular.curv_radius_mm / 1e3, ... % Convert radius to meters
        [parameters.transducer.annular.elem_id_mm; parameters.transducer.annular.elem_od_mm] / 1e3, ... % Element dimensions in meters
        repmat(profile_sim.velocity, 1, parameters.transducer.annular.elem_n), ... % Velocity array
        parameters.transducer.annular.elem_phase_rad, ... % Source phases [radians]
        parameters.transducer.freq_hz, ... % Source frequency [Hz]
        parameters.medium_properties.water.sound_speed, ... % Sound speed in water [m/s]
        parameters.medium_properties.water.density, ... % Water density [kg/m^3]
        (axial_position - 0.5) * 1e-3); % Axial positions (adjusted, in meters)
    
    % Convert pressure to intensities [W/cm^2]
    i_axial_oneil = p_axial_oneil .^ 2 / (2 * parameters.medium_properties.water.sound_speed * parameters.medium_properties.water.density) * 1e-4;

    % Plot intensity along the beam axis
    figure('Position', [10, 10, 900, 500]);
    plot(axial_position, i_axial_oneil, ...
        'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'O''Neil Analytical Solution');
    hold on;
    plot(axial_position, profile_sim.axial_intensity, ...
        '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Inital Simulated Intensity');
    plot(profile_target.axial_distance_bowl, profile_target.axial_intensity, ...
        'LineWidth', 2, 'Color', [1 0 0], 'DisplayName', 'Desired Profile');
    if isfield(parameters.transducer(1), 'focal_distance_bowl')
        xline(parameters.transducer(1).focal_distance_bowl, '--', ...
            'LineWidth', 1.2, 'DisplayName', 'Expected Focal Distance (mm from bowl)', 'Color', [1 0 0]);
    end
    if isfield(parameters.transducer(1), 'focal_distance_offset')
        xline(parameters.transducer(1).focal_distance_offset, '--', ...
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
    fig_path = fullfile(parameters.io.outputs_folder, sprintf('Initial_Simulation_F_%.2f_I_%.2f_%s.png', ...
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