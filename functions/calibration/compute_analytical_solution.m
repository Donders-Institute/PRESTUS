function [profile_analytical, simulated_analytical_scaling] = compute_analytical_solution(parameters, profile_sim, profile_target)
% COMPUTE_ANALYTICAL_SOLUTION  Compute the analytical forward-model solution and compare with simulation
%
% Evaluates the configured forward model (O'Neil or Rayleigh), converts to
% intensity, and plots the result alongside the simulated and target profiles.
% Returns the analytical profile and a scaling factor for aligning simulation
% to the analytical solution.
%
% The forward model is selected by parameters.calibration.forward_model:
%   'oneil'    (default) — O'Neil closed-form via focusedAnnulusONeil
%   'rayleigh'           — Rayleigh–Sommerfeld via rayleigh_axial_intensity
%
% Use as:
%   [profile_analytical, simulated_analytical_scaling] = ...
%       compute_analytical_solution(parameters, profile_sim, profile_target)
%
% Input:
%   parameters     - PRESTUS config; uses transducer.annular geometry,
%                    medium_properties.water, calibration.desired_focal_distance_ep [mm],
%                    calibration.desired_intensity [W/cm²], calibration.equipment_name,
%                    calibration.forward_model ('oneil'|'rayleigh', optional)
%   profile_sim    - struct with axial_intensity [W/cm²], axial_distance_bowl [mm], velocity [m/s]
%   profile_target - struct with axial_intensity and axial_distance_bowl [mm]
%
% Output:
%   profile_analytical          - struct with axial_intensity [W/cm²] and axial_distance_bowl [mm]
%   simulated_analytical_scaling - scaling factor aligning simulated to analytical peak intensity
%
% See also: CALIBRATION_TRANSDUCER, FIT_VELOCITY_TO_INTENSITY, RECOMPUTE_ANALYTICAL_SOLUTION,
%           EXTRACT_ANALYTICAL_PROFILE, RAYLEIGH_AXIAL_INTENSITY

arguments
    parameters   (1,1) struct
    profile_sim  (1,1) struct
    profile_target (1,1) struct
end

    if isfield(parameters.calibration, 'forward_model')
        forward_model = parameters.calibration.forward_model;
    else
        forward_model = 'oneil';
    end

    axial_position = profile_sim.axial_distance_bowl;

    % Flip the desired profile
    profile_target.axial_intensity = profile_target.axial_intensity';

    if strcmp(forward_model, 'rayleigh')
        i_axial = rayleigh_axial_intensity( ...
            parameters.transducer.annular.elem_phase_rad, profile_sim.velocity, ...
            parameters, axial_position);
    elseif strcmp(forward_model, 'oneil')
        p_axial = focusedAnnulusONeil( ...
            parameters.transducer.annular.curv_radius_mm / 1e3, ...
            [parameters.transducer.annular.elem_id_mm; parameters.transducer.annular.elem_od_mm] / 1e3, ...
            repmat(profile_sim.velocity, 1, parameters.transducer.annular.elem_n), ...
            parameters.transducer.annular.elem_phase_rad, ...
            parameters.transducer.freq_hz, ...
            parameters.medium_properties.water.sound_speed, ...
            parameters.medium_properties.water.density, ...
            (axial_position - 0.5) * 1e-3);
        i_axial = p_axial .^ 2 / (2 * parameters.medium_properties.water.sound_speed * parameters.medium_properties.water.density) * 1e-4;
    else
        error('compute_analytical_solution: unknown forward_model ''%s''; expected ''oneil'' or ''rayleigh''.', forward_model);
    end

    % Plot intensity along the beam axis
    figure('Position', [10, 10, 900, 500]);
    plot(axial_position, i_axial, ...
        'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', sprintf('Analytical Solution (%s)', forward_model));
    hold on;
    plot(axial_position, profile_sim.axial_intensity, ...
        '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Initial Simulated Intensity');
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

    img_folder = fullfile(parameters.io.outputs_folder, 'img_calibration');
    if ~exist(img_folder, 'dir'); mkdir(img_folder); end
    fig_path = fullfile(img_folder, sprintf('Initial_Simulation_F_%.2f_I_%.2f_%s.png', ...
        parameters.calibration.desired_focal_distance_ep, parameters.calibration.desired_intensity, ...
        parameters.calibration.equipment_name));
    saveas(gcf, fig_path);
    close(gcf);

    [~, max_idx] = max(i_axial);
    fprintf('Estimated distance to maximum intensity: %.2f mm\n', axial_position(max_idx));

    simulated_analytical_scaling = max(profile_sim.axial_intensity(:)) / max(i_axial(:));

    profile_analytical.axial_intensity     = i_axial;
    profile_analytical.axial_distance_bowl = axial_position;

end
