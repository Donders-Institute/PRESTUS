function profile_oneil_opt = recompute_oneil_solution(parameters, profile_oneil, profile_target, opt_phases, opt_velocity)
% RECOMPUTE_ONEIL_SOLUTION  Recompute and plot O'Neil solution with optimised phases and velocity
%
% Evaluates the focusedAnnulusONeil model with the optimised parameters,
% plots the result alongside the original analytical and target profiles,
% and saves the figure.
%
% Use as:
%   profile_oneil_opt = recompute_oneil_solution(parameters, profile_oneil, profile_target, opt_phases, opt_velocity)
%
% Input:
%   parameters     - PRESTUS config; uses transducer.annular geometry,
%                    medium_properties.water, calibration.desired_focal_distance_ep [mm],
%                    calibration.desired_intensity [W/cm²], calibration.equipment_name
%   profile_oneil  - struct with axial_distance_bowl [mm] and axial_intensity [W/cm²]
%   profile_target - struct with axial_distance_bowl [mm] and axial_intensity [W/cm²]
%   opt_phases     - optimised phases per element [rad]
%   opt_velocity   - optimised particle velocity [m/s]
%
% Output:
%   profile_oneil_opt - struct with axial_intensity [W/cm²] and axial_distance_bowl [mm]
%
% See also: COMPUTE_ONEIL_SOLUTION, PLOT_OPT_SIM_RESULTS, CALIBRATION_TRANSDUCER

arguments
    parameters     (1,1) struct
    profile_oneil  (1,1) struct
    profile_target (1,1) struct
    opt_phases     (1,:) {mustBeNumeric}
    opt_velocity   (1,1) {mustBeNumeric}
end

    % Get axial position
    axial_position = profile_oneil.axial_distance_bowl;

    % Compute optimized analytical pressure profile
    p_axial_oneil_opt = focusedAnnulusONeil(...
        parameters.transducer.annular.curv_radius_mm / 1e3, ...
        [parameters.transducer.annular.elem_id_mm; parameters.transducer.annular.elem_od_mm] / 1e3, ...
        repmat(opt_velocity, 1, parameters.transducer.annular.elem_n), ...
        [opt_phases], ...
        parameters.transducer.freq_hz, ...
        parameters.medium_properties.water.sound_speed, ...
        parameters.medium_properties.water.density, ...
        (axial_position - 0.5) * 1e-3);
    
    % Convert pressure to intensity
    i_axial_oneil = profile_oneil.axial_intensity;
    i_axial_oneil_opt = p_axial_oneil_opt .^ 2 / (2 * parameters.medium_properties.water.sound_speed * parameters.medium_properties.water.density) * 1e-4;

    % Plot comparison of profiles
    figure('Position', [10, 10, 900, 500]);
    plot(axial_position, i_axial_oneil, ...
        'LineWidth', 1, 'Color', [0.5 0.5 0.5], 'DisplayName', 'Original Profile (Analytical)');
    hold on;
    plot(axial_position, i_axial_oneil_opt, ...
        'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', 'Optimized Profile (Analytical)');
    plot(profile_target.axial_distance_bowl, profile_target.axial_intensity, ...
        'LineWidth', 2, 'Color', [1 0 0], 'DisplayName', 'Desired Profile');
    hold off;
    xline(parameters.transducer(1).focal_distance_bowl, '--', 'DisplayName', 'Focal Point wrt Mid-Bowl');
    yline(parameters.calibration.desired_intensity, '--', 'DisplayName', 'Desired Intensity');

    xlabel('Distance wrt Mid-Bowl of Transducer [mm]');
    ylabel('Intensity [W/cm^2]');
    legend('Location', 'NorthEast');
    legend('boxoff')
    title('Pressure Along the Beam Axis');
    ylim([0 inf]);
    xlim([-5 inf]);
    
    % Save the profile comparison figure
    fig_path = fullfile(parameters.io.outputs_folder, ...
        strcat('Recalculated_oneil_at_F_', num2str(parameters.calibration.desired_focal_distance_ep), ...
        '_at_I_', num2str(parameters.calibration.desired_intensity), ...
        '_', parameters.calibration.equipment_name, '.png'));
    saveas(gcf, fig_path);
    close(gcf);

    % exclude the near-field if requested
    i_axial_oneil_opt_summary = i_axial_oneil_opt;
    axial_position_adj = axial_position;
    if parameters.calibration.skip_front_peak_mm ~=0
        i_remove = axial_position_adj <= parameters.calibration.skip_front_peak_mm;
        axial_position_adj(i_remove) = [];
        i_axial_oneil_opt_summary(i_remove) = [];
    end

    fprintf('Estimated distance to the point of maximum intensity: %.2f mm\n', ...
        axial_position_adj(i_axial_oneil_opt_summary == max(i_axial_oneil_opt_summary)))
    try
        fprintf('Estimated distance to the center of half-maximum range: %.2f mm\n', ...
            get_flhm_center_position(axial_position_adj, i_axial_oneil_opt_summary))
    catch
        warning('Could not estimated distance to the center of half-maximum range');
    end

    %% Collect in structure

    profile_oneil_opt.axial_intensity = i_axial_oneil_opt;
    profile_oneil_opt.axial_distance_bowl = axial_position;

    
end