function plot_opt_sim_results(parameters, profile_target, profile_oneil, profile_oneil_opt, profile_sim, profile_sim_opt, min_err)
% PLOT_OPT_SIM_RESULTS  Plot and save comparison of calibration optimisation profiles
%
% Generates a figure comparing the target, original analytical, optimised
% analytical, original simulated, and optimised simulated axial intensity
% profiles. Saves the figure to calibration.path_output.
%
% Use as:
%   plot_opt_sim_results(parameters, profile_target, profile_oneil, profile_oneil_opt, ...
%                        profile_sim, profile_sim_opt, min_err)
%
% Input:
%   parameters       - PRESTUS config; uses calibration.path_output,
%                      calibration.desired_focal_distance_ep [mm],
%                      calibration.desired_intensity [W/cm²], calibration.equipment_name
%   profile_target   - struct with axial_distance_bowl [mm] and axial_intensity [W/cm²]
%   profile_oneil    - struct with axial_intensity from original O'Neil solution [W/cm²]
%   profile_oneil_opt- struct with axial_intensity from optimised O'Neil solution [W/cm²]
%   profile_sim      - struct with axial_intensity from original simulation [W/cm²]
%   profile_sim_opt  - struct with axial_intensity from optimised simulation [W/cm²]
%   min_err          - minimum optimisation error
%
% See also: CALIBRATION_TRANSDUCER, PERFORM_GLOBAL_SEARCH, RECOMPUTE_ONEIL_SOLUTION

arguments
    parameters        (1,1) struct
    profile_target    (1,1) struct
    profile_oneil     (1,1) struct
    profile_oneil_opt (1,1) struct
    profile_sim       (1,1) struct
    profile_sim_opt   (1,1) struct
    min_err           (1,1) {mustBeNumeric}
end
    
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