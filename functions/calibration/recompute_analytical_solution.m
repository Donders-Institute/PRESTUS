function profile_analytical_opt = recompute_analytical_solution(parameters, profile_analytical, profile_target, opt_phases, opt_velocity)
% RECOMPUTE_ANALYTICAL_SOLUTION  Recompute and plot analytical solution with optimised phases and velocity
%
% Evaluates the configured forward model with the optimised parameters,
% plots the result alongside the original analytical and target profiles,
% and saves the figure.
%
% The forward model is selected by parameters.calibration.forward_model:
%   'oneil'    (default) — O'Neil closed-form via focusedAnnulusONeil
%   'rayleigh'           — Rayleigh–Sommerfeld via rayleigh_axial_intensity
%
% Use as:
%   profile_analytical_opt = recompute_analytical_solution(parameters, profile_analytical, ...
%                                profile_target, opt_phases, opt_velocity)
%
% Input:
%   parameters          - PRESTUS config; uses transducer.annular geometry,
%                         medium_properties.water, calibration.desired_focal_distance_ep [mm],
%                         calibration.desired_intensity [W/cm²], calibration.equipment_name,
%                         calibration.forward_model ('oneil'|'rayleigh', optional)
%   profile_analytical  - struct with axial_distance_bowl [mm] and axial_intensity [W/cm²]
%   profile_target      - struct with axial_distance_bowl [mm] and axial_intensity [W/cm²]
%   opt_phases          - optimised phases per element [rad]
%   opt_velocity        - optimised particle velocity [m/s]
%
% Output:
%   profile_analytical_opt - struct with axial_intensity [W/cm²] and axial_distance_bowl [mm]
%
% See also: COMPUTE_ANALYTICAL_SOLUTION, PLOT_OPT_SIM_RESULTS, CALIBRATION_TRANSDUCER,
%           RAYLEIGH_AXIAL_INTENSITY

arguments
    parameters          (1,1) struct
    profile_analytical  (1,1) struct
    profile_target      (1,1) struct
    opt_phases          (1,:) {mustBeNumeric}
    opt_velocity        (1,1) {mustBeNumeric}
end

    if isfield(parameters.calibration, 'forward_model')
        forward_model = parameters.calibration.forward_model;
    else
        forward_model = 'oneil';
    end

    axial_position = profile_analytical.axial_distance_bowl;

    % Evaluate forward model with optimised parameters
    if strcmp(forward_model, 'rayleigh')
        % Temporarily set phases so rayleigh_axial_intensity picks them up via parameters
        params_opt = parameters;
        params_opt.transducer.annular.elem_phase_rad = opt_phases;
        i_axial_opt = rayleigh_axial_intensity(opt_phases, opt_velocity, params_opt, axial_position);
    elseif strcmp(forward_model, 'oneil')
        p_axial_opt = focusedAnnulusONeil( ...
            parameters.transducer.annular.curv_radius_mm / 1e3, ...
            [parameters.transducer.annular.elem_id_mm; parameters.transducer.annular.elem_od_mm] / 1e3, ...
            repmat(opt_velocity, 1, parameters.transducer.annular.elem_n), ...
            [opt_phases], ...
            parameters.transducer.freq_hz, ...
            parameters.medium_properties.water.sound_speed, ...
            parameters.medium_properties.water.density, ...
            (axial_position - 0.5) * 1e-3);
        i_axial_opt = p_axial_opt .^ 2 / (2 * parameters.medium_properties.water.sound_speed * parameters.medium_properties.water.density) * 1e-4;
    else
        error('recompute_analytical_solution: unknown forward_model ''%s''; expected ''oneil'' or ''rayleigh''.', forward_model);
    end

    % Plot comparison
    figure('Position', [10, 10, 900, 500]);
    plot(axial_position, profile_analytical.axial_intensity, ...
        'LineWidth', 1, 'Color', [0.5 0.5 0.5], 'DisplayName', sprintf('Original (%s)', forward_model));
    hold on;
    plot(axial_position, i_axial_opt, ...
        'LineWidth', 2, 'Color', [0 0 0], 'DisplayName', sprintf('Optimised (%s)', forward_model));
    plot(profile_target.axial_distance_bowl, profile_target.axial_intensity, ...
        'LineWidth', 2, 'Color', [1 0 0], 'DisplayName', 'Desired Profile');
    hold off;
    xline(parameters.transducer(1).focal_distance_bowl, '--', 'DisplayName', 'Focal Point wrt Mid-Bowl');
    yline(parameters.calibration.desired_intensity, '--', 'DisplayName', 'Desired Intensity');
    xlabel('Distance wrt Mid-Bowl of Transducer [mm]');
    ylabel('Intensity [W/cm^2]');
    legend('Location', 'NorthEast');
    legend('boxoff');
    title('Pressure Along the Beam Axis');
    ylim([0 inf]);
    xlim([-5 inf]);

    img_folder = fullfile(parameters.io.outputs_folder, 'img_calibration');
    if ~exist(img_folder, 'dir'); mkdir(img_folder); end
    fig_path = fullfile(img_folder, ...
        strcat('Recalculated_analytical_at_F_', num2str(parameters.calibration.desired_focal_distance_ep), ...
        '_at_I_', num2str(parameters.calibration.desired_intensity), ...
        '_', parameters.calibration.equipment_name, '.png'));
    saveas(gcf, fig_path);
    close(gcf);

    % Exclude near-field if requested
    i_axial_opt_summary = i_axial_opt;
    axial_position_adj  = axial_position;
    if parameters.calibration.skip_front_peak_mm ~= 0
        i_remove = axial_position_adj <= parameters.calibration.skip_front_peak_mm;
        axial_position_adj(i_remove)      = [];
        i_axial_opt_summary(i_remove) = [];
    end

    fprintf('Estimated distance to the point of maximum intensity: %.2f mm\n', ...
        axial_position_adj(i_axial_opt_summary == max(i_axial_opt_summary)));
    try
        fprintf('Estimated distance to the center of half-maximum range: %.2f mm\n', ...
            get_flhm_center_position(axial_position_adj, i_axial_opt_summary));
    catch
        warn('Could not estimate distance to the center of half-maximum range');
    end

    profile_analytical_opt.axial_intensity     = i_axial_opt;
    profile_analytical_opt.axial_distance_bowl = axial_position;

end
