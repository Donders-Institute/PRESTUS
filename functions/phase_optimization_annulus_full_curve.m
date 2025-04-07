function [error, ax1, ax2, h] = phase_optimization_annulus_full_curve(phase, parameters, velocity, axial_position, desired_intensity_curve, plot_results, opt_limits, weights)

% PHASE_OPTIMIZATION_ANNULUS_FULL_CURVE Optimizes the phase profile for a focused annular transducer.
%
% This function computes the acoustic intensity profile along the axial direction 
% for a focused annular transducer using O'Neil's model. It compares the computed 
% profile to a desired intensity curve and calculates an error metric. Optionally, 
% it visualizes the results and the cost function used for optimization.
%
% Input:
%   phase                  - Array specifying the phase shifts for each transducer element.
%   parameters             - Struct containing transducer and medium properties.
%   velocity               - Array specifying particle velocities for each transducer element.
%   axial_position         - Array specifying axial positions (0 = transducer surface).
%   desired_intensity_curve- Array specifying the desired intensity profile along the axial direction.
%   plot_results           - Boolean flag to enable/disable visualization of results (default: 0).
%   opt_limits             - [1x2] array specifying limits for optimization (default: [1, max(axial_position)]).
%   weights                - Array specifying weights for the cost function (default: Gaussian weights around FLHM center).
%
% Output:
%   error                  - Scalar value representing the mean squared error between computed and desired profiles.
%   ax1                    - Handle to the first subplot (intensity profiles and weights).
%   ax2                    - Handle to the second subplot (error visualization).
%   h                      - Handle to the figure containing plots.

    arguments
        phase
        parameters
        velocity
        axial_position % 0 is the transducer surface
        desired_intensity_curve
        plot_results = 0
        opt_limits = [1, max(axial_position)]
        weights = 0
    end

    %% Restrict axial positions to optimization limits
    limit_ind = (axial_position >= opt_limits(1) & axial_position <= opt_limits(2));

    %% Compute acoustic intensity profile using O'Neil's model
    % Focused annulus model based on transducer geometry and medium properties
    p_axial_oneil = focusedAnnulusONeil(parameters.transducer.curv_radius_mm / 1e3, ...
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm] / 1e3, ...
        repmat(velocity, 1, parameters.transducer.n_elements), ...
        [0 phase], parameters.transducer.source_freq_hz, ...
        parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, ...
        (axial_position - 0.5) * 1e-3);

    % Convert pressure profile to intensity profile
    i_axial_oneil = p_axial_oneil.^2 / (2 * parameters.medium.water.sound_speed * parameters.medium.water.density) * 1e-4;

    %% Generate weights if not provided
    if weights == 0
        % Find FLHM center and generate Gaussian weights around it
        [flhm_center, flhm_center_index] = get_flhm_center_position(axial_position, desired_intensity_curve);
        weights = normpdf(axial_position, axial_position(flhm_center_index) + 0.5, axial_position(flhm_center_index) / 3);
    end

    % Normalize weights to sum to 1
    weights = weights / sum(weights);

    %% Calculate error metric
    % Compute weighted squared error between computed and desired profiles
    error_v = (i_axial_oneil - desired_intensity_curve).^2 .* weights;
    error_v = error_v(limit_ind); % Restrict error calculation to optimization limits
    error = mean(error_v); % Mean squared error

    %% Visualization of results (optional)
    if plot_results
        h = figure;

        % Plot computed and desired profiles along with weights
        ax1 = subplot(1, 2, 1);
        hold on;
        plot(axial_position, i_axial_oneil);
        plot(axial_position, desired_intensity_curve);
        yyaxis right;
        plot(axial_position, weights);
        hold off;
        legend(["fitted profile", "real profile", "cost function", "error"]);

        % Plot error values within optimization limits
        ax2 = subplot(1, 2, 2);
        plot(axial_position(limit_ind), error_v,'-o');
        legend(["error"]);
    end

end