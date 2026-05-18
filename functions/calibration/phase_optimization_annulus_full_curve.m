function [error, ax1, ax2, h] = phase_optimization_annulus_full_curve(phase, parameters, velocity, axial_position, desired_intensity_curve, plot_results, opt_limits, weights)

% PHASE_OPTIMIZATION_ANNULUS_FULL_CURVE  Compute intensity profile error for full-curve phase optimisation
%
% Evaluates the O'Neil axial intensity model, restricts comparison to the
% optimisation range, applies optional Gaussian weighting around the FLHM
% centre, and returns the weighted mean squared error. Optionally plots the
% profiles and cost function.
%
% Use as:
%   [error, ax1, ax2, h] = phase_optimization_annulus_full_curve(phase, parameters, velocity, ...
%       axial_position, desired_intensity_curve, plot_results, opt_limits, weights)
%
% Input:
%   phase                   - phase shifts per transducer element [rad]
%   parameters              - PRESTUS config with transducer.annular geometry
%                             and medium_properties.water
%   velocity                - particle velocities per element [m/s]
%   axial_position          - axial positions (0 = transducer surface) [mm]
%   desired_intensity_curve - desired intensity profile along the axial direction [W/cm²]
%   plot_results            - enable visualisation (default: 0)
%   opt_limits              - [1x2] distance limits for error computation [mm]
%                             (default: [1, max(axial_position)])
%   weights                 - weights for the cost function; 0 = Gaussian around FLHM centre
%                             (default: 0)
%
% Output:
%   error - weighted mean squared error between computed and desired profiles
%   ax1   - handle to the intensity and weights subplot
%   ax2   - handle to the error visualisation subplot
%   h     - handle to the figure
%
% See also: PERFORM_GLOBAL_SEARCH, PHASE_OPTIMIZATION_ANNULUS

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

    %% Compute acoustic intensity profile
    use_rayleigh = isfield(parameters, 'calibration') && ...
                   isfield(parameters.calibration, 'forward_model') && ...
                   strcmp(parameters.calibration.forward_model, 'rayleigh');

    if use_rayleigh
        i_axial_oneil = rayleigh_axial_intensity(phase, velocity, parameters, axial_position);
    else
        p_axial_oneil = focusedAnnulusONeil(parameters.transducer.annular.curv_radius_mm / 1e3, ...
            [parameters.transducer.annular.elem_id_mm; parameters.transducer.annular.elem_od_mm] / 1e3, ...
            repmat(velocity, 1, parameters.transducer.annular.elem_n), ...
            [phase], parameters.transducer.freq_hz, ...
            parameters.medium_properties.water.sound_speed, ...
            parameters.medium_properties.water.density, ...
            (axial_position - 0.5) * 1e-3);
        i_axial_oneil = p_axial_oneil.^2 / (2 * parameters.medium_properties.water.sound_speed * parameters.medium_properties.water.density) * 1e-4;
    end

    %% Generate weights if not provided

    if weights == 0
        % UNIFORM: equal weight everywhere (optimize entire profile)
        weights = ones(size(axial_position));
    elseif weights >= 1
        % FWHM GAUSSIAN: weights controls narrowness
        [~, flhm_center_index] = get_flhm_center_position(axial_position, desired_intensity_curve);
        center_pos = axial_position(flhm_center_index);
        sigma = center_pos / weights;  % 1=wide FWHM, 10=narrow peak
        weights = normpdf(axial_position, center_pos, sigma);
    end
    % Weights always sum to one.
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
        plot(axial_position, i_axial_oneil, 'Color', [0 0 0]);
        plot(axial_position, desired_intensity_curve, 'Color', [1 0 0]);
        yyaxis right;
        plot(axial_position, weights);
        hold off;
        legend(["fitted profile", "real profile", "cost function", "error"], 'Location', 'South');
        legend('boxoff')

        % Plot error values within optimization limits
        ax2 = subplot(1, 2, 2);
        plot(axial_position(limit_ind), error_v,'-o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0]);
        legend(["error"]);
    end

end