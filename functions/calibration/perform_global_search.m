function [opt_phases, opt_velocity, min_err, opt_params_raw] = perform_global_search(parameters, profile_target, velocity)
% PERFORM_GLOBAL_SEARCH  Global optimisation of transducer phases and particle velocity
%
% Runs a global search (surrogate or genetic algorithm) to minimise the
% difference between the O'Neil analytical profile and the target empirical
% profile across the optimisation distance range.
%
% Use as:
%   [opt_phases, opt_velocity, min_err] = perform_global_search(parameters, profile_target, velocity)
%   [opt_phases, opt_velocity, min_err, opt_params_raw] = perform_global_search(...)
%
% Input:
%   parameters     - PRESTUS config with transducer.annular geometry and
%                    calibration.opt_limits [mm], calibration.opt_weights,
%                    calibration.opt_seed, calibration.opt_upper_velocity [m/s]
%   profile_target - struct with axial_distance_bowl [mm] and axial_intensity [W/cm²]
%   velocity       - initial particle velocity estimate [m/s]
%
% Output:
%   opt_phases      - optimised phases per element [rad]
%   opt_velocity    - optimised particle velocity [m/s]
%   min_err         - minimum error achieved during optimisation
%   opt_params_raw  - raw optimisation parameter vector before phase expansion:
%                     'linear':   [phase_start_rad, phase_step_rad_per_elem, velocity_m_s]
%                     'monotonic':[phase_start_rad, delta_1..delta_{N-1}_rad, velocity_m_s]
%                     otherwise:  [phase_1..phase_N_rad, velocity_m_s]
%
% See also: FIT_VELOCITY_TO_INTENSITY, PHASE_OPTIMIZATION_ANNULUS_FULL_CURVE,
%           CALIBRATION_TRANSDUCER

arguments
    parameters     (1,1) struct
    profile_target (1,1) struct
    velocity       (1,1) {mustBeNumeric}
end
    
    if ~isfield(parameters.calibration, 'opt_limits') || isempty(parameters.calibration.opt_limits)
        % By default set to min and max distance with non-NAN intensity
        min_available = min(profile_target.axial_distance_bowl(~isnan(profile_target.axial_intensity)));
        max_available = max(profile_target.axial_distance_bowl(~isnan(profile_target.axial_intensity)));
        opt_limits = [min_available, max_available];

    else
        opt_limits = parameters.calibration.opt_limits;
    end

    if ~isfield(parameters.calibration, 'opt_weights') 
        weights = 0; % uniform weighting
    else
        weights = parameters.calibration.opt_weights;
    end

    % Phase precession mode (calibration.opt_phase_precession):
    %   false      - independent per-element phases (default)
    %   'linear'   - strict linear ramp: phase[i] = phase_start + i * phase_step
    %                (2 free phase parameters)
    %   'monotonic'- ordered phases via cumulative increments: phase[i] = phase[i-1] + delta[i]
    %                with delta[i] >= 0 (N free phase parameters, but monotonically constrained)
    if isfield(parameters.calibration, 'opt_phase_precession')
        precession_mode = parameters.calibration.opt_phase_precession;
    else
        precession_mode = false;
    end

    n_elem = parameters.transducer.annular.elem_n;

    % Helper: reconstruct full phase vector from optimisation parameters
    switch precession_mode
        case 'linear'
            % [phase_start, phase_step, velocity]
            phases_from_params = @(p) mod(p(1) + (0:n_elem-1) * p(2), 2*pi);
        case 'monotonic'
            % [phase_start, delta_1, ..., delta_{N-1}, velocity]
            % delta_i >= 0 ensures monotonic ordering; mod wraps to [0, 2pi]
            phases_from_params = @(p) mod(p(1) + [0, cumsum(p(2:n_elem))], 2*pi);
        otherwise
            phases_from_params = @(p) p(1:end-1);
    end

    optimize_phases = @(p) phase_optimization_annulus_full_curve(...
        phases_from_params(p), ...
        parameters, ...
        p(end), ...
        profile_target.axial_distance_bowl, ...
        profile_target.axial_intensity, ...
        0, ...
        opt_limits, ...
        weights);

    % Set a random seed for reproducibility.
    if isfield(parameters.calibration, 'opt_seed')
        rng(parameters.calibration.opt_seed, 'twister');
    end

    % Define initial guess, bounds, and options for the optimization problem.
    if ~isfield(parameters.calibration, 'opt_upper_velocity') || isempty(parameters.calibration.opt_upper_velocity)
        parameters.calibration.opt_upper_velocity = 0.2; % set default for upper velocity to 20 mm/s;
    end

    switch precession_mode
        case 'linear'
            % [phase_start (rad), phase_step (rad/elem), velocity (m/s)]
            initial_guess = [rand() * 2*pi, (rand()-0.5) * 2*pi, velocity];
            lower_bounds  = [0,    -2*pi, 0.001];
            upper_bounds  = [2*pi,  2*pi, parameters.calibration.opt_upper_velocity];
        case 'monotonic'
            % [phase_start, delta_1, ..., delta_{N-1}, velocity]
            initial_guess = [rand() * 2*pi, rand(1, n_elem-1) * 2*pi/(n_elem-1), velocity];
            lower_bounds  = [0,    zeros(1, n_elem-1), 0.001];
            upper_bounds  = [2*pi, 2*pi * ones(1, n_elem-1), parameters.calibration.opt_upper_velocity];
        otherwise
            initial_guess = [randi(360, [1, n_elem]) / 180 * pi, velocity];
            lower_bounds  = [zeros(1, n_elem), 0.001];
            upper_bounds  = [2*pi * ones(1, n_elem), parameters.calibration.opt_upper_velocity];
    end

    if ~isfield(parameters.calibration, 'opt_method') || strcmp(parameters.calibration.opt_method, 'FEXminimize')
        % by default use FEXminimize
        func = optimize_phases;
        options = setoptimoptions('popsize', 5000, 'FinDiffType', 'central', 'TolCon', 1e-8);

        % Perform the global search to find optimal phases and velocity.
        [opt_phases_and_velocity, min_err] = minimize(...
            func, ...
            initial_guess, ...
            [],...
            [],...
            [],...
            [],...
            lower_bounds, ...
            upper_bounds, ...
            [], ...
            options);
    elseif strcmp(parameters.calibration.opt_method, 'GlobalSearch')
        % The GlobalSearch functionality is part of MATLAB's Global Optimization Toolbox.
        % Initialize a GlobalSearch object to perform the optimization.
        gs = GlobalSearch;
        problem = createOptimProblem('fmincon', ...
            'x0', initial_guess, ...
            'objective', optimize_phases, ...
            'lb', lower_bounds, ...
            'ub', upper_bounds, ...
            'options', optimoptions('fmincon', 'OptimalityTolerance', 1e-8));

        % Perform the global search to find optimal phases and velocity.
        [opt_phases_and_velocity, min_err] = run(gs, problem);
    end

    % Extract optimized phases and velocity from the result.
    opt_phases      = phases_from_params(opt_phases_and_velocity);
    opt_velocity    = opt_phases_and_velocity(end);
    opt_params_raw  = opt_phases_and_velocity;

    switch precession_mode
        case 'linear'
            fprintf('Precession params (linear): start=%.2f deg, step=%.2f deg/elem\n', ...
                opt_phases_and_velocity(1)/pi*180, opt_phases_and_velocity(2)/pi*180);
        case 'monotonic'
            deltas = opt_phases_and_velocity(2:n_elem);
            fprintf('Precession params (monotonic): start=%.2f deg, increments=%s deg\n', ...
                opt_phases_and_velocity(1)/pi*180, mat2str(round(deltas/pi*180)));
        otherwise
            % unconstrained: phases reported below via standard output
    end

    % Plot and evaluate the optimization result.
    phase_optimization_annulus_full_curve(...
        opt_phases, ... % Optimized phases
        parameters, ... % Simulation and transducer parameters
        opt_velocity, ... % Optimized velocity
        profile_target.axial_distance_bowl, ... % Distance vector
        profile_target.axial_intensity,... % Desired intensity profile
        1, ... % Enable plotting
        opt_limits, ...
        weights);

    % Save the figure
    img_folder = fullfile(parameters.io.outputs_folder, 'img_calibration');
    if ~exist(img_folder, 'dir'); mkdir(img_folder); end
    fig_path = fullfile(img_folder, sprintf('GlobalSearch.png'));
    saveas(gcf, fig_path);
    close(gcf);

    % The left plot above shows the real and the fitted profiles along with the 
    % cost function used for fitting, while the right shows the error in fitting (the 
    % squared difference between the real and the fitted profile weighted by the cost 
    % function). 

    % Display optimization results.
    fprintf('Optimal phases: %s deg; velocity: %.4f m/s; optimization error: %.4f \n', ...
        mat2str(round(opt_phases / pi * 180)), opt_velocity, min_err);
    
end