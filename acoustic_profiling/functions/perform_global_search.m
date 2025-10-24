function [opt_phases, opt_velocity, min_err] = perform_global_search(parameters, axial_position, adjusted_profile_focus, velocity)
    % Perform a global search to optimize transducer phases and particle velocity.
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - axial_position: Distance vector for the intensity profile [mm from transducer bowl].
    % - adjusted_profile_focus: Adjusted desired intensity profile [W/cm^2].
    % - velocity: Initial particle velocity estimate [m/s].
    %
    % Returns:
    % - opt_phases: Optimized phases for each transducer element [rad].
    % - opt_velocity: Optimized particle velocity [m/s].
    % - min_err: Minimum error achieved during optimization.
    
    if ~isfield(parameters.calibration, 'opt_limits')
        opt_limits = [1, max(axial_position)];
    else
        opt_limits = parameters.calibration.opt_limits;
    end

    if ~isfield(parameters.calibration, 'weights') 
        weights = 1;
    else
        weights = parameters.calibration.weights;
    end

    % Define the objective function for optimization.
    % This function evaluates the error based on the current phases and velocity.
    optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve(...
        phases_and_velocity(1:end-1), ... % Phases for transducer elements
        parameters, ... % Simulation and transducer parameters
        phases_and_velocity(end), ... % Particle velocity
        axial_position, ... % Distance vector
        adjusted_profile_focus,... % Desired intensity profile
        0, ... % Disable plotting
        opt_limits, ...
        weights);
    
    % Set a random seed for reproducibility.
    if isfield(parameters.calibration, 'seed') 
        rng(parameters.calibration.seed, 'twister');
    end

    % Define initial guess, bounds, and options for the optimization problem.
    initial_guess = [randi(360, [1, parameters.transducer.n_elements]) / 180 * pi, velocity];
    lower_bounds = zeros(1, parameters.transducer.n_elements+1); % Lower bounds: [0 rad, 0 m/s]
    upper_bounds = [2 * pi * ones(1, parameters.transducer.n_elements), 0.2]; % Upper bounds: [2pi rad, 0.2 m/s]

    if ~isfield(parameters.calibration, 'optmethod') || strcmp(parameters.calibration.optmethod, 'FEXminimize')
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
    elseif strcmp(parameters.calibration.optmethod, 'GlobalSearch')
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
    opt_phases = opt_phases_and_velocity(1:end-1);
    opt_velocity = opt_phases_and_velocity(end);

    % Plot and evaluate the optimization result.
    phase_optimization_annulus_full_curve(...
        opt_phases, ... % Optimized phases
        parameters, ... % Simulation and transducer parameters
        opt_velocity, ... % Optimized velocity
        axial_position, ... % Distance vector
        adjusted_profile_focus,... % Desired intensity profile
        1, ... % Enable plotting
        opt_limits, ...
        weights);

    % The left plot above shows the real and the fitted profiles along with the 
    % cost function used for fitting, while the right shows the error in fitting (the 
    % squared difference between the real and the fitted profile weighted by the cost 
    % function). 

    % Display optimization results.
    fprintf('Optimal phases: %s deg; velocity: %.4f m/s; optimization error: %.4f \n', ...
        mat2str(round(opt_phases / pi * 180)), opt_velocity, min_err);
    
end