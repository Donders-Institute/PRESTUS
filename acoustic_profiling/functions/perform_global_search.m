function [opt_phases, opt_velocity, min_err] = perform_global_search(parameters, dist_exit_plane, adjusted_profile_focus, velocity)
    % Perform a global search to optimize transducer phases and particle velocity.
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - dist_exit_plane: Distance vector for the desired intensity profile [mm].
    % - adjusted_profile_focus: Adjusted desired intensity profile [W/cm^2].
    % - velocity: Initial particle velocity estimate [m/s].
    %
    % Returns:
    % - opt_phases: Optimized phases for each transducer element [rad].
    % - opt_velocity: Optimized particle velocity [m/s].
    % - min_err: Minimum error achieved during optimization.

    % Initialize a GlobalSearch object to perform the optimization.
    gs = GlobalSearch;
    
    % Define the objective function for optimization.
    % This function evaluates the error based on the current phases and velocity.
    % Define the objective function for optimization.
    % This function evaluates the error based on the current phases and velocity.
    optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve(...
        phases_and_velocity(1:(parameters.transducer.n_elements - 1)), ... % Phases for transducer elements
        parameters, ... % Simulation and transducer parameters
        phases_and_velocity(parameters.transducer.n_elements), ... % Particle velocity
        dist_exit_plane, ... % Distance vector
        adjusted_profile_focus); % Desired intensity profile
    
    % Set a random seed for reproducibility.
    rng(195, 'twister');

    % Define initial guess, bounds, and options for the optimization problem.
    initial_guess = [randi(360, [1, parameters.transducer.n_elements - 1]) / 180 * pi, velocity];
    lower_bounds = zeros(1, parameters.transducer.n_elements); % Lower bounds: [0 rad, 0 m/s]
    upper_bounds = [2 * pi * ones(1, parameters.transducer.n_elements - 1), 0.2]; % Upper bounds: [2pi rad, 0.2 m/s]
    
    problem = createOptimProblem('fmincon', ...
        'x0', initial_guess, ...
        'objective', optimize_phases, ...
        'lb', lower_bounds, ...
        'ub', upper_bounds, ...
        'options', optimoptions('fmincon', 'OptimalityTolerance', 1e-8));
    
    % Perform the global search to find optimal phases and velocity.
    [opt_phases_and_velocity, min_err] = run(gs, problem);

    % Extract optimized phases and velocity from the result.
    opt_phases = opt_phases_and_velocity(1:(parameters.transducer.n_elements - 1));
    opt_velocity = opt_phases_and_velocity(parameters.transducer.n_elements);

    % Plot and evaluate the optimization result.
    phase_optimization_annulus_full_curve(...
        opt_phases, ... % Optimized phases
        parameters, ... % Simulation and transducer parameters
        opt_velocity, ... % Optimized velocity
        dist_exit_plane, ... % Distance vector
        adjusted_profile_focus, ... % Desired intensity profile
        1); % Enable plotting

    % Display optimization results.
    fprintf('Optimal phases: %s deg; velocity: %.4f m/s; optimization error: %.4f \n', ...
        mat2str(round(opt_phases / pi * 180)), opt_velocity, min_err);
    
end