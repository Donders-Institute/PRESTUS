function [opt_phases, opt_velocity, min_err] = perform_global_search(parameters, profile_target, velocity)
    % Perform a global search to optimize transducer phases and particle velocity.
    %
    % Arguments:
    % - parameters: Structure containing simulation and transducer parameters.
    % - profile_target.axial_distance_bowl: Distance vector for the intensity profile [mm from transducer bowl].
    % - profile_target.axial_intensity: Adjusted desired intensity profile [W/cm^2].
    % - velocity: Initial particle velocity estimate [m/s].
    %
    % Returns:
    % - opt_phases: Optimized phases for each transducer element [rad].
    % - opt_velocity: Optimized particle velocity [m/s].
    % - min_err: Minimum error achieved during optimization.
    
    if ~isfield(parameters.calibration, 'opt_limits')
        % By default set to min and max distance with non-NAN intensity
        min_available = min(profile_target.axial_distance_bowl(~isnan(profile_target.axial_intensity)));
        max_available = max(profile_target.axial_distance_bowl(~isnan(profile_target.axial_intensity)));
        opt_limits = [min_available, max_available];

    else
        opt_limits = parameters.calibration.opt_limits;
    end

    if ~isfield(parameters.calibration, 'weights') 
        weights = 0; % uniform weighting
    else
        weights = parameters.calibration.weights;
    end

    % Define the objective function for optimization.
    % This function evaluates the error based on the current phases and velocity.
    optimize_phases = @(phases_and_velocity) phase_optimization_annulus_full_curve(...
        phases_and_velocity(1:end-1), ... % Phases for transducer elements
        parameters, ... % Simulation and transducer parameters
        phases_and_velocity(end), ... % Particle velocity
        profile_target.axial_distance_bowl, ... % Distance vector
        profile_target.axial_intensity,... % Desired intensity profile
        0, ... % Disable plotting
        opt_limits, ...
        weights);
    
    % Set a random seed for reproducibility.
    if isfield(parameters.calibration, 'seed') 
        rng(parameters.calibration.seed, 'twister');
    end

    % Define initial guess, bounds, and options for the optimization problem.
    if ~isfield(parameters.calibration, 'gs_upper_velocity') || isempty(parameters.calibration.gs_upper_velocity)
        parameters.calibration.gs_upper_velocity = 0.2; % set default for upper velocity to 20 mm/s;
    end
    initial_guess = [randi(360, [1, parameters.transducer.n_elements]) / 180 * pi, velocity];
    lower_bounds = [zeros(1, parameters.transducer.n_elements), 0.001]; % Lower bounds: [0 rad, 1 mm/s]
    upper_bounds = [2 * pi * ones(1, parameters.transducer.n_elements), parameters.calibration.gs_upper_velocity]; % Upper bounds: [2pi rad, 200 mm/s]

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
        profile_target.axial_distance_bowl, ... % Distance vector
        profile_target.axial_intensity,... % Desired intensity profile
        1, ... % Enable plotting
        opt_limits, ...
        weights);

    % Save the figure
    fig_path = fullfile(parameters.outputs_folder, sprintf('GlobalSearch.png'));
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