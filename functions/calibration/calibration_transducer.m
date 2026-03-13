function [opt_source_amp, opt_source_phase_deg, opt_source_phase_rad] = calibration_transducer(...
    profile_empirical,...
    equipment_name, ...
    desired_intensity, ...
    desired_focal_distance_ep, ...
    parameters, ...
    sim_id)

% calibration_transducer - Simulate and optimize ultrasonic transducer acoustic profile
%
% Inputs:
%   profile_empirical
%       axial_intensity     - Measured or theoretical intensity profile along beam axis
%       axial_distance_bowl - Distance from transducer reference point (mm)
%   equipment_name          - Name/identifier for the transducer or equipment
%   desired_intensity       - Target focal intensity (W/cm^2)
%   desired_focal_distance_ep - Desired focal distance from transducer exit plane (mm)
%   parameters              - Structure with simulation parameters and paths
%   parameters.calibration  - Structure with calibration settings
%   sim_id                  - Numeric ID of the subject (or simulation)
%
% Outputs:
%   opt_source_amp          | Optimized amplitude
%   opt_source_phase_deg    | Optimized phases (degrees)
%   opt_source_phase_rad    | Optimized phases (radians)
%
% Description:
%   This function scales the input intensity profile to the desired value,
%   runs an initial simulation using the specified submission method,
%   computes analytical and simulated intensity profiles,
%   performs optimization of transducer element phases and amplitudes 
%   to match the desired acoustic profile,
%   reruns the simulation with optimized parameters,
%   visualizes results, and saves optimized values.

    %% Attach more information to parameters

    parameters.calibration.equipment_name = equipment_name;
    parameters.calibration.desired_intensity = desired_intensity;
    parameters.calibration.desired_focal_distance_ep = desired_focal_distance_ep;

    % Scale the requested profile to the desired intensity
    [profile_target, ~] = scale_real_intensity_profile(...
        parameters, ...
        profile_empirical, ...
        desired_intensity);

    % We continue with the scaled empirical profile
    clear profile_empirical;

    %% Initial water simulation

    disp('Run initial simulation...')

    % Set up simulation parameters

    % Run all simulations in calibration folder (default)
    if parameters.calibration.save_in_calibration_folder
        parameters.data_path = parameters.calibration.path_output;
        parameters.seg_path = parameters.calibration.path_output;
        parameters.sim_path = parameters.calibration.path_output;
    end
    disp(['Saving free-water calibration in ', parameters.sim_path]);

    % if a subfolder is requested, move outputs to subfolders
    if parameters.subject_subfolder
        parameters.outputs_folder = sprintf('%s/sub-%03d', parameters.sim_path, sim_id);
    else
        parameters.outputs_folder = sprintf('%s', parameters.sim_path);
    end

    % Copy calibration settings to relevant entries in simulation config
    sim_param = parameters;
    % Force water medium
    sim_param.simulation_medium = 'water';
    % Force save result matrices
    sim_param.savemat = 1;
    % Overwrite transducer kwavearray modeling (if specified)
    if isfield(parameters.calibration, 'force_kwavearray') && ...
            parameters.calibration.force_kwavearray == 1
        sim_param.use_kwavearray = 1; % force to run with kwavearray setup
    end
    % Convert from default 3D to 2D axisymmetric simulation (if requested)
    if isfield(parameters.calibration, 'axisymmetric2D') && ...
            parameters.calibration.axisymmetric2D == 1
        parameters.n_sim_dims = 2;
        parameters.axisymmetric = 1;
        if numel(parameters.default_grid_dims)==3
            parameters.default_grid_dims(2) = [];
        end
    end
    % Force deactivate interactive mode
    sim_param.interactive = 0;
    
    % Run the simulation based on the submission method
    sim_param.hpc_wait_for_completion = true;
    prestus_pipeline_start(sim_id, sim_param);

    %% Load initial results

    initial_res = load(sprintf('%s/sub-%03d_water_results%s.mat', ...
        sim_param.outputs_folder, sim_id, sim_param.results_filename_affix));

    initial_params = initial_res.acoustic_info.parameters;
    initial_params.calibration.prefix = 'Initial_';
    
    %% Extract simulated intensity along the focal axis

    [profile_sim] = extract_simulated_profile(initial_res, initial_params);
    
    %% Optimization
    
    % Compute analytical O'Neil solution and scaling factor to simulated intensity
    [profile_oneil, simulated_analytical_scaling] = ...
        compute_oneil_solution(...
        initial_params, ...
        profile_sim, ...
        profile_target);

    % Optimize phases [rad] and source amplitude to match real profile
    [opt_phases, opt_velocity, min_err] = ...
        perform_global_search(...
        initial_params, ...
        profile_target, ...
        profile_sim.velocity);

    % Fit velocity to match desired peak intensity exactly
    if ~isfield(parameters.calibration, 'fit_velocity_to_intensity') || ...
            parameters.calibration.fit_velocity_to_intensity
        [opt_velocity, ~, ~] = fit_velocity_to_intensity(...
            initial_params, profile_oneil, opt_phases, opt_velocity, ...
            parameters.calibration.desired_intensity, simulated_analytical_scaling);
    end

    % Recalculate analytical solution with optimized phases and velocity
    profile_oneil_opt = ...
        recompute_oneil_solution(...
        initial_params, ...
        profile_oneil, ...
        profile_target, ...
        opt_phases, ...
        opt_velocity);
    
    % Calculate optimized source amplitude
    opt_source_amp = round(opt_velocity / profile_sim.velocity * ...
        initial_params.transducer.source_amp / ...
        simulated_analytical_scaling);

    % Collect phases
    opt_source_phase_rad = opt_phases;
    opt_source_phase_deg = opt_phases/pi*180;

    fprintf('The optimized source_amp = %i\n', opt_source_amp(1));

    %% Rerun water simulation with optimized phases and source amplitude

    opt_param = sim_param;
    opt_param.transducer.source_amp = opt_source_amp;
    opt_param.transducer.source_phase_rad = opt_source_phase_rad;
    opt_param.transducer.source_phase_deg = opt_source_phase_deg;
    opt_param.results_filename_affix = '_optimized';

    sim_param.hpc_wait_for_completion = true;
    prestus_pipeline_start(sim_id, sim_param);

    %% Load optimized simulation results    
    opt_res = load(sprintf('%s/sub-%03d_water_results%s.mat', ...
        opt_param.outputs_folder, sim_id, opt_param.results_filename_affix));

    opt_params = opt_res.acoustic_info.parameters;
    opt_params.calibration.prefix = 'Opt_';
    
    %% Extract simulated intensity along the focal axis
    [profile_sim_opt] = extract_simulated_profile(opt_res, opt_params);

    % Plot optimized simulation results
    plot_opt_sim_results(...
        opt_param, ...
        profile_target, ...
        profile_oneil, ...
        profile_oneil_opt, ...
        profile_sim, ...
        profile_sim_opt, ...
        min_err)

    % Save optimized values
    save_optimized_values(...
        opt_param);
    
end
