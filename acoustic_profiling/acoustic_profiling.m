function acoustic_profiling(...
    profile_empirical,...
    equipment_name, ...
    desired_intensity, ...
    parameters, ...
    sim_id)

% acoustic_profiling - Simulate and optimize ultrasonic transducer acoustic profile
%
% Inputs:
%   profile_empirical
%       profile_focus        - Measured or theoretical intensity profile along beam axis
%       dist_from_tran       - Distance from transducer reference point (mm)
%       focus_wrt_exit_plane - Distance related to transducer exit plane (mm)
%   equipment_name       - Name/identifier for the transducer or equipment
%   desired_intensity    - Target focal intensity (W/cm^2)
%   parameters           - Structure with simulation parameters and paths
%   .calibration.path_output_profiles  - Directory path for saving outputs
%   .calibration.filename_calibrated_CSV - Virtual path for saving optimized parameter files
%   sim_id               - Numeric ID of the subject (or simulation)
%
% Outputs:
%   None (results and plots are saved to disk)
%
% Description:
%   This function scales the input intensity profile to the desired value,
%   runs an initial simulation using the specified submission method,
%   computes analytical and simulated pressure profiles,
%   performs optimization of transducer element phases and amplitudes 
%   to match the desired acoustic profile,
%   reruns the simulation with optimized parameters,
%   visualizes results, and saves optimized values.

    %% Scale the profile to desired intensity

    profile_opt.adjusted_profile_focus = scale_real_intensity_profile(...
        parameters, ...
        desired_intensity, ...
        profile_empirical.profile_focus);

    %% Initial water simulation

    disp('Run initial simulation...')

    % Set up simulation parameters

    % Run all simulations in calibration folder (default)
    if parameters.calibration.save_in_calibration_folder
        parameters.data_path = parameters.calibration.path_output;
        parameters.subject_subfolder = 0;
        parameters.sim_path = parameters.calibration.path_output;
    end

    % Copy calibration settings to relevant entries in simulation config
    sim_param = parameters;
    sim_param.submit_medium = parameters.calibration.submit_medium;

    % Manage the submission setup
    sim_param.simulation_medium = 'water';
    if strcmp(sim_param.submit_medium, 'matlab') == true
        sim_param.code_type = 'matlab_cpu';
        sim_param.using_donders_hpc = 0;
    else
        sim_param.overwrite_files = 'always';
        sim_param.interactive = 0;
    end
    
    % Run the simulation based on the submission method
    switch sim_param.submit_medium
        case 'qsub'
            single_subject_pipeline_with_qsub(sim_id, sim_param, true);
        case 'slurm'
            single_subject_pipeline_with_slurm(sim_id, sim_param, true);
        case 'matlab'
            single_subject_pipeline(sim_id, sim_param);
        otherwise
            error('Submit medium does not correspond to available options.');
    end

    % Load initial results
    outputs_folder = sprintf('%s', sim_param.data_path);
    initial_res = load(sprintf('%s/sub-%03d_water_results%s.mat', ...
        outputs_folder, sim_id, sim_param.results_filename_affix),...
        'sensor_data','parameters');
    
    % Get maximum pressure
    p_max = gather(initial_res.sensor_data.p_max_all); % transform from GPU array to normal array
    
    % Plot 2D intensity map
    figure;
    imagesc((1:size(p_max, 1)) * initial_res.parameters.grid_step_mm, ...
        (1:size(p_max, 3)) * initial_res.parameters.grid_step_mm, ...
        squeeze(p_max(:, initial_res.parameters.transducer.pos_grid(2), :))')
    axis image;
    colormap(getColorMap);
    xlabel('Lateral Position [mm]');
    ylabel('Axial Position [mm]');
    axis image;
    cb = colorbar;
    title('Pressure for the focal plane')
    
    % Save the intensity map
    fig_path = fullfile(parameters.calibration.path_output_profiles, ...
        strcat('Intensity_map_2D_at_F_', num2str(profile_empirical.focus_wrt_exit_plane), ...
        '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    saveas(gcf, fig_path);

    %% Optimization

    % Extract simulated pressure along the focal axis
    i_x = initial_res.parameters.transducer.pos_grid(1);
    i_y = initial_res.parameters.transducer.pos_grid(2);
    pred_axial_pressure = squeeze(p_max(i_x, i_y,:));
    
    % Compute O'Neil solution and related parameters
    [p_axial_oneil, simulated_grid_adj_factor, velocity, axial_position] = ...
        compute_oneil_solution(...
        initial_res.parameters, ...
        pred_axial_pressure, ...
        profile_empirical.dist_from_tran, ...
        profile_opt.adjusted_profile_focus, ...
        profile_empirical.focus_wrt_exit_plane, ...
        desired_intensity, ...
        equipment_name);
    
    % Optimize phases and source amplitude to match real profile
    [opt_phases, opt_velocity, min_err] = ...
        perform_global_search(...
        initial_res.parameters, ...
        profile_empirical.dist_from_tran, ...
        profile_opt.adjusted_profile_focus, ...
        velocity);

    % Recalculate analytical solution with optimized phases and velocity
    p_axial_oneil_opt = ...
        recalculate_analytical_sol(...
        initial_res.parameters, ...
        p_axial_oneil, ...
        opt_phases, ...
        opt_velocity, ...
        profile_empirical.dist_from_tran, ...
        profile_opt.adjusted_profile_focus, ...
        axial_position, ...
        profile_empirical.focus_wrt_exit_plane, ...
        desired_intensity, ...
        equipment_name);
    
    % Calculate optimized source amplitude
    opt_source_amp = round(opt_velocity / velocity * ...
        initial_res.parameters.transducer.source_amp / ...
        simulated_grid_adj_factor);

    fprintf('The optimized source_amp = %i\n', opt_source_amp(1));

    %% Rerun water simulation with optimized phases and source amplitude

    opt_param = sim_param;
    opt_param.transducer.source_amp = opt_source_amp;
    opt_param.transducer.source_phase_rad = [0 opt_phases];
    opt_param.transducer.source_phase_deg = [0 opt_phases]/pi*180;
    opt_param.results_filename_affix = '_optimized';

    switch sim_param.submit_medium
        case 'qsub'
            single_subject_pipeline_with_qsub(sim_id, opt_param, true);
        case 'slurm'
            single_subject_pipeline_with_slurm(sim_id, opt_param, true);
        case 'matlab'
            single_subject_pipeline(sim_id, opt_param);
        otherwise
            error('Submit medium does not correspond to available options.');
    end

    % Plot optimized simulation results
    plot_opt_sim_results(...
        opt_param, ...
        outputs_folder, ...
        sim_id, ...
        axial_position, ...
        profile_empirical.dist_from_tran, ...
        profile_opt.adjusted_profile_focus, ...
        p_axial_oneil_opt, ...
        p_axial_oneil, ...
        profile_empirical.focus_wrt_exit_plane, ...
        desired_intensity, ...
        equipment_name, ...
        min_err)

    % Save optimized values
    save_optimized_values(...
        parameters, ...
        profile_empirical.focus_wrt_exit_plane, ...
        desired_intensity, ...
        opt_param, ...
        equipment_name);

end