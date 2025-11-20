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
%       focus_wrt_exit_plane - Desired focal distance from transducer exit plane (mm)
%   equipment_name       - Name/identifier for the transducer or equipment
%   desired_intensity    - Target focal intensity (W/cm^2)
%   parameters           - Structure with simulation parameters and paths
%   parameters.calibration - Structure with calibration settings
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

    %% Initial water simulation

    disp('Run initial simulation...')

    % Set up simulation parameters

    % Run all simulations in calibration folder (default)
    if parameters.calibration.save_in_calibration_folder
        parameters.data_path = parameters.calibration.path_output;
        parameters.seg_path = parameters.calibration.path_output;
        parameters.sim_path = parameters.calibration.path_output;
    end
    disp(['Saving free-water calibration in ', parameters.sim_path{1}]);

    % if a subfolder is requested, move outputs to subfolders
    if parameters.subject_subfolder
        parameters.outputs_folder = sprintf('%s/sub-%03d', parameters.sim_path, sim_id);
    else
        parameters.outputs_folder = sprintf('%s', parameters.sim_path);
    end

    % Copy calibration settings to relevant entries in simulation config
    sim_param = parameters;
    sim_param.submit_medium = parameters.calibration.submit_medium;

    % Manage the submission setup
    sim_param.simulation_medium = 'water';
    sim_param.save_mat = 0; % do not save full results
    if isfield(parameters.calibration, 'force_kwavearray') && ...
            parameters.calibration.force_kwavearray == 1
        sim_param.use_kwavearray = 1; % force to run with kwavearray setup
    end
    if strcmp(sim_param.submit_medium, 'matlab') == true
        sim_param.code_type = 'matlab_cpu';
        sim_param.using_donders_hpc = 0;
    else
        sim_param.overwrite_files = 'never';
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

    %% Load initial results

    initial_res = load(sprintf('%s/sub-%03d_water_results%s.mat', ...
        sim_param.outputs_folder, sim_id, sim_param.results_filename_affix),...
        'sensor_data','parameters');
    
    %% Get maximum pressure

    p_max = gather(initial_res.sensor_data.p_max_all); % transform from GPU array to normal array
    
    %% Plot 2D intensity map
    
    figure;
    p_distance = (1:size(p_max, 1)) * initial_res.parameters.grid_step_mm;
    p_width = (1:size(p_max, sim_param.n_sim_dims)) * initial_res.parameters.grid_step_mm;
    if sim_param.n_sim_dims == 2
        p_axialprofile = squeeze(p_max(:, :))';
    elseif sim_param.n_sim_dims == 3
        p_axialprofile = squeeze(p_max(:, initial_res.parameters.transducer.pos_grid(2), :))';
    end
    imagesc(p_distance, p_width, p_axialprofile);
    clear p_distance p_width p_axialprofile;
    axis image;
    colormap(getColorMap);
    xlabel('Lateral Position [mm]');
    ylabel('Axial Position [mm]');
    axis image;
    cb = colorbar;
    title('Pressure for the focal plane')
    
    % Save the intensity map
    fig_path = fullfile(sim_param.outputs_folder, ...
        strcat('Initial_Intensity_map_2D_at_F_', num2str(profile_empirical.focus_wrt_exit_plane), ...
        '_at_I_', num2str(desired_intensity), '_', equipment_name, '.png'));
    saveas(gcf, fig_path);
    close(gcf);

    %% Optimization

    % Extract simulated pressure along the focal axis
    if sim_param.n_sim_dims == 2
        i_x = initial_res.parameters.transducer.pos_grid(1);
        pred_axial_pressure = squeeze(p_max(i_x,:));
        clear i_x;
    elseif sim_param.n_sim_dims == 3
        i_x = initial_res.parameters.transducer.pos_grid(1);
        i_y = initial_res.parameters.transducer.pos_grid(2);
        pred_axial_pressure = squeeze(p_max(i_x, i_y,:));
        clear i_x i_y;
    end

    % Scale the profile to the desired intensity
    profile_opt.adjusted_profile_focus = scale_real_intensity_profile(...
        initial_res.parameters, ...
        desired_intensity, ...
        profile_empirical.profile_focus);
    
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

    % Optimize phases [rad] and source amplitude to match real profile
    [opt_phases, opt_velocity, min_err] = ...
        perform_global_search(...
        initial_res.parameters, ...
        profile_empirical.dist_from_tran, ...
        profile_opt.adjusted_profile_focus', ...
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
    opt_param.transducer.source_phase_rad = opt_phases;
    opt_param.transducer.source_phase_deg = opt_phases/pi*180;
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
        opt_param, ...
        profile_empirical.focus_wrt_exit_plane, ...
        desired_intensity, ...
        equipment_name);

end