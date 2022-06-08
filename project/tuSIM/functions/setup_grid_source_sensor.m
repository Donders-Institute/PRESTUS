function [kgrid, source, sensor, source_labels] = setup_grid_source_sensor(parameters, max_sound_speed, trans_pos_final, focus_pos_final)
    % create simulation grid
    if parameters.n_sim_dims == 3
        kgrid = kWaveGrid(parameters.grid_dims(1), parameters.grid_step_m, ...
                      parameters.grid_dims(2), parameters.grid_step_m, ...
                      parameters.grid_dims(3), parameters.grid_step_m);
    else
        kgrid = kWaveGrid(parameters.grid_dims(1), parameters.grid_step_m, ...
                      parameters.grid_dims(2), parameters.grid_step_m);
    end
        
	% calculate the time step using an integer number of points per period
    points_per_wavelength = max_sound_speed / (parameters.transducer.source_freq_hz * parameters.grid_step_m);    % points per wavelength
    cfl = 0.3;                                                % CFL number (kwave default)
    points_per_period = ceil(points_per_wavelength / cfl);    % points per period
    wave_period   = 1 / parameters.transducer.source_freq_hz;                   % period [s]
    grid_time_step = (wave_period / points_per_period)/2;     % time step [s]

    % calculate the number of time steps to reach steady state
    t_end = sqrt(kgrid.x_size.^2  + kgrid.z_size.^2 + kgrid.y_size.^2) / max_sound_speed;  % [s]
    simulation_time_points = round(t_end / grid_time_step);

    % create the time array
    kgrid.setTime(simulation_time_points, grid_time_step);

    % create source
    parameters.kwave_source_filename  = fullfile(parameters.data_path, sprintf('sub-%03d_%s_kwave_source%s.mat', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    if confirm_overwriting(parameters.kwave_source_filename, parameters)
        [source, source_labels, ~] = setup_source(parameters, kgrid, trans_pos_final, focus_pos_final);
        save(parameters.kwave_source_filename, 'source', 'source_labels','-v7.3');
    else
        load(parameters.kwave_source_filename);
    end    

    % create sensor
    sensor = struct();
    sensor.mask = ones(parameters.grid_dims);
    sensor.record = {'p_max_all','p_final'};

    % record the last 3 cycles in steady state
    num_periods = 3;
    time_points_to_record = round(num_periods * wave_period / kgrid.dt);
    sensor.record_start_index = simulation_time_points - time_points_to_record + 1;

end