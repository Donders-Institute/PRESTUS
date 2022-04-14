function [kgrid, source, sensor, source_labels] = setup_grid_source_sensor(parameters, max_sound_speed, trans_pos_final, focus_pos_final)
    % create simulation grid
    kgrid = kWaveGrid(parameters.grid_dims(1), parameters.grid_step_m, ...
                      parameters.grid_dims(2), parameters.grid_step_m, ...
                      parameters.grid_dims(3), parameters.grid_step_m);

    source = struct();

    [source.p_mask, source_labels] = transducer_setup(parameters.transducer, trans_pos_final, focus_pos_final, ...
                                                            parameters.grid_dims, parameters.grid_step_mm);
    
	% calculate the time step using an integer number of points per period
    points_per_wavelength = max_sound_speed / (parameters.transducer.source_freq_hz * parameters.grid_step_m);    % points per wavelength
    cfl = 0.3;                                                % cfl number????
    points_per_period = ceil(points_per_wavelength / cfl);    % points per period
    wave_period   = 1 / parameters.transducer.source_freq_hz;                   % period [s]
    grid_time_step = (wave_period / points_per_period)/2;     % time step [s]

    % calculate the number of time steps to reach steady state
    t_end = sqrt(kgrid.x_size.^2  + kgrid.z_size.^2 + kgrid.y_size.^2) / max_sound_speed;  % [s]
    simulation_time_points = round(t_end / grid_time_step);

    % create the time array
    kgrid.setTime(simulation_time_points, grid_time_step);

    % define the input signal for each element
    cw_signal = createCWSignals(kgrid.t_array, parameters.transducer.source_freq_hz, parameters.transducer.source_amp, parameters.transducer.source_phase_rad);

    % ensure each source point has an assigned time series of the source signal
    p_mask_source_p = source_labels;
    p_mask_source_p(p_mask_source_p(:,:,:) == 0) = [];
    p_mask_source_p = reshape(p_mask_source_p,[],1);
    source.p = zeros(length(p_mask_source_p),length(cw_signal));

    for ii = 1 : length(p_mask_source_p)
        source.p(ii, :) = cw_signal(p_mask_source_p(ii), :);
    end

    % create sensor
    sensor = struct();
    sensor.mask = ones(parameters.grid_dims);
    sensor.record = {'p_max_all'};

    % record the last 3 cycles in steady state
    num_periods = 3;
    time_points_to_record = round(num_periods * wave_period / kgrid.dt);
    sensor.record_start_index = simulation_time_points - time_points_to_record + 1;

end