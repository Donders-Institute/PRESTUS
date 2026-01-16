function [kgrid, source, sensor, source_labels] = source_sensor_setup(parameters, max_sound_speed, trans_pos_final, focus_pos_final, grid_time_step)
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                  Set up the transducer and sensor                 %
    %                                                                   %
    % This function sets up the transducer in the grid, the timeperiod  %
    % during which simulations will take place and the sensor that      %
    % records the pressure-levels in the grid.                          %     
    % Time axis is set from the (first) transducer source frequency.    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % Creates a simulation grid in 3 or 2 dimensions
    if parameters.n_sim_dims == 3
        kgrid = kWaveGrid(parameters.grid_dims(1), parameters.grid_step_m, ...
                      parameters.grid_dims(2), parameters.grid_step_m, ...
                      parameters.grid_dims(3), parameters.grid_step_m);
    elseif parameters.n_sim_dims == 2
        kgrid = kWaveGrid(parameters.grid_dims(1), parameters.grid_step_m, ...
                      parameters.grid_dims(2), parameters.grid_step_m);
    end

    % Backward-compatible access to (first) transducer
    tx = parameters.transducer(1);
    if numel(parameters.transducer)>1
        warning('Individual source frequencies per transducer are not supported yet. Using %i Hz for grid time axis.', tx.source_freq_hz)
    end

    wave_period = 1 / tx.source_freq_hz;    % period [s]
    
    % Check the number of input arguments
    % As a default the time step is based on the default CFL number of 0.3.
    % If grid_time_step is given, it means that the transducer and sensor
    % has to be set up again with a smaller time step due to simulation
    % instability.
    if nargin < 5
        % Calculate the time step using an integer number of points per period
        points_per_wavelength = max_sound_speed /(tx.source_freq_hz * parameters.grid_step_m);
        % PPW: Spatial samples per wavelength at source freq; ensures dx resolves waves (target ≥3).
        cfl = 0.3;                                                          
        % Courant-Friedrichs-Lewy: Fraction of dx/c for dt; k-Wave default.
        points_per_period = ceil(points_per_wavelength / cfl);              
        % Temporal samples per wave period; ceil enforces integer ≥ PPW/CFL.
        grid_time_step = (wave_period / points_per_period)/2;               
        % dt: Half the standard dt, conservative; wave_period=1/source_freq.    % time step [s]  
     end

    % Calculate the number of time steps to reach steady state
    t_end = sqrt(kgrid.x_size.^2 + kgrid.z_size.^2 + kgrid.y_size.^2) / max_sound_speed;    % [s]
    simulation_time_points = round(t_end / grid_time_step);

    % Create the time array
    kgrid.setTime(simulation_time_points, grid_time_step);

    % Create source (transducer with focuspoint)
    parameters.kwave_source_filename  = fullfile(parameters.output_dir, sprintf('sub-%03d_%s_kwave_source%s.mat', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    if confirm_overwriting(parameters.kwave_source_filename, parameters)

        % build per-transducer positions for geometry when not using kWaveArray
        if numel(parameters.transducer)>1 && parameters.use_kWaveArray == 0
            nT = numel(parameters.transducer);
            trans_pos = zeros(nT, parameters.n_sim_dims);
            focus_pos = zeros(nT, parameters.n_sim_dims);

            for ti = 1:nT
                tr = parameters.transducer(ti);

                % transducer position in grid
                if isfield(tr, 'pos_grid') && ~isempty(tr.pos_grid)
                    pos = tr.pos_grid;
                    if size(pos,1) > size(pos,2)
                        pos = pos';
                    end
                    trans_pos(ti,:) = pos;
                else
                    trans_pos(ti,:) = trans_pos_final;
                end

                % focus position in grid
                if isfield(tr, 'focus_pos_grid') && ~isempty(tr.focus_pos_grid)
                    fpos = tr.focus_pos_grid;
                    if size(fpos,1) > size(fpos,2)
                        fpos = fpos';
                    end
                    focus_pos(ti,:) = fpos;
                else
                    focus_pos(ti,:) = focus_pos_final;
                end
            end
        else
            % legacy / single-transducer / kWaveArray case
            trans_pos = trans_pos_final;
            focus_pos = focus_pos_final;
        end
        [source, source_labels, ~] = source_create(parameters, kgrid, trans_pos, focus_pos);
        save(parameters.kwave_source_filename, 'source', 'source_labels','-v7.3');
    else
        load(parameters.kwave_source_filename);
    end

    % Creates a sensor that records the the maximum and final pressure
    % values in every point of the grid
    sensor = struct();
    sensor.mask = ones(parameters.grid_dims);
    sensor.record = {'p_max_all','p_final'};

    % Record the last 3 cycles in steady state (when sonic waves have traversed the entire medium)
    num_periods = 3;
    time_points_to_record = round(num_periods * wave_period / kgrid.dt);
    sensor.record_start_index = simulation_time_points - time_points_to_record + 1;
    % alternative (EM): sensor.record_start_index = simulation_time_points - 3*points_per_period;

end
