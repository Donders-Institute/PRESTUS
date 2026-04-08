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
    if numel(parameters.grid.dims) == 3
        kgrid = kWaveGrid(parameters.grid.dims(1), parameters.grid.resolution_mm/1e3, ...
                      parameters.grid.dims(2), parameters.grid.resolution_mm/1e3, ...
                      parameters.grid.dims(3), parameters.grid.resolution_mm/1e3);
    elseif numel(parameters.grid.dims) == 2
        kgrid = kWaveGrid(parameters.grid.dims(1), parameters.grid.resolution_mm/1e3, ...
                      parameters.grid.dims(2), parameters.grid.resolution_mm/1e3);
    end

    % Backward-compatible access to (first) transducer
    tx = parameters.transducer(1);
    wave_period = 1 / tx.freq_hz;    % period [s]
    
    % Check the number of input arguments
    % As a default the time step is based on the default CFL number of 0.3.
    % If grid_time_step is given, it means that the transducer and sensor
    % has to be set up again with a smaller time step due to simulation
    % instability.
    if nargin < 5
        % Calculate the time step using an integer number of points per period
        % PPW: Spatial samples per wavelength at source freq; ensures dx resolves waves (target ≥3).
        if ~isfield(parameters.grid, 'source_ppw') || isempty(parameters.grid.source_ppw)
            points_per_wavelength = max_sound_speed /(tx.freq_hz * parameters.grid.resolution_mm/1e3);
        else
            points_per_wavelength = parameters.grid.source_ppw;
        end
        % Courant-Friedrichs-Lewy: Fraction of dx/c for dt; k-Wave default.
        if ~isfield(parameters.grid, 'source_cfl') || isempty(parameters.grid.source_cfl)
            cfl = 0.3;
        else
            cfl = parameters.grid.source_cfl;
        end
        % Temporal samples per wave period.
        points_per_period = ceil(points_per_wavelength / cfl);              
        % Calculate time step dt  
        grid_time_step = (wave_period / points_per_period);
        % Print parameter summary
        fprintf('Calculating the time step of %.1d based on PPW %d, CFL= %.1f and PPP %.1d.\n', ...
            grid_time_step, round(points_per_wavelength), cfl, points_per_period);
    else
        fprintf('Using the time step of %.1d based on approx. max. stable time step.\n', ...
            grid_time_step);
    end

    % Calculate the number of time steps to reach steady state
    t_end = sqrt(kgrid.x_size.^2 + kgrid.z_size.^2 + kgrid.y_size.^2) / max_sound_speed;    % [s]
    simulation_time_points = round(t_end / grid_time_step);

    % Create the time array
    kgrid.setTime(simulation_time_points, grid_time_step);

    % Create source (transducer with focuspoint).
    % Use io.preproc_affix when set (uncertainty pipeline points simulation
    % variants at the stage-1 cache with affix ''), otherwise fall back to
    % io.output_affix.
    if isfield(parameters.io, 'preproc_affix')
        source_affix = parameters.io.preproc_affix;
    else
        source_affix = parameters.io.output_affix;
    end
    parameters.io.kwave_source_filename  = fullfile(parameters.io.output_dir, ...
        sprintf('sub-%03d_%s_kwave_source%s.mat', ...
        parameters.subject_id, parameters.simulation.medium, source_affix));
    if confirm_overwriting(parameters.io.kwave_source_filename, parameters)

        % build per-transducer positions for geometry when not using kWaveArray
        if numel(parameters.transducer)>1 && parameters.grid.use_kWaveArray == 0
            nT = numel(parameters.transducer);
            trans_pos = zeros(nT, numel(parameters.grid.dims));
            focus_pos = zeros(nT, numel(parameters.grid.dims));

            for ti = 1:nT
                tr = parameters.transducer(ti);

                % transducer position in grid
                if isfield(tr, 'trans_pos') && ~isempty(tr.trans_pos)
                    pos = tr.trans_pos;
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
        % Create the source matrix
        [source, source_labels, ~] = source_create(parameters, kgrid, trans_pos, focus_pos);
        % Save the source matrix (controlled by io.save_source_matrices,
        % falling back to io.save_matrices, defaulting to true)
        if should_save_output(parameters.io, 'save_source_matrices')
            save(parameters.io.kwave_source_filename, 'source', 'source_labels','-v7.3');
        else
            disp('Not saving kwave source matrix ...')
        end
    else
        load(parameters.io.kwave_source_filename);
    end

    % Creates a sensor that records the the maximum and final pressure
    % values in every point of the grid
    sensor = struct();
    sensor.mask = ones(parameters.grid.dims);
    sensor.record = {'p_max_all','p_final'};

    % Record the last 3 cycles in steady state (when sonic waves have traversed the entire medium)
    num_periods = 3;
    time_points_to_record = round(num_periods * wave_period / kgrid.dt);
    sensor.record_start_index = simulation_time_points - time_points_to_record + 1;
    % alternative (EM): sensor.record_start_index = simulation_time_points - 3*points_per_period;

end
