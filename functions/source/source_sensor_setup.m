function [kgrid, source, sensor, source_labels] = source_sensor_setup(parameters, max_sound_speed, trans_pos_final, focus_pos_final, grid_time_step, min_sound_speed)
% SOURCE_SENSOR_SETUP  Create kWaveGrid, CW source, and full-grid sensor for a TUS simulation
%
% Builds the k-Wave computational grid (kWaveGrid) in 2-D or 3-D from
% parameters.grid.dims and resolution_mm. The time axis is chosen from
% the first transducer frequency: dt is derived from points-per-wavelength
% (PPW) and the CFL number (grid.source_cfl, default 0.3) so that the
% temporal sampling is an integer divisor of the wave period. A PPW
% check against grid.min_ppw (default 6) warns when the grid is too
% coarse. If grid_time_step is supplied (non-empty), that value is used
% directly (retry after instability). The sensor records p_max_all and
% p_final over the last 3 wave periods to capture steady state.
%
% Use as:
%   [kgrid, source, sensor, source_labels] = ...
%       source_sensor_setup(parameters, max_sound_speed, trans_pos_final, focus_pos_final)
%   [kgrid, source, sensor, source_labels] = ...
%       source_sensor_setup(parameters, max_sound_speed, trans_pos_final, focus_pos_final, grid_time_step)
%   [kgrid, source, sensor, source_labels] = ...
%       source_sensor_setup(parameters, max_sound_speed, trans_pos_final, focus_pos_final, grid_time_step, min_sound_speed)
%
% Input:
%   parameters      - PRESTUS config; must contain grid.dims, grid.resolution_mm [mm],
%                     grid.source_cfl, grid.min_ppw, grid.use_kWaveArray, io.dir_cache,
%                     and transducer(1).freq_hz [Hz]
%   max_sound_speed - maximum sound speed across all media [m/s]
%   trans_pos_final - transducer position in grid indices
%   focus_pos_final - focus position in grid indices
%   grid_time_step  - override dt [s] (optional, default: auto)
%   min_sound_speed - minimum sound speed for PPW check [m/s] (optional, default: max_sound_speed)
%
% Output:
%   kgrid         - kWaveGrid with time axis set
%   source        - struct with source.p_mask and source.p
%   sensor        - struct with sensor.mask, sensor.record, sensor.record_start_index
%   source_labels - grid of integer element labels (0 = inactive)
%
% See also: SOURCE_CREATE, GRID_TISSUE_SETUP, GRID_TRANSDUCER_LOCATION

arguments
    parameters      (1,1) struct
    max_sound_speed (1,1) {mustBeNumeric, mustBePositive}
    trans_pos_final (1,:) {mustBeNumeric}
    focus_pos_final (1,:) {mustBeNumeric}
    grid_time_step  {mustBeNumeric} = []
    min_sound_speed (1,1) {mustBeNumeric} = []
end

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
    % As a default the time step is based on the CFL number (grid.source_cfl,
    % default 0.15). If grid_time_step is given (non-empty), the function was
    % called a second time with a smaller time step due to simulation instability.
    if nargin < 5 || isempty(grid_time_step)
        % Calculate the time step using an integer number of points per period
        % PPW: Spatial samples per wavelength at source freq; ensures dx resolves waves (minimum grid.min_ppw, default 6).
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
        % Validate PPW against minimum requirement using the slowest medium
        % speed (shortest wavelength = worst-case spatial sampling).
        % max_sound_speed is intentionally NOT used here — it gives the
        % best-case PPW and would miss under-sampling in slow tissues.
        if isfield(parameters.grid, 'min_ppw') && ~isempty(parameters.grid.min_ppw)
            min_ppw = parameters.grid.min_ppw;
        else
            min_ppw = 6;
        end
        if nargin >= 6 && ~isempty(min_sound_speed)
            c_check = min_sound_speed;
        else
            c_check = max_sound_speed;  % fallback: conservative (overestimates PPW)
        end
        ppw_actual = c_check / (tx.freq_hz * parameters.grid.resolution_mm / 1e3);
        if ppw_actual < min_ppw
            dx_required_mm = c_check / (tx.freq_hz * min_ppw) * 1e3;
            prev_state = warning('off', 'backtrace');
            warn(['PPW %.1f < minimum %.0f at %.0f Hz with %.2f mm resolution ' ...
                     '(based on min. sound speed %.0f m/s). ' ...
                     'Reduce grid.resolution_mm to ≤ %.2f mm.'], ...
                ppw_actual, min_ppw, tx.freq_hz, parameters.grid.resolution_mm, ...
                c_check, dx_required_mm);
            warn(prev_state);
        end
        % Print parameter summary
        fprintf('Calculating the time step of %.1d based on PPW %.1f, CFL= %.2f and PPP %d.\n', ...
            grid_time_step, points_per_wavelength, cfl, points_per_period);
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
    % Always use output_affix so each simulation variant (default, liberal,
    % conservative) computes its own kwave_source.mat with the correct time
    % axis for that variant's medium speed.  preproc_affix is intentionally
    % NOT used here — it only applies to head-preprocessing file lookups in
    % preproc_head.m.
    source_affix = parameters.io.output_affix;
    parameters.io.filename_kwave_source  = fullfile(parameters.io.dir_cache, ...
        sprintf('sub-%03d_%s_kwave_source%s.mat', ...
        parameters.subject_id, parameters.simulation.medium, source_affix));
    if confirm_overwriting(parameters.io.filename_kwave_source, parameters)

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
            save(parameters.io.filename_kwave_source, 'source', 'source_labels','-v7.3');
        else
            disp('Not saving kwave source matrix ...')
        end
    else
        load(parameters.io.filename_kwave_source);
    end

    % Creates a sensor that records the the maximum and final pressure
    % values in every point of the grid
    sensor = struct();
    sensor.mask = ones(parameters.grid.dims);
    sensor.record = {'p_max_all','p_final','p'};

    % Record the last 3 cycles in steady state (when sonic waves have traversed the entire medium)
    num_periods = 3;
    time_points_to_record = round(num_periods * wave_period / kgrid.dt);
    sensor.record_start_index = simulation_time_points - time_points_to_record + 1;
    % alternative (EM): sensor.record_start_index = simulation_time_points - 3*points_per_period;

end
