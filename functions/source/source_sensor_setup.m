function [kgrid, source, sensor, source_labels] = source_sensor_setup(parameters, max_sound_speed, trans_pos_final, focus_pos_final, grid_time_step, min_sound_speed)
% SOURCE_SENSOR_SETUP  Create kWaveGrid, CW source, and full-grid sensor for a TUS simulation
%
% Builds the k-Wave computational grid (kWaveGrid) in 2-D or 3-D from
% parameters.grid.dims and resolution_mm.
%
% --- Time axis ---
%
% dt  is derived from two quantities:
%   PPW  = max_sound_speed / (freq_hz * dx)   [spatial samples per wavelength]
%   CFL  = grid.source_cfl  (default 0.15)    [Courant-Friedrichs-Lewy number]
%   PPP  = ceil(PPW / CFL)                    [temporal samples per period]
%   dt   = wave_period / PPP                   [time step, s]
%
% max_sound_speed sets PPW (and therefore dt): using the fastest medium
% gives the smallest dt, which is numerically stable everywhere.
% Override dt by passing a non-empty grid_time_step (used on a stability
% retry after checkStability flags instability).
%
% A PPW check is run against grid.min_ppw (default 6) using min_sound_speed
% (shortest wavelength, worst-case spatial sampling). A warning is issued
% when the grid is too coarse.
%
% t_end is the simulation duration:
%   t_end = grid_diagonal / min_sound_speed
%
% min_sound_speed (the slowest medium, typically water at ~1500 m/s) is
% used here — NOT max_sound_speed. Using max_sound_speed would shorten the
% simulation when bone is present (e.g. skull at 2800 m/s cuts t_end by
% ~1.9×), causing the wavefront traveling through slow tissue to not reach
% steady state at and beyond the focal peak, producing an artefactual
% post-peak pressure drop.
%
% The sensor records p_max_all and p_final over the last 3 wave periods
% to capture steady state.
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
%   min_sound_speed - minimum sound speed across all media [m/s] (optional, default: max_sound_speed).
%                     Used for (1) the PPW spatial-sampling check and (2) computing t_end so
%                     the simulation runs long enough for waves in slow media to reach steady state.
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
    min_sound_speed {mustBeNumeric} = []
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

    % The time step (dt) is constrained by the highest frequency present: a shorter
    % wave period requires finer temporal sampling to satisfy the CFL condition and
    % avoid staircasing the waveform. Using the lowest frequency would under-sample
    % higher-frequency transducers. The PPW check is also run against the highest
    % frequency (shortest wavelength = worst-case spatial sampling).
    freqs_all  = [parameters.transducer.freq_hz];
    freq_hz_max = max(freqs_all);
    tx = parameters.transducer(freqs_all == freq_hz_max);
    tx = tx(1);   % pick first if multiple transducers share the max frequency
    if numel(unique(freqs_all)) > 1
        fprintf(['Multi-frequency setup: grid time step and PPW sized for highest frequency ' ...
                 '%g Hz (transducer frequencies: %s Hz).\n'], freq_hz_max, num2str(freqs_all));
    end
    wave_period = 1 / tx.freq_hz;    % period [s] — based on highest frequency
    
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
        % Courant-Friedrichs-Lewy: Fraction of dx/c for dt; matches config_default.yaml.
        if ~isfield(parameters.grid, 'source_cfl') || isempty(parameters.grid.source_cfl)
            % jwave time-domain solvers (PSTD, FD) use a Stokes (y=2) absorber whose
            % stability margin depends on alpha_coeff.  High-attenuation skull tissue
            % (>~15 dB/cm/MHz) tightens the effective stability limit beyond the pure
            % acoustic CFL bound of 1/(pi*sqrt(ndim)).  A conservative default of 0.10
            % keeps the absorption metric well within the stable regime; users can
            % override with grid.source_cfl if the medium is low-absorption.
            is_jwave_td = isfield(parameters, 'acoustic_backend') && ...
                startsWith(lower(parameters.acoustic_backend), 'jwave') && ...
                (~isfield(parameters, 'jwave') || ~isfield(parameters.jwave, 'solver') || ...
                 ~strcmpi(parameters.jwave.solver, 'helmholtz'));
            if is_jwave_td
                cfl = 0.10;
                fprintf(['[source_sensor_setup] jwave time-domain backend detected: ' ...
                         'using conservative CFL=%.2f (override with grid.source_cfl).\n'], cfl);
            else
                cfl = 0.15;
            end
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
            msg = sprintf(['PPW %.1f < minimum %.0f at %.0f Hz with %.2f mm resolution ' ...
                     '(based on min. sound speed %.0f m/s). ' ...
                     'Reduce grid.resolution_mm to \x2264 %.2f mm.'], ...
                ppw_actual, min_ppw, tx.freq_hz, parameters.grid.resolution_mm, ...
                c_check, dx_required_mm);
            % jwave's Stokes absorption amplifies aliased high-k modes at PPW < 6,
            % causing divergence (NaN). Treat as hard error for jwave backends.
            is_jwave = isfield(parameters, 'acoustic_backend') && ...
                       startsWith(parameters.acoustic_backend, 'jwave');
            if is_jwave
                error('source_sensor_setup:ppwTooLow', '%s\n%s', msg, ...
                    'jwave PSTD/FD with absorption diverges (NaN) at PPW < min_ppw. Reduce resolution.');
            else
                prev_state = warning('off', 'backtrace');
                warn(msg);
                warn(prev_state);
            end
        end
        % Print parameter summary
        fprintf('Calculating the time step of %.1d based on PPW %.1f, CFL= %.2f and PPP %d.\n', ...
            grid_time_step, points_per_wavelength, cfl, points_per_period);
    else
        fprintf('Using the time step of %.1d based on approx. max. stable time step.\n', ...
            grid_time_step);
    end

    % Calculate the number of time steps to reach steady state.
    % Use min_sound_speed so that the simulation runs long enough for the
    % wave to traverse the full grid through the slowest medium (e.g. water
    % at 1500 m/s).  Using max_sound_speed underestimates t_end when bone is
    % present, cutting the simulation short before the steady-state field
    % develops beyond the focal peak.
    if nargin >= 6 && ~isempty(min_sound_speed)
        t_end_speed = min_sound_speed;
    else
        t_end_speed = max_sound_speed;
    end
    t_end = sqrt(kgrid.x_size.^2 + kgrid.z_size.^2 + kgrid.y_size.^2) / t_end_speed;    % [s]

    % jwave steady-state safety margin.
    % jwave's PML sits INSIDE the declared grid (k-Wave's sits outside), so
    % run_jwave_sim.py's auto_pad_for_pml silently EXPANDS the domain actually
    % propagated over whenever the transducer source lands within pml_size
    % voxels of the FOV boundary -- which can still happen even after
    % acoustic_fov_crop.m's own PML pre-expansion (observed empirically: a
    % source landed exactly at the pre-expanded boundary, 1 voxel from the
    % edge, triggering a further ~pml_size-voxel Python-side pad). t_end above
    % has no visibility into that later, Python-side expansion, so without
    % this margin the "last num_periods" recording window can start BEFORE
    % the wavefront has even crossed the domain that is actually simulated
    % (confirmed via run_jwave_sim.py --debug: margin of -5.23 periods, i.e.
    % negative, on a benchmark run with this exact source-near-PML geometry).
    % Add worst-case margin (pml on both ends of one axis) so recording start
    % stays safely past the true one-way transit regardless of whether the
    % Python-side pad actually triggers.
    if isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'jwave')
        if isfield(parameters.grid, 'pml_size') && isnumeric(parameters.grid.pml_size)
            pml_voxels = max(parameters.grid.pml_size(:));
        else
            pml_voxels = 40;   % worst case: resolve_jwave_pml_size.m clips 'auto' to [10, 40]
        end
        pml_margin_s = 2 * pml_voxels * parameters.grid.resolution_mm * 1e-3 / t_end_speed;
        fprintf('[source_sensor_setup] jwave steady-state safety margin: +%.2f periods (t_end +%.2f us) for possible PML auto-padding\n', ...
            pml_margin_s / wave_period, pml_margin_s * 1e6);
        t_end = t_end + pml_margin_s;
    end

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

        % build per-transducer positions for geometry (both kWaveArray and non-kWaveArray paths)
        if numel(parameters.transducer) > 1
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
                if isfield(tr, 'focus_pos') && ~isempty(tr.focus_pos)
                    fpos = tr.focus_pos;
                    if size(fpos,1) > size(fpos,2)
                        fpos = fpos';
                    end
                    focus_pos(ti,:) = fpos;
                else
                    focus_pos(ti,:) = focus_pos_final;
                end
            end
        else
            % single-transducer case
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
    sensor.record = {'p_max_all','p_min_all','p_final'};

    % Record the last 3 cycles in steady state (when sonic waves have traversed the entire medium)
    num_periods = 3;
    time_points_to_record = round(num_periods * wave_period / kgrid.dt);
    sensor.record_start_index = simulation_time_points - time_points_to_record + 1;
    % alternative (EM): sensor.record_start_index = simulation_time_points - 3*points_per_period;

    fprintf('[source_sensor_setup] ---- simulation grid summary ----\n');
    fprintf('[source_sensor_setup] Grid dims         : %s\n', num2str(parameters.grid.dims));
    fprintf('[source_sensor_setup] Resolution        : %.4f mm\n', parameters.grid.resolution_mm);
    fprintf('[source_sensor_setup] Frequency         : %.2f kHz\n', tx.freq_hz / 1e3);
    if ~isempty(min_sound_speed)
        ppw_report = min_sound_speed / (tx.freq_hz * parameters.grid.resolution_mm / 1e3);
    else
        ppw_report = max_sound_speed / (tx.freq_hz * parameters.grid.resolution_mm / 1e3);
    end
    fprintf('[source_sensor_setup] PPW (min speed)   : %.2f\n', ppw_report);
    fprintf('[source_sensor_setup] dt                : %.4e s\n', kgrid.dt);
    fprintf('[source_sensor_setup] Nt                : %d\n', simulation_time_points);
    fprintf('[source_sensor_setup] t_end             : %.4f ms\n', simulation_time_points * kgrid.dt * 1e3);
    fprintf('[source_sensor_setup] Record start idx  : %d  (steady-state window: last %d periods)\n', sensor.record_start_index, num_periods);
    fprintf('[source_sensor_setup] ------------------------------------\n');

end
