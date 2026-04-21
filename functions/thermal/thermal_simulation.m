function [thermal_diff_obj, time_status_seq, T_max, T_focal, T_cur, CEM43_max, CEM43_focal, CEM43_cur, timeseries, CEM43_iso_max, CEM43_iso_focal, CEM43_iso_cur] = ...
    thermal_simulation(...
    parameters, sensor_data, kgrid, kwave_medium, sensor, source, transf, medium_masks)
% THERMAL_SIMULATION  Simulate thermal effects of transcranial ultrasound using k-Wave
%
% Converts the acoustic pressure field (p_max_all from sensor_data) into a
% volumetric heat source Q via Q = alpha_np * p^2 / (rho * c). Runs
% kWaveDiffusion over the pulse protocol defined by thermal_parameters,
% looping through PT ON, PT OFF, PTRI-OFF, and post-PTRD cool-off phases.
% At each recorded step, temperature (T) and two CEM43 variants (k-Wave
% native; ISO 13 threshold) are accumulated into focal-plane snapshots,
% running maxima, and per-layer timeseries. GPU arrays are used when
% simulation.code_type is 'matlab_gpu' or 'cpp_gpu'.
%
% Use as:
%   [thermal_diff_obj, time_status_seq, T_max, T_focal, T_cur, ...
%    CEM43_max, CEM43_focal, CEM43_cur, timeseries, ...
%    CEM43_iso_max, CEM43_iso_focal, CEM43_iso_cur] = ...
%       thermal_simulation(parameters, sensor_data, kgrid, kwave_medium, sensor, source, transf, medium_masks)
%
% Input:
%   parameters   - PRESTUS config; must contain thermal, timing, grid.dims,
%                  transducer(1).trans_pos, transducer(1).freq_hz [Hz],
%                  pct.enabled, simulation.code_type, simulation.precision,
%                  simulation.interactive, io.adopted_cem43 (optional)
%   sensor_data  - struct from acoustic simulation; must contain p_max_all [Pa]
%   kgrid        - kWaveGrid object
%   kwave_medium - struct with fields: alpha_coeff, alpha_power,
%                  density [kg/m^3], sound_speed [m/s],
%                  thermal_conductivity [W/(m K)], specific_heat [J/(kg K)],
%                  perfusion_coeff [1/s], absorption_fraction, temp_0 [°C]
%   sensor       - struct with sensor.mask (used to define focal axis recording window)
%   source       - struct; Q is overwritten internally from p_max_all
%   transf       - affine forward transform from T1 to simulation grid
%                  (used only when io.adopted_cem43 is specified)
%   medium_masks - layer label map for per-tissue timeseries
%
% Output:
%   thermal_diff_obj  - kWaveDiffusion object after simulation
%   time_status_seq   - struct array with fields time [s], step, status, recorded
%   T_max             - maximum temperature over all time steps [°C]
%   T_focal           - temperature in focal plane [x, z, time] [°C]
%   T_cur             - temperature at the final recorded step [°C]
%   CEM43_max         - running max CEM43 (k-Wave method) [min]
%   CEM43_focal       - CEM43 in focal plane [x, z, time] [min]
%   CEM43_cur         - CEM43 at the final recorded step [min]
%   timeseries        - struct from THERMAL_UPDATE_TIMESERIES
%   CEM43_iso_max     - running max ISO CEM43 [min]
%   CEM43_iso_focal   - ISO CEM43 in focal plane [x, z, time] [min]
%   CEM43_iso_cur     - ISO CEM43 at the final recorded step [min]
%
% See also: THERMAL_PARAMETERS, THERMAL_PLOT_PROTOCOL, THERMAL_UPDATE_TIMESERIES,
%   THERMAL_ANALYSIS, MEDIUM_SETUP

arguments
    parameters   (1,1) struct
    sensor_data  (1,1) struct
    kgrid        (1,1)
    kwave_medium (1,1) struct
    sensor       (1,1) struct
    source       (1,1) struct
    transf
    medium_masks {mustBeNumeric}
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Notes:                                                            %
%                                                                   %
% Instead of modelling each oscillation for measuring heating       %
% effects, k-wave divides each duty cycle up into segments called   %
% 'sim_time_steps' with the total of one duty cycle being           %
% represented as the 'on_off_step_duration'.                        %
% The stable state of the acoustic simulations (what is seen in the %
% figures and nifti files) is used as the input for the temperature %
% simulations.                                                      %
%                                                                   %
% For a detailed explanation on how to correctly configure your     %
% thermal parameters, see doc_thermal_simulations.                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Clear the sensor mask
sensor.mask = zeros(size(sensor.mask));

% Transducer position in grid coordinates
tr = parameters.transducer(1);

% Place sensor along the focal axis 
heating_window_dims = ones(2,3);
for i = 1:2
    heating_window_dims(:,i) = [...
        max(1, -parameters.thermal.sensor_xy_halfsize + tr.trans_pos(i)), ...
        min(parameters.grid.dims(i), parameters.thermal.sensor_xy_halfsize + tr.trans_pos(i))];
end

% Define a simulation window
if length(size(parameters.grid.dims))>2
    heating_window_dims(2,3) = parameters.grid.dims(3);
end
sensor.mask(heating_window_dims(1,1):heating_window_dims(2,1), heating_window_dims(1,2):heating_window_dims(2,2), :) = 1;

% [(pseudo-)CT] if 'k-plan' mapping is used for density:
% skull bone density is fixed for thermal simulation
% (note that downstream metrics such as sound speed contribute to pressure, but are not input to kWaveDiffusion)
% https://dispatch.k-plan.io/static/docs/simulation-pipeline.html
if parameters.pct.enabled ==1 && ...
    (~isfield(parameters, 'pct') || ~isfield(parameters.pct, 'mapping_density') || strcmp(parameters.pct.mapping_density, 'k-plan'))
    skull_i = find(ismember(fieldnames(parameters.layers), {'skull', 'skull_cortical', 'skull_trabecular'}));
    % [k-Plan] use fixed density for thermal simulation
    kwave_medium.density(ismember(medium_masks,skull_i)) = 1850;
    clear skull_i;
end

% convert attenuation into absorption coefficients
kwave_medium.alpha_coeff = ...
    kwave_medium.alpha_coeff .* kwave_medium.absorption_fraction;

% convert the absorption coefficients to nepers/m (!)
% see also fitPowerLawParamsMulti.m
w = 2*pi*tr.freq_hz; % Frequency [rad/s]
a0_np = db2neper(kwave_medium.alpha_coeff, kwave_medium.alpha_power); % [Nepers/((rad/s)^y m)]
alpha_np = a0_np.*w.^kwave_medium.alpha_power;
clear w a0_np;

% alternative simplified conversion dB to Nepers
% alpha_np = (100 * kwave_medium.alpha_coeff .* (tr.freq_hz/10^6)^kwave_medium.alpha_power)/8.686;

% save absorption coefficient [Np] for debugging
if contains(parameters.simulation.medium, {'layered', 'phantom'}) && parameters.simulation.debug == 1
    try
        filename_absorption = fullfile(parameters.io.output_dir, 'debug', sprintf('matrix_absorption'));
        niftiwrite(alpha_np, filename_absorption, 'Compressed',true);
    catch
        warning("Error with saving absorption debug image...")
    end
end

% Get the maximum pressure (in Pa) and calculate Q, the volume rate of heat deposition
p = gather(abs(sensor_data.p_max_all));
source.Q = (alpha_np .* p.^2) ./ (kwave_medium.density .* kwave_medium.sound_speed); % Heat delivered to the system (W/m3)
source.T0 = kwave_medium.temp_0; %parameters.thermal.temp_0; % Initial temperature distribution
clear alpha_np p;

% remove unexpected fields from kWave_medium
kwave_medium = rmfield(kwave_medium, 'temp_0');
kwave_medium = rmfield(kwave_medium, 'absorption_fraction');

% if perfusion was specified, implement it by specifying blood temperature
if ~all(isnan(kwave_medium.perfusion_coeff))
    kwave_medium.blood_ambient_temperature = 37; % [degC]
else
    kwave_medium = rmfield(kwave_medium, 'perfusion_coeff');
end

% create kWaveDiffusion object
if isfield(parameters.thermal,'record_t_at_every_step') && ~parameters.thermal.record_t_at_every_step 
    sensor = [];
end

% Check if DataCast is supported
try
    use_datacast = any(regexp(fileread(which('kWaveDiffusion')), '''DataCast'''));
    if ~use_datacast
        warning('MATLAB GPU support requested but kWaveDiffusion lacks ''DataCast'' (kWave ≥1.4.1 required).');
    end
catch
    warning('Cannot verify kWaveDiffusion DataCast support; assuming GPU-compatible version.');
    use_datacast = true;
end

% Set precision and enable GPU mode (if requested)
if use_datacast == true
    if strcmp(parameters.simulation.code_type, 'matlab_gpu') || strcmp(parameters.simulation.code_type, 'cpp_gpu')
        datacast = ['gpuArray-', char(parameters.simulation.precision)];
    else
        datacast = char(parameters.simulation.precision);
    end
end

% Build final input args
thermal_args = {'PlotSim', boolean(parameters.simulation.interactive)};
if use_datacast
    thermal_args = [thermal_args, {'DataCast', datacast}];
end

% Run simulation
thermal_diff_obj = kWaveDiffusion(...
    kgrid, ...
    kwave_medium, ...
    source, ...
    sensor, ...
    thermal_args{:});

% initialize field temperature
if strcmp(parameters.simulation.code_type, 'matlab_gpu') || strcmp(parameters.simulation.code_type, 'cpp_gpu')
    thermal_diff_obj.T = gpuArray(thermal_diff_obj.T);
end
T_max = thermal_diff_obj.T;

% initialize field for cem43
if isfield(parameters.io, 'adopted_cem43') && ~isempty(parameters.io.adopted_cem43) && isfile(parameters.io.adopted_cem43)
    cumulative_heat_image = niftiread(parameters.io.adopted_cem43);
    fprintf('\nAdopting CEM43 heatmap %s from previous simulation\n', parameters.io.adopted_cem43)
    thermal_diff_obj.cem43 = double(tformarray(cumulative_heat_image, ...
        maketform("affine", transf), ...
        makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], size(medium_masks), [], 0));
else
    thermal_diff_obj.cem43 = zeros(size(thermal_diff_obj.T));
end
if strcmp(parameters.simulation.code_type, 'matlab_gpu') || strcmp(parameters.simulation.code_type, 'cpp_gpu')
    thermal_diff_obj.cem43 = gpuArray(thermal_diff_obj.cem43);
end
CEM43_max = thermal_diff_obj.cem43;

% initialize field for cem43 (iso variant)
tmp_obj.cem43_iso = zeros(size(thermal_diff_obj.T));
CEM43_iso_max = tmp_obj.cem43_iso;

% Convert the thermal parameters for the simulation
params_thermal = thermal_parameters(parameters);

% Visualize the thermal parameters
thermal_plot_protocol(params_thermal, parameters);

% Total timepoints estimation
% Records: 1 per ON, 1 per OFF (inside PT), 1 per PTRI-off block, 1 per post-PTRI step
n_ptri_reps = params_thermal.n_ptri_reps;
n_pulses_pt = params_thermal.n_pulses_per_pt;
snapshots_per_ptri = n_pulses_pt * 2 + (params_thermal.ptri_off_steps_n > 0);
total_timepoints = 1 + n_ptri_reps * snapshots_per_ptri + params_thermal.post_ptri_steps_n;

% Resize focal arrays
if ndims(T_max) == 3
    T_focal         = NaN([size(T_max,[1,3]) total_timepoints]);
    CEM43_focal     = NaN([size(T_max,[1,3]) total_timepoints]);
    CEM43_iso_focal = NaN([size(T_max,[1,3]) total_timepoints]);
    T_focal(:,:,1)         = squeeze(T_max(:,tr.trans_pos(2),:));
    CEM43_focal(:,:,1)     = squeeze(CEM43_max(:,tr.trans_pos(2),:));
    CEM43_iso_focal(:,:,1) = squeeze(CEM43_iso_max(:,tr.trans_pos(2),:));
elseif ndims(T_max) == 2
    T_focal         = NaN([size(T_max) total_timepoints]);
    CEM43_focal     = NaN([size(CEM43_max) total_timepoints]);
    CEM43_iso_focal = NaN([size(CEM43_iso_max) total_timepoints]);
    T_focal(:,:,1)         = T_max;
    CEM43_focal(:,:,1)     = CEM43_max;
    CEM43_iso_focal(:,:,1) = CEM43_iso_max;
end

% setup structure for layer-specific timerseries
timeseries = thermal_update_timeseries(parameters, medium_masks);

cur_timepoint = 1;
time_status_seq = struct('status', {'off'}, 'time', {0}, 'step', {0}, 'recorded', {1});

fprintf('Starting %d PTRI reps (PTRD=%.1fs), each with %d pulses (PTD=%.2fs)\n', ...
    n_ptri_reps, params_thermal.ptrd, n_pulses_pt, params_thermal.ptd);

% === Outer loop - PTRI (Pulse Train Repetitions) ===
for rep_i = 1:n_ptri_reps
    fprintf('Pulse Train Repetition %d/%d\n', rep_i, n_ptri_reps);
    
    % === 1. Inner loop: Pulse Trains (multiple pulses) ===
    for pulse_i = 1:n_pulses_pt
        fprintf('  Pulse train %d/%d\n', pulse_i, n_pulses_pt);
        
        % PULSE ON
        thermal_diff_obj.Q = source.Q;
        tmp_status = 'on';
        
        thermal_diff_obj.takeTimeStep(params_thermal.pt_on_steps_n, params_thermal.pt_on_steps_dur);
        
        % CEM43 iso update
        tmp_obj.cem43_iso = tmp_obj.cem43_iso + params_thermal.pt_on_steps_n*params_thermal.pt_on_steps_dur ./ 60 .* ...
            (0 .* (thermal_diff_obj.T < 39 & T_max < 43) + ...
             0.25 .* (thermal_diff_obj.T >= 39 & thermal_diff_obj.T < 43 & T_max < 43) + ...
             0.5 .* (thermal_diff_obj.T >= 43 | T_max >= 43)).^(43 - thermal_diff_obj.T);
        tmp_obj.cem43_iso(thermal_diff_obj.T >= 57 | T_max >= 57) = Inf;
        
        % Record status
        time_status_seq = [time_status_seq, struct(...
            'time', num2cell(max([time_status_seq(:).time]) + params_thermal.pt_on_steps_n*params_thermal.pt_on_steps_dur), ...
            'step', num2cell(max([time_status_seq(:).step]) + params_thermal.pt_on_steps_n), ...
            'status', tmp_status, 'recorded', 1)];
        
        % Update focal / max
        T_cur         = thermal_diff_obj.T;
        CEM43_cur     = thermal_diff_obj.cem43;
        CEM43_iso_cur = tmp_obj.cem43_iso;
        cur_timepoint = cur_timepoint + 1;
        if ndims(T_max) == 3
            T_focal(:,:,cur_timepoint)         = squeeze(T_cur(:,tr.trans_pos(2),:));
            CEM43_focal(:,:,cur_timepoint)     = squeeze(CEM43_cur(:,tr.trans_pos(2),:));
            CEM43_iso_focal(:,:,cur_timepoint) = squeeze(CEM43_iso_cur(:,tr.trans_pos(2),:));
        else
            T_focal(:,:,cur_timepoint)         = T_cur;
            CEM43_focal(:,:,cur_timepoint)     = CEM43_cur;
            CEM43_iso_focal(:,:,cur_timepoint) = CEM43_iso_cur;
        end
        T_max         = max(T_max,         T_cur);
        CEM43_max     = max(CEM43_max,     CEM43_cur);
        CEM43_iso_max = max(CEM43_iso_max, CEM43_iso_cur);
        timeseries = thermal_update_timeseries(parameters, medium_masks, timeseries, T_cur, CEM43_cur, CEM43_iso_cur);
        
        % PULSE OFF (within PT)
        if params_thermal.pt_off_steps_n > 0
            thermal_diff_obj.Q(:,:,:) = 0;
            thermal_diff_obj.takeTimeStep(params_thermal.pt_off_steps_n, params_thermal.pt_off_steps_dur);
            
            % CEM43 iso update (reuse your OFF logic, with new params)
            tmp_obj.cem43_iso = tmp_obj.cem43_iso + params_thermal.pt_off_steps_n*params_thermal.pt_off_steps_dur ./ 60 .* ...
                (0 .* (thermal_diff_obj.T < 39 & T_max < 43) + ...
                 0.25 .* (thermal_diff_obj.T >= 39 & thermal_diff_obj.T < 43 & T_max < 43) + ...
                 0.5 .* (thermal_diff_obj.T >= 43 | T_max >= 43)).^(43 - thermal_diff_obj.T);
            tmp_obj.cem43_iso(thermal_diff_obj.T >= 57 | T_max >= 57) = Inf;
            
            % Record status
            time_status_seq = [time_status_seq, struct(...
                'time', num2cell(max([time_status_seq(:).time]) + params_thermal.pt_off_steps_n*params_thermal.pt_off_steps_dur), ...
                'step', num2cell(max([time_status_seq(:).step]) + params_thermal.pt_off_steps_n), ...
                'status', 'off', 'recorded', 1)];
            
            % Update focal / max
            T_cur         = thermal_diff_obj.T;
            CEM43_cur     = thermal_diff_obj.cem43;
            CEM43_iso_cur = tmp_obj.cem43_iso;
            cur_timepoint = cur_timepoint + 1;
            if ndims(T_max) == 3
                T_focal(:,:,cur_timepoint)         = squeeze(T_cur(:,tr.trans_pos(2),:));
                CEM43_focal(:,:,cur_timepoint)     = squeeze(CEM43_cur(:,tr.trans_pos(2),:));
                CEM43_iso_focal(:,:,cur_timepoint) = squeeze(CEM43_iso_cur(:,tr.trans_pos(2),:));
            else
                T_focal(:,:,cur_timepoint)         = T_cur;
                CEM43_focal(:,:,cur_timepoint)     = CEM43_cur;
                CEM43_iso_focal(:,:,cur_timepoint) = CEM43_iso_cur;
            end
            T_max         = max(T_max,         T_cur);
            CEM43_max     = max(CEM43_max,     CEM43_cur);
            CEM43_iso_max = max(CEM43_iso_max, CEM43_iso_cur);
            timeseries = thermal_update_timeseries(parameters, medium_masks, timeseries, T_cur, CEM43_cur, CEM43_iso_cur);

        end  % end pulse OFF
    end  % end pulse loop
    
    % === 2. PTRI-OFF period (coarse timestep) ===
    if params_thermal.ptri_off_steps_n > 0
        thermal_diff_obj.Q(:,:,:) = 0;
        for i_step = 1:params_thermal.ptri_off_steps_n
            thermal_diff_obj.takeTimeStep(1, params_thermal.ptri_off_step_dur);
            
            % CEM43 iso
            tmp_obj.cem43_iso = tmp_obj.cem43_iso + 1*params_thermal.ptri_off_step_dur ./ 60 .* ...
                (0 .* (thermal_diff_obj.T < 39 & T_max < 43) + ...
                 0.25 .* (thermal_diff_obj.T >= 39 & thermal_diff_obj.T < 43 & T_max < 43) + ...
                 0.5 .* (thermal_diff_obj.T >= 43 | T_max >= 43)).^(43 - thermal_diff_obj.T);
            tmp_obj.cem43_iso(thermal_diff_obj.T >= 57 | T_max >= 57) = Inf;
            
            % Record time parameters
            time_status_seq = [time_status_seq struct(...
                'time', num2cell(max([time_status_seq(:).time]) + params_thermal.ptri_off_step_dur), ...
                'step', num2cell(max([time_status_seq(:).step]) + 1), ...
                'status', 'off', 'recorded', 1)];
            
            % Update output matrices
            T_cur         = thermal_diff_obj.T;
            CEM43_cur     = thermal_diff_obj.cem43;
            CEM43_iso_cur = tmp_obj.cem43_iso;
            cur_timepoint = cur_timepoint + 1;
            if ndims(T_max) == 3
                T_focal(:,:,cur_timepoint)         = squeeze(T_cur(:,tr.trans_pos(2),:));
                CEM43_focal(:,:,cur_timepoint)     = squeeze(CEM43_cur(:,tr.trans_pos(2),:));
                CEM43_iso_focal(:,:,cur_timepoint) = squeeze(CEM43_iso_cur(:,tr.trans_pos(2),:));
            else
                T_focal(:,:,cur_timepoint)         = T_cur;
                CEM43_focal(:,:,cur_timepoint)     = CEM43_cur;
                CEM43_iso_focal(:,:,cur_timepoint) = CEM43_iso_cur;
            end
            T_max         = max(T_max,         T_cur);
            CEM43_max     = max(CEM43_max,     CEM43_cur);
            CEM43_iso_max = max(CEM43_iso_max, CEM43_iso_cur);
            timeseries = thermal_update_timeseries(parameters, medium_masks, timeseries, T_cur, CEM43_cur, CEM43_iso_cur);

        end
    end
end  % end PTRI rep loop

% === 3. Post whole-PTRD cooloff (coarse timestep) ===
if params_thermal.post_ptri_steps_n > 0
    fprintf('Post-PTRD cooloff: %d steps\n', params_thermal.post_ptri_steps_n);
    thermal_diff_obj.Q(:,:,:) = 0;
    for i_step = 1:params_thermal.post_ptri_steps_n
        thermal_diff_obj.takeTimeStep(1, params_thermal.post_ptri_step_dur);
        
        % CEM43 iso update (identical to PTRI-off above)
        tmp_obj.cem43_iso = tmp_obj.cem43_iso + 1*params_thermal.post_ptri_step_dur ./ 60 .* ...
            (0 .* (thermal_diff_obj.T < 39 & T_max < 43) + ...
             0.25 .* (thermal_diff_obj.T >= 39 & thermal_diff_obj.T < 43 & T_max < 43) + ...
             0.5 .* (thermal_diff_obj.T >= 43 | T_max >= 43)).^(43 - thermal_diff_obj.T);
        tmp_obj.cem43_iso(thermal_diff_obj.T >= 57 | T_max >= 57) = Inf;
        
        % Record status (identical)
        time_status_seq = [time_status_seq struct(...
            'time', num2cell(max([time_status_seq(:).time]) + params_thermal.post_ptri_step_dur), ...
            'step', num2cell(max([time_status_seq(:).step]) + 1), ...
            'status', 'off', 'recorded', 1)];
        
        % Update focal / max
        T_cur         = thermal_diff_obj.T;
        CEM43_cur     = thermal_diff_obj.cem43;
        CEM43_iso_cur = tmp_obj.cem43_iso;
        cur_timepoint = cur_timepoint + 1;
        if ndims(T_max) == 3
            T_focal(:,:,cur_timepoint)         = squeeze(T_cur(:,tr.trans_pos(2),:));
            CEM43_focal(:,:,cur_timepoint)     = squeeze(CEM43_cur(:,tr.trans_pos(2),:));
            CEM43_iso_focal(:,:,cur_timepoint) = squeeze(CEM43_iso_cur(:,tr.trans_pos(2),:));
        else
            T_focal(:,:,cur_timepoint)         = T_cur;
            CEM43_focal(:,:,cur_timepoint)     = CEM43_cur;
            CEM43_iso_focal(:,:,cur_timepoint) = CEM43_iso_cur;
        end
        T_max         = max(T_max,         T_cur);
        CEM43_max     = max(CEM43_max,     CEM43_cur);
        CEM43_iso_max = max(CEM43_iso_max, CEM43_iso_cur);
        timeseries = thermal_update_timeseries(parameters, medium_masks, timeseries, T_cur, CEM43_cur, CEM43_iso_cur);

    end
end

% Trim unused timepoints (if too many were originally required)
T_focal         = T_focal(:,:,1:cur_timepoint);
CEM43_focal     = CEM43_focal(:,:,1:cur_timepoint);
CEM43_iso_focal = CEM43_iso_focal(:,:,1:cur_timepoint);

% Apply gather (if variables are GPU arrays)
T_max         = gather(T_max);
T_focal       = gather(T_focal);
T_cur         = gather(T_cur);
CEM43_max     = gather(CEM43_max);
CEM43_focal   = gather(CEM43_focal);
CEM43_cur     = gather(CEM43_cur);
CEM43_iso_max   = gather(CEM43_iso_max);
CEM43_iso_focal = gather(CEM43_iso_focal);
CEM43_iso_cur   = gather(CEM43_iso_cur);

fprintf('Thermal simulation complete. Recorded %d timepoints.\n', cur_timepoint);

end