function [thermal_diff_obj, time_status_seq, T_max, T_focal, CEM43_max, CEM43_focal] = ...
    run_heating_simulations(sensor_data, kgrid, kwave_medium, sensor, source, parameters, trans_pos, final_transformation_matrix, medium_masks)

% RUN_HEATING_SIMULATIONS Simulates thermal effects of ultrasound using k-Wave.
%
% This function runs heating simulations based on acoustic simulation results 
% from k-Wave. It models the temperature rise and thermal dose (CEM43) caused 
% by ultrasound transducers over multiple trials and pulses.
%
% Input:
%   sensor_data  - Struct containing acoustic simulation results (e.g., pressure fields).
%   kgrid        - Struct representing the k-Wave computational grid (e.g., `kWaveGrid`).
%   kwave_medium - Struct containing medium properties (e.g., density, sound speed).
%   sensor       - Struct defining the sensor mask for thermal simulations.
%   source       - Struct defining the heat deposition source for thermal simulations.
%   parameters   - Struct containing simulation parameters (e.g., thermal properties).
%   trans_pos    - [1x3] array specifying the transducer position in grid coordinates.
%
% Output:
%   thermal_diff_obj - kWaveDiffusion object used for thermal simulations.
%   time_status_seq  - Struct array tracking the status of each time step (e.g., on/off).
%   T_max            - Maximum temperature reached during the simulation.
%   T_focal          - Temperature at the focal plane over time. [x,z,time]
%   CEM43_max        - Maximum cumulative equivalent minutes at 43°C (CEM43) during the simulation.
%   CEM43_focal      - CEM43 values at the focal plane over time. [x,z,time]

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

% convert the absorption coefficients to nepers/m (!)
% see also fitPowerLawParamsMulti.m
w = 2*pi*parameters.transducer.source_freq_hz; % Frequency [rad/s]
a0_np = db2neper(kwave_medium.alpha_coeff, kwave_medium.alpha_power); % [Nepers/((rad/s)^y m)]
alpha_np = a0_np.*w.^kwave_medium.alpha_power;
clear w a0_np;

% alternative simplified conversion dB to Nepers
% alpha_np = (100 * kwave_medium.alpha_coeff .* (parameters.transducer.source_freq_hz/10^6)^kwave_medium.alpha_power)/8.686;

% save absorption coefficient [Np] for debugging
if (contains(parameters.simulation_medium, 'skull') || ...
        contains(parameters.simulation_medium, 'layered') || ...
        contains(parameters.simulation_medium, 'phantom'))
    if ~exist(fullfile(parameters.output_dir, 'debug')); mkdir(parameters.output_dir, 'debug'); end
    try
        filename_absorption = fullfile(parameters.output_dir, 'debug', sprintf('matrix_absorption'));
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

if strcmp(parameters.code_type, 'matlab_gpu') || strcmp(parameters.code_type, 'cuda')
    datacast = 'gpuArray-double';
else
    datacast = 'off';
end    

thermal_diff_obj = kWaveDiffusion(...
    kgrid, ...
    kwave_medium, ...
    source, ...
    sensor, ...
    'PlotSim', boolean(parameters.interactive), ...
    'DataCast', datacast);

% initialize field temperature
if strcmp(parameters.code_type, 'matlab_gpu') || strcmp(parameters.code_type, 'cuda')
    thermal_diff_obj.T = gpuArray(thermal_diff_obj.T);
end    
T_max = thermal_diff_obj.T;

% initialize field for cem43
if isfield(parameters, 'adopted_cumulative_heat')
    cumulative_heat_image = niftiread(parameters.adopted_cumulative_heat);
    thermal_diff_obj.cem43 = double(tformarray(cumulative_heat_image, ...
        maketform("affine", final_transformation_matrix), ...
        makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], size(medium_masks), [], 0));
else
    thermal_diff_obj.cem43 = zeros(size(thermal_diff_obj.T));
end
if strcmp(parameters.code_type, 'matlab_gpu') || strcmp(parameters.code_type, 'cuda')
    thermal_diff_obj.cem43 = gpuArray(thermal_diff_obj.cem43);
end
CEM43_max = thermal_diff_obj.cem43;

% initialize field for cem43 (iso variant)
tmp_obj.cem43_iso = zeros(size(thermal_diff_obj.T));

[on_steps_n,  on_steps_dur, off_steps_n, off_steps_dur, ...
    post_stim_step_n, post_stim_step_dur] = ...
    check_thermal_parameters(parameters);

time_status_seq = struct('status', {'off'}, 'time', {0}, 'step', {0}, 'recorded', {1});

% Set final parameters for simulation
total_timepoints = 1+parameters.thermal.n_trials*(double(on_steps_n>0)+1);
cur_timepoint = 1;

% Set results for the focal axis (3D; for 2D assume that focal axis is provided)
if ndims(T_max) == 3
    T_focal = NaN([size(T_max,[1,3]) total_timepoints]);
    T_focal(:,:,cur_timepoint) = squeeze(T_max(:,trans_pos(2),:));
    CEM43_focal = NaN([size(T_max,[1,3]),total_timepoints]);
    CEM43_focal(:,:,cur_timepoint) = squeeze(CEM43_max(:,trans_pos(2),:));
elseif ndims(T_max) == 2
    T_focal = NaN([size(T_max,[1,2]) total_timepoints]);
    T_focal(:,:,cur_timepoint) = squeeze(T_max);
    CEM43_focal = NaN([size(T_max,[1,2]),total_timepoints]);
    CEM43_focal(:,:,cur_timepoint) = squeeze(CEM43_max);
end

% Loop over trials, so that the heat can accumulate
for trial_i = 1:parameters.thermal.n_trials
    fprintf('Trial %i\n', trial_i)
    
    is_break_period = 0;
    if isfield(parameters, 'start_break_trials')
        for i_break = 1:length(parameters.start_break_trials)
            if trial_i >= parameters.start_break_trials(i_break) && trial_i <= parameters.stop_break_trials(i_break)
                is_break_period = 1;
            end
        end
    end
    
    % Calculate the heat accumulation within a trial
    for pulse_i = 1:1
        fprintf('Pulse %i\n', pulse_i)
        
        if is_break_period == 0
          thermal_diff_obj.Q = source.Q;
          tmp_status = 'on';
        else
          thermal_diff_obj.Q(:,:,:) = 0;
          tmp_status = 'off';
          disp('Break active');
        end
        
        thermal_diff_obj.takeTimeStep(on_steps_n, on_steps_dur);
        % calculate iso-variant of CEM43 (not updated within kWaveDiffusion)
        tmp_obj.cem43_iso = tmp_obj.cem43_iso + on_steps_n*on_steps_dur ./ 60 .* ...
            (0 .* (thermal_diff_obj.T < 39 & T_max < 43) + ...
            0.25 .* (thermal_diff_obj.T >= 39 & thermal_diff_obj.T < 43 & T_max < 43) + ...
            0.5 .* (thermal_diff_obj.T >= 43 | T_max >= 43)).^(43 - thermal_diff_obj.T);
        tmp_obj.cem43_iso(thermal_diff_obj.T >= 57 | T_max >= 57) = Inf;

        time_status_seq = [time_status_seq, ...
          struct('time', num2cell(max([time_status_seq(:).time]) + (on_steps_n)*on_steps_dur), ...
                 'step', num2cell(max([time_status_seq(:).step]) + (on_steps_n)), ...
                 'status', tmp_status,...
                 'recorded',1)];
        curT = thermal_diff_obj.T;
        if isfield(parameters.thermal, 'cem43_iso') && parameters.thermal.cem43_iso == 1
            curCEM43 = tmp_obj.cem43_iso;
        else
            curCEM43 = thermal_diff_obj.cem43;
        end
        cur_timepoint = cur_timepoint+1;
        
        if ndims(T_max) == 3
            T_focal(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
            CEM43_focal(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));
        elseif ndims(T_max) == 2
            T_focal(:,:,cur_timepoint) = squeeze(curT);
            CEM43_focal(:,:,cur_timepoint) = squeeze(curCEM43);
        end
        
        % update T_max and CEM43_max where applicable
        T_max = max(T_max, curT);
        CEM43_max = max(CEM43_max, curCEM43);
        
        % Pulse off 
        if off_steps_n>0
            thermal_diff_obj.Q(:,:,:) = 0;
            thermal_diff_obj.takeTimeStep(off_steps_n, off_steps_dur);
            % calculate iso-variant of CEM43 (not updated within kWaveDiffusion)
            tmp_obj.cem43_iso = tmp_obj.cem43_iso + off_steps_n*off_steps_dur ./ 60 .* ...
                (0 .* (thermal_diff_obj.T < 39 & T_max < 43) + ...
                0.25 .* (thermal_diff_obj.T >= 39 & thermal_diff_obj.T < 43 & T_max < 43) + ...
                0.5 .* (thermal_diff_obj.T >= 43 | T_max >= 43)).^(43 - thermal_diff_obj.T);
            tmp_obj.cem43_iso(thermal_diff_obj.T >= 57 | T_max >= 57) = Inf;
            time_status_seq = [time_status_seq,...
              struct('time', num2cell(max([time_status_seq(:).time]) + (off_steps_n)*off_steps_dur), ...
              'step', num2cell(max([time_status_seq(:).step]) + (off_steps_n)),...
              'status', "off",...
              'recorded',1)];
            cur_timepoint = cur_timepoint+1;
            curT = thermal_diff_obj.T;
            if isfield(parameters.thermal, 'cem43_iso') && parameters.thermal.cem43_iso == 1
                curCEM43 = tmp_obj.cem43_iso;
            else
                curCEM43 = thermal_diff_obj.cem43;
            end
            if ndims(T_max) == 3
                T_focal(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
                CEM43_focal(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));
            elseif ndims(T_max) == 2
                T_focal(:,:,cur_timepoint) = squeeze(curT);
                CEM43_focal(:,:,cur_timepoint) = squeeze(curCEM43);
            end
        end
    end
end

% Post-stimulation period
if post_stim_step_n > 0 
    thermal_diff_obj.Q(:,:,:) = 0;
    thermal_diff_obj.takeTimeStep(post_stim_step_n, post_stim_step_dur);
    % calculate iso-variant of CEM43 (not updated within kWaveDiffusion)
    tmp_obj.cem43_iso = tmp_obj.cem43_iso + post_stim_step_n*post_stim_step_dur ./ 60 .* ...
        (0 .* (thermal_diff_obj.T < 39 & T_max < 43) + ...
        0.25 .* (thermal_diff_obj.T >= 39 & thermal_diff_obj.T < 43 & T_max < 43) + ...
        0.5 .* (thermal_diff_obj.T >= 43 | T_max >= 43)).^(43 - thermal_diff_obj.T);
    tmp_obj.cem43_iso(thermal_diff_obj.T >= 57 | T_max >= 57) = Inf;
    time_status_seq = [time_status_seq struct(...
      'time', num2cell(max([time_status_seq(:).time]) + (post_stim_step_n)*post_stim_step_dur), ...
      'step', num2cell(max([time_status_seq(:).step]) + (post_stim_step_n)), ...
      'status', "off",...
      'recorded',1)];
    cur_timepoint = cur_timepoint+1;
    curT = thermal_diff_obj.T;
    if isfield(parameters.thermal, 'cem43_iso') && parameters.thermal.cem43_iso == 1
        curCEM43 = tmp_obj.cem43_iso;
    else
        curCEM43 = thermal_diff_obj.cem43;
    end
    if ndims(T_max) == 3
        T_focal(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
        CEM43_focal(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));
    elseif ndims(T_max) == 2
        T_focal(:,:,cur_timepoint) = squeeze(curT);
        CEM43_focal(:,:,cur_timepoint) = squeeze(curCEM43);
    end
end

end