function [thermal_diff_obj, time_status_seq, maxT,focal_planeT, maxCEM43, CEM43] = run_heating_simulations(sensor_data, kgrid, kwave_medium, sensor, source, parameters, trans_pos)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                    Runs the heating simulations                   %
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
% thermal parameters, see thermal simulations getting started.      %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Convert the absorption coefficients to nepers/m (!)
% for the following, see also fitPowerLawParamsMulti.m
% define frequency in rad/s
w = 2*pi*parameters.transducer.source_freq_hz;
% convert absorption to Nepers/((rad/s)^y m)
a0_np = db2neper(kwave_medium.alpha_coeff, kwave_medium.alpha_power);
alpha_np = a0_np.*w.^kwave_medium.alpha_power;
clear w a0_np;

% alternative simplified conversion dB to Nepers
% alpha_np = (100 * kwave_medium.alpha_coeff .* (parameters.transducer.source_freq_hz/10^6)^kwave_medium.alpha_power)/8.686;

% Get the maximum pressure (in Pa) and calculate Q, the volume rate of heat deposition
p = gather(abs(sensor_data.p_max_all));
source.Q = (alpha_np .* p.^2) ./ (kwave_medium.density .* kwave_medium.sound_speed); % Heat delivered to the system (W/m3)
source.T0 = kwave_medium.temp_0; %parameters.thermal.temp_0; % Initial temperature distribution
clear alpha_np p;

% split temp_0 off from kWave_medium
kwave_medium = rmfield(kwave_medium, 'temp_0');

% create kWaveDiffusion object
if isfield(parameters.thermal,'record_t_at_every_step') && ~parameters.thermal.record_t_at_every_step 
    sensor = [];
end
thermal_diff_obj = kWaveDiffusion(kgrid, kwave_medium, source, sensor, 'PlotSim', boolean(parameters.interactive));

% New field temperature
thermal_diff_obj.T = gpuArray(thermal_diff_obj.T);
maxT = thermal_diff_obj.T;

% new field for cem43 (replicate above)
thermal_diff_obj.cem43 = gpuArray(zeros(size(thermal_diff_obj.T)));
maxCEM43 = thermal_diff_obj.cem43;

[~, on_off_repetitions, on_steps_n,  on_steps_dur, off_steps_n, off_steps_dur, post_stim_steps_n, post_stim_time_step_dur] = ...
    check_thermal_parameters(parameters);

time_status_seq = struct('status', {'off'}, 'time', {0}, 'step', {0}, 'recorded', {1});

% Set final parameters for simulation
total_timepoints = 1+parameters.thermal.n_trials*(on_off_repetitions*(1+double(off_steps_n>0))+1);
cur_timepoint = 1;
focal_planeT = gpuArray(NaN([size(maxT,[1,3]) total_timepoints]));
focal_planeT(:,:,cur_timepoint) = squeeze(maxT(:,trans_pos(2),:));

CEM43 = gpuArray(NaN([size(maxT,[1,3]),total_timepoints]));
CEM43(:,:,cur_timepoint) = squeeze(maxCEM43(:,trans_pos(2),:));

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
  for pulse_i = 1:on_off_repetitions
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
      
      time_status_seq = [time_status_seq, ...
          struct('time', num2cell(max([time_status_seq(:).time]) + (1:on_steps_n)*on_steps_dur), ...
                 'step', num2cell(max([time_status_seq(:).step]) + (1:on_steps_n)), ...
                 'status', tmp_status,...
                 'recorded',1)];
      curT = thermal_diff_obj.T;
      curCEM43 = thermal_diff_obj.cem43;
      cur_timepoint = cur_timepoint+1;
      focal_planeT(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
      CEM43(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));

      % update maxT and maxCEM43 where applicable
      maxT = max(maxT, curT);
      maxCEM43 = max(maxCEM43, curCEM43);

      % Pulse off 
      if off_steps_n>0
          thermal_diff_obj.Q(:,:,:) = 0;
          thermal_diff_obj.takeTimeStep(off_steps_n, off_steps_dur);
          time_status_seq = [time_status_seq,...
              struct('time', num2cell(max([time_status_seq(:).time]) + (1:off_steps_n)*off_steps_dur), ...
              'step', num2cell(max([time_status_seq(:).step]) + (1:off_steps_n)),...
              'status', "off",...
              'recorded',1)];
          cur_timepoint = cur_timepoint+1;
          curT = thermal_diff_obj.T;
          curCEM43 = thermal_diff_obj.cem43;
          focal_planeT(:,:,cur_timepoint) = ...
              squeeze(curT(:,trans_pos(2),:));
          CEM43(:,:,cur_timepoint) = ...
              squeeze(curCEM43(:,trans_pos(2),:));
      end
  end
end

% Post-stimulation period
if post_stim_steps_n > 0 
  thermal_diff_obj.Q(:,:,:) = 0;
  thermal_diff_obj.takeTimeStep(post_stim_steps_n, post_stim_time_step_dur);
  time_status_seq = [time_status_seq struct(...
      'time', num2cell(max([time_status_seq(:).time]) + (1:off_steps_n)*off_steps_dur), ...
      'step', num2cell(max([time_status_seq(:).step]) + (1:off_steps_n)), ...
      'status', "off",...
      'recorded',1)];
  cur_timepoint = cur_timepoint+1;
  curT = thermal_diff_obj.T;
  curCEM43 = thermal_diff_obj.cem43;
  focal_planeT(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
  CEM43(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));
end

end