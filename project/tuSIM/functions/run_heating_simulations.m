function [thermal_diff_obj, time_status_seq, maxT,focal_planeT] = run_heating_simulations(sensor_data, kgrid, kwave_medium, sensor, source, parameters, trans_pos)

% =========================================================================
% CALCULATE HEATING
% =========================================================================

% convert the absorption coefficient to nepers/m
alpha_np = db2neper(kwave_medium.alpha_coeff, kwave_medium.alpha_power) * ...
    (2 * pi * parameters.transducer.source_freq_hz).^kwave_medium.alpha_power;

% extract the pressure amplitude at each position
% p = extractAmpPhase(sensor_data.p_max_all, 1/kgrid.dt, parameters.transducer.source_freq_hz);
% 
% % reshape the data, and calculate the volume rate of heat deposition
% p = reshape(p, kgrid.Nx, kgrid.Ny);
p = sensor_data.p_max_all;
source.Q = alpha_np .* p.^2 ./ (kwave_medium.density .* kwave_medium.sound_speed);
source.T0 = parameters.thermal.temp_0;

% create kWaveDiffusion object
if isfield(parameters.thermal,'record_t_at_every_step') && ~parameters.thermal.record_t_at_every_step 
    sensor = [];
end
thermal_diff_obj = kWaveDiffusion(kgrid, kwave_medium, source, sensor, 'PlotSim', boolean(parameters.interactive));

thermal_diff_obj.T = gpuArray(thermal_diff_obj.T);
maxT = thermal_diff_obj.T;


% % 
% parameters.thermal.duty_cycle = 1; % share of the stimulation duration during which the stimulation is on
% parameters.thermal.sim_time_steps = 0.1; % [s] simulation time steps during the stimulation period
% 
% parameters.thermal.stim_duration = 0.6; % [s] stimulation duration within a trial
% 
% parameters.thermal.iti = 0.6; % interval between the trials, from the start of one trial to the start of another [s]
% parameters.thermal.n_trials = 12; % number of trials to simulate; the total simulated duration is then n_trials*iti seconds
% %parameters.thermal.nostim_period_time_steps = 0.1; % [s] simulation time steps during the period within a trial after the stimulation 
if ~isfield(parameters.thermal,'on_off_step_duration')
    if parameters.thermal.duty_cycle > 0 && parameters.thermal.duty_cycle < 1
        on_off_step_duration = parameters.thermal.sim_time_steps/min([parameters.thermal.duty_cycle 1-parameters.thermal.duty_cycle ]); % one on+off cycle duration
    else
        on_off_step_duration = parameters.thermal.sim_time_steps;
    end
else 
    on_off_step_duration = parameters.thermal.on_off_step_duration;
end
fprintf('Duration of 1 repetition of the on+off cycle: %.3f s\n', on_off_step_duration)

on_off_repetitions = parameters.thermal.stim_duration/on_off_step_duration;

on_off_repetitions = round_if_integer(on_off_repetitions, sprintf('The total stimulation duration (%.3f s) should be divisible by the on+off cycle duration (%.3f s)\n',parameters.thermal.stim_duration, on_off_step_duration));

on_steps_n = round_if_integer(on_off_step_duration*parameters.thermal.duty_cycle/parameters.thermal.sim_time_steps, 'The number of ''on'' steps within the on+off cycle should be integer');
off_steps_n = round_if_integer(on_off_step_duration*(1-parameters.thermal.duty_cycle)/parameters.thermal.sim_time_steps, 'The number of ''off'' steps within the on+off cycle should be integer');

if ~isfield(parameters.thermal,'post_stim_time_step_dur')
    parameters.thermal.post_stim_time_step_dur = parameters.thermal.sim_time_steps;
end

post_stim_period = parameters.thermal.iti - parameters.thermal.stim_duration;
post_stim_steps_n = post_stim_period / parameters.thermal.post_stim_time_step_dur;
post_stim_steps_n = round_if_integer(post_stim_steps_n, 'Number of simulation steps must be integer');

if  isfield(parameters.thermal,'equal_steps') && parameters.thermal.equal_steps == 0
    on_steps_dur = on_steps_n * parameters.thermal.sim_time_steps;
    off_steps_dur = off_steps_n * parameters.thermal.sim_time_steps;
    on_steps_n  = 1;
    off_steps_n = 1;
else
    on_steps_dur = parameters.thermal.sim_time_steps;
    off_steps_dur = parameters.thermal.sim_time_steps;
end

time_status_seq = struct('status', {'off'}, 'time', {0}, 'step', {0}, 'recorded', {1});

total_timepoints = 1+parameters.thermal.n_trials*(on_off_repetitions*(1+double(off_steps_n>0))+1);
cur_timepoint = 1;
focal_planeT = gpuArray(zeros([size(maxT,[1,3]) total_timepoints]));
focal_planeT(:,:,cur_timepoint) = squeeze(maxT(:,trans_pos(2),:));
% loop over trials
for trial_i = 1:parameters.thermal.n_trials
  fprintf('Trial %i\n', trial_i)
  % stimulation period
  for pulse_i = 1:on_off_repetitions
      fprintf('Pulse %i\n', pulse_i)
      % pulse on
      thermal_diff_obj.Q = source.Q;
      thermal_diff_obj.takeTimeStep(on_steps_n, on_steps_dur);
      time_status_seq = [time_status_seq struct('time', num2cell(max([time_status_seq(:).time]) + (1:on_steps_n)*on_steps_dur), ...
                                                'step', num2cell(max([time_status_seq(:).step]) + (1:on_steps_n)), 'status', "on",'recorded',0)];
      curT = thermal_diff_obj.T;
      cur_timepoint = cur_timepoint+1;
      time_status_seq(end).recorded = 1;
      focal_planeT(:,:,cur_timepoint) = cat(3, squeeze(curT(:,trans_pos(2),:)));
      if any(maxT>curT,'all')
          maxT = curT;
      end
      % pulse off 
      if off_steps_n>0
          thermal_diff_obj.Q = 0;
          thermal_diff_obj.takeTimeStep(off_steps_n, off_steps_dur);
          time_status_seq = [time_status_seq struct('time', num2cell(max([time_status_seq(:).time]) + (1:off_steps_n)*off_steps_dur), 'step', num2cell(max([time_status_seq(:).step]) + (1:off_steps_n)), 'status', "off",'recorded',0)];
          cur_timepoint = cur_timepoint+1;
          time_status_seq(end).recorded = 1;
          curT = thermal_diff_obj.T;

          focal_planeT(:,:,cur_timepoint) = cat(3, squeeze(curT(:,trans_pos(2),:)));

      end
  end
  
  % post-stimulation period
  if post_stim_steps_n > 0 
      thermal_diff_obj.Q = 0;
      thermal_diff_obj.takeTimeStep(post_stim_steps_n, parameters.thermal.post_stim_time_step_dur);
      time_status_seq = [time_status_seq struct('time', num2cell(max([time_status_seq(:).time]) + (1:post_stim_steps_n)*parameters.thermal.post_stim_time_step_dur), 'step', num2cell(max([time_status_seq(:).step]) + (1:post_stim_steps_n)), 'status', "off",'recorded',0)];
      time_status_seq(end).recorded = 1;
      cur_timepoint = cur_timepoint+1;
      curT = thermal_diff_obj.T;

      focal_planeT(:,:,cur_timepoint) = cat(3, squeeze(curT(:,trans_pos(2),:)));
      

  end
end


end