function [thermal_diff_obj, time_status_seq, maxT,focal_planeT] = run_heating_simulations(sensor_data, kgrid, kwave_medium, sensor, source, parameters, trans_pos)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                    Runs the heating simulations                   %
%                                                                   %
% Instead of modelling each oscillation for measuring heating       %
% effects, k-wave divides each duty cycle up into segments called   %
% 'sim_time_steps' with the total of one duty cycle being           %
% represented as the 'on_off_step_duration'.                        %
% The stable state of the acoustic simulations (what is seen in the %
% figures) is used as the input for the temperature simulations.    %
%                                                                   %
% Some notes:                                                       %
% As can be seen in the explanations in lines 51 to 53, these two   %
% values can be any value between 0 and 1 but they have to meet the %
% following constraints:                                            %
% on_off_step_duration * duty_cycle / stim_time_steps = integer     %
% on_off_step_duration * (1-duty_cycle) / stim_time_steps = integer %
% cycle_duration / onf_off_step_duration = integer                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Convert the absorption coefficient to nepers/m
alpha_np = db2neper(kwave_medium.alpha_coeff, kwave_medium.alpha_power) * ...
    (2 * pi * parameters.transducer.source_freq_hz).^kwave_medium.alpha_power;

% % Extract the pressure amplitude at each position
% p = extractAmpPhase(sensor_data.p_max_all, 1/kgrid.dt, parameters.transducer.source_freq_hz);
% 
% % Reshape the data, and calculate the volume rate of heat deposition
% p = reshape(p, kgrid.Nx, kgrid.Ny);

p = sensor_data.p_max_all;
source.Q = alpha_np .* p.^2 ./ (kwave_medium.density .* kwave_medium.sound_speed); % Heat delivered to the system
source.T0 = parameters.thermal.temp_0; % Initial temperature distribution

% create kWaveDiffusion object
if isfield(parameters.thermal,'record_t_at_every_step') && ~parameters.thermal.record_t_at_every_step 
    sensor = [];
end
thermal_diff_obj = kWaveDiffusion(kgrid, kwave_medium, source, sensor, 'PlotSim', boolean(parameters.interactive));

% New field temperature
thermal_diff_obj.T = gpuArray(thermal_diff_obj.T);
maxT = thermal_diff_obj.T;

[~, on_off_repetitions, on_steps_n,  on_steps_dur, off_steps_n, off_steps_dur, post_stim_steps_n, post_stim_time_step_dur] = check_thermal_parameters(parameters);

time_status_seq = struct('status', {'off'}, 'time', {0}, 'step', {0}, 'recorded', {1});

% Set up some last parameters for simulation
total_timepoints = 1+parameters.thermal.n_trials*(on_off_repetitions*(1+double(off_steps_n>0))+1);
cur_timepoint = 1;
focal_planeT = gpuArray(zeros([size(maxT,[1,3]) total_timepoints]));
focal_planeT(:,:,cur_timepoint) = squeeze(maxT(:,trans_pos(2),:));

% Loop over trials, so that the heat can accumulate
for trial_i = 1:parameters.thermal.n_trials
  fprintf('Trial %i\n', trial_i)
  % Calculate the heat accumulation within a trial
  for pulse_i = 1:on_off_repetitions
      fprintf('Pulse %i\n', pulse_i)
      
      % Pulse on
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
      
      % Pulse off 
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
  
  % Post-stimulation period
  if post_stim_steps_n > 0 
      thermal_diff_obj.Q = 0;
      thermal_diff_obj.takeTimeStep(post_stim_steps_n, post_stim_time_step_dur);
      time_status_seq = [time_status_seq struct('time', num2cell(max([time_status_seq(:).time]) + (1:post_stim_steps_n)*parameters.thermal.post_stim_time_step_dur), 'step', num2cell(max([time_status_seq(:).step]) + (1:post_stim_steps_n)), 'status', "off",'recorded',0)];
      time_status_seq(end).recorded = 1;
      cur_timepoint = cur_timepoint+1;
      curT = thermal_diff_obj.T;

      focal_planeT(:,:,cur_timepoint) = cat(3, squeeze(curT(:,trans_pos(2),:)));
      
  end
  
  
end

end