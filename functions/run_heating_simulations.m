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

% new field for cem43 (replicate above)
thermal_diff_obj.cem43 = gpuArray(zeros(size(thermal_diff_obj.T)));
maxCEM43 = thermal_diff_obj.cem43;

[~, on_off_repetitions, on_steps_n,  on_steps_dur, off_steps_n, off_steps_dur, post_stim_steps_n, post_stim_time_step_dur] = check_thermal_parameters(parameters);

time_status_seq = struct('status', {'off'}, 'time', {0}, 'step', {0}, 'recorded', {1});

% Set up some last parameters for simulation
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
        for i = 1:length(parameters.start_break_trials)
            if trial_i >= parameters.start_break_trials(i) && trial_i <= parameters.stop_break_trials(i)
                is_break_period = 1;
            end
        end
    end
    
  % Calculate the heat accumulation within a trial
  for pulse_i = 1:on_off_repetitions
      fprintf('Pulse %i\n', pulse_i)

      if is_break_period == 0
          thermal_diff_obj.Q = source.Q;
      else
          thermal_diff_obj.Q(:,:,:) = 0;   
          disp('Break active');
      end

      thermal_diff_obj.takeTimeStep(on_steps_n, on_steps_dur);
      time_status_seq = [time_status_seq struct('time', num2cell(max([time_status_seq(:).time]) + (1:on_steps_n)*on_steps_dur), ...
                                                'step', num2cell(max([time_status_seq(:).step]) + (1:on_steps_n)), 'status', "on",'recorded',0)];
      curT = thermal_diff_obj.T;
      curCEM43 = thermal_diff_obj.cem43;
      cur_timepoint = cur_timepoint+1;
      time_status_seq(end).recorded = 1;
      focal_planeT(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
      CEM43(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));

      % JQK: shouldn't this comparison be the other way around?
      if any(maxT>curT,'all')
          maxT = curT;
      end

      if any(maxCEM43>curCEM43, 'all')
          maxCEM43 = curCEM43;
      end

      % Pulse off 
      if off_steps_n>0
          thermal_diff_obj.Q(:,:,:) = 0;
          thermal_diff_obj.takeTimeStep(off_steps_n, off_steps_dur);
          time_status_seq = [time_status_seq struct('time', num2cell(max([time_status_seq(:).time]) + (1:off_steps_n)*off_steps_dur), 'step', num2cell(max([time_status_seq(:).step]) + (1:off_steps_n)), 'status', "off",'recorded',0)];
          cur_timepoint = cur_timepoint+1;
          time_status_seq(end).recorded = 1;
          curT = thermal_diff_obj.T;
          curCEM43 = thermal_diff_obj.cem43;

          focal_planeT(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
          CEM43(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));
      end
  end
  
end
% Post-stimulation period
if post_stim_steps_n > 0 
  thermal_diff_obj.Q = 0;
  thermal_diff_obj.takeTimeStep(post_stim_steps_n, post_stim_time_step_dur);
  time_status_seq = [time_status_seq struct('time', num2cell(max([time_status_seq(:).time]) + (1:off_steps_n)*off_steps_dur), 'step', num2cell(max([time_status_seq(:).step]) + (1:off_steps_n)), 'status', "off",'recorded',0)];
  time_status_seq(end).recorded = 1;
  cur_timepoint = cur_timepoint+1;
  curT = thermal_diff_obj.T;
  curCEM43 = thermal_diff_obj.cem43;
  focal_planeT(:,:,cur_timepoint) = squeeze(curT(:,trans_pos(2),:));
  CEM43(:,:,cur_timepoint) = squeeze(curCEM43(:,trans_pos(2),:));
                end

end