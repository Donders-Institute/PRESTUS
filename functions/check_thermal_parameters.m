function [on_off_step_duration, on_off_repetitions, on_steps_n,  on_steps_dur, off_steps_n, off_steps_dur, post_stim_steps_n, post_stim_time_step_dur] = check_thermal_parameters(parameters)
arguments
    parameters struct
end

disp('Checking thermal parameters...');

% Calculates the duration of a duty cycle in seconds
% If no on_off_step_duration is specified, a duty cycle must be
if ~isfield(parameters.thermal,'on_off_step_duration')
    if parameters.thermal.duty_cycle > 0 && parameters.thermal.duty_cycle < 1
        on_off_step_duration = parameters.thermal.sim_time_steps/min([parameters.thermal.duty_cycle 1-parameters.thermal.duty_cycle ]); % one on+off cycle duration
    else
        on_off_step_duration = parameters.thermal.sim_time_steps;
    end
else
    on_off_step_duration = parameters.thermal.on_off_step_duration;
end

% Amount of steps to be simulated
on_off_repetitions = parameters.thermal.stim_duration/on_off_step_duration;

fprintf('Duration of 1 repetition of the on+off cycle: %.3f s; repeated for %.2f times.\n', on_off_step_duration, on_off_repetitions)

% Shows an argument that must be met for k-wave to accept the parameters of
% 'on_off_step_duration' and 'sim_time_step'
on_off_repetitions = round_if_integer(on_off_repetitions, sprintf('The total stimulation duration (%.3f s) should be divisible by the on+off cycle duration (%.3f s)\n',parameters.thermal.stim_duration, on_off_step_duration));

if  isfield(parameters.thermal,'equal_steps') && parameters.thermal.equal_steps == 1

    on_steps_n = on_off_step_duration*parameters.thermal.duty_cycle/parameters.thermal.sim_time_steps;
    off_steps_n = on_off_step_duration*(1-parameters.thermal.duty_cycle)/parameters.thermal.sim_time_steps;
    
    fprintf('1 on+off cycle contains %.2f on steps and %.2f off steps assuming equal steps of %.2f s each.\n', on_steps_n, off_steps_n, parameters.thermal.sim_time_steps)
    
    on_steps_n = round_if_integer(on_steps_n, 'The number of ''on'' steps within the on+off cycle should be integer');
    off_steps_n = round_if_integer(off_steps_n, 'The number of ''off'' steps within the on+off cycle should be integer');

end

% Allows for a different step duration after stimulation
if ~isfield(parameters.thermal,'post_stim_time_step_dur')
    post_stim_time_step_dur = parameters.thermal.sim_time_steps;
else
    post_stim_time_step_dur = parameters.thermal.post_stim_time_step_dur;
end

% Defines the number of post stimulation steps, sets the post_stim_period to 0 if a continuous protocol is used
if isfield(parameters.thermal,'continuous_protocol') && parameters.thermal.continuous_protocol == 1
    post_stim_period = 0;
    post_stim_steps_n = 0;
else
    post_stim_period = parameters.thermal.iti - parameters.thermal.stim_duration;
    post_stim_steps_n = post_stim_period / post_stim_time_step_dur;
end

fprintf('Post-stimulation off period is %.2f s, consists of %.0f simulation steps, each taking %.2f s.\n', ...
    post_stim_period, post_stim_steps_n, post_stim_time_step_dur)

post_stim_steps_n = round_if_integer(post_stim_steps_n, 'Number of simulation steps must be integer');

% This can be used if the number for 'sim_time_steps' is hard to find,
% but determining the 'sim_time_steps' is still ideal
if  isfield(parameters.thermal,'equal_steps') && parameters.thermal.equal_steps == 0
    on_steps_dur = on_off_step_duration*parameters.thermal.duty_cycle;
    off_steps_dur = on_off_step_duration*(1-parameters.thermal.duty_cycle);
    on_steps_n = 1;
    off_steps_n = 1;
    fprintf('Each cycle = 1 on and 1 off step with durations of %.3f s & %.3f s respectively.\n', ...
        on_steps_dur, off_steps_dur)
else
    on_steps_dur = parameters.thermal.sim_time_steps;
    off_steps_dur = parameters.thermal.sim_time_steps;
    fprintf('Each cycle = %.0f on  and %.0f off steps assuming equal steps of %.3f s each.\n', ...
        on_steps_n, off_steps_n, parameters.thermal.sim_time_steps)
end

end