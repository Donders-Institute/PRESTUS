function [on_off_step_duration, on_off_repetitions, on_steps_n,  on_steps_dur, off_steps_n, off_steps_dur, post_stim_steps_n, post_stim_time_step_dur] = check_thermal_parameters(parameters)
    
% CHECK_THERMAL_PARAMETERS Validates and computes thermal simulation parameters.
%
% This function checks and computes various parameters required for thermal 
% simulations. It ensures compatibility with simulation constraints, such as 
% step durations and duty cycles, and calculates the number of steps for 
% stimulation and post-stimulation periods.
%
% Input:
%   parameters - Struct containing thermal simulation parameters (e.g., duty cycle, 
%                step durations, stimulation duration).
%
% Output:
%   on_off_step_duration    - Duration of one on+off cycle (in seconds).
%   on_off_repetitions      - Number of repetitions of the on+off cycle.
%   on_steps_n              - Number of "on" steps in one cycle.
%   on_steps_dur            - Duration of the "on" step (in seconds).
%   off_steps_n             - Number of "off" steps in one cycle.
%   off_steps_dur           - Duration of the "off" step (in seconds).
%   post_stim_steps_n       - Number of post-stimulation steps.
%   post_stim_time_step_dur - Duration of each post-stimulation step (in seconds).

    arguments
        parameters struct
    end

    disp('Checking thermal parameters...');

    % Calculate the duration of a duty cycle if not explicitly provided
    if ~isfield(parameters.thermal,'on_off_step_duration')
        if parameters.thermal.duty_cycle > 0 && parameters.thermal.duty_cycle < 1
            on_off_step_duration = parameters.thermal.sim_time_steps / ...
                min([parameters.thermal.duty_cycle 1 - parameters.thermal.duty_cycle]);
        else
            on_off_step_duration = parameters.thermal.sim_time_steps;
        end
    else
        on_off_step_duration = parameters.thermal.on_off_step_duration;
    end

    % Determine the number of repetitions for the on+off cycle
    on_off_repetitions = parameters.thermal.stim_duration / on_off_step_duration;
    fprintf('Duration of 1 repetition of the on+off cycle: %.3f s; repeated for %.2f times.\n', ...
        on_off_step_duration, on_off_repetitions);

    % Ensure the total stimulation duration is divisible by the cycle duration
    on_off_repetitions = round_if_integer(on_off_repetitions, ...
        sprintf('The total stimulation duration (%.3f s) should be divisible by the on+off cycle duration (%.3f s)\n', ...
        parameters.thermal.stim_duration, on_off_step_duration));

    % Calculate "on" and "off" steps if equal step durations are specified
    if isfield(parameters.thermal,'equal_steps') && parameters.thermal.equal_steps == 1
        on_steps_n = on_off_step_duration * parameters.thermal.duty_cycle / parameters.thermal.sim_time_steps;
        off_steps_n = on_off_step_duration * (1 - parameters.thermal.duty_cycle) / parameters.thermal.sim_time_steps;

        fprintf('1 on+off cycle contains %.2f on steps and %.2f off steps assuming equal steps of %.2f s each.\n', ...
            on_steps_n, off_steps_n, parameters.thermal.sim_time_steps);

        % Ensure "on" and "off" steps are integers
        on_steps_n = round_if_integer(on_steps_n, 'The number of ''on'' steps within the on+off cycle should be integer');
        off_steps_n = round_if_integer(off_steps_n, 'The number of ''off'' steps within the on+off cycle should be integer');
    end

    % Set or validate post-stimulation step duration
    if ~isfield(parameters.thermal,'post_stim_time_step_dur')
        post_stim_time_step_dur = parameters.thermal.sim_time_steps;
    else
        post_stim_time_step_dur = parameters.thermal.post_stim_time_step_dur;
    end

    % Calculate post-stimulation period and number of steps
    if isfield(parameters.thermal,'continuous_protocol') && parameters.thermal.continuous_protocol == 1
        post_stim_period = 0;
        post_stim_steps_n = 0;
    else
        post_stim_period = parameters.thermal.iti - parameters.thermal.stim_duration;
        post_stim_steps_n = post_stim_period / post_stim_time_step_dur;
    end

    fprintf('Post-stimulation off period is %.2f s, consists of %.0f simulation steps, each taking %.2f s.\n', ...
        post_stim_period, post_stim_steps_n, post_stim_time_step_dur);

    % Ensure the number of post-stimulation steps is an integer
    post_stim_steps_n = round_if_integer(post_stim_steps_n, 'Number of simulation steps must be integer');

    % Handle cases where unequal step durations are allowed
    if isfield(parameters.thermal,'equal_steps') && parameters.thermal.equal_steps == 0
        on_steps_dur = on_off_step_duration * parameters.thermal.duty_cycle;
        off_steps_dur = on_off_step_duration * (1 - parameters.thermal.duty_cycle);
        on_steps_n = 1;
        off_steps_n = 1;

        fprintf('Each cycle = 1 on and 1 off step with durations of %.3f s & %.3f s respectively.\n', ...
            on_steps_dur, off_steps_dur);
    else
        % Default to equal step durations if not specified otherwise
        on_steps_dur = parameters.thermal.sim_time_steps;
        off_steps_dur = parameters.thermal.sim_time_steps;

        fprintf('Each cycle = %.0f on and %.0f off steps assuming equal steps of %.3f s each.\n', ...
            on_steps_n, off_steps_n, parameters.thermal.sim_time_steps);
    end

end
