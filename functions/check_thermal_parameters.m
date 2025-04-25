function [on_steps_n,  on_steps_dur, ...
    off_steps_n, off_steps_dur, post_stim_step_n, post_stim_step_dur] = ...
    check_thermal_parameters(parameters)
    
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
%   on_steps_n              - Number of "on" setps in one pulse repetition interval.
%   on_steps_dur            - Duration of the "on" pulse (in seconds).
%   off_steps_n             - Number of "off" steps in one pulse repetition interval.
%   off_steps_dur           - Duration of the "off" step (in seconds).
%   post_stim_step_n        - Number of post-stimulation steps.
%   post_stim_step_dur      - Duration of each post-stimulation step (in seconds).

    arguments
        parameters struct
    end

    disp('Checking thermal parameters...');

    % Establish duration of an on-off pulse repetition interval (PRI)
    if ~isfield(parameters.thermal,'pri_duration')
        error("Duration of pulse repetition interval is not specified.");
    else
        pri_duration = parameters.thermal.pri_duration;
    end

    % Set or validate post-stimulation step duration
    if ~isfield(parameters.thermal,'post_stim_dur')
        disp("No post-stimulation time requested...")
        post_stim_step_dur = 0;
        post_stim_step_n = 0;
    else
        post_stim_dur = parameters.thermal.post_stim_dur;
        if isfield(parameters.thermal, 'post_time_steps')
            post_stim_step_dur = parameters.thermal.post_time_steps;
        else
            post_stim_step_dur = parameters.thermal.sim_time_steps;
        end
        post_stim_step_n = post_stim_dur/post_stim_step_dur;
        % Ensure the number of post-stimulation steps is an integer
        post_stim_step_n = round_if_integer(post_stim_step_n, 'Number of simulation steps must be integer');
        fprintf('Post-stimulation off period is %.2f s, consists of %.0f simulation steps, each taking %.2f s.\n', ...
            post_stim_dur, post_stim_step_n, post_stim_step_dur);
    end

    % Calculate number of "on" and "off" steps
    if isfield(parameters.thermal,'equal_steps') && parameters.thermal.equal_steps == 0
        % Variable step duration between on & off [fixed step count = 1]
        on_steps_dur = pri_duration * parameters.thermal.duty_cycle;
        off_steps_dur = pri_duration * (1 - parameters.thermal.duty_cycle);
        on_steps_n = 1;
        off_steps_n = 1;

        fprintf('Each cycle = 1 on and 1 off step with durations of %.3f s & %.3f s respectively.\n', ...
            on_steps_dur, off_steps_dur);
    else
        % Variable step number between on & off [fixed step duration: sim_time_steps]
        on_steps_dur = parameters.thermal.sim_time_steps;
        off_steps_dur = parameters.thermal.sim_time_steps;
        on_steps_n = pri_duration * parameters.thermal.duty_cycle / on_steps_dur;
        on_steps_n = round_if_integer(on_steps_n, 'The number of ''on'' steps within the on+off cycle should be integer');
        off_steps_n = pri_duration * (1 - parameters.thermal.duty_cycle) / off_steps_dur;
        off_steps_n = round_if_integer(off_steps_n, 'The number of ''off'' steps within the on+off cycle should be integer');
    
        fprintf('1 on+off cycle contains %.2f on steps and %.2f off steps assuming equal steps of %.2f s each.\n', ...
            on_steps_n, off_steps_n, parameters.thermal.sim_time_steps);
    end

end
