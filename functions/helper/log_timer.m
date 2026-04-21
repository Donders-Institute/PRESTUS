function varargout = log_timer(action, label, varargin)
% LOG_TIMER  Track elapsed time, peak RAM, and drive space for labelled code blocks
%
% Persistent timers accumulate per-label state, so multiple independent
% blocks can be timed concurrently. A background MATLAB timer samples RAM
% every 0.5 s to capture peak usage.
%
% Use as:
%   log_timer('start', 'label')
%   log_timer('start', 'label', '/path/to/monitor')
%   result = log_timer('stop', 'label')
%
% Input:
%   action   - 'start' | 'stop'
%   label    - string identifier for the timed block; must match between start/stop
%   varargin - (start only) optional path to monitor for free disk space
%
% Output (stop only):
%   varargout{1} - struct with fields:
%                    .elapsed_s    — wall-clock time [s]
%                    .peak_ram_gb  — peak RAM usage during block [GB]
%                    .free_disk_gb — free disk space at stop (if path provided) [GB]
%
% Note:
%   Uses MATLAB persistent variables and timer() objects. If MATLAB exits
%   unexpectedly during a timed block, background samplers may not be cleaned
%   up. Call log_timer('stop', label) explicitly to release timers.
%
% See also: PATH_LOG_SETUP

persistent TIMER_STATES COUNTER MONITOR_PATH SAMPLER_TIMERS

if isempty(TIMER_STATES), TIMER_STATES = containers.Map(); end
if isempty(COUNTER), COUNTER = containers.Map(); end
if isempty(SAMPLER_TIMERS), SAMPLER_TIMERS = containers.Map(); end
if isempty(MONITOR_PATH), MONITOR_PATH = pwd; end

switch action
    case 'start'
        counter_val = 0; if COUNTER.isKey(label), counter_val = COUNTER(label); end
        COUNTER(label) = counter_val + 1;
        
        MONITOR_PATH_LOCAL = MONITOR_PATH;
        if nargin>2 && exist(varargin{1},'dir'), MONITOR_PATH_LOCAL = varargin{1}; end
        
        state.start_time = tic; state.start_ram = log_process_memory();
        state.start_used_disk = log_disk_space(MONITOR_PATH_LOCAL);
        state.counter = COUNTER(label); state.monitor_path = MONITOR_PATH_LOCAL;
        state.peak_ram_gb = state.start_ram;  % Track in state
        TIMER_STATES(label) = state;
        
        fprintf('▶ %s [#%d] 💿%s\n', label, state.counter, state.monitor_path);
        
        % Start sampling resource use (every 0.5 s) to update max. RAM requirement
        sampler = timer('Period', 0.5, 'ExecutionMode', 'fixedRate', ...
                       'TimerFcn', @(t,e)sample_peak_closure(label));
        start(sampler);
        SAMPLER_TIMERS(label) = sampler;
        
    case 'stop'
        % CLEANUP
        if SAMPLER_TIMERS.isKey(label)
            sampler = SAMPLER_TIMERS(label);
            if isvalid(sampler), stop(sampler); delete(sampler); end
            SAMPLER_TIMERS.remove(label);
        end
        
        state = TIMER_STATES(label);
        time_sec = toc(state.start_time);
        end_ram = log_process_memory();
        peak_ram = state.peak_ram_gb;  % From state!
        
        ram_delta = end_ram - state.start_ram; % This would be the RAM between onset and offset of function (unused).
        ram_peak_delta = peak_ram - state.start_ram; % This captures the change in RAM use during this processing stage.
        
        [used_now_bytes, ~] = log_disk_space(state.monitor_path);
        used_delta_bytes = used_now_bytes - state.start_used_disk;
        used_delta_gb = used_delta_bytes / 1e9; used_delta_mb = used_delta_bytes / 1e6;
        
        fprintf('⏱ %-20s %.2fs | ΔRAM %.2fGB | PEAKRAM %.2fGB | ΔHDD %.2fGB (%.0fMB) [#%d]\n', ...
            label, time_sec, ram_peak_delta, state.peak_ram_gb, used_delta_gb, used_delta_mb, state.counter);

        TIMER_STATES.remove(label);
        
    case 'set_monitor', MONITOR_PATH = label;
    otherwise, error('Invalid action');
end


function [used_bytes, free_bytes] = log_disk_space(monitor_path) % fetch disk space of specified directory
    if isunix
        cmd_used = sprintf('df --block-size=1 --output=used "%s" | tail -n1', monitor_path);
        cmd_free = sprintf('df --block-size=1 --output=avail "%s" | tail -n1', monitor_path);
        [~, out_used] = system(cmd_used);
        [~, out_free] = system(cmd_free);
        used_bytes = str2double(strtrim(out_used));
        free_bytes = str2double(strtrim(out_free));
    elseif ispc
        drive = monitor_path(1);
        driveLetter = char(drive);
        driveLetter = driveLetter(1);
        [~, out_used] = system(sprintf('powershell -c "(Get-PSDrive ''%c'').Used"', driveLetter));
        [~, out_free] = system(sprintf('powershell -c "(Get-PSDrive ''%c'').Free"', driveLetter));
        used_bytes = str2double(strtrim(out_used));
        free_bytes = str2double(strtrim(out_free));
    else
        used_bytes = NaN; free_bytes = NaN;
    end
end

function ram_gb = log_process_memory() % fetch memory use of a process
    if ispc
        try m = memory; ram_gb = m.MemUsedMATLAB/1e9; catch, ram_gb = NaN; end
    else
        [~, out] = system(sprintf('ps -o rss= -p %d', feature('getpid')));
        ram_gb = str2double(strtrim(out)) / 1024^2;
    end
end

function sample_peak_closure(label)
    % TIMER_STATES is a persistent variable in the parent scope. If the
    % parent function was cleared mid-run (e.g. due to an error), it will
    % no longer exist when this timer callback fires. Guard defensively and
    % self-stop the timer to avoid repeated "Reference to a cleared variable"
    % warnings.
    try
        if ~TIMER_STATES.isKey(label)
            return   % Timer already stopped or entry removed — nothing to do.
        end
        state = TIMER_STATES(label);
        current_ram = log_process_memory();
        if current_ram > state.peak_ram_gb
            state.peak_ram_gb = current_ram;
            TIMER_STATES(label) = state;
        end
    catch
        % Parent scope was cleared; stop and delete this timer so it doesn't
        % keep firing.
        t_list = timerfindall('TimerFcn', @(t,e)sample_peak_closure(label));
        for k = 1:numel(t_list)
            if isvalid(t_list(k)), stop(t_list(k)); delete(t_list(k)); end
        end
    end
end
end