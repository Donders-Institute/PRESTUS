function varargout = log_timer(action, label, varargin)
% LOG_TIMER - Track USED disk space Δ (GB/MB) per label independently
% Usage:
% log_timer('start', 'func1', '/path/to/monitor');
% log_timer('start', 'func2');  % Independent!
% log_timer('stop', 'func2');
% log_timer('stop', 'func1');

persistent TIMER_STATES COUNTER MONITOR_PATH

if isempty(TIMER_STATES)
    TIMER_STATES = containers.Map();
    COUNTER = containers.Map();
end

switch action
    case 'start'
        if ~TIMER_STATES.isKey(label)
            COUNTER(label) = 0;
        end
        COUNTER(label) = COUNTER(label) + 1;
        
        % Parse monitor path (3rd arg optional)
        monitor_path = pwd;
        if length(varargin) > 0 && exist(varargin{1}, 'dir')
            monitor_path = varargin{1};
        elseif ~isempty(MONITOR_PATH)
            monitor_path = MONITOR_PATH;
        end
        
        state.start_time = tic;
        state.start_ram = log_process_memory();
        state.start_used_disk = log_disk_space(monitor_path);
        state.counter = COUNTER(label);
        state.monitor_path = monitor_path;
        TIMER_STATES(label) = state;
        
        fprintf('▶️ %s [#%d] 💿%s\n', label, state.counter, monitor_path);
        
    case 'stop'
        if ~TIMER_STATES.isKey(label)
            fprintf('⚠️ log_timer stop %s: No active timer\n', label);
            return;
        end
        
        state = TIMER_STATES(label);
        time_sec = toc(state.start_time);
        ram_delta = log_process_memory() - state.start_ram;
        [used_now_bytes, ~] = log_disk_space(state.monitor_path);
        used_delta_bytes = used_now_bytes - state.start_used_disk;
        used_delta_gb = used_delta_bytes / 1e9;
        used_delta_mb = used_delta_bytes / 1e6;
        
        fprintf('⏱️ %-20s %.2fs | ΔRAM %.2fGB | ΔUSED %.2fGB (%.0fMB) [#%d]\n', ...
            label, time_sec, ram_delta, used_delta_gb, used_delta_mb, state.counter);
        
        % Remove state
        TIMER_STATES.remove(label);
        
    case 'set_monitor'
        MONITOR_PATH = label;  % Use label as path here
        if ~exist(MONITOR_PATH, 'dir'), mkdir(MONITOR_PATH); end
        fprintf('💿 Monitoring USED space: %s\n', MONITOR_PATH);
        
    otherwise
        error('log_timer:Action', '''start'', ''stop'', or ''set_monitor''');
end
end

%% HELPERS (unchanged)
function [used_bytes, free_bytes] = log_disk_space(monitor_path)
    if isunix
        cmd_used = sprintf('df --block-size=1 --output=used "%s" | tail -n1', monitor_path);
        cmd_free = sprintf('df --block-size=1 --output=avail "%s" | tail -n1', monitor_path);
        [~, out_used] = system(cmd_used);
        [~, out_free] = system(cmd_free);
        used_bytes = str2double(strtrim(out_used));
        free_bytes = str2double(strtrim(out_free));
    elseif ispc
        drive = monitor_path(1);
        [~, out_used] = system(sprintf('powershell -c "(Get-PSDrive ''%c'').Used"', drive));
        [~, out_free] = system(sprintf('powershell -c "(Get-PSDrive ''%c'').Free"', drive));
        used_bytes = str2double(strtrim(out_used));
        free_bytes = str2double(strtrim(out_free));
    else
        used_bytes = NaN; free_bytes = NaN;
    end
end

function ram_gb = log_process_memory()
    if ispc
        try, m = memory; ram_gb = m.MemUsedMATLAB/1e9; catch, ram_gb = NaN; end
    else
        [~, out] = system(sprintf('ps -o rss= -p %d', feature('getpid')));
        ram_gb = str2double(strtrim(out)) / 1024^2;
    end
end
