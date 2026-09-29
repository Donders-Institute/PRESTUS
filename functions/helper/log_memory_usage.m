function mem = log_memory_usage(tag)
% LOG_MEMORY_USAGE  Print this MATLAB process's RSS and the enclosing cgroup's usage.
%
% Diagnostic for out-of-memory kills on SLURM: the OOM handler acts on the
% job's cgroup, which also contains child processes launched via
% run_command (e.g. the BabelBrain Python solver). Logging both lets a
% crashed run show how the limit was split between MATLAB and its child.
%
% Use as:
%   log_memory_usage('before BabelBrain launch')
%
% Input:
%   tag - (1,:) char, label printed with the line
%
% Output:
%   mem - struct with fields (GB, NaN when unavailable):
%         .rss_gb, .hwm_gb           - this process: current / peak resident set
%         .cg_current_gb, .cg_peak_gb, .cg_limit_gb - cgroup usage / peak / limit
%
% Linux only (/proc, /sys/fs/cgroup); prints nothing and returns NaNs elsewhere.

mem = struct('rss_gb', NaN, 'hwm_gb', NaN, ...
             'cg_current_gb', NaN, 'cg_peak_gb', NaN, 'cg_limit_gb', NaN);
if ~isunix || ismac || ~isfile('/proc/self/status')
    return
end

try
    status_txt = fileread('/proc/self/status');
    mem.rss_gb = proc_status_kb(status_txt, 'VmRSS') / 1024^2;
    mem.hwm_gb = proc_status_kb(status_txt, 'VmHWM') / 1024^2;
    [mem.cg_current_gb, mem.cg_peak_gb, mem.cg_limit_gb] = cgroup_memory();
catch
    % Diagnostics must never break a simulation.
end

fprintf(['[memory] %s: MATLAB rss=%.1f GB (peak %.1f GB) | ' ...
         'cgroup current=%.1f GB, peak=%.1f GB, limit=%.1f GB\n'], ...
    tag, mem.rss_gb, mem.hwm_gb, mem.cg_current_gb, mem.cg_peak_gb, mem.cg_limit_gb);
end

%% ------------------------------------------------------------------------
function kb = proc_status_kb(status_txt, key)
tok = regexp(status_txt, [key ':\s+(\d+)\s*kB'], 'tokens', 'once');
if isempty(tok), kb = NaN; else, kb = str2double(tok{1}); end
end

function [cur_gb, peak_gb, limit_gb] = cgroup_memory()
% Supports cgroup v2 (memory.current/peak/max) and v1 (memory controller).
% The limit may sit on an ancestor (SLURM job vs step/task cgroup), so walk up
% to the first finite one.
cur_gb = NaN; peak_gb = NaN; limit_gb = NaN;
cg_txt = fileread('/proc/self/cgroup');
% MATLAB regexp has no \b; delimit the controller name by ':' or ','.
v2 = regexp(cg_txt, '^0::(\S*)$', 'tokens', 'once', 'lineanchors');
v1 = regexp(cg_txt, '^\d+:([^:\n]*,)?memory(,[^:\n]*)?:(\S*)$', 'tokens', 'once', 'lineanchors');
if ~isempty(v1), v1 = v1(end); end
if ~isempty(v1) && isfolder('/sys/fs/cgroup/memory')
    root = '/sys/fs/cgroup/memory'; rel = v1{1};
    f_cur = 'memory.usage_in_bytes'; f_peak = 'memory.max_usage_in_bytes'; f_lim = 'memory.limit_in_bytes';
elseif ~isempty(v2)
    root = '/sys/fs/cgroup'; rel = v2{1};
    f_cur = 'memory.current'; f_peak = 'memory.peak'; f_lim = 'memory.max';
else
    return
end
d = fullfile(root, rel);
cur_gb  = read_bytes_gb(fullfile(d, f_cur));
peak_gb = read_bytes_gb(fullfile(d, f_peak));
while true
    limit_gb = read_bytes_gb(fullfile(d, f_lim));
    if isfinite(limit_gb) || strcmp(d, root) || numel(d) <= numel(root)
        break
    end
    d = fileparts(d);
end
end

function gb = read_bytes_gb(path)
gb = NaN;
if ~isfile(path), return, end
v = str2double(strtrim(fileread(path)));
if isfinite(v) && v < 2^60   % v1 reports ~2^63 for "unlimited"; v2 reports "max" (NaN)
    gb = v / 1024^3;
end
end
