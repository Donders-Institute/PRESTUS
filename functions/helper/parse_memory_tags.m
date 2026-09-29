function mem = parse_memory_tags(text)
% PARSE_MEMORY_TAGS  Extract peak memory figures from [memory]/[memory_bb] log lines.
%
% Two sources write these lines into the console output, which lands in the
% pipeline log (diary), the SLURM output log, and the BabelBrain run logs:
%   [memory]    - MATLAB, from log_memory_usage:
%                 "<tag>: MATLAB rss=X GB (peak Y GB) | cgroup current=X GB, peak=Y GB, limit=Z GB"
%   [memory_bb] - the BabelBrain Python child, from run_babelbrain_case.py:
%                 "<tag>: python rss=X GB (peak Y GB) | cgroup current=X GB, peak=Y GB"
%
% Every figure is a high-water mark, so each output is the maximum over all
% matching lines. Parsing the same lines twice (e.g. pipeline log + SLURM log)
% therefore changes nothing, and several BabelBrain runs in one job (acoustic
% + water-only) reduce to the overall peak. The cgroup peak is the kernel's
% continuous high-water mark of the job cgroup (MATLAB + child), so the
% newest line of either tag carries the job peak up to that point; the dense
% [memory_bb] lines are what make it current during the BabelBrain run.
%
% Use as:
%   mem = parse_memory_tags(fileread(log_file))
%
% Input:
%   text - (1,:) char or string, log text (may be empty)
%
% Output:
%   mem - struct with fields (GB, NaN when no matching line / value):
%         .bb_peak_gb      - peak RSS of the BabelBrain Python process
%         .matlab_peak_gb  - peak RSS of the MATLAB process
%         .cgroup_peak_gb  - peak usage of the job's memory cgroup
%         .cgroup_limit_gb - cgroup memory limit (from [memory] lines)
%         .n_bb_lines, .n_matlab_lines - number of lines matched per tag
%
% See also: LOG_MEMORY_USAGE, HPC_JOB_REPORT

mem = struct('bb_peak_gb', NaN, 'matlab_peak_gb', NaN, ...
             'cgroup_peak_gb', NaN, 'cgroup_limit_gb', NaN, ...
             'n_bb_lines', 0, 'n_matlab_lines', 0);
text = char(text);
if isempty(text)
    return
end

num = '(\S+)';   % "12.3", "NaN" (MATLAB) or "nan" (Python); str2double maps both to NaN
bb = regexp(text, ['\[memory_bb\][^\n]*?python rss=\S+ GB \(peak ' num ' GB\)' ...
                   '[^\n]*?cgroup current=\S+ GB, peak=' num ' GB'], 'tokens');
ml = regexp(text, ['\[memory\][^\n]*?MATLAB rss=\S+ GB \(peak ' num ' GB\)' ...
                   '[^\n]*?cgroup current=\S+ GB, peak=' num ' GB, limit=' num ' GB'], 'tokens');

mem.n_bb_lines     = numel(bb);
mem.n_matlab_lines = numel(ml);
bb = token_matrix(bb, 2);
ml = token_matrix(ml, 3);

mem.bb_peak_gb      = nanmax_or_nan(bb(:, 1));
mem.matlab_peak_gb  = nanmax_or_nan(ml(:, 1));
mem.cgroup_peak_gb  = nanmax_or_nan([bb(:, 2); ml(:, 2)]);
mem.cgroup_limit_gb = nanmax_or_nan(ml(:, 3));
end

%% ------------------------------------------------------------------------
function m = token_matrix(tok, n_cols)
m = NaN(numel(tok), n_cols);
for i = 1:numel(tok)
    m(i, :) = str2double(tok{i});
end
end

function v = nanmax_or_nan(x)
x = x(isfinite(x));
if isempty(x), v = NaN; else, v = max(x); end
end
