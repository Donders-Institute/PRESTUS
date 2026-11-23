function prestus_progress(stage, k, n, detail)
% PRESTUS_PROGRESS  Print a progress line that the Nexus reads from the run log
%
% Prints one line of the form
%   [progress] <stage> <k>/<n> @<unix time in s> <detail>
% to the console (the MATLAB log locally, the SLURM output log on a
% cluster). The Nexus parses these lines into its progress bar and estimates
% the time left from the time stamps. k = 0 marks the start of a solver (its
% detail names the solver); k = n marks the end.
%
% Lines for 0 < k < n are rate-limited per stage: printed only when at least
% 1 % of n has passed since the last line and at least 2 s have elapsed, so a
% long loop adds at most about 100 lines to the log.
%
% Use as:
%   prestus_progress('acoustic', 0, 1, 'solver matlab_gpu')
%   prestus_progress('thermal', seg_i, n_seg, @() sprintf('Tmax %.2f C', max(T(:))))
%
% Input:
%   stage  - short label without spaces ('acoustic', 'thermal', ...)
%   k, n   - steps done and total steps
%   detail - optional text, or a function handle returning text; a handle is
%            evaluated only when the line is printed (for values that are
%            costly to compute, such as a maximum over the grid)
%
% See also: LIVE_VIEW, LOG_TIMER

    persistent LAST
    if isempty(LAST); LAST = struct(); end
    if nargin < 4; detail = ''; end
    key = matlab.lang.makeValidName(char(stage));

    now_s = posixtime(datetime('now', 'TimeZone', 'UTC'));
    if k > 0 && k < n && isfield(LAST, key) && LAST.(key).n == n && k >= LAST.(key).k
        last = LAST.(key);
        if (k - last.k) < n / 100 || (now_s - last.t) < 2
            return
        end
    end
    LAST.(key) = struct('k', k, 'n', n, 't', now_s);

    if isa(detail, 'function_handle')
        try
            detail = detail();
        catch
            detail = '';     % progress must never break the run
        end
    end
    fprintf('[progress] %s %d/%d @%.1f %s\n', char(stage), round(k), round(n), now_s, char(detail));
end
