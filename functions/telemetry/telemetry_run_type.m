function rt = telemetry_run_type()
% TELEMETRY_RUN_TYPE  Resolve whether this run is 'technical', 'development' or 'deployment'.
%
% Controlled by the PRESTUS_RUN_TYPE environment variable:
%   'technical'   - a test/CI/smoke run, reported separately from real usage
%   'development' - a run while developing PRESTUS (set it in your shell or
%                   startup.m on a development machine), filtered out of usage
%                   statistics; telemetry also sends prestus_branch and
%                   prestus_dirty, which mark runs from a modified checkout
%   anything else (including unset) - 'deployment', an actual simulation run
%
% HPC_SUBMIT_JOB writes the run type into the job scripts it submits, so a
% job reports the run type of the session that submitted it.
%
% This lets pass/fail rates from technical runs (e.g. the test suite, see
% RUN_ALL_TESTS) be tracked independently of real deployment usage, so one
% does not skew the other's success-rate statistics.
%
% Use as:
%   rt = telemetry_run_type()
%
% Output:
%   rt - char; 'technical', 'development' or 'deployment'
%
% See also: TRACK_USAGE, RUN_ALL_TESTS, HPC_SUBMIT_JOB

    rt = 'deployment';
    try
        v = lower(strtrim(getenv('PRESTUS_RUN_TYPE')));
        if any(strcmp(v, {'technical', 'development'}))
            rt = v;
        end
    catch
        rt = 'deployment';
    end
end
