function rt = telemetry_run_type()
% TELEMETRY_RUN_TYPE  Resolve whether this run is 'technical' or 'deployment'.
%
% Controlled by the PRESTUS_RUN_TYPE environment variable:
%   'technical'  - a test/CI/smoke run, reported separately from real usage
%   anything else (including unset) - 'deployment', an actual simulation run
%
% This lets pass/fail rates from technical runs (e.g. the test suite, see
% RUN_ALL_TESTS) be tracked independently of real deployment usage, so one
% does not skew the other's success-rate statistics.
%
% Use as:
%   rt = telemetry_run_type()
%
% Output:
%   rt - char; 'technical' or 'deployment'
%
% See also: TRACK_USAGE, RUN_ALL_TESTS

    rt = 'deployment';
    try
        if strcmpi(strtrim(getenv('PRESTUS_RUN_TYPE')), 'technical')
            rt = 'technical';
        end
    catch
        rt = 'deployment';
    end
end
