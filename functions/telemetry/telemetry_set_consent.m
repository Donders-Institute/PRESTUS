function telemetry_set_consent(choice)
% TELEMETRY_SET_CONSENT  Record which usage statistics this user shares.
%
% Call this once before running the pipeline — from any MATLAB session or
% script — or any time you want to change your choice.
%
% Use as:
%   telemetry_set_consent(true)                       % share all levels
%   telemetry_set_consent({'basic', 'performance'})   % share these levels
%   telemetry_set_consent(false)                      % share nothing
%
% Levels (TELEMETRY_LEVEL_INFO): 'basic' (usage and errors), 'performance',
% 'outcomes' (summary outcomes, rounded and anonymous). An empty list shares nothing.
% The decision is stored in ~/.prestus/telemetry.json and is read by every
% subsequent pipeline run. To revert to the undecided state, call
% telemetry_setup_reset().
%
% See also: TELEMETRY_SETUP, TELEMETRY_SETUP_RESET, TELEMETRY_LEVEL_INFO

    arguments
        choice
    end

    names = {telemetry_level_info().name};
    if (islogical(choice) || isnumeric(choice)) && isscalar(choice)
        levels = names(1:numel(names) * logical(choice));
    else
        requested = lower(strtrim(cellstr(string(choice))));
        bad = setdiff(requested, names);
        if ~isempty(bad)
            error('telemetry_set_consent:level', 'Unknown level ''%s''; use %s.', ...
                  bad{1}, strjoin(names, ', '));
        end
        levels = names(ismember(names, requested));
    end
    telemetry_write_consent(levels);

    if isempty(levels)
        fprintf('Usage statistics: nothing shared. Change with telemetry_set_consent.\n');
    else
        fprintf('Usage statistics: sharing %s. Change with telemetry_set_consent.\n', ...
                strjoin(levels, ', '));
    end
end
