function [mode, levels] = telemetry_env_override()
% TELEMETRY_ENV_OVERRIDE  Read the PRESTUS_TELEMETRY environment variable.
%
% PRESTUS_TELEMETRY sets the telemetry decision for the current process
% only. It suppresses every prompt and is never written to
% ~/.prestus/telemetry.json. 'off' takes precedence over a recorded
% decision; sharing applies only when no decision is recorded, so it never
% overrides a recorded opt-out:
%   'off', '0', 'false' - opted out: no prompt, nothing sent
%   'on',  '1', 'true'  - all levels shared, no prompt
%   a comma-separated list of levels, e.g. 'basic,performance' - those
%                         levels (TELEMETRY_LEVEL_INFO), no prompt
%   unset or empty      - no override; the recorded decision (or the
%                         prompt, see TELEMETRY_SETUP) applies
% Values are case-insensitive and surrounding whitespace is ignored. Any
% other value raises the warning 'prestus:telemetry:invalidEnv' (once per
% distinct value per MATLAB session) and is treated as unset.
%
% Use as:
%   [mode, levels] = telemetry_env_override()
%
% Output:
%   mode   - char; 'off', 'on', or '' (no override)
%   levels - cellstr of the levels 'on' shares (all for on/1/true); empty
%            unless mode is 'on'
%
% See also: TELEMETRY_SETUP, TELEMETRY_LEVELS, TRACK_USAGE, TELEMETRY_SET_CONSENT

    persistent warned_values
    if isempty(warned_values)
        warned_values = {};
    end

    mode = '';
    levels = {};
    try
        val = lower(strtrim(getenv('PRESTUS_TELEMETRY')));
    catch
        return
    end
    if isempty(val)
        return
    end

    names = {telemetry_level_info().name};
    switch val
        case {'off', '0', 'false'}
            mode = 'off';
            return
        case {'on', '1', 'true'}
            mode = 'on';
            levels = names;
            return
    end
    parts = strtrim(strsplit(val, ','));
    if all(ismember(parts, names))
        mode = 'on';
        levels = names(ismember(names, parts));
        return
    end
    if ~ismember(val, warned_values)
        warned_values{end+1} = val;
        warning('prestus:telemetry:invalidEnv', ...
            ['PRESTUS_TELEMETRY=''%s'' is not recognised (use off/0/false, on/1/true, ' ...
             'or levels such as basic,performance,outcomes); ignoring it.'], val);
    end
end
