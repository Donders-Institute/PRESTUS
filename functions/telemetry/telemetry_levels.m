function [levels, notice] = telemetry_levels()
% TELEMETRY_LEVELS  The usage-statistics levels this process may send.
%
% Decision order (as TELEMETRY_SETUP):
%   1. PRESTUS_TELEMETRY=off          -> nothing
%   2. a recorded decision (~/.prestus/telemetry.json) -> its levels; a
%      recorded opt-out is never overridden by the environment
%   3. PRESTUS_TELEMETRY=on or a list -> those levels (on: all), this
%      process only
%   4. otherwise                      -> nothing
% A record written before levels existed (notice 1 or 2, opt_in true) keeps
% what was agreed to then: 'basic' and 'performance', plus the skull
% property medians ('skull') under notice 2. Unreadable records send
% nothing.
%
% Use as:
%   [levels, notice] = telemetry_levels()
%
% Output:
%   levels - cellstr, a subset of {'basic', 'performance', 'outcomes'} (or
%            'skull' for a notice-2 record); empty: send nothing
%   notice - notice version the decision was made under
%
% See also: TELEMETRY_LEVEL_INFO, TELEMETRY_ENV_OVERRIDE, TRACK_USAGE

    levels = {};
    notice = 1;
    [env, env_levels] = telemetry_env_override();
    if strcmp(env, 'off'), return, end
    cfg_file = fullfile(prefdir_prestus(), 'telemetry.json');
    if isfile(cfg_file)
        try
            cfg = jsondecode(fileread(cfg_file));
        catch
            return
        end
        if ~isstruct(cfg) || ~isfield(cfg, 'opt_in') || ~isequal(cfg.opt_in, true)
            return
        end
        if isfield(cfg, 'notice') && isnumeric(cfg.notice) && isscalar(cfg.notice)
            notice = double(cfg.notice);
        end
        if isfield(cfg, 'levels')
            levels = known(cfg.levels);
        else
            levels = {'basic', 'performance'};
            if notice >= 2
                levels{end+1} = 'skull';
            end
        end
        return
    end
    if strcmp(env, 'on')
        levels = env_levels;
        notice = telemetry_notice_version();
    end
end

% -------------------------------------------------------------------------
function out = known(v)
% The known level names in v (cellstr, string array, or one char), in prompt order.
    try
        v = cellstr(string(v));
    catch
        v = {};
    end
    names = {telemetry_level_info().name};
    out = names(ismember(names, lower(strtrim(v))));
end
