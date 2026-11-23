function info = telemetry_level_info()
% TELEMETRY_LEVEL_INFO  The usage-statistics levels a user can share, in prompt order.
%
% Each level is a group of fields (see doc_telemetry.md). Consent records
% the levels a user ticked; "Share all" ticks all of them, and unticking all
% of them shares nothing. Every record carries the same small envelope
% (event, run and install ID, PRESTUS version, run type, status, and the
% consent itself) whichever levels are shared.
%   basic        - which features, platforms and versions are used, and
%                  whether runs succeed (error identifier, failing stage)
%   performance  - hardware, requested resources, grid size, stage times and
%                  peak memory; these train the Nexus run-time estimates
%   outcomes     - eight rounded summary statistics of the simulated field:
%                  brain intensity and its ratio to free water, brain MI,
%                  targeting error, brain and skull temperature rise, duty
%                  cycle and sonication time; plus the modelled skull
%                  property medians. Never maps, coordinates or target names.
% The Nexus shows the same levels (prestus_nexus/telemetry.py, LEVELS).
%
% Use as:
%   info = telemetry_level_info()
%
% Output:
%   info - struct array with fields name, label, description
%
% See also: TELEMETRY_LEVELS, TELEMETRY_SETUP, TELEMETRY_SET_CONSENT

    info = struct( ...
        'name',  {'basic', 'performance', 'outcomes'}, ...
        'label', {'Usage and errors', 'Performance', 'Summary outcomes (anonymous)'}, ...
        'description', { ...
            'Which features, platforms, and versions are used, and whether runs succeed', ...
            'Hardware, grid size, run times, and memory: improves run-time estimates', ...
            ['Rounded summary statistics: brain intensity, temperature rise, targeting error, ' ...
             'and skull properties (no maps, coordinates, or targets)']});
end
