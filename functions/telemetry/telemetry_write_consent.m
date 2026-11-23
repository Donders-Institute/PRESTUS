function cfg_file = telemetry_write_consent(levels)
% TELEMETRY_WRITE_CONSENT  Write the consent record ~/.prestus/telemetry.json.
%
% The record holds only:
%   opt_in       - true when any level is shared
%   levels       - the shared levels (TELEMETRY_LEVEL_INFO names)
%   decided_on   - ISO-8601 date (UTC, no time)
%   prestus_ver  - PRESTUS git hash at decision time
%   notice       - notice version the decision was made under
% The Nexus writes the same record (prestus_nexus/telemetry.py).
%
% Use as:
%   cfg_file = telemetry_write_consent({'basic', 'outcomes'})
%   cfg_file = telemetry_write_consent({})       % share nothing
%
% See also: TELEMETRY_SET_CONSENT, TELEMETRY_SETUP, TELEMETRY_LEVELS

    cfg_file = fullfile(prefdir_prestus(), 'telemetry.json');
    cfg             = struct();
    cfg.opt_in      = ~isempty(levels);
    cfg.levels      = {};
    cfg.decided_on  = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd'));
    cfg.prestus_ver = prestus_version();
    cfg.notice      = telemetry_notice_version();
    txt = jsonencode(cfg);
    % jsonencode writes a 1-element cell as a plain string; always write a list
    levels = cellstr(levels);
    txt = strrep(txt, '"levels":[]', ['"levels":' list_json(levels)]);

    fid = fopen(cfg_file, 'w');
    fprintf(fid, '%s\n', txt);
    fclose(fid);
end

% -------------------------------------------------------------------------
function s = list_json(c)
    if isempty(c)
        s = '[]';
    else
        s = ['["' strjoin(c, '","') '"]'];
    end
end
