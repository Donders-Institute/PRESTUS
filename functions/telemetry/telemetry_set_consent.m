function telemetry_set_consent(opt_in)
% TELEMETRY_SET_CONSENT  Record a telemetry opt-in or opt-out decision.
%
% Call this once before running the pipeline — from any MATLAB session or
% script — or any time you want to change your preference.
%
% Use as:
%   telemetry_set_consent(true)   % opt in
%   telemetry_set_consent(false)  % opt out
%
% The decision is stored in ~/.prestus/telemetry.json and is read by
% every subsequent pipeline run. To revert to the undecided state, call
% telemetry_setup_reset().
%
% See also: TELEMETRY_SETUP, TELEMETRY_SETUP_RESET

    arguments
        opt_in (1,1) logical
    end

    cfg_file = fullfile(prefdir_prestus(), 'telemetry.json');

    cfg             = struct();
    cfg.opt_in      = opt_in;
    cfg.decided_on  = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd'));
    cfg.prestus_ver = prestus_version();

    fid = fopen(cfg_file, 'w');
    fprintf(fid, '%s\n', jsonencode(cfg));
    fclose(fid);

    if opt_in
        fprintf('Telemetry enabled. Opt out any time with telemetry_set_consent(false).\n');
    else
        fprintf('Telemetry disabled. Opt in any time with telemetry_set_consent(true).\n');
    end
end
