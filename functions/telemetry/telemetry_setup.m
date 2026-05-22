function telemetry_setup()
% TELEMETRY_SETUP  Show the opt-in prompt for anonymous usage statistics.
%
% Called at the start of prestus_pipeline. Skipped if the user has already
% made an explicit decision (i.e. ~/.prestus/telemetry.json exists).
%
% In interactive sessions the prompt is shown every run until the user
% explicitly answers 'y' or 'n'. No decision is recorded and no data is
% ever sent without an explicit answer.
%
% In non-interactive sessions (HPC batch jobs) the full notice is printed
% to stdout on every run but no decision is recorded and no data is sent.
% To opt in from an HPC environment, either:
%   1. Run PRESTUS once interactively and answer 'y', or
%   2. Create ~/.prestus/telemetry.json with {"opt_in":true} manually.
%
% The consent record written to ~/.prestus/telemetry.json contains only:
%   opt_in         - true/false
%   decided_on     - ISO-8601 date (UTC, no time)
%   prestus_ver    - PRESTUS git hash at decision time
%
% Users can withdraw consent at any time by editing that file and setting
% opt_in to false, or by running:
%   telemetry_setup_reset()   % deletes the file so the prompt reappears
%
% See also: TRACK_USAGE, TELEMETRY_SETUP_RESET

    cfg_file = fullfile(prefdir_prestus(), 'telemetry.json');
    if isfile(cfg_file)
        return   % explicit decision already recorded
    end

    fprintf('\n========================================\n');
    fprintf('PRESTUS – TELEMETRY\n');
    fprintf('========================================\n');
    fprintf([...
        'PRESTUS can optionally collect anonymous usage statistics to help\n', ...
        'the developers understand which features are used and on which\n', ...
        'platforms. This helps us prioritise development and find bugs.\n\n', ...
        'What IS collected (examples):\n', ...
        '  - PRESTUS version, MATLAB version, k-Wave version, OS platform\n', ...
        '  - Which pipeline modules are enabled\n', ...
        '  - Transducer type, frequency, and element count\n', ...
        '  - Simulation medium, layer names, whether pCT is used\n', ...
        '  - Run duration, success/failure, and error type (no message text)\n', ...
        '  - A random local ID (not linked to you or your machine)\n\n', ...
        'What is NEVER collected:\n', ...
        '  - Subject IDs, file paths, coordinates, or any free-text\n', ...
        '  - IP addresses or hostnames\n', ...
        '  - Acoustic pressure values, simulation results, or tissue property values\n\n', ...
        'Full details: https://github.com/Donders-Institute/PRESTUS/blob/main/documentation/doc_telemetry.md\n\n', ...
        'To silence this message, answer below or edit:\n', ...
        '  %s\n\n'], cfg_file);

    if ~feature('ShowFigureWindows')
        % Non-interactive (HPC batch job): print notice but do not record a
        % decision and do not send data. The message will reappear on every
        % run until the user opts in or out interactively.
        fprintf(['Non-interactive session: no data will be sent.\n', ...
                 'To opt in, run PRESTUS once interactively, or create:\n', ...
                 '  %s\n', ...
                 'containing: {"opt_in":true,"decided_on":"YYYY-MM-DD","prestus_ver":""}\n\n'], ...
                cfg_file);
        return   % no decision recorded – prompt will reappear next run
    end

    reply = input('Allow anonymous usage statistics? [y/n]: ', 's');
    reply = lower(strtrim(reply));

    if ~ismember(reply, {'y', 'n'})
        fprintf('No answer recorded – you will be asked again on the next run.\n\n');
        return   % no decision recorded – prompt will reappear next run
    end

    opted = strcmp(reply, 'y');

    cfg             = struct();
    cfg.opt_in      = opted;
    cfg.decided_on  = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd'));
    cfg.prestus_ver = prestus_version();

    fid = fopen(cfg_file, 'w');
    fprintf(fid, '%s\n', jsonencode(cfg));
    fclose(fid);

    if opted
        fprintf('Thank you! Telemetry enabled. Opt out any time by setting opt_in=false in:\n  %s\n\n', cfg_file);
    else
        fprintf('Understood. No data will be sent. You can opt in later by setting opt_in=true in:\n  %s\n\n', cfg_file);
    end
end
