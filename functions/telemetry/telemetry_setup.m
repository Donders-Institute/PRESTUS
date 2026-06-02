function telemetry_setup()
% TELEMETRY_SETUP  Show the opt-in prompt for anonymous usage statistics.
%
% Called at the start of prestus_pipeline. Skipped if the user has already
% made an explicit decision (i.e. ~/.prestus/telemetry.json exists).
%
% The prompt blocks until the user answers 'y' or 'n'. To avoid being
% prompted at all (e.g. before an HPC job), call telemetry_set_consent()
% once beforehand to record the decision.
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
        'Answer below. This prompt will loop until you choose.\n', ...
        'To opt out later, set opt_in=false in:\n', ...
        '  %s\n\n'], cfg_file);

    reply = '';
    while ~ismember(reply, {'y', 'n'})
        reply = lower(strtrim(input('Allow anonymous usage statistics? [y/n]: ', 's')));
        if ~ismember(reply, {'y', 'n'})
            fprintf('Please answer y or n.\n');
        end
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
