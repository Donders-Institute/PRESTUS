function telemetry_setup(parent_fig)
% TELEMETRY_SETUP  Show the opt-in prompt for anonymous usage statistics.
%
% Called at the start of prestus_pipeline, and (with a figure handle) at
% GUI startup. Skipped if the user has already made an explicit decision
% (i.e. ~/.prestus/telemetry.json exists) or PRESTUS_TELEMETRY is set.
% Never waits for input in a non-interactive session.
%
% Decision order:
%   1. PRESTUS_TELEMETRY set to off/0/false (see TELEMETRY_ENV_OVERRIDE)
%      -> nothing is sent from this process, whatever was recorded; one
%      line is printed, nothing is prompted or recorded.
%   2. Consent file (~/.prestus/telemetry.json) exists -> the recorded
%      decision is used, no prompt. A recorded opt-out is never overridden:
%      PRESTUS_TELEMETRY=on cannot turn sending back on.
%   3. PRESTUS_TELEMETRY set to on/1/true and no decision recorded -> sent
%      from this process; one line is printed, nothing is recorded.
%   4. parent_fig given (called from the GUI) -> a modal dialog is shown,
%      parented to that figure.
%   5. Interactive console session -> the y/n prompt below is shown until
%      the user answers (blocks until a decision is made).
%   6. Otherwise (non-interactive: deployed app, SLURM/PBS job,
%      `matlab -batch` with or without a display, or no desktop and no
%      display; see TELEMETRY_IS_NONINTERACTIVE) -> only the notice is
%      printed; nothing is recorded or sent, and nothing waits for input.
%
% To record a decision without any prompt (e.g. before an HPC job), call
% telemetry_set_consent().
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
% Use as:
%   telemetry_setup()
%   telemetry_setup(fig)   % GUI context: show a dialog parented to fig
%
% Input:
%   parent_fig - (optional) a uifigure handle. When given, the opt-in
%                decision (if not already made) is always asked via a
%                modal dialog rather than console input(), since a GUI
%                session's Command Window prompt is easy to miss and the
%                simulation itself may run on a background worker with no
%                console at all.
%
% See also: TRACK_USAGE, TELEMETRY_SETUP_RESET, TELEMETRY_ENV_OVERRIDE,
%           TELEMETRY_IS_NONINTERACTIVE

    arguments
        parent_fig = []
    end

    env = telemetry_env_override();
    if strcmp(env, 'off')
        fprintf('Telemetry off for this process (PRESTUS_TELEMETRY): no data sent.\n');
        return
    end

    cfg_file = fullfile(prefdir_prestus(), 'telemetry.json');
    if isfile(cfg_file)
        % Explicit decision already recorded; an opt-out is never overridden.
        if strcmp(env, 'on')
            fprintf(['PRESTUS_TELEMETRY=on ignored: the decision recorded in %s applies ' ...
                     '(change it with telemetry_set_consent).\n'], cfg_file);
        end
        return
    end

    if strcmp(env, 'on')
        fprintf(['Telemetry on for this process (PRESTUS_TELEMETRY): usage statistics ' ...
                 'are sent; no decision recorded.\n']);
        return
    end

    if ~isempty(parent_fig)
        prompt_via_dialog(cfg_file, parent_fig);
        return
    end

    fprintf('\n========================================\n');
    fprintf('PRESTUS – TELEMETRY\n');
    fprintf('========================================\n');
    fprintf('%s\n', notice_body(cfg_file));

    if ~telemetry_is_noninteractive()
        reply = '';
        while ~ismember(reply, {'y', 'n'})
            reply = lower(strtrim(input('Allow anonymous usage statistics? [y/n]: ', 's')));
            if ~ismember(reply, {'y', 'n'})
                fprintf('Please answer y or n.\n');
            end
        end
        write_decision(cfg_file, strcmp(reply, 'y'));
        return
    end

    fprintf(['Non-interactive session: no decision recorded, no data sent.\n' ...
             'Run telemetry_set_consent(true|false) once, or set PRESTUS_TELEMETRY=off|on,\n' ...
             'to stop this notice.\n\n']);
end

% -------------------------------------------------------------------------
function prompt_via_dialog(cfg_file, parent_fig)
% Show a modal Allow/Decline/Ask-me-later dialog parented to the GUI figure
% in place of console input(). Used for GUI sessions only; a session with
% no figure never gets a dialog (see TELEMETRY_IS_NONINTERACTIVE).

    % uiconfirm renders a cell array as one line per element; a plain char
    % array with embedded \n is not guaranteed to wrap the same way.
    lines = strsplit(notice_body(cfg_file), '\n');
    selection = uiconfirm(parent_fig, lines, 'PRESTUS - Telemetry', ...
        'Options',       {'Allow', 'Decline', 'Ask me later'}, ...
        'DefaultOption', 3, ...
        'CancelOption',  3, ...
        'Icon',          'question');

    switch selection
        case 'Allow'
            write_decision(cfg_file, true);
        case 'Decline'
            write_decision(cfg_file, false);
        otherwise
            % "Ask me later" (or dialog dismissed): no decision recorded,
            % the prompt reappears on the next call.
    end
end

% -------------------------------------------------------------------------
function write_decision(cfg_file, opted)
    cfg             = struct();
    cfg.opt_in      = opted;
    cfg.decided_on  = char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd'));
    cfg.prestus_ver = prestus_version();
    cfg.notice      = telemetry_notice_version();

    fid = fopen(cfg_file, 'w');
    fprintf(fid, '%s\n', jsonencode(cfg));
    fclose(fid);

    if opted
        fprintf('Thank you! Telemetry enabled. Opt out any time by setting opt_in=false in:\n  %s\n\n', cfg_file);
    else
        fprintf('Understood. No data will be sent. You can opt in later by setting opt_in=true in:\n  %s\n\n', cfg_file);
    end
end

% -------------------------------------------------------------------------
function txt = notice_body(cfg_file)
    txt = sprintf([...
        'PRESTUS can optionally collect anonymous usage statistics to help\n', ...
        'the developers understand which features are used and on which\n', ...
        'platforms. This helps us prioritise development and find bugs.\n\n', ...
        'What IS collected (examples):\n', ...
        '  - PRESTUS version, MATLAB version, k-Wave version, OS platform\n', ...
        '  - CPU model, core count, RAM, GPU model, installed toolboxes, HPC name\n', ...
        '  - Which pipeline modules are enabled\n', ...
        '  - Transducer type, frequency, and element count\n', ...
        '  - Simulation medium, layer names, whether pCT is used\n', ...
        '  - Run start, duration, success/failure, error type and failing stage\n', ...
        '    (no message text)\n', ...
        '  - Requested job resources, queue wait, and solver settings\n', ...
        '  - Simulation grid size, time steps, stage durations and peak RAM\n', ...
        '  - Median sound speed, density and attenuation of the modelled skull\n', ...
        '    (all skull voxels, and on the beam path)\n', ...
        '  - A random local ID (not linked to you or your machine)\n\n', ...
        'What is NEVER collected:\n', ...
        '  - Subject IDs, file paths, coordinates, or any free-text\n', ...
        '  - IP addresses or hostnames\n', ...
        '  - Acoustic pressure values, simulation results, or tissue property maps\n\n', ...
        'Full details: https://github.com/Donders-Institute/PRESTUS/blob/main/doc/doc_telemetry.md\n\n', ...
        'Your answer is stored in:\n', ...
        '  %s\n', ...
        'To opt out later, set opt_in=false in that file, or run telemetry_set_consent(false).'], ...
        cfg_file);
end
