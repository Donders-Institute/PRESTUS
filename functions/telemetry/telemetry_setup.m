function telemetry_setup(parent_fig)
% TELEMETRY_SETUP  Show the opt-in prompt for anonymous usage statistics.
%
% Called at the start of prestus_pipeline, and (with a figure handle) at
% GUI startup. Skipped if the user has already made an explicit decision
% (i.e. ~/.prestus/telemetry.json exists).
%
% Decision order:
%   1. Consent file (~/.prestus/telemetry.json) exists -> the recorded
%      decision is used, no prompt.
%   2. parent_fig given (called from the GUI) -> a modal dialog is shown,
%      parented to that figure.
%   3. Interactive console session -> the y/n prompt below is shown until
%      the user answers (blocks until a decision is made).
%   4. matlab -batch run locally with a display attached, i.e. not a
%      deployed app or a SLURM/PBS job (see TELEMETRY_HAS_DISPLAY,
%      TELEMETRY_IS_HEADLESS_CLUSTER) -> a modal dialog is shown, since
%      there is no console to answer input() on.
%   5. Otherwise (deployed app, SLURM/PBS job, or batch with no display at
%      all) -> only the notice is printed; nothing is recorded or sent.
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
% See also: TRACK_USAGE, TELEMETRY_SETUP_RESET, TELEMETRY_HAS_DISPLAY

    arguments
        parent_fig = []
    end

    cfg_file = fullfile(prefdir_prestus(), 'telemetry.json');
    if isfile(cfg_file)
        return   % explicit decision already recorded
    end

    if ~isempty(parent_fig)
        prompt_via_dialog(cfg_file, parent_fig);
        return
    end

    fprintf('\n========================================\n');
    fprintf('PRESTUS – TELEMETRY\n');
    fprintf('========================================\n');
    fprintf('%s\n', notice_body(cfg_file));

    if is_interactive_session()
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

    if ~telemetry_is_headless_cluster() && telemetry_has_display()
        % matlab -batch run locally with a display: no console to answer
        % input() on, but a dialog can still be shown.
        prompt_via_dialog(cfg_file, []);
        return
    end

    fprintf(['Non-interactive session: no decision recorded, no data sent.\n' ...
             'Run telemetry_set_consent(true|false) once to stop this notice.\n\n']);
end

% -------------------------------------------------------------------------
function tf = is_interactive_session()
    tf = ~telemetry_is_headless_cluster() && ~batchStartupOptionUsed;
end

% -------------------------------------------------------------------------
function prompt_via_dialog(cfg_file, parent_fig)
% Show a modal Allow/Decline/Ask-me-later dialog in place of console
% input(). Used for GUI sessions (parent_fig given) and for local -batch
% sessions with a display but no interactive stdin (parent_fig empty, a
% throwaway invisible uifigure is created to host the dialog).
    owns_fig = isempty(parent_fig);
    if owns_fig
        parent_fig = uifigure('Visible', 'off', 'Name', 'PRESTUS');
    end
    cleanup_fig = onCleanup(@() delete_if_owned(owns_fig, parent_fig)); %#ok<NASGU>

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
function delete_if_owned(owns_fig, fig)
    if owns_fig && isvalid(fig)
        delete(fig);
    end
end

% -------------------------------------------------------------------------
function write_decision(cfg_file, opted)
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
        '  - Run duration, success/failure, and error type (no message text)\n', ...
        '  - A random local ID (not linked to you or your machine)\n\n', ...
        'What is NEVER collected:\n', ...
        '  - Subject IDs, file paths, coordinates, or any free-text\n', ...
        '  - IP addresses or hostnames\n', ...
        '  - Acoustic pressure values, simulation results, or tissue property values\n\n', ...
        'Full details: https://github.com/Donders-Institute/PRESTUS/blob/main/doc/doc_telemetry.md\n\n', ...
        'Your answer is stored in:\n', ...
        '  %s\n', ...
        'To opt out later, set opt_in=false in that file, or run telemetry_set_consent(false).'], ...
        cfg_file);
end
