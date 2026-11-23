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
%      levels are used, no prompt. A recorded opt-out is never overridden:
%      PRESTUS_TELEMETRY cannot turn sending back on.
%   3. PRESTUS_TELEMETRY set to on or a list of levels, and no decision
%      recorded -> those levels are sent from this process; one line is
%      printed, nothing is recorded.
%   4. parent_fig given (called from the GUI) -> a modal dialog is shown,
%      parented to that figure.
%   5. Interactive console session -> the prompt below asks which levels
%      to share (Enter: all) until the user answers.
%   6. Otherwise (non-interactive: deployed app, SLURM/PBS job,
%      `matlab -batch` with or without a display, or no desktop and no
%      display; see TELEMETRY_IS_NONINTERACTIVE) -> only the notice is
%      printed; nothing is recorded or sent, and nothing waits for input.
%
% Both prompts offer "Share all (recommended)" with the levels
% (TELEMETRY_LEVEL_INFO) indented below it, all ticked; unticking every
% level shares nothing. To record a choice without any prompt (e.g. before
% an HPC job), call telemetry_set_consent(). The record is written by
% TELEMETRY_WRITE_CONSENT. Users can change it at any time with
% telemetry_set_consent, or forget it with telemetry_setup_reset().
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
% See also: TRACK_USAGE, TELEMETRY_LEVEL_INFO, TELEMETRY_SETUP_RESET, TELEMETRY_ENV_OVERRIDE,
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
            fprintf(['PRESTUS_TELEMETRY ignored: the decision recorded in %s applies ' ...
                     '(change it with telemetry_set_consent).\n'], cfg_file);
        end
        % An opt-in given under an earlier notice keeps what it agreed to; the
        % levels are offered once, in an interactive session, as a new choice.
        if ~allowed_before_levels(cfg_file) || (isempty(parent_fig) && telemetry_is_noninteractive())
            return
        end
        fprintf(['\nUsage statistics now come in levels, with a new one: summary ' ...
                 'outcomes (anonymous). You keep sharing what you allowed until you choose.\n']);
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
    fprintf('PRESTUS – USAGE STATISTICS\n');
    fprintf('========================================\n');
    fprintf('%s\n', notice_body(cfg_file));

    if ~telemetry_is_noninteractive()
        info = telemetry_level_info();
        fprintf('  [x] Share all (recommended)\n');
        for k = 1:numel(info)
            fprintf('      [x] %d  %s: %s\n', k, info(k).label, info(k).description);
        end
        fprintf('%s\n\n', never_sent());
        while true
            reply = strtrim(input(sprintf(['Press Enter to share all, or type the numbers ' ...
                'to share (e.g. 1 2); 0 unticks all: ']), 's'));
            [levels, ok] = parse_reply(reply, info);
            if ok, break, end
            fprintf('Please press Enter, or type numbers from 0 to %d.\n', numel(info));
        end
        record(levels);
        return
    end

    fprintf('%s\n\n', never_sent());
    fprintf(['Non-interactive session: no decision recorded, no data sent.\n' ...
             'Run telemetry_set_consent(true|false|levels) once, or set PRESTUS_TELEMETRY,\n' ...
             'to stop this notice.\n\n']);
end

% -------------------------------------------------------------------------
function tf = allowed_before_levels(cfg_file)
% An opt-in recorded under a notice before levels existed (notice < 3).
    tf = false;
    if ~isfile(cfg_file)
        return
    end
    try
        cfg = jsondecode(fileread(cfg_file));
    catch
        % Unreadable or corrupt consent file: not an earlier opt-in; nothing is sent and nothing is asked (TELEMETRY_LEVELS).
        return
    end
    tf = isstruct(cfg) && isfield(cfg, 'opt_in') && isequal(cfg.opt_in, true) && ...
         ~isfield(cfg, 'levels') && (~isfield(cfg, 'notice') || cfg.notice < 3);
end

% -------------------------------------------------------------------------
function [levels, ok] = parse_reply(reply, info)
% '' -> all levels; '0' -> none; '1 3' or '1,3' -> those levels.
    names = {info.name};
    ok = true;
    if isempty(reply)
        levels = names;
        return
    end
    nums = str2double(regexp(reply, '\d+', 'match'));
    if isempty(nums) || ~isempty(regexprep(reply, '[\d\s,]', '')) || any(nums > numel(names))
        levels = {};
        ok = false;
        return
    end
    levels = names(ismember(1:numel(names), nums));
end

% -------------------------------------------------------------------------
function prompt_via_dialog(cfg_file, parent_fig)
% A modal dialog parented to the GUI figure in place of console input():
% "Share all (recommended)" in bold with the levels indented below it, all
% ticked; Save records the ticked levels (none ticked: nothing is shared),
% Ask me later (or closing) records nothing.
    info = telemetry_level_info();
    n = numel(info);
    h = 236 + 46 * n;
    d = uifigure('Name', 'PRESTUS - Usage statistics', 'WindowStyle', 'modal', ...
                 'Position', [0 0 560 h], 'Resize', 'off');
    movegui(d, 'center');
    if ~isempty(parent_fig) && isvalid(parent_fig)
        d.Position(1:2) = parent_fig.Position(1:2) + (parent_fig.Position(3:4) - d.Position(3:4)) / 2;
    end
    uilabel(d, 'Position', [20 h-60 520 44], 'WordWrap', 'on', 'Text', ...
        ['Share anonymous usage statistics with the PRESTUS developers? Each run sends ' ...
         'one record; each level below adds its fields.']);
    all_box = uicheckbox(d, 'Text', 'Share all (recommended)', 'Value', true, ...
                         'FontWeight', 'bold', 'Position', [20 h-92 400 22]);
    boxes = gobjects(1, n);
    for k = 1:n
        y = h - 92 - 46 * k;
        boxes(k) = uicheckbox(d, 'Text', info(k).label, 'Value', true, ...
                              'Position', [44 y+18 480 22]);
        uilabel(d, 'Text', info(k).description, 'FontColor', [0.4 0.4 0.4], ...
                'Position', [66 y 470 20]);
    end
    uilabel(d, 'Position', [20 54 520 80], 'WordWrap', 'on', 'Text', strrep(never_sent(), newline, ' '));
    all_box.ValueChangedFcn = @(src, ~) set(boxes, 'Value', src.Value);
    for k = 1:n
        boxes(k).ValueChangedFcn = @(~, ~) set(all_box, 'Value', all([boxes.Value]));
    end
    uibutton(d, 'Text', 'Ask me later', 'Position', [300 16 110 28], ...
             'ButtonPushedFcn', @(~, ~) uiresume(d));
    uibutton(d, 'Text', 'Save', 'Position', [420 16 110 28], ...
             'ButtonPushedFcn', @(~, ~) save_and_close());
    decided = false;
    picked = {};
    uiwait(d);
    if isvalid(d), delete(d), end
    if decided
        record(picked);
    end     % else Ask me later, or closed: no decision recorded

    function save_and_close()
        picked = {info([boxes.Value]).name};
        decided = true;
        uiresume(d);
    end
end

% -------------------------------------------------------------------------
function record(levels)
    cfg_file = telemetry_write_consent(levels);
    if isempty(levels)
        fprintf('Understood: nothing is shared. Change with telemetry_set_consent; your choice is in:\n  %s\n\n', cfg_file);
    else
        fprintf('Thank you! Sharing %s. Change with telemetry_set_consent; your choice is in:\n  %s\n\n', ...
                strjoin(levels, ', '), cfg_file);
    end
end

% -------------------------------------------------------------------------
function txt = notice_body(cfg_file)
    txt = sprintf([...
        'PRESTUS can share anonymous usage statistics with its developers: which\n', ...
        'features are used, how long runs take, and what they simulate. This\n', ...
        'guides development, improves run-time estimates, and lets us publish\n', ...
        'statistics across sites (for example, skull transmission).\n\n', ...
        'Every record holds: a random ID (not linked to you or your machine),\n', ...
        'the PRESTUS version, the run type, and whether the run succeeded.\n', ...
        'Each level you share adds (details in doc_telemetry.md):\n', ...
        '  1 Usage and errors: versions, platform, modules, medium, transducer\n', ...
        '    model and frequency, error type and failing stage (no message text)\n', ...
        '  2 Performance: CPU, GPU, RAM, requested resources, grid size, stage\n', ...
        '    times, and peak memory\n', ...
        '  3 Summary outcomes (anonymous), rounded: brain intensity and its ratio to\n', ...
        '    free water, brain mechanical index, targeting error, brain and skull\n', ...
        '    temperature rise with duty cycle and duration, and median skull sound\n', ...
        '    speed, density, and attenuation\n\n', ...
        'Full details: https://github.com/Donders-Institute/PRESTUS/blob/main/doc/doc_telemetry.md\n', ...
        'Your choice is stored in %s'], cfg_file);
end

% -------------------------------------------------------------------------
function txt = never_sent()
    txt = sprintf(['Never stored in the usage statistics: subject IDs, file paths, coordinates,\n' ...
                   'target names, free text, IP addresses, host names, or any map or image of\n' ...
                   'the simulation. Like any web service, the hosting provider (Supabase)\n' ...
                   'keeps a request log with the sender''s IP address and approximate location\n' ...
                   'for 1 day; it is not linked to the statistics.']);
end
