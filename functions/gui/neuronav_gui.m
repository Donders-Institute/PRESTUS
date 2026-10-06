function fig = neuronav_gui(init)
%NEURONAV_GUI  Stand-alone GUI for Localite marker ingestion (neuronav Stage 1).
%
%   Front-end to NEURONAV_INGEST_MARKERS. Works in two steps:
%
%     1. Detect series — runs the ingest file selection and series
%        statistics into a temporary folder and shows every detected series
%        with its QC statistics (pulse count, position/rotation deviation),
%        a 3-D plot of all pulse positions (one dot per pulse), per-pulse
%        deviation traces, the raw marker files found (and which one was
%        used), and the log.
%     2. Write outputs — calls NEURONAV_INGEST_MARKERS again with a
%        target_map built from the series table (included rows, target
%        names, colours) and writes
%          <output root>/<sub>/<ses>/InstrumentMarker_<sub>_<ses>.xml / .json
%
%   Coil table: leave it empty for the single-coil path (newest
%   non-duplicate TriggerMarkers_Coil0 file, GUMMarkers fallback; one row
%   per detected series). Add one row per coil for multi-coil setups; the
%   ingest then reads one file per coil and yields one series per coil
%   (parameters.neuronav.coil_map). A coil that records several blocks
%   takes comma-separated target names (one coil_map sequence each).
%
%   Batch window (▦ Batch…): scan, label map and NEURONAV_INGEST_BATCH over
%   many sessions. Its ◎ MNI view… window maps the summary CSV to MNI space
%   (NEURONAV_SUMMARY_TO_MNI; needs path.seg with the SimNIBS toMNI warps)
%   and plots transducer, focus and beam axis of every subject/session over
%   an MNI head outline, with tick filters per subject, session and target.
%
%   Usage:
%     neuronav_gui()
%     neuronav_gui(parameters)  % pre-fill from a PRESTUS config struct
%     neuronav_gui('config.yaml')
%     fig = neuronav_gui(...)    % return the figure handle
%     p = fig.UserData.settings()          % the parameters detect/batch use
%     fig.UserData.save_settings('x.yaml') % save them (as ⬇ Save settings…)
%
%   Fields read from a config (all optional):
%     path.localite_raw (falls back to path.localite), path.localite_post,
%     placement.localite.loc_folder / voxel_size / expected_segment_length /
%     session, io.overwrite_files, neuronav.coil_map (.coil, .target or
%     .sequences(:).target)
%
%   Requirements: MATLAB R2023b, PRESTUS functions on the MATLAB path.
%
% See also: NEURONAV_INGEST_MARKERS, NEURONAV_INGEST_BATCH,
%           NEURONAV_SUMMARY_TO_MNI, NEURONAV_BUILD_TARGET_MAP, PRESTUS_GUI

    if nargin < 1, init = []; end

    if exist('neuronav_ingest_markers', 'file') ~= 2
        error('neuronav_gui:path', ...
            'PRESTUS functions are not on the MATLAB path (neuronav_ingest_markers not found).');
    end

%% ── State ─────────────────────────────────────────────────────────────────

st = sty();

S = struct();
S.results     = {};      % per-series stats from the detect step
S.markertype  = '';
S.sources     = {};      % raw file name(s) used by the detect step
S.candidates  = [];      % named planned targets (neuronav_read_label_candidates)
S.suggest     = [];      % per-series naming result (neuronav_suggest_target_names)
S.D           = [];      % detection result (neuronav_detect_series)
S.N           = [];      % naming result (neuronav_name_series)
S.base        = struct();% loaded config; GUI fields override it (see collect_settings)
S.target_labels = [];    % label map (neuronav.target_labels), edited in the Batch window
S.batch_write = 'all';   % neuronav.batch_write
B = [];                  % batch window handles (open_batch)
MV = [];                 % MNI view window handles and data (open_mni)
S.coil_ids    = [];      % coil index per series (coil mode) / [] (single)
S.snapshot    = [];      % input settings at detect time
S.xml_path    = '';
S.json_path   = '';

%% ── Figure ────────────────────────────────────────────────────────────────

fig = uifigure( ...
    'Name',     'PRESTUS · Localite ingest', ...
    'Tag',      'neuronav_gui', ...
    'Position', [80 60 1220 780], ...
    'Color',    st.bg_fig, ...
    'Visible',  'off');

root = uigridlayout(fig, [3 1], ...
    'RowHeight', {64, 3, '1x'}, 'ColumnWidth', {'1x'}, ...
    'Padding', [0 0 0 0], 'RowSpacing', 0, 'BackgroundColor', st.bg_fig);

% Header bar: logo (white background, so the bar is light) + title
hdr = uigridlayout(root, [2 3], ...
    'RowHeight', {'1x', '1x'}, 'ColumnWidth', {86, '1x', 'fit'}, ...
    'Padding', [12 4 12 4], 'RowSpacing', 0, 'ColumnSpacing', 12, ...
    'BackgroundColor', [1 1 1]);
logo_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'report', 'assets', 'logo_PRESTUS.png');
if isfile(logo_file)
    im = uiimage(hdr, 'ImageSource', logo_file, 'ScaleMethod', 'fit', ...
        'Tooltip', 'PRESTUS');
    place(im, [1 2], 1);
end
o = uilabel(hdr, 'Text', 'Localite marker ingestion', 'FontName', st.font, ...
    'FontSize', st.fs_xl, 'FontWeight', 'bold', 'FontColor', st.text, ...
    'VerticalAlignment', 'bottom');
place(o, 1, 2);
o = uilabel(hdr, 'Text', ['neuronav Stage 1 — raw TriggerMarkers / GUMMarkers → ' ...
    'InstrumentMarker XML + QC sidecar'], 'FontName', st.font, 'FontSize', st.fs_sm, ...
    'FontColor', st.text_sub, 'VerticalAlignment', 'top');
place(o, 2, 2);
gl_hb = uigridlayout(hdr, [3 3], 'RowHeight', {'1x', 28, '1x'}, ...
    'ColumnWidth', {'fit', 'fit', 'fit'}, 'Padding', [0 0 0 0], 'ColumnSpacing', 6, ...
    'BackgroundColor', [1 1 1]);
place(gl_hb, [1 2], 3);
b = uibutton(gl_hb, 'Text', '⬆ Load settings…', 'FontName', st.font, ...
    'FontSize', st.fs_sm, 'BackgroundColor', st.btn_lt, ...
    'Tooltip', 'Read path / placement.localite / neuronav settings from a PRESTUS YAML config', ...
    'ButtonPushedFcn', @(~,~) cb_load_config());
place(b, 2, 1);
b = uibutton(gl_hb, 'Text', '⬇ Save settings…', 'FontName', st.font, ...
    'FontSize', st.fs_sm, 'BackgroundColor', st.btn_lt, 'Tag', 'btn_save_settings', ...
    'Tooltip', 'Write these settings as a YAML config, usable by load_parameters in scripts', ...
    'ButtonPushedFcn', @(~,~) cb_save_config());
place(b, 2, 2);
b = uibutton(gl_hb, 'Text', '▦ Batch…', 'FontName', st.font, 'FontWeight', 'bold', ...
    'FontSize', st.fs_sm, 'BackgroundColor', st.accent, 'FontColor', [1 1 1], 'Tag', 'btn_batch', ...
    'Tooltip', 'Scan all subjects/sessions, map target labels and ingest many sessions', ...
    'ButtonPushedFcn', @(~,~) open_batch());
place(b, 2, 3);
uipanel(root, 'BackgroundColor', st.accent, 'BorderType', 'none');

body = uigridlayout(root, [1 2], ...
    'ColumnWidth', {400, '1x'}, 'Padding', [0 0 0 0], 'ColumnSpacing', 0, ...
    'BackgroundColor', st.bg_fig);

%% ── Left panel: inputs ───────────────────────────────────────────────────

left = uigridlayout(body, [2 1], 'RowHeight', {'1x', 'fit'}, ...
    'Padding', [0 0 0 0], 'RowSpacing', 0, 'BackgroundColor', st.bg_panel);
L = uigridlayout(left, [20 4], ...
    'RowHeight', [repmat({26}, 1, 18), {110}, {26}], ...
    'ColumnWidth', {130, '1x', 'fit', 34}, ...
    'Padding', [16 14 12 14], 'RowSpacing', 4, 'ColumnSpacing', 6, ...
    'BackgroundColor', st.bg_panel, 'Scrollable', 'on');

sec(L, 1, 'Input');
lbl(L, 2, 'Raw Localite root');
h.raw = uieditfield(L, 'text', 'FontName', st.font, 'FontSize', st.fs, ...
    'Tooltip', 'path.localite_raw — contains <sub>/<ses>/<loc folder>/…', ...
    'ValueChangedFcn', @(~,~) cb_rescan_subjects());
place(h.raw, 2, [2 3]);
brw(L, 2, h.raw, @cb_rescan_subjects);

lbl(L, 3, 'Output root');
h.out = uieditfield(L, 'text', 'FontName', st.font, 'FontSize', st.fs, ...
    'Tooltip', 'path.localite_post — outputs go to <root>/<sub>/<ses>/');
place(h.out, 3, [2 3]);
brw(L, 3, h.out, []);

lbl(L, 4, 'Subject');
h.sub = uidropdown(L, 'Items', {}, 'Editable', 'on', 'FontName', st.font, ...
    'FontSize', st.fs, 'Placeholder', 'sub-001', ...
    'ValueChangedFcn', @(~,~) cb_rescan_sessions());
place(h.sub, 4, [2 3]);
b = uibutton(L, 'Text', '↻', 'FontSize', st.fs, 'Tooltip', 'Rescan raw root', ...
    'BackgroundColor', st.btn_lt, 'ButtonPushedFcn', @(~,~) cb_rescan_subjects());
place(b, 4, 4);

lbl(L, 5, 'Session');
h.ses = uidropdown(L, 'Items', {}, 'Editable', 'on', 'FontName', st.font, ...
    'FontSize', st.fs, 'Placeholder', 'ses-01', ...
    'Tooltip', 'ses-XX, or a number (2 → ses-02)', ...
    'ValueChangedFcn', @(~,~) cb_inputs_changed());
place(h.ses, 5, [2 3]);

lbl(L, 6, 'Localite folder');
h.loc = uieditfield(L, 'text', 'Value', 'loc', 'FontName', st.font, 'FontSize', st.fs, ...
    'Placeholder', 'none', ...
    'Tooltip', 'placement.localite.loc_folder — sub-folder of <sub>/<ses>/; empty = export directly under <ses>/', ...
    'ValueChangedFcn', @(~,~) cb_inputs_changed());
place(h.loc, 6, [2 3]);

sec(L, 7, 'Series detection');
lbl(L, 8, 'Voxel size');
h.vox = uieditfield(L, 'numeric', 'Value', 0.9, 'Limits', [0 Inf], ...
    'LowerLimitInclusive', 'off', 'FontName', st.font, 'FontSize', st.fs, ...
    'Tooltip', 'placement.localite.voxel_size — only affects position_dev_vox in the returned stats');
place(h.vox, 8, 2);
unit(L, 8, 'mm');

lbl(L, 9, 'Pulses per series');
h.seglen = uieditfield(L, 'text', 'Value', '80', 'FontName', st.font, 'FontSize', st.fs, ...
    'Placeholder', 'accept all', ...
    'Tooltip', 'placement.localite.expected_segment_length — blank accepts every series');
place(h.seglen, 9, 2);
note(L, 10, 'TriggerMarkers series whose pulse count differs from this by more than 20 are dropped. Leave blank to keep all.');
L.RowHeight{10} = 30;

h.overwrite = uicheckbox(L, 'Text', 'Overwrite existing outputs', 'Value', true, ...
    'FontName', st.font, 'FontSize', st.fs, ...
    'Tooltip', 'Unchecked = io.overwrite_files ''never'' (skip if XML and JSON exist)');
place(h.overwrite, 11, [2 4]);

sec(L, 12, 'Target names');
h.suggest = uicheckbox(L, 'Text', 'Suggest from planned targets', 'Value', true, ...
    'FontName', st.font, 'FontSize', st.fs, ...
    'Tooltip', ['Name each series after the EntryTarget / GUMMarkers target nearest ' ...
                'its beam axis (neuronav_suggest_target_names)']);
place(h.suggest, 13, [2 4]);
lbl(L, 14, 'Axis tolerance');
h.tol = uieditfield(L, 'numeric', 'Value', 10, 'Limits', [0 Inf], ...
    'FontName', st.font, 'FontSize', st.fs, ...
    'Tooltip', 'Max. distance between a planned target and the beam axis');
place(h.tol, 14, 2);
unit(L, 14, 'mm');
lbl(L, 15, 'Target depth');
gl_depth = uigridlayout(L, [1 3], 'ColumnWidth', {'1x', 'fit', '1x'}, ...
    'Padding', [0 0 0 0], 'ColumnSpacing', 4, 'BackgroundColor', st.bg_panel);
place(gl_depth, 15, 2);
h.depth_min = uieditfield(gl_depth, 'numeric', 'Value', 0, 'Limits', [0 Inf], ...
    'FontName', st.font, 'FontSize', st.fs, 'Tooltip', 'Min. depth along the beam axis');
uilabel(gl_depth, 'Text', '–', 'FontName', st.font, 'FontSize', st.fs);
h.depth_max = uieditfield(gl_depth, 'numeric', 'Value', 150, 'Limits', [0 Inf], ...
    'FontName', st.font, 'FontSize', st.fs, 'Tooltip', 'Max. depth along the beam axis');
unit(L, 15, 'mm');
b = uibutton(L, 'Text', '↻ Re-apply to series', 'FontName', st.font, 'FontSize', st.fs_sm, ...
    'BackgroundColor', st.btn_lt, 'Tag', 'btn_suggest', ...
    'Tooltip', 'Re-run the naming on the detected series (overwrites names)', ...
    'ButtonPushedFcn', @(~,~) cb_suggest());
place(b, 16, [2 4]);

sec(L, 17, 'Coils (optional)');
note(L, 18, 'Empty = single-coil path (Coil0, GUMMarkers fallback). One row per coil = multi-coil path; comma-separate names for several blocks per coil.');
L.RowHeight{18} = 30;
h.coils = uitable(L, 'Data', cell(0, 2), ...
    'ColumnName', {'Coil', 'Target name'}, ...
    'ColumnEditable', [true true], 'ColumnFormat', {'numeric', 'char'}, ...
    'ColumnWidth', {60, 'auto'}, 'RowName', {}, 'FontName', st.font, 'FontSize', st.fs_sm, ...
    'CellEditCallback', @(~,~) cb_inputs_changed());
place(h.coils, 19, [1 4]);

gl_cb = uigridlayout(L, [1 2], 'Padding', [0 0 0 0], 'ColumnSpacing', 6, ...
    'BackgroundColor', st.bg_panel);
place(gl_cb, 20, [1 4]);
uibutton(gl_cb, 'Text', '+ Add coil', 'FontName', st.font, 'FontSize', st.fs_sm, ...
    'BackgroundColor', st.btn_lt, 'ButtonPushedFcn', @(~,~) cb_add_coil());
uibutton(gl_cb, 'Text', '− Remove last', 'FontName', st.font, 'FontSize', st.fs_sm, ...
    'BackgroundColor', st.btn_lt, 'ButtonPushedFcn', @(~,~) cb_remove_coil());

% Footer: actions + status, always visible below the scrolling inputs.
F = uigridlayout(left, [3 1], 'RowHeight', {30, 30, 44}, ...
    'Padding', [16 10 12 10], 'RowSpacing', 6, 'BackgroundColor', st.bg_fig);

h.btn_detect = uibutton(F, 'Text', '▶  1. Detect series', 'Tag', 'btn_detect', ...
    'FontName', st.font, 'FontSize', st.fs, 'FontWeight', 'bold', ...
    'BackgroundColor', st.accent, 'FontColor', [1 1 1], ...
    'ButtonPushedFcn', @(~,~) cb_detect());
place(h.btn_detect, 1, 1);

h.btn_write = uibutton(F, 'Text', '✓  2. Write outputs', 'Tag', 'btn_write', ...
    'FontName', st.font, 'FontSize', st.fs, 'FontWeight', 'bold', ...
    'BackgroundColor', st.accent_dk, 'FontColor', [1 1 1], 'Enable', 'off', ...
    'ButtonPushedFcn', @(~,~) cb_write());
place(h.btn_write, 2, 1);

h.status = uilabel(F, 'Text', '', 'Tag', 'lbl_status', 'FontName', st.font, ...
    'FontSize', st.fs_sm, 'FontColor', st.text_sub, 'WordWrap', 'on', ...
    'VerticalAlignment', 'top');
place(h.status, 3, 1);

%% ── Right panel: series table + result tabs ──────────────────────────────

R = uigridlayout(body, [3 1], ...
    'RowHeight', {22, 190, '1x'}, 'Padding', [12 12 12 12], 'RowSpacing', 6, ...
    'BackgroundColor', st.bg_fig);

uilabel(R, 'Text', ['Detected series — tick the ones to write and check the names. ' ...
    'Deviations: distance from the series mean position (mm), angle from the mean orientation (rad).'], ...
    'FontName', st.font, 'FontSize', st.fs_sm, 'FontColor', st.text_sub);

series_cols = {'Write', '#', 'Coil', 'Target name', 'Name from', 'Color', 'Pulses', ...
    'Pos (mm)', 'Pos SD', 'Rot (rad)', 'Rot SD', 'Flag'};
h.series = uitable(R, 'Data', cell(0, numel(series_cols)), 'Tag', 'tbl_series', ...
    'ColumnName', series_cols, 'RowName', {}, ...
    'ColumnEditable', [true false false true false true false false false false false false], ...
    'ColumnFormat', {'logical', 'numeric', 'char', 'char', 'char', 'char', 'numeric', ...
                     'numeric', 'numeric', 'numeric', 'numeric', 'char'}, ...
    'ColumnWidth', {46, 32, 40, 120, 150, 66, 52, 62, 54, 62, 54, 'auto'}, ...
    'FontName', st.font, 'FontSize', st.fs_sm, ...
    'CellEditCallback', @cb_series_edited);

tg = uitabgroup(R);
t_pos = uitab(tg, 'Title', 'Positions (3-D)');
t_dev = uitab(tg, 'Title', 'Deviation per pulse');
t_raw = uitab(tg, 'Title', 'Raw marker files');
t_out = uitab(tg, 'Title', 'Output');
t_log = uitab(tg, 'Title', 'Log');

g = uigridlayout(t_pos, [1 1], 'Padding', [4 4 4 4]);
h.ax_pos = uiaxes(g);
g = uigridlayout(t_dev, [1 1], 'Padding', [4 4 4 4]);
h.ax_dev = uiaxes(g);

g = uigridlayout(t_raw, [1 1], 'Padding', [4 4 4 4]);
h.raw_files = uitable(g, 'Data', cell(0, 5), 'Tag', 'tbl_raw', ...
    'ColumnName', {'Type', 'File', 'Folder (relative)', 'Size (kB)', 'Status'}, ...
    'RowName', {}, 'ColumnWidth', {110, 'auto', 'auto', 70, 150}, ...
    'FontName', st.font, 'FontSize', st.fs_sm);

g = uigridlayout(t_out, [2 3], 'RowHeight', {28, '1x'}, ...
    'ColumnWidth', {'1x', 'fit', 'fit'}, 'Padding', [6 6 6 6]);
h.out_path = uilabel(g, 'Text', 'No outputs written yet.', 'FontName', st.font, ...
    'FontSize', st.fs_sm, 'FontColor', st.text_sub);
h.btn_folder = uibutton(g, 'Text', 'Open output folder', 'Enable', 'off', ...
    'FontName', st.font, 'FontSize', st.fs_sm, 'ButtonPushedFcn', @(~,~) cb_open_folder());
h.btn_xml = uibutton(g, 'Text', 'Show XML', 'Enable', 'off', ...
    'FontName', st.font, 'FontSize', st.fs_sm, 'ButtonPushedFcn', @(~,~) cb_show_xml());
h.out_text = uitextarea(g, 'Editable', 'off', 'FontName', 'Menlo', 'FontSize', st.fs_sm, ...
    'Tag', 'txt_output');
place(h.out_text, 2, [1 3]);

g = uigridlayout(t_log, [1 1], 'Padding', [4 4 4 4]);
h.log = uitextarea(g, 'Editable', 'off', 'FontName', 'Menlo', 'FontSize', st.fs_sm, ...
    'Tag', 'txt_log');

reset_plots();

%% ── Initial values, then show ────────────────────────────────────────────

if ischar(init) || isstring(init)
    apply_config(load_yaml(char(init)));
elseif isstruct(init) && ~isempty(init)
    apply_config(init);
end

% Programmatic access (scripts/tests): settings as used by detect/batch, and save
fig.UserData = struct('settings', @() collect_settings(), ...
    'save_settings', @(file) cb_save_config(file));

fig.Visible = 'on';
drawnow;

%% ════════════════════════════════════════════════════════════════════════
%%  CALLBACKS
%% ════════════════════════════════════════════════════════════════════════

    function cb_rescan_subjects()
        raw_root = strtrim(h.raw.Value);
        subs = list_dirs(raw_root, 'sub-*');
        cur  = dd_value(h.sub);
        h.sub.Items = subs;
        if any(strcmp(subs, cur))
            h.sub.Value = cur;
        elseif ~isempty(subs)
            h.sub.Value = subs{1};
        else
            h.sub.Value = cur;
        end
        cb_rescan_sessions();
    end

    function cb_rescan_sessions()
        sub_dir = fullfile(strtrim(h.raw.Value), dd_value(h.sub));
        sess = list_dirs(sub_dir, 'ses-*');
        cur  = dd_value(h.ses);
        h.ses.Items = sess;
        if any(strcmp(sess, cur))
            h.ses.Value = cur;
        elseif ~isempty(sess)
            h.ses.Value = sess{1};
        else
            h.ses.Value = cur;
        end
        cb_inputs_changed();
    end

    function cb_inputs_changed()
        % Any input change invalidates a previous detect step.
        if ~isempty(S.snapshot)
            S.snapshot = [];
            h.btn_write.Enable = 'off';
            set_status('Inputs changed — run “Detect series” again before writing.', st.text_sub);
        end
        refresh_raw_files();
    end

    function cb_add_coil()
        d = h.coils.Data;
        if isempty(d)
            next_coil = 0;
        else
            next_coil = max(cell2mat(d(:, 1))) + 1;
        end
        h.coils.Data = [d; {next_coil, ''}];
        cb_inputs_changed();
    end

    function cb_remove_coil()
        d = h.coils.Data;
        if ~isempty(d)
            h.coils.Data = d(1:end-1, :);
            cb_inputs_changed();
        end
    end

    function cb_load_config()
        [f, d] = uigetfile({'*.yaml;*.yml', 'PRESTUS config (*.yaml)'}, 'Load settings from config');
        figure(fig);
        if isequal(f, 0), return; end
        try
            apply_config(load_yaml(fullfile(d, f)));
            set_status(sprintf('Loaded settings from %s', f), st.success);
        catch ME
            uialert(fig, ME.message, 'Load error', 'Icon', 'error');
        end
    end

    function cb_save_config(file)
        % Save the settings that define ingestion and naming, as a YAML
        % config: neuronav_ingest_batch(load_parameters(file)) reproduces them.
        if nargin < 1
            [f, d] = uiputfile({'*.yaml', 'PRESTUS config (*.yaml)'}, 'Save settings', 'neuronav_settings.yaml');
            figure(fig);
            if isequal(f, 0), return; end
            file = fullfile(d, f);
        end
        [p, msg] = collect_settings();
        if ~isempty(msg), set_status(['✗ ' msg], st.error); return; end
        out = struct();
        for key = {'localite_raw', 'localite_post', 'localite_t1', 'localite_t1_pattern'}
            if isfield(p.path, key{1}) && ~isempty(p.path.(key{1})), out.path.(key{1}) = p.path.(key{1}); end
        end
        lc = p.placement.localite;
        out.placement.localite = struct('loc_folder', lc.loc_folder, ...
            'expected_segment_length', lc.expected_segment_length, 'voxel_size', lc.voxel_size);
        if isfield(lc, 'session_folder_pattern')
            out.placement.localite.session_folder_pattern = lc.session_folder_pattern;
        end
        out.io.overwrite_files = p.io.overwrite_files;
        nv = p.neuronav;
        if isempty(nv.coil_map), nv.coil_map = []; end
        if isempty(nv.target_labels), nv.target_labels = []; end
        out.neuronav = struct('coil_map', {nv.coil_map}, 'suggest_names', nv.suggest_names, ...
            'axis_tolerance_mm', nv.axis_tolerance_mm, 'target_depth_mm', nv.target_depth_mm, ...
            'target_labels', {nv.target_labels}, 'batch_write', nv.batch_write);
        try
            yaml.dumpFile(file, out);
            set_status(sprintf('✓ Settings saved to %s', file), st.success);
        catch ME
            uialert(fig, ME.message, 'Save error', 'Icon', 'error');
        end
    end

    function cb_detect()
        [p, sub_id, ses_id, msg] = collect_inputs();
        if ~isempty(msg)
            set_status(['✗ ' msg], st.error); return;
        end

        busy(true, '⟳ Detecting series…');
        [out, txt, ME] = run_captured(@() neuronav_detect_series(p, sub_id, ses_id), 1);
        busy(false);

        S.D = []; S.N = []; S.results = {}; S.markertype = ''; S.sources = {};
        S.snapshot = []; S.candidates = []; S.suggest = [];
        h.btn_write.Enable = 'off';
        head = sprintf('── Detect  %s %s  (%s)', sub_id, ses_label(ses_id), datestr(now, 'HH:MM:SS')); %#ok<TNOW1,DATST>

        if ~isempty(ME)
            log_append(head, txt);
            show_series(); refresh_raw_files(); plot_series();
            set_status(['✗ ' ME.message], st.error);
            log_append('', ME.getReport('basic'));
            return;
        end

        D = out{1};
        log_append(head, [D.log newline txt]);
        S.D = D;
        S.results    = D.results;
        S.markertype = D.markertype;
        S.sources    = unique(D.sources(~cellfun(@isempty, D.sources)), 'stable');
        S.coil_ids   = D.coil_ids;
        S.candidates = D.candidates;

        show_series();
        apply_naming();
        refresh_raw_files();
        plot_series();

        n_ok = sum(~cellfun(@isempty, S.results));
        if n_ok == 0
            set_status('✗ No usable series found — see Log and Raw marker files.', st.error);
            return;
        end
        S.snapshot = snapshot();
        h.btn_write.Enable = 'on';
        set_status(sprintf('✓ %d series detected (%s, %s). Check series and names, then write.', ...
            n_ok, S.markertype, strjoin(S.sources, ', ')), st.success);
    end

    function cb_write()
        if isempty(S.snapshot) || ~isequal(S.snapshot, snapshot())
            set_status('Inputs changed — run “Detect series” again before writing.', st.error);
            h.btn_write.Enable = 'off';
            return;
        end
        [p, sub_id, ses_id, msg] = collect_inputs();
        if ~isempty(msg), set_status(['✗ ' msg], st.error); return; end
        if isempty(p.path.localite_post)
            set_status('✗ Output root is required for writing.', st.error); return;
        end

        [target_map, msg] = collect_target_map();
        if ~isempty(msg), set_status(['✗ ' msg], st.error); return; end

        busy(true, '⟳ Writing outputs…');
        [out, txt, ME] = run_captured(@() neuronav_ingest_markers(p, sub_id, ses_id, target_map), 3);
        log_append(sprintf('── Write  %s %s  (%s)', sub_id, ses_label(ses_id), datestr(now, 'HH:MM:SS')), txt); %#ok<TNOW1,DATST>
        busy(false);

        if ~isempty(ME)
            set_status(['✗ ' ME.message], st.error);
            log_append('', ME.getReport('basic'));
            return;
        end

        [xml_path, json_path, res] = out{:};
        if isempty(xml_path)
            set_status('✗ Ingestion wrote nothing — see Log.', st.error); return;
        end
        S.xml_path = xml_path; S.json_path = json_path;
        h.out_path.Text = xml_path;
        h.btn_folder.Enable = 'on'; h.btn_xml.Enable = 'on';
        h.out_text.Value = splitlines(fileread(json_path));
        tg.SelectedTab = t_out;
        if isempty(res)
            set_status(sprintf('Outputs already exist and overwrite is off — kept existing files in %s', ...
                fileparts(xml_path)), st.warn);
        else
            set_status(sprintf('✓ Wrote %d target(s) to %s', numel(target_map), fileparts(xml_path)), ...
                st.success);
        end
    end

    function cb_suggest()
        if isempty(S.D)
            set_status('Detect series first.', st.text_sub); return;
        end
        apply_naming();
        plot_series();
    end

    function apply_naming()
        % Same rule as neuronav_ingest_batch: coil table > heuristic (label
        % map applied) > series_<k>; see neuronav_name_series.
        p = collect_settings();
        N = neuronav_name_series(S.D, p);
        S.N = N;
        S.suggest = N;
        d = h.series.Data;
        for k = 1:numel(N)
            if isempty(S.results{k}), continue; end
            g = N(k);
            d{k, 1} = g.tracked;
            d{k, 4} = g.name;
            d{k, 5} = describe_name(g);
        end
        h.series.Data = d;
        log_append('', sprintf('Names: %d from the coil table, %d from planned targets, %d default (%d planned target(s)).', ...
            sum(strcmp({N.rule}, 'coil_map')), sum(strcmp({N.rule}, 'heuristic')), ...
            sum(strcmp({N.rule}, 'default') & [N.tracked]), numel(S.candidates)));
    end

    function t = describe_name(g)
        lab = g.label;
        if ~isempty(lab) && ~strcmp(lab, regexprep(g.name, '_(rep\d+|target)$', ''))
            lab = sprintf('"%s" → %s', lab, regexprep(g.name, '_(rep\d+|target)$', ''));
        elseif ~isempty(lab)
            lab = sprintf('"%s"', lab);
        end
        switch g.rule
            case 'coil_map'
                t = 'coil table';
                if ~isempty(g.label), t = sprintf('coil table (nearest planned: %s)', g.label); end
            case 'heuristic'
                if isnan(g.dist_mm)
                    t = sprintf('GUMMarkers pair %s', lab);
                else
                    t = sprintf('%s %s · %.1f mm off axis, %.0f mm deep', ...
                        strtok(g.source), lab, g.dist_mm, g.depth_mm);
                end
            otherwise
                if ~g.tracked
                    t = 'untracked';
                elseif ~h.suggest.Value
                    t = 'default';
                elseif isempty(S.candidates)
                    t = 'default (no planned targets)';
                elseif ~isnan(g.dist_mm)
                    t = sprintf('default (nearest %.1f mm off axis)', g.dist_mm);
                else
                    t = 'default (no match)';
                end
        end
    end

    function cb_series_edited(~, evt)
        % Keep empty (dropped-out) series un-ticked; replot on colour edits.
        r = evt.Indices(1); c = evt.Indices(2);
        if c == 1 && evt.NewData && isempty(S.results{r})
            h.series.Data{r, 1} = false;
        end
        if c == 6 || c == 1
            plot_series();
        end
    end

    function cb_open_folder()
        d = fileparts(S.xml_path);
        if ~isfolder(d), return; end
        if ismac
            system(['open ' shell_quote(d)]);
        elseif ispc
            winopen(d);
        else
            system(['xdg-open ' shell_quote(d) ' &']);
        end
    end

    function cb_show_xml()
        if isfile(S.xml_path)
            h.out_text.Value = splitlines(fileread(S.xml_path));
        end
    end

%% ════════════════════════════════════════════════════════════════════════
%%  INPUTS / OUTPUTS
%% ════════════════════════════════════════════════════════════════════════

    function [p, msg] = collect_settings()
        % The parameter struct a script would load from the study YAML: the
        % loaded config (S.base) with the GUI fields on top.
        p = S.base; msg = '';
        raw_root = strtrim(h.raw.Value);
        if isempty(raw_root) || ~isfolder(raw_root)
            msg = 'Raw Localite root is missing or not a folder.';
        end

        seg = strtrim(h.seglen.Value);
        seg_len = [];
        if ~isempty(seg)
            seg_len = str2double(seg);
            if isnan(seg_len) || seg_len <= 0
                msg = 'Pulses per series must be a positive number or blank.';
                seg_len = [];
            end
        end

        p.path.localite_raw  = raw_root;
        p.path.localite_post = strtrim(h.out.Value);
        p.placement.localite.loc_folder = strtrim(h.loc.Value);
        p.placement.localite.voxel_size = h.vox.Value;
        p.placement.localite.expected_segment_length = seg_len;
        if h.overwrite.Value
            p.io.overwrite_files = 'always';
        else
            p.io.overwrite_files = 'never';
        end
        p.neuronav.suggest_names     = logical(h.suggest.Value);
        p.neuronav.axis_tolerance_mm = h.tol.Value;
        p.neuronav.target_depth_mm   = [h.depth_min.Value, h.depth_max.Value];
        p.neuronav.target_labels     = S.target_labels;
        p.neuronav.batch_write       = S.batch_write;
        p.neuronav.coil_map          = [];

        d = h.coils.Data;
        if ~isempty(d)
            coils = d(:, 1);
            if any(cellfun(@(v) ~isnumeric(v) || isempty(v) || isnan(v) || v < 0 || v ~= round(v), coils))
                msg = 'Coil indices must be non-negative integers.'; return;
            end
            coils = cell2mat(coils);
            if numel(unique(coils)) < numel(coils)
                msg = 'Each coil may appear only once in the coil table.'; return;
            end
            % Target cell: one name, or comma-separated names for a coil that
            % records several blocks (one coil_map sequence each); blank = none.
            coil_map = struct('coil', {}, 'sequences', {});
            for k = 1:numel(coils)
                names = strtrim(strsplit(char(d{k, 2}), ','));
                names = names(~cellfun(@isempty, names));
                if isempty(names), names = {sprintf('coil%d', coils(k))}; end
                coil_map(k).coil = coils(k);
                coil_map(k).sequences = struct('target', names);
            end
            p.neuronav.coil_map = coil_map;
        end
    end

    function [p, sub_id, ses_id, msg] = collect_inputs()
        [p, msg] = collect_settings();
        sub_id = dd_value(h.sub);
        ses_id = parse_session(dd_value(h.ses));
        if ~isempty(msg), return; end
        if isempty(sub_id), msg = 'Subject is required.'; return; end
        if isempty(ses_id), msg = 'Session is required.'; return; end
    end

    function [target_map, msg] = collect_target_map()
        target_map = []; msg = '';
        d = h.series.Data;
        use = cell2mat(d(:, 1)) & ~cellfun(@isempty, S.results(:));
        if ~any(use)
            msg = 'Tick at least one series to write.'; return;
        end
        rows  = find(use);
        for k = 1:numel(rows)
            r = rows(k);
            if isempty(S.coil_ids), tid = 0; else, tid = S.coil_ids(r); end
            entry = struct('series_index', r, 'transducer_id', tid, ...
                'target_name', strtrim(char(d{r, 4})), 'color', strtrim(char(d{r, 6})));
            if isempty(target_map), target_map = entry; else, target_map(end+1) = entry; end %#ok<AGROW>
        end
        msg = neuronav_check_target_map(target_map);
        if ~isempty(msg), target_map = []; end
    end

    function show_series()
        n = numel(S.results);
        d = cell(n, 12);
        for k = 1:n
            r = S.results{k};
            if isempty(S.coil_ids)
                coil = '0';
                if ~strcmp(S.markertype, 'TriggerMarkers'), coil = '—'; end
            else
                coil = sprintf('%d', S.coil_ids(k));
            end
            if isempty(r)
                d(k, :) = {false, k, coil, '', 'no data', '', 0, NaN, NaN, NaN, NaN, 'no data'};
            else
                d(k, :) = {true, k, coil, sprintf('series_%d', k), 'default', '', r.N_pulses, ...
                    r.position_dev_mm.mean, r.position_dev_mm.std, ...
                    r.rotation_dev_rad.mean, r.rotation_dev_rad.std, neuronav_series_qc_flag(r)};
            end
        end
        h.series.Data = d;
    end

    function refresh_raw_files()
        base = fullfile(strtrim(h.raw.Value), dd_value(h.sub), ...
            ses_label(parse_session(dd_value(h.ses))), strtrim(h.loc.Value));
        rows = cell(0, 5);
        if isfolder(base)
            trig = dir(fullfile(base, '**', 'TMSTrigger', 'TriggerMarkers_Coil*.xml'));
            gum  = dir(fullfile(base, '**', 'GUMMarkers*.xml'));
            et   = dir(fullfile(base, '**', 'EntryTarget*.xml'));
            files = [trig; gum; et];
            types = [repmat({'TriggerMarkers'}, numel(trig), 1); ...
                     repmat({'GUMMarkers'}, numel(gum), 1); repmat({'EntryTarget'}, numel(et), 1)];
            label_src = {};
            if ~isempty(S.candidates), label_src = unique({S.candidates.source}); end
            for k = 1:numel(files)
                f = files(k);
                rel = strrep(f.folder, [base filesep], '');
                if strcmp(f.folder, base), rel = '.'; end
                status = '';
                if any(strcmp(S.sources, f.name))
                    status = '✓ used';
                elseif strcmp(types{k}, 'TriggerMarkers') && f.bytes <= 10000
                    status = 'ignored (≤ 10 kB)';
                end
                if any(strcmp(label_src, [types{k} ' ' f.name]))
                    status = strtrim([status '  target labels']);
                end
                rows(end+1, :) = {types{k}, f.name, rel, round(f.bytes / 1024, 1), status}; %#ok<AGROW>
            end
        end
        h.raw_files.Data = rows;
    end

    function apply_config(p)
        S.base = p;
        nv = get_nested(p, {'neuronav'}, struct());
        if isfield(nv, 'suggest_names') && ~isempty(nv.suggest_names), h.suggest.Value = logical(nv.suggest_names); end
        if isfield(nv, 'axis_tolerance_mm') && ~isempty(nv.axis_tolerance_mm), h.tol.Value = double(nv.axis_tolerance_mm); end
        if isfield(nv, 'target_depth_mm') && numel(nv.target_depth_mm) == 2
            h.depth_min.Value = double(nv.target_depth_mm(1));
            h.depth_max.Value = double(nv.target_depth_mm(2));
        end
        if isfield(nv, 'target_labels'), S.target_labels = nv.target_labels; end
        if isfield(nv, 'batch_write') && ~isempty(nv.batch_write), S.batch_write = char(nv.batch_write); end

        raw = get_nested(p, {'path', 'localite_raw'}, '');
        if isempty(raw), raw = get_nested(p, {'path', 'localite'}, ''); end
        if ~isempty(raw), h.raw.Value = char(raw); end
        post = get_nested(p, {'path', 'localite_post'}, '');
        if ~isempty(post), h.out.Value = char(post); end

        lc = get_nested(p, {'placement', 'localite'}, struct());
        if isfield(lc, 'loc_folder') && ~isempty(lc.loc_folder), h.loc.Value = char(lc.loc_folder); end
        if isfield(lc, 'voxel_size') && ~isempty(lc.voxel_size), h.vox.Value = double(lc.voxel_size); end
        if isfield(lc, 'expected_segment_length')
            h.seglen.Value = num2str(lc.expected_segment_length);
        end
        ow = get_nested(p, {'io', 'overwrite_files'}, '');
        if ~isempty(ow), h.overwrite.Value = ~strcmpi(char(ow), 'never'); end

        cm = get_nested(p, {'neuronav', 'coil_map'}, []);
        if isstruct(cm), cm = num2cell(cm); end
        if iscell(cm) && ~isempty(cm)
            d = cell(numel(cm), 2);
            for k = 1:numel(cm)
                e = cm{k};
                d{k, 1} = double(e.coil);
                d{k, 2} = '';
                if isfield(e, 'target') && ~isempty(e.target)
                    d{k, 2} = char(e.target);
                elseif isfield(e, 'sequences') && ~isempty(e.sequences)
                    seqs = e.sequences;
                    if isstruct(seqs), seqs = num2cell(seqs); end
                    d{k, 2} = strjoin(cellfun(@(q) char(q.target), seqs, 'UniformOutput', false), ', ');
                end
            end
            h.coils.Data = d;
        end

        cb_rescan_subjects();
        ses = get_nested(p, {'placement', 'localite', 'session'}, []);
        if ~isempty(ses)
            label = ses_label(parse_session(ses));
            if ~any(strcmp(h.ses.Items, label)), h.ses.Items = [h.ses.Items, {label}]; end
            h.ses.Value = label;
            cb_inputs_changed();
        end
    end

    function s = snapshot()
        s = {strtrim(h.raw.Value), dd_value(h.sub), ses_label(parse_session(dd_value(h.ses))), ...
             strtrim(h.loc.Value), h.vox.Value, strtrim(h.seglen.Value), h.coils.Data};
    end

%% ════════════════════════════════════════════════════════════════════════
%%  BATCH WINDOW
%% ════════════════════════════════════════════════════════════════════════

    function open_batch()
        if ~isempty(B) && isvalid(B.fig), figure(B.fig); return; end
        B.fig = uifigure('Name', 'PRESTUS · Localite batch', 'Tag', 'neuronav_gui_batch', ...
            'Position', [140 80 1180 780], 'Color', st.bg_fig);
        G = uigridlayout(B.fig, [6 1], 'RowHeight', {44, 30, '1x', 30, 34, '1x'}, ...
            'Padding', [12 10 12 12], 'RowSpacing', 6, 'BackgroundColor', st.bg_fig);
        uilabel(G, 'Text', ['Uses the settings of the main window (paths, Localite folder, pulses per ' ...
            'series, coil table, target naming). Tick sessions, check the label map, run. The same ' ...
            'run as a script: neuronav_ingest_batch(load_parameters(<saved settings>), selection).'], ...
            'WordWrap', 'on', 'FontName', st.font, 'FontSize', st.fs_sm, 'FontColor', st.text_sub);

        T = uigridlayout(G, [1 7], 'ColumnWidth', {'fit', 'fit', 'fit', 'fit', '1x', 'fit', 'fit'}, ...
            'Padding', [0 0 0 0], 'ColumnSpacing', 6, 'BackgroundColor', st.bg_fig);
        uibutton(T, 'Text', '↻ Scan sessions', 'FontName', st.font, 'FontWeight', 'bold', ...
            'BackgroundColor', st.accent, 'FontColor', [1 1 1], 'Tag', 'btn_scan', ...
            'ButtonPushedFcn', @(~,~) cb_scan());
        uibutton(T, 'Text', 'All', 'FontName', st.font, 'BackgroundColor', st.btn_lt, ...
            'Tag', 'btn_all', 'ButtonPushedFcn', @(~,~) select_sessions('all'));
        uibutton(T, 'Text', 'None', 'FontName', st.font, 'BackgroundColor', st.btn_lt, ...
            'ButtonPushedFcn', @(~,~) select_sessions('none'));
        uibutton(T, 'Text', 'Not yet ingested', 'FontName', st.font, 'BackgroundColor', st.btn_lt, ...
            'ButtonPushedFcn', @(~,~) select_sessions('new'));
        B.count = uilabel(T, 'Text', 'Scan to list the sessions under the raw root.', ...
            'FontName', st.font, 'FontSize', st.fs_sm, 'FontColor', st.text_sub);
        uibutton(T, 'Text', 'Group labels by stem', 'FontName', st.font, 'BackgroundColor', st.btn_lt, ...
            'Tag', 'btn_stem', 'Tooltip', 'Unmapped labels: target = label without a trailing _<number> (mapped labels are kept)', ...
            'ButtonPushedFcn', @(~,~) set_label_targets('stem'));
        uibutton(T, 'Text', 'Use labels as found', 'FontName', st.font, 'BackgroundColor', st.btn_lt, ...
            'ButtonPushedFcn', @(~,~) set_label_targets('label'));

        M = uigridlayout(G, [1 2], 'ColumnWidth', {'1.4x', '1x'}, 'Padding', [0 0 0 0], ...
            'ColumnSpacing', 10, 'BackgroundColor', st.bg_fig);
        B.sessions = uitable(M, 'Data', cell(0, 9), 'Tag', 'tbl_batch_sessions', ...
            'ColumnName', {'Include', 'Subject', 'Session', 'Trigger', 'GUM', 'EntryTarget', ...
                           'Planned labels', 'Ingested', 'Status'}, ...
            'ColumnEditable', [true false(1, 8)], ...
            'ColumnFormat', {'logical', 'char', 'char', 'numeric', 'numeric', 'numeric', 'char', 'logical', 'char'}, ...
            'ColumnWidth', {54, 80, 60, 54, 44, 74, 'auto', 60, 'auto'}, 'RowName', {}, ...
            'FontName', st.font, 'FontSize', st.fs_sm);
        B.labels = uitable(M, 'Data', cell(0, 5), 'Tag', 'tbl_batch_labels', ...
            'ColumnName', {'Label (as recorded)', 'Sessions', 'Subjects', 'Stem', 'Target name'}, ...
            'ColumnEditable', [false false false false true], ...
            'ColumnWidth', {'auto', 60, 60, 'auto', 'auto'}, 'RowName', {}, ...
            'FontName', st.font, 'FontSize', st.fs_sm, ...
            'CellEditCallback', @(~,~) labels_to_map());

        R2 = uigridlayout(G, [1 6], 'ColumnWidth', {'fit', 110, 'fit', '1x', 'fit', 'fit'}, ...
            'Padding', [0 0 0 0], 'ColumnSpacing', 8, 'BackgroundColor', st.bg_fig);
        uilabel(R2, 'Text', 'Write', 'FontName', st.font, 'FontSize', st.fs);
        B.write = uidropdown(R2, 'Items', {'all', 'named'}, 'Value', S.batch_write, ...
            'FontName', st.font, 'FontSize', st.fs, 'Tag', 'dd_batch_write', ...
            'Tooltip', 'all: every tracked series · named: only series named by the coil table or a planned target', ...
            'ValueChangedFcn', @(dd, ~) set_batch_write(dd.Value));
        B.run = uibutton(R2, 'Text', '▶  Run batch', 'FontName', st.font, 'FontWeight', 'bold', ...
            'BackgroundColor', st.accent_dk, 'FontColor', [1 1 1], 'Tag', 'btn_run_batch', ...
            'ButtonPushedFcn', @(~,~) cb_run_batch());
        B.progress = uilabel(R2, 'Text', '', 'FontName', st.font, 'FontSize', st.fs_sm, ...
            'FontColor', st.text_sub, 'Tag', 'lbl_batch_progress');
        B.open = uibutton(R2, 'Text', 'Open summary CSV', 'FontName', st.font, 'Enable', 'off', ...
            'ButtonPushedFcn', @(~,~) open_path(B.summary_file));
        uibutton(R2, 'Text', '◎ MNI view…', 'FontName', st.font, 'BackgroundColor', st.btn_lt, ...
            'Tag', 'btn_mni_view', ...
            'Tooltip', 'Show the written targets of all subjects and sessions in MNI space (summary CSV + SimNIBS warps)', ...
            'ButtonPushedFcn', @(~,~) open_mni());
        uilabel(G, 'Text', 'Written targets (summary CSV also has Localite pose, beam direction and transducer/focus coordinates)', ...
            'FontName', st.font, 'FontSize', st.fs_sm, 'FontColor', st.text_sub);
        B.results = uitable(G, 'Data', cell(0, 10), 'Tag', 'tbl_batch_results', ...
            'ColumnName', {'Subject', 'Session', 'Target', 'Rule', 'Label', 'Off axis (mm)', ...
                           'Depth (mm)', 'Pulses', 'Pos dev (mm)', 'Status'}, ...
            'ColumnWidth', {80, 60, 'auto', 70, 'auto', 80, 70, 54, 80, 'auto'}, 'RowName', {}, ...
            'FontName', st.font, 'FontSize', st.fs_sm);
        B.summary_file = '';
        B.scan = [];
        drawnow;
        if ~isempty(strtrim(h.raw.Value)) && isfolder(strtrim(h.raw.Value)), cb_scan(); end
    end

    function cb_scan()
        [p, msg] = collect_settings();
        if ~isempty(msg), B.count.Text = ['✗ ' msg]; B.count.FontColor = st.error; return; end
        B.count.Text = '⟳ Scanning…'; B.count.FontColor = st.accent; drawnow;
        [out, txt, ME] = run_captured(@() neuronav_scan_sessions(p), 2);
        log_append('── Batch scan', txt);
        if ~isempty(ME), B.count.Text = ['✗ ' ME.message]; B.count.FontColor = st.error; return; end
        [sess, labs] = out{:};
        B.scan = sess;
        n = height(sess);
        d = cell(n, 9);
        for i = 1:n
            usable = sess.n_trigger(i) + sess.n_gum(i) > 0;
            d(i, :) = {usable, sess.sub{i}, sess.ses{i}, sess.n_trigger(i), sess.n_gum(i), ...
                sess.n_entrytarget(i), sess.labels{i}, logical(sess.ingested(i)), ...
                ternary(usable, '', 'no marker files')};
        end
        B.sessions.Data = d;
        L = cell(height(labs), 5);
        for i = 1:height(labs)
            L(i, :) = {labs.label{i}, labs.n_sessions(i), labs.n_subjects(i), labs.stem{i}, labs.target{i}};
        end
        B.labels.Data = L;
        B.count.Text = sprintf('%d session(s) in %d subject(s); %d distinct planned-target label(s).', ...
            n, numel(unique(sess.sub)), height(labs));
        B.count.FontColor = st.text_sub;
    end

    function select_sessions(mode)
        d = B.sessions.Data;
        if isempty(d), return; end
        usable = cell2mat(d(:, 4)) + cell2mat(d(:, 5)) > 0;
        switch mode
            case 'all',  d(:, 1) = num2cell(usable);
            case 'none', d(:, 1) = {false};
            case 'new',  d(:, 1) = num2cell(usable & ~cell2mat(d(:, 8)));
        end
        B.sessions.Data = d;
    end

    function set_label_targets(mode)
        L = B.labels.Data;
        if isempty(L), return; end
        if strcmp(mode, 'stem')
            % Only labels not mapped yet: a curated entry is kept
            unmapped = strcmp(L(:, 5), L(:, 1));
            L(unmapped, 5) = L(unmapped, 4);
        else
            L(:, 5) = L(:, 1);
        end
        B.labels.Data = L;
        labels_to_map();
    end

    function labels_to_map()
        % Label table -> neuronav.target_labels (only labels whose target
        % differs from the label itself; an empty map keeps labels as found).
        L = B.labels.Data;
        S.target_labels = struct('target', {}, 'labels', {});
        if isempty(L), return; end
        L(:, 5) = cellfun(@(v) strtrim(char(v)), L(:, 5), 'UniformOutput', false);
        blank = cellfun(@isempty, L(:, 5));
        L(blank, 5) = L(blank, 1);
        B.labels.Data = L;
        changed = ~strcmp(L(:, 1), L(:, 5));
        [u, ~, idx] = unique(L(changed, 5), 'stable');
        lab = L(changed, 1);
        for k = 1:numel(u)
            S.target_labels(k).target = u{k};
            S.target_labels(k).labels = lab(idx == k)';
        end
        if ~isempty(S.D), apply_naming(); plot_series(); end
    end

    function set_batch_write(v)
        S.batch_write = v;
    end

    function cb_run_batch()
        d = B.sessions.Data;
        if isempty(d), B.progress.Text = 'Scan sessions first.'; return; end
        sel = find(cell2mat(d(:, 1)));
        if isempty(sel), B.progress.Text = 'Tick at least one session.'; return; end
        [p, msg] = collect_settings();
        if isempty(msg) && isempty(p.path.localite_post), msg = 'Set the output root in the main window.'; end
        if ~isempty(msg), B.progress.Text = ['✗ ' msg]; B.progress.FontColor = st.error; return; end
        selection = d(sel, 2:3);
        B.run.Enable = 'off'; B.progress.FontColor = st.accent;
        B.progress.Text = sprintf('⟳ 0 / %d sessions…', numel(sel)); drawnow;
        [out, txt, ME] = run_captured(@() neuronav_ingest_batch(p, selection, ...
            'progress', @(k, n, sub, ses, status) batch_progress(sel(k), k, n, status)), 1);
        B.run.Enable = 'on';
        log_append(sprintf('── Batch  %d session(s)  (%s)', numel(sel), datestr(now, 'HH:MM:SS')), txt); %#ok<TNOW1,DATST>
        if ~isempty(ME)
            B.progress.Text = ['✗ ' ME.message]; B.progress.FontColor = st.error; return;
        end
        T = out{1};
        B.results.Data = [T.sub, T.ses, T.target, T.rule, T.label, num2cell(T.axis_dist_mm), ...
            num2cell(T.depth_mm), num2cell(T.n_pulses), num2cell(T.pos_dev_mm), T.status];
        B.summary_file = fullfile(p.path.localite_post, 'neuronav_ingest_summary.csv');
        B.open.Enable = ternary(isfile(B.summary_file), 'on', 'off');
        n_written = sum(~cellfun(@isempty, T.target));
        B.progress.Text = sprintf('✓ %d target(s) from %d session(s). Summary: %s', ...
            n_written, numel(sel), B.summary_file);
        B.progress.FontColor = st.success;
        cb_scan_keep_selection(sel);
    end

    function batch_progress(row, k, n, status)
        d = B.sessions.Data;
        d{row, 9} = status;
        B.sessions.Data = d;
        B.progress.Text = sprintf('⟳ %d / %d sessions…', k, n);
        drawnow limitrate;
    end

    function cb_scan_keep_selection(sel)
        % Refresh the 'Ingested' column without losing ticks and statuses.
        d_old = B.sessions.Data;
        cb_scan();
        d = B.sessions.Data;
        if size(d, 1) == size(d_old, 1)
            d(:, 1) = {false}; d(sel, 1) = {true};
            d(:, 9) = d_old(:, 9);
            B.sessions.Data = d;
        end
    end

    function open_path(f)
        if isempty(f) || ~isfile(f), return; end
        if ismac
            system(['open ' shell_quote(f)]);
        elseif ispc
            winopen(f);
        else
            system(['xdg-open ' shell_quote(f) ' &']);
        end
    end

%% ════════════════════════════════════════════════════════════════════════
%%  MNI VIEW WINDOW
%% ════════════════════════════════════════════════════════════════════════

    function open_mni()
        if ~isempty(MV) && isvalid(MV.fig), figure(MV.fig); return; end
        MV = struct();
        MV.fig = uifigure('Name', 'PRESTUS · Localite targets in MNI space', 'Tag', 'neuronav_gui_mni', ...
            'Position', [200 90 1200 780], 'Color', st.bg_fig);
        G = uigridlayout(MV.fig, [1 2], 'ColumnWidth', {360, '1x'}, 'Padding', [0 0 0 0], ...
            'ColumnSpacing', 0, 'BackgroundColor', st.bg_fig);
        C = uigridlayout(G, [16 4], ...
            'RowHeight', {26, 26, 26, 26, 44, 30, 44, 26, 26, 26, 26, '1x', '1x', '1x', 30, 1}, ...
            'ColumnWidth', {100, '1x', 'fit', 34}, 'Padding', [14 12 10 12], ...
            'RowSpacing', 4, 'ColumnSpacing', 6, 'BackgroundColor', st.bg_panel);

        sec(C, 1, 'Input');
        lbl(C, 2, 'Summary CSV');
        MV.csv = uieditfield(C, 'text', 'FontName', st.font, 'FontSize', st.fs, 'Tag', 'ed_mni_csv', ...
            'Value', fullfile(strtrim(h.out.Value), 'neuronav_ingest_summary.csv'), ...
            'Tooltip', 'neuronav_ingest_summary.csv written by the batch run (all ingested sessions)');
        place(MV.csv, 2, [2 3]);
        mni_browse(C, 2, MV.csv, 'file', {'*.csv', 'Summary CSV (*.csv)'});
        lbl(C, 3, 'SimNIBS root');
        MV.seg = uieditfield(C, 'text', 'FontName', st.font, 'FontSize', st.fs, 'Tag', 'ed_mni_seg', ...
            'Value', char(get_nested(S.base, {'path', 'seg'}, '')), 'Placeholder', 'path.seg', ...
            'Tooltip', 'path.seg — contains m2m_<sub>/final_tissues.nii.gz and m2m_<sub>/toMNI/Conform2MNI_nonl.nii.gz');
        place(MV.seg, 3, [2 3]);
        mni_browse(C, 3, MV.seg, 'dir', []);
        lbl(C, 4, 'Head template');
        MV.tpl = uieditfield(C, 'text', 'FontName', st.font, 'FontSize', st.fs, 'Tag', 'ed_mni_template', ...
            'Value', neuronav_find_mni_template(S.base), 'Placeholder', 'first subject''s final_tissues_MNI', ...
            'Tooltip', 'MNI-space NIfTI for the head outline (default: SimNIBS MNI152_T1_1mm.nii.gz)', ...
            'ValueChangedFcn', @(~,~) plot_mni());
        place(MV.tpl, 4, [2 3]);
        mni_browse(C, 4, MV.tpl, 'file', {'*.nii;*.nii.gz', 'NIfTI (*.nii, *.nii.gz)'});
        note(C, 5, ['Planning-image coordinates are mapped with each subject''s SimNIBS nonlinear ' ...
            'warp (neuronav_summary_to_mni). The planning T1 is found as in the main window.']);
        b = uibutton(C, 'Text', '▶  Load & map to MNI', 'FontName', st.font, 'FontWeight', 'bold', ...
            'BackgroundColor', st.accent, 'FontColor', [1 1 1], 'Tag', 'btn_mni_map', ...
            'ButtonPushedFcn', @(~,~) cb_mni_map());
        place(b, 6, [1 4]);
        MV.status = uilabel(C, 'Text', 'Load the summary CSV to map its targets.', 'FontName', st.font, ...
            'FontSize', st.fs_sm, 'FontColor', st.text_sub, 'WordWrap', 'on', ...
            'VerticalAlignment', 'top', 'Tag', 'lbl_mni_status');
        place(MV.status, 7, [1 4]);

        sec(C, 8, 'Display');
        g = uigridlayout(C, [1 3], 'Padding', [0 0 0 0], 'ColumnSpacing', 4, 'BackgroundColor', st.bg_panel);
        place(g, 9, [1 4]);
        MV.show_trans = mni_check(g, 'Transducer', true, 'Bowl rear centre (filled circle)');
        MV.show_focus = mni_check(g, 'Focus', true, 'Acoustic focus (diamond)');
        MV.show_beam  = mni_check(g, 'Beam axis', true, 'Line from transducer to focus');
        g = uigridlayout(C, [1 3], 'Padding', [0 0 0 0], 'ColumnSpacing', 4, 'BackgroundColor', st.bg_panel);
        place(g, 10, [1 4]);
        MV.show_pose = mni_check(g, 'Tracker pose', false, 'Mean Localite instrument position (triangle)');
        MV.show_head = mni_check(g, 'Head outline', true, 'Isosurface of the head template');
        lbl(C, 11, 'Colour by');
        MV.color_by = uidropdown(C, 'Items', {'target', 'subject', 'session'}, 'Value', 'target', ...
            'FontName', st.font, 'FontSize', st.fs, 'Tag', 'dd_mni_color', ...
            'ValueChangedFcn', @(~,~) plot_mni());
        place(MV.color_by, 11, [2 3]);

        MV.subs    = mni_filter_table(C, 12, 'Subject', 'tbl_mni_subjects');
        MV.sess    = mni_filter_table(C, 13, 'Session', 'tbl_mni_sessions');
        MV.targets = mni_filter_table(C, 14, 'Target',  'tbl_mni_targets');

        g = uigridlayout(C, [1 3], 'ColumnWidth', {'1x', '1x', '1x'}, 'Padding', [0 0 0 0], ...
            'ColumnSpacing', 6, 'BackgroundColor', st.bg_panel);
        place(g, 15, [1 4]);
        uibutton(g, 'Text', 'Tick all', 'FontName', st.font, 'FontSize', st.fs_sm, ...
            'BackgroundColor', st.btn_lt, 'ButtonPushedFcn', @(~,~) mni_tick_all(true));
        uibutton(g, 'Text', 'Untick all', 'FontName', st.font, 'FontSize', st.fs_sm, ...
            'BackgroundColor', st.btn_lt, 'ButtonPushedFcn', @(~,~) mni_tick_all(false));
        MV.save = uibutton(g, 'Text', '⬇ Save CSV…', 'FontName', st.font, 'FontSize', st.fs_sm, ...
            'BackgroundColor', st.btn_lt, 'Enable', 'off', 'Tag', 'btn_mni_save', ...
            'Tooltip', 'Write the summary with its MNI columns', ...
            'ButtonPushedFcn', @(~,~) cb_mni_save());

        g = uigridlayout(G, [1 1], 'Padding', [8 8 8 8], 'BackgroundColor', st.bg_fig);
        MV.ax = uiaxes(g, 'Tag', 'ax_mni');
        MV.T = [];
        MV.head = struct('key', '', 'faces', [], 'vertices', []);
        MV.view = [];
        plot_mni();
    end

    function cb = mni_check(parent, text, value, tip)
        cb = uicheckbox(parent, 'Text', text, 'Value', value, 'FontName', st.font, ...
            'FontSize', st.fs_sm, 'Tooltip', tip, 'ValueChangedFcn', @(~,~) plot_mni());
    end

    function t = mni_filter_table(parent, row, what, tag)
        t = uitable(parent, 'Data', cell(0, 3), 'Tag', tag, ...
            'ColumnName', {'Show', what, 'Targets'}, 'ColumnEditable', [true false false], ...
            'ColumnFormat', {'logical', 'char', 'numeric'}, 'ColumnWidth', {44, 'auto', 56}, ...
            'RowName', {}, 'FontName', st.font, 'FontSize', st.fs_sm, ...
            'CellEditCallback', @(~,~) plot_mni());
        place(t, row, [1 4]);
    end

    function mni_browse(parent, row, target, kind, filter)
        o = uibutton(parent, 'Text', '···', 'FontName', st.font, 'FontSize', st.fs, ...
            'BackgroundColor', st.btn_lt, 'ButtonPushedFcn', @(~,~) pick());
        place(o, row, 4);
        function pick()
            start = strtrim(target.Value);
            if strcmp(kind, 'dir')
                if ~isfolder(start), start = pwd; end
                d = uigetdir(start, 'Select folder');
                figure(MV.fig);
                if isequal(d, 0), return; end
                target.Value = d;
            else
                if ~isfolder(start), start = fileparts(start); end
                if ~isfolder(start), start = pwd; end
                [f, d] = uigetfile(filter, 'Select file', [start filesep]);
                figure(MV.fig);
                if isequal(f, 0), return; end
                target.Value = fullfile(d, f);
                if target == MV.tpl, plot_mni(); end
            end
        end
    end

    function mni_tick_all(v)
        for t = [MV.subs, MV.sess, MV.targets]
            d = t.Data;
            if ~isempty(d), d(:, 1) = {v}; t.Data = d; end
        end
        plot_mni();
    end

    function cb_mni_map()
        csv = strtrim(MV.csv.Value);
        if ~isfile(csv), mni_status(['✗ No summary CSV: ' csv], st.error); return; end
        [p, msg] = collect_settings();
        if ~isempty(msg), mni_status(['✗ ' msg ' (the planning T1 is looked up there)'], st.error); return; end
        p.path.seg = strtrim(MV.seg.Value);
        if ~isfolder(p.path.seg), mni_status('✗ Set the SimNIBS root (path.seg).', st.error); return; end
        mni_status('⟳ Mapping to MNI…', st.accent); drawnow;
        [out, txt, ME] = run_captured(@() neuronav_summary_to_mni(csv, p, 'progress', ...
            @(k, n, sub, ~) mni_status(sprintf('⟳ %d / %d subjects (%s)…', k, n, sub), st.accent)), 1);
        log_append(sprintf('── MNI mapping  %s  (%s)', csv, datestr(now, 'HH:MM:SS')), txt); %#ok<TNOW1,DATST>
        if ~isempty(ME), mni_status(['✗ ' ME.message], st.error); return; end
        T = out{1};
        T = T(~cellfun(@isempty, T.target), :);   % sessions without written targets
        MV.T = T;
        MV.save.Enable = ternary(isempty(T), 'off', 'on');
        mni_fill_table(MV.subs, T.sub);
        mni_fill_table(MV.sess, T.ses);
        mni_fill_table(MV.targets, T.target);

        ok = cellfun(@isempty, T.mni_status);
        [bad, ib] = unique(T.sub(~ok), 'stable');
        why = T.mni_status(~ok); why = why(ib);
        text = sprintf('✓ %d target(s) of %d subject(s) mapped.', sum(ok), numel(unique(T.sub(ok))));
        if ~isempty(bad)
            text = sprintf('%s Not mapped: %s', text, strjoin(strcat(bad, {' ('}, why, {')'}), '; '));
        end
        mni_status(text, ternary(isempty(bad), st.success, st.warn));
        MV.view = [];
        plot_mni();
    end

    function mni_status(text, color)
        MV.status.Text = text; MV.status.FontColor = color;
        drawnow limitrate;
    end

    function mni_fill_table(t, names)
        % One row per distinct name; keep the tick of names already listed.
        [u, ~, idx] = unique(names);
        d = [num2cell(true(numel(u), 1)), u(:), num2cell(accumarray(idx(:), 1, [numel(u) 1]))];
        old = t.Data;
        for i = 1:size(d, 1)
            j = find(strcmp(old(:, 2), d{i, 2}), 1);
            if ~isempty(j), d{i, 1} = old{j, 1}; end
        end
        t.Data = d;
    end

    function cb_mni_save()
        [f, d] = uiputfile({'*.csv', 'CSV (*.csv)'}, 'Save MNI coordinates', ...
            fullfile(fileparts(strtrim(MV.csv.Value)), 'neuronav_ingest_summary_mni.csv'));
        figure(MV.fig);
        if isequal(f, 0), return; end
        try
            writetable(MV.T, fullfile(d, f));
            mni_status(sprintf('✓ Saved %s', fullfile(d, f)), st.success);
        catch ME
            mni_status(['✗ ' ME.message], st.error);
        end
    end

    function plot_mni()
        ax = MV.ax;
        if ~isempty(ax.Children) && ~isempty(MV.T)
            [az, el] = view(ax); MV.view = [az el];   % keep the user's rotation
        end
        cla(ax); legend(ax, 'off');
        hold(ax, 'on');
        if MV.show_head.Value
            mni_head();
            if ~isempty(MV.head.faces)
                patch(ax, 'Faces', MV.head.faces, 'Vertices', MV.head.vertices, ...
                    'FaceColor', [0.86 0.78 0.72], 'EdgeColor', 'none', 'FaceAlpha', 0.18, ...
                    'FaceLighting', 'gouraud', 'HandleVisibility', 'off', 'PickableParts', 'none');
            end
        end
        T = MV.T;
        n_shown = 0;
        hl = gobjects(0); leg = {};
        if ~isempty(T)
            use = cellfun(@isempty, T.mni_status) & ...
                ticked(MV.subs, T.sub) & ticked(MV.sess, T.ses) & ticked(MV.targets, T.target);
            key = T.(ternary(strcmp(MV.color_by.Value, 'target'), 'target', ...
                ternary(strcmp(MV.color_by.Value, 'subject'), 'sub', 'ses')));
            groups = unique(key(use));
            palette = lines(max(numel(groups), 1));
            if numel(groups) > 7, palette = turbo(numel(groups)); end
            for k = 1:numel(groups)
                r = use & strcmp(key, groups{k});
                c = palette(k, :);
                tr = [T.trans_mni_x(r), T.trans_mni_y(r), T.trans_mni_z(r)];
                fo = [T.focus_mni_x(r), T.focus_mni_y(r), T.focus_mni_z(r)];
                po = [T.pose_mni_x(r),  T.pose_mni_y(r),  T.pose_mni_z(r)];
                tips = {T.sub(r), T.ses(r), T.target(r)};
                hk = gobjects(0);
                if MV.show_beam.Value
                    seg = permute(cat(3, tr, fo, nan(size(tr))), [3 1 2]);
                    seg = reshape(seg, [], 3);
                    plot3(ax, seg(:, 1), seg(:, 2), seg(:, 3), '-', 'Color', [c 0.6], ...
                        'LineWidth', 0.8, 'HandleVisibility', 'off');
                end
                if MV.show_trans.Value, hk(end+1) = mni_scatter(ax, tr, 'o', 36, c, true, tips); end %#ok<AGROW>
                if MV.show_focus.Value, hk(end+1) = mni_scatter(ax, fo, 'd', 36, c, true, tips); end %#ok<AGROW>
                if MV.show_pose.Value,  hk(end+1) = mni_scatter(ax, po, '^', 30, c, false, tips); end %#ok<AGROW>
                if ~isempty(hk)
                    hl(end+1) = hk(1); %#ok<AGROW>
                    leg{end+1} = sprintf('%s (%d)', groups{k}, sum(r)); %#ok<AGROW>
                    set(hk(2:end), 'HandleVisibility', 'off');
                end
            end
            n_shown = sum(use);
        end
        hold(ax, 'off');
        axis(ax, 'equal'); grid(ax, 'on');
        xlabel(ax, 'x MNI (mm, +R)'); ylabel(ax, 'y MNI (mm, +A)'); zlabel(ax, 'z MNI (mm, +S)');
        if isempty(MV.view), view(ax, [-135 25]); else, view(ax, MV.view); end
        if MV.show_head.Value && ~isempty(MV.head.faces), camlight(ax, 'headlight'); end
        if ~isempty(hl)
            legend(ax, hl, leg, 'Interpreter', 'none', 'Location', 'northeastoutside');
        end
        title(ax, sprintf('Targets in MNI space (%d shown)', n_shown));
        subtitle(ax, '● transducer · ◆ focus · ▲ tracker pose · line = beam axis', ...
            'FontSize', 9, 'Color', st.text_sub);
    end

    function mni_head()
        % Head outline from the template (or the first mapped subject's
        % final_tissues_MNI), cached per source file.
        f = strtrim(MV.tpl.Value);
        if isempty(f) && ~isempty(MV.T)
            ok = find(cellfun(@isempty, MV.T.mni_status), 1);
            if ~isempty(ok)
                f = fullfile(strtrim(MV.seg.Value), sprintf('m2m_%s', MV.T.sub{ok}), ...
                    'toMNI', 'final_tissues_MNI.nii.gz');
            end
        end
        if strcmp(f, MV.head.key), return; end
        MV.head = struct('key', f, 'faces', [], 'vertices', []);
        if isempty(f) || ~isfile(f), return; end
        try
            [MV.head.faces, MV.head.vertices] = neuronav_head_surface(f);
        catch ME
            mni_status(['Head outline: ' ME.message], st.warn);
        end
    end

    function hs = mni_scatter(ax, xyz, marker, sz, c, filled, tips)
        if filled
            hs = scatter3(ax, xyz(:, 1), xyz(:, 2), xyz(:, 3), sz, c, marker, 'filled', ...
                'MarkerEdgeColor', c * 0.6);
        else
            hs = scatter3(ax, xyz(:, 1), xyz(:, 2), xyz(:, 3), sz, c, marker, 'LineWidth', 1.2);
        end
        try
            rows = hs.DataTipTemplate.DataTipRows;
            hs.DataTipTemplate.DataTipRows = [dataTipTextRow('Subject', tips{1}); ...
                dataTipTextRow('Session', tips{2}); dataTipTextRow('Target', tips{3}); rows(:)];
            hs.DataTipTemplate.Interpreter = 'none';
        catch
        end
    end

%% ════════════════════════════════════════════════════════════════════════
%%  PLOTS
%% ════════════════════════════════════════════════════════════════════════

    function reset_plots()
        for ax = [h.ax_pos, h.ax_dev]
            cla(ax); legend(ax, 'off');
        end
        title(h.ax_pos, 'Pulse positions per series'); subtitle(h.ax_pos, '');
        xlabel(h.ax_pos, 'x (mm)'); ylabel(h.ax_pos, 'y (mm)'); zlabel(h.ax_pos, 'z (mm)');
        title(h.ax_dev, 'Distance of each pulse from its series mean');
        xlabel(h.ax_dev, 'Pulse within series');
        ylabel(h.ax_dev, 'Distance (mm)');
    end

    function plot_series()
        reset_plots();
        d = h.series.Data;
        n = numel(S.results);
        if n == 0, return; end
        palette = lines(max(n, 1));
        hold(h.ax_pos, 'on'); hold(h.ax_dev, 'on');
        leg = {}; hp = gobjects(0); hd = gobjects(0);
        for k = 1:n
            r = S.results{k};
            if isempty(r), continue; end
            c = hex_to_rgb(d{k, 6}, palette(k, :));
            P = reshape(r.matrix4d_all(:, 1:3, 4), [], 3);
            % One small dot per pulse; the mean is an open ring so the
            % pulse cloud (typically a few mm) stays visible inside it.
            hp(end+1) = scatter3(h.ax_pos, P(:, 1), P(:, 2), P(:, 3), 6, c, 'filled', ...
                'MarkerFaceAlpha', 0.7); %#ok<AGROW>
            M = reshape(r.matrix4d_mean, 4, 4);
            scatter3(h.ax_pos, M(1, 4), M(2, 4), M(3, 4), 110, c * 0.6, 'o', ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
            if ~strcmp(S.markertype, 'GUMMarkers')
                % Column 1 = beam direction (coil centre → head), as used by
                % localite_matrix_to_positions. GUMMarkers matrices are
                % synthesised with an arbitrary column 1, so no arrow there.
                v = M(1:3, 1) * 30;
                quiver3(h.ax_pos, M(1, 4), M(2, 4), M(3, 4), v(1), v(2), v(3), 0, ...
                    'Color', c, 'LineWidth', 1.5, 'MaxHeadSize', 0.4, 'HandleVisibility', 'off');
            end
            hd(end+1) = plot(h.ax_dev, 1:numel(r.position_dev_mm.all), r.position_dev_mm.all, ...
                '.-', 'Color', c, 'LineWidth', 1.2, 'MarkerSize', 9); %#ok<AGROW>
            leg{end+1} = sprintf('%d · %s', k, char(d{k, 4})); %#ok<AGROW>
        end
        if ~isempty(S.candidates)
            xyz = vertcat(S.candidates.xyz);
            scatter3(h.ax_pos, xyz(:, 1), xyz(:, 2), xyz(:, 3), 40, [0.2 0.2 0.2], 'x', ...
                'LineWidth', 1.5, 'HandleVisibility', 'off');
            text(h.ax_pos, xyz(:, 1), xyz(:, 2), xyz(:, 3), ...
                strcat({'  '}, {S.candidates.description}'), 'FontSize', 9, ...
                'Color', [0.3 0.3 0.3], 'Interpreter', 'none');
        end
        for k = 1:numel(S.suggest)
            g = S.suggest(k);
            if isempty(S.results{k}) || isempty(g.name) || isnan(g.depth_mm), continue; end
            M = reshape(S.results{k}.matrix4d_mean, 4, 4);
            a = M(1:3, 4); e = a + g.depth_mm * M(1:3, 1) / norm(M(1:3, 1));
            plot3(h.ax_pos, [a(1) e(1)], [a(2) e(2)], [a(3) e(3)], ':', ...
                'Color', hex_to_rgb(d{k, 6}, palette(k, :)), 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end
        hold(h.ax_pos, 'off'); hold(h.ax_dev, 'off');
        axis(h.ax_pos, 'equal'); grid(h.ax_pos, 'on'); view(h.ax_pos, 3);
        grid(h.ax_dev, 'on');
        if ~isempty(hp)
            legend(h.ax_pos, hp, leg, 'Interpreter', 'none', 'Location', 'northeastoutside');
            legend(h.ax_dev, hd, leg, 'Interpreter', 'none', 'Location', 'northeastoutside');
        end
        if ~strcmp(S.markertype, 'GUMMarkers')
            subtitle(h.ax_pos, ['dots = pulses · ring = mean · arrow = beam direction (30 mm) · ' ...
                '× = planned target · dotted = beam axis to matched target'], ...
                'FontSize', 9, 'Color', st.text_sub);
        end
    end

%% ════════════════════════════════════════════════════════════════════════
%%  UI HELPERS
%% ════════════════════════════════════════════════════════════════════════

    function busy(on, text)
        if on
            h.btn_detect.Enable = 'off'; h.btn_write.Enable = 'off';
            set_status(text, st.accent);
        else
            h.btn_detect.Enable = 'on';
            if ~isempty(S.snapshot), h.btn_write.Enable = 'on'; end
        end
        drawnow;
    end

    function set_status(text, color)
        h.status.Text = text;
        h.status.FontColor = color;
    end

    function log_append(header, txt)
        lines_ = {};
        if ~isempty(header), lines_{end+1} = header; end
        if ~isempty(txt)
            t = splitlines(strtrim(txt));
            lines_ = [lines_, t(:)'];
        end
        v = h.log.Value;
        if isscalar(v) && isempty(v{1}), v = {}; end
        h.log.Value = [v(:); lines_(:)];
        try, scroll(h.log, 'bottom'); catch, end
    end

    function place(obj, row, col)
        obj.Layout.Row = row; obj.Layout.Column = col;
    end

    function lbl(parent, row, text)
        o = uilabel(parent, 'Text', text, 'FontName', st.font, 'FontSize', st.fs, ...
            'FontColor', st.text, 'HorizontalAlignment', 'right');
        place(o, row, 1);
    end

    function unit(parent, row, text)
        o = uilabel(parent, 'Text', text, 'FontName', st.font, 'FontSize', st.fs_sm, ...
            'FontColor', st.text_sub);
        place(o, row, 3);
    end

    function note(parent, row, text)
        o = uilabel(parent, 'Text', text, 'FontName', st.font, 'FontSize', st.fs_sm, ...
            'FontColor', st.text_sub, 'WordWrap', 'on');
        place(o, row, [1 4]);
    end

    function sec(parent, row, text)
        o = uilabel(parent, 'Text', ['  ' upper(text)], 'FontName', st.font, ...
            'FontSize', st.fs_sm, 'FontWeight', 'bold', 'FontColor', st.accent, ...
            'BackgroundColor', st.accent_lt);
        place(o, row, [1 4]);
    end

    function sp(parent, row)
        o = uilabel(parent, 'Text', '');
        place(o, row, 1);
    end

    function brw(parent, row, target, after_fcn)
        o = uibutton(parent, 'Text', '···', 'FontName', st.font, 'FontSize', st.fs, ...
            'BackgroundColor', st.btn_lt, 'ButtonPushedFcn', @(~,~) browse(target, after_fcn));
        place(o, row, 4);
    end

    function browse(target, after_fcn)
        start = strtrim(target.Value);
        if ~isfolder(start), start = pwd; end
        d = uigetdir(start, 'Select folder');
        figure(fig);
        if isequal(d, 0), return; end
        target.Value = d;
        if ~isempty(after_fcn), after_fcn(); end
    end

end % neuronav_gui


%% ════════════════════════════════════════════════════════════════════════
%%  LOCAL FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════

function [out, txt, ME] = run_captured(fcn, nout)
% Run fcn with nout outputs, capturing command-window text and warnings.
    out = cell(1, nout); ME = [];
    bt = warning('off', 'backtrace');
    restore = onCleanup(@() warning(bt));
    % try/catch inside the evalc string keeps the text printed before an error.
    txt = evalc('try, [out{1:nout}] = fcn(); catch ME, end');
    txt = regexprep(txt, '</?a[^>]*>', '');   % hyperlinks in warning text
end

function s = sty()
% Palette shared with prestus_gui (derived from the PRESTUS logo).
    s.bg_fig    = [0.945 0.945 0.947];
    s.bg_panel  = [1.000 1.000 1.000];
    s.bg_dark   = [0.122 0.122 0.125];
    s.accent    = [0.043 0.482 0.600];
    s.accent_lt = [0.851 0.933 0.945];
    s.accent_dk = [0.027 0.322 0.404];
    s.btn_lt    = [0.88 0.91 0.96];
    s.error     = [0.663 0.114 0.145];
    s.warn      = [0.702 0.420 0.000];
    s.success   = [0.180 0.631 0.239];
    s.text      = [0.082 0.082 0.090];
    s.text_sub  = [0.400 0.400 0.420];
    s.font      = 'Helvetica Neue';
    s.fs        = 13;
    s.fs_sm     = 11;
    s.fs_xl     = 20;
end

function v = dd_value(dd)
% Editable dropdown value as trimmed char (an empty Items list yields {}).
    v = strtrim(char(dd.Value));
end

function names = list_dirs(parent, pattern)
    names = {};
    if isempty(parent) || ~isfolder(parent), return; end
    d = dir(fullfile(parent, pattern));
    d = d([d.isdir]);
    names = sort({d.name});
end

function ses = parse_session(v)
% '2' / 2 → 2 (ingest formats it as ses-02); 'ses-02' stays a string.
    if isnumeric(v)
        ses = v; return;
    end
    v = strtrim(char(v));
    if ~isempty(regexp(v, '^\d+$', 'once'))
        ses = str2double(v);
    else
        ses = v;
    end
end

function label = ses_label(ses)
    if isnumeric(ses) && ~isempty(ses)
        label = sprintf('ses-%02d', double(ses));
    else
        label = char(ses);
    end
end

function tf = ticked(t, names)
% Rows whose name is ticked in a Show/name/count filter table (empty = all).
    d = t.Data;
    if isempty(d), tf = true(size(names)); return; end
    tf = ismember(names, d(cell2mat(d(:, 1)), 2));
end




function rgb = hex_to_rgb(hex, fallback)
    hex = strtrim(char(hex));
    if isempty(regexp(hex, '^#[0-9a-fA-F]{6}$', 'once'))
        rgb = fallback; return;
    end
    rgb = [hex2dec(hex(2:3)), hex2dec(hex(4:5)), hex2dec(hex(6:7))] / 255;
end

function val = get_nested(s, fields, default)
    val = s;
    for i = 1:numel(fields)
        if isstruct(val) && isscalar(val) && isfield(val, fields{i})
            val = val.(fields{i});
        else
            val = default; return;
        end
    end
end

function p = load_yaml(f)
    try
        p = yaml.loadFile(f, 'ConvertToArray', true);
    catch
        % Lists of maps with differing keys (e.g. coil_map entries mixing
        % .target and .sequences) cannot be concatenated; keep them as cells.
        p = yaml.loadFile(f, 'ConvertToArray', false);
    end
end

function remove_dir(d)
    if isfolder(d)
        try, rmdir(d, 's'); catch, end
    end
end

function q = shell_quote(s)
    q = ['''' strrep(s, '''', '''\''''') ''''];
end

function v = ternary(c, a, b)
    if c, v = a; else, v = b; end
end
