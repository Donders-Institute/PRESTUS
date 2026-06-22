function report_path = generate_group_report(parameters, subject_list)
% GENERATE_GROUP_REPORT  Generate a self-contained group-level HTML report
%
% Aggregates per-subject PRESTUS outputs (CSV tables and PNG images) across
% all subjects sharing the current simulation.medium and io.output_affix,
% and assembles a single portable HTML file with:
%   - a cross-subject Exposure Dashboard (min / median / mean / max / SD,
%     colour-coded against ITRUSST consensus limits)
%   - a subject roster table
%   - group-level box plots of acoustic and (if available) thermal metrics
%   - one collapsible card per subject embedding their key images
%   - an interactive filter bar (JavaScript) that toggles subjects in/out
%     and recomputes dashboard stats client-side
%
% Mirrors the look-and-feel of GENERATE_SIMULATION_REPORT and reuses its
% HTML utilities (html_utils, table2html, css_styles_base, risk_color,
% get_risk_limits).
%
% Use as:
%   report_path = generate_group_report(parameters)
%   report_path = generate_group_report(parameters, subject_list)
%
% Input:
%   parameters   - PRESTUS parameters struct (path.sim, simulation.medium,
%                  and io.output_affix are required)
%   subject_list - (optional) numeric array of subject IDs. If omitted or
%                  empty, DISCOVER_GROUP_SUBJECTS is used to auto-scan
%                  path.sim for subjects sharing the medium / affix.
%
% Output:
%   report_path - path to the generated HTML file, or '' if generation
%                 failed or no subjects were found.
%
% Naming convention:
%   <path.sim>/group_<medium>_report<affix>.html
%
% See also: GENERATE_SIMULATION_REPORT, DISCOVER_GROUP_SUBJECTS,
%           PRESTUS_GROUP_REPORT_START, HTML_UTILS, TABLE2HTML

arguments
    parameters   (1,1) struct
    subject_list       = []
end

report_path = '';

try
    %% Resolve subject list (auto-discover if not provided)
    if isempty(subject_list)
        subject_list = discover_group_subjects(parameters);
    end
    subject_list = unique(subject_list(:).');  % dedup + row vector

    if isempty(subject_list)
        fprintf('generate_group_report: no subjects found for medium=%s, affix=%s in %s — skipping.\n', ...
            parameters.simulation.medium, get_affix(parameters), parameters.path.sim);
        return
    end

    %% Determine output path
    medium = parameters.simulation.medium;
    is_layered = contains(medium, {'layered', 'phantom'});
    affix = get_affix(parameters);

    report_filename = sprintf('group_%s_report%s.html', medium, affix);
    report_path = fullfile(parameters.path.sim, report_filename);

    %% Aggregate per-subject data
    [group_table, subjects_meta] = aggregate_subjects(parameters, subject_list, medium, affix);

    if isempty(group_table)
        fprintf('generate_group_report: no per-subject CSVs could be read — skipping.\n');
        report_path = '';
        return
    end

    %% Render group-level box plots (best-effort)
    group_plots_dir = fullfile(parameters.path.sim, 'group_plots');
    if ~exist(group_plots_dir, 'dir'); mkdir(group_plots_dir); end

    plot_path = @(tag) fullfile(group_plots_dir, ...
        sprintf('group_%s_%s%s.png', medium, tag, affix));

    acoustic_specs = struct('path', {}, 'caption', {}, 'alt', {});
    acoustic_specs(end+1) = struct( ...
        'path',    render_group_boxplot(group_table, get_intensity_plot_cols(is_layered), ...
                       'Spatial-peak pulse-average intensity (W/cm^2)', plot_path('intensity')), ...
        'caption', 'Spatial-peak pulse-average intensity (Isppa) per location, in W/cm².', ...
        'alt',     'Intensity box plot across subjects');
    if is_layered
        acoustic_specs(end+1) = struct( ...
            'path',    render_group_boxplot(group_table, get_mi_plot_cols(), ...
                           'Mechanical Index', plot_path('mi')), ...
            'caption', 'Mechanical Index per tissue (unitless; NSR limit 1.9).', ...
            'alt',     'Mechanical Index box plot across subjects');
    end

    thermal_specs = struct('path', {}, 'caption', {}, 'alt', {});
    if is_layered && isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims
        thermal_specs(end+1) = struct( ...
            'path',    render_group_boxplot(group_table, get_maxT_plot_cols(), ...
                           'Maximum temperature (^\circC)', plot_path('maxT')), ...
            'caption', 'Maximum tissue temperature reached during the protocol, in °C.', ...
            'alt',     'Maximum temperature box plot across subjects');
        thermal_specs(end+1) = struct( ...
            'path',    render_group_boxplot(group_table, get_riseT_plot_cols(), ...
                           'Temperature rise above baseline (^\circC)', plot_path('riseT')), ...
            'caption', 'Temperature rise above the 37 °C baseline, in °C.', ...
            'alt',     'Temperature rise box plot across subjects');
        thermal_specs(end+1) = struct( ...
            'path',    render_group_boxplot(group_table, get_cem43_plot_cols(), ...
                           'CEM43 thermal dose (eq. min)', plot_path('cem43')), ...
            'caption', 'CEM43 thermal dose per tissue, in equivalent minutes at 43 °C (includes iso variants).', ...
            'alt',     'CEM43 box plot across subjects');
    end

    %% Build HTML
    limits = get_risk_limits(is_layered);
    n_subjects = numel(subject_list);

    html_parts = {};
    html_parts{end+1} = '<!DOCTYPE html>';
    html_parts{end+1} = '<html lang="en">';
    html_parts{end+1} = '<head>';
    html_parts{end+1} = '<meta charset="UTF-8">';
    html_parts{end+1} = '<meta name="viewport" content="width=device-width, initial-scale=1.0">';
    html_parts{end+1} = sprintf('<title>PRESTUS Group Report — %s%s (N=%d)</title>', ...
        medium, affix, n_subjects);
    html_parts{end+1} = '<style>';
    html_parts{end+1} = css_styles_base();
    html_parts{end+1} = group_css_extra();
    html_parts{end+1} = '</style>';
    html_parts{end+1} = '</head>';
    html_parts{end+1} = '<body>';

    % Table of contents (sticky)
    try
        html_parts{end+1} = build_toc();
    catch ME
        html_parts{end+1} = html_utils.section_error('Table of Contents', ME);
    end

    % Section 1: Header
    try
        html_parts{end+1} = build_group_header(parameters, subject_list, medium, affix, is_layered);
    catch ME
        html_parts{end+1} = html_utils.section_error('Header', ME);
    end

    % Section 2: Filter bar (subject checkboxes)
    try
        html_parts{end+1} = build_filter_bar(subjects_meta);
    catch ME
        html_parts{end+1} = html_utils.section_error('Filter', ME);
    end

    % Section 3: Group Exposure Dashboard
    try
        html_parts{end+1} = build_group_dashboard(subjects_meta, limits, is_layered);
    catch ME
        html_parts{end+1} = html_utils.section_error('Group Exposure Dashboard', ME);
    end

    % Section 4: Subject roster table
    try
        html_parts{end+1} = html_utils.collapsible('Subject Roster', ...
            build_subject_roster(group_table, subjects_meta, parameters, medium, affix, limits, is_layered), ...
            true, 'roster');
    catch ME
        html_parts{end+1} = html_utils.section_error('Subject Roster', ME);
    end

    % Section 5: Group acoustic plots
    try
        html_parts{end+1} = html_utils.collapsible('Group Acoustic Summary', ...
            build_group_plot_section(acoustic_specs, 'No acoustic columns available for plotting.'), ...
            true, 'group-acoustic');
    catch ME
        html_parts{end+1} = html_utils.section_error('Group Acoustic Summary', ME);
    end

    % Section 6: Group thermal plots (conditional)
    if is_layered && isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims
        try
            html_parts{end+1} = html_utils.collapsible('Group Thermal Summary', ...
                build_group_plot_section(thermal_specs, 'No thermal columns available for plotting.'), ...
                true, 'group-thermal');
        catch ME
            html_parts{end+1} = html_utils.section_error('Group Thermal Summary', ME);
        end
    end

    % Section 7: Per-subject cards
    try
        html_parts{end+1} = html_utils.collapsible('Per-Subject Details', ...
            build_subject_cards(subjects_meta, parameters, medium, affix, is_layered), ...
            true, 'subjects');
    catch ME
        html_parts{end+1} = html_utils.section_error('Per-Subject Details', ME);
    end

    % Section 8: Configuration summary
    try
        html_parts{end+1} = html_utils.collapsible('Configuration Summary', ...
            build_group_config_summary(parameters), false, 'config');
    catch ME
        html_parts{end+1} = html_utils.section_error('Configuration Summary', ME);
    end

    % Footer
    html_parts{end+1} = '<footer>';
    html_parts{end+1} = sprintf('<p>Generated by PRESTUS (group report) — %s</p>', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    html_parts{end+1} = '</footer>';

    % Lightbox overlay
    html_parts{end+1} = html_utils.lightbox();

    % Embedded JSON payload + filter JS
    try
        html_parts{end+1} = build_filter_script(subjects_meta, limits);
    catch ME
        html_parts{end+1} = html_utils.section_error('Filter Script', ME);
    end

    html_parts{end+1} = '</body>';
    html_parts{end+1} = '</html>';

    %% Write file
    fid = fopen(report_path, 'w', 'n', 'UTF-8');
    if fid == -1
        warning('generate_group_report:fileOpen', 'Cannot open %s for writing.', report_path);
        report_path = '';
        return
    end
    fprintf(fid, '%s\n', html_parts{:});
    fclose(fid);

    fprintf('Group HTML report saved to: %s\n', report_path);

catch ME
    warning('generate_group_report:failed', ...
        'Group report generation failed: %s\n%s', ME.message, getReport(ME, 'extended'));
    report_path = '';
end

end


%% =========================================================================
%  AGGREGATION
%  =========================================================================

function [group_table, subjects_meta] = aggregate_subjects(parameters, subject_list, medium, affix)
% Read each subject's CSV, vertically concatenate, and build a per-subject
% struct array with the metrics needed for the dashboard and JS filter.

    group_table = table();
    subjects_meta = struct('subject_id', {}, 'csv_path', {}, 'row_index', {}, ...
                           'sub_dir', {}, 'img_dir', {}, 'report_path', {}, ...
                           'metrics', {});

    for k = 1:numel(subject_list)
        id = subject_list(k);
        sub_dir  = fullfile(parameters.path.sim, sprintf('sub-%03d', id));
        img_dir  = fullfile(sub_dir, 'img');
        csv_path = fullfile(sub_dir, sprintf('sub-%03d_%s%s.csv', id, medium, affix));
        rpt_path = fullfile(sub_dir, sprintf('sub-%03d_%s_report%s.html', id, medium, affix));

        if ~isfile(csv_path)
            continue
        end

        try
            t = readtable(csv_path, 'VariableNamingRule', 'preserve');
        catch
            continue
        end
        if isempty(t)
            continue
        end

        % Always pick the last row of each per-subject CSV (matches the
        % per-subject report's csv_value() convention).
        row = t(end, :);
        if ~ismember('subject_id', row.Properties.VariableNames)
            row = addvars(row, id, 'Before', 1, 'NewVariableNames', {'subject_id'});
        end

        % Append to combined table (outer-join semantics — union columns,
        % missing values become NaN/empty). If a per-subject row has an
        % incompatible column type (rare — schemas should match), skip the
        % table append for this subject but keep its struct so it still
        % appears in subjects_meta / the JSON payload.
        try
            if isempty(group_table)
                group_table = row;
            else
                group_table = outerjoin_tables(group_table, row);
            end
        catch ME
            warning('generate_group_report:tableAppend', ...
                ['sub-%03d: CSV row could not be merged into group table (%s); ' ...
                 'subject still included in JSON / dashboard.'], id, ME.message);
        end

        % Build per-subject metrics struct for JSON payload
        metrics = struct('subject_id', id);
        cols = row.Properties.VariableNames;
        for c = 1:numel(cols)
            v = row{1, c};
            if isnumeric(v) && isscalar(v)
                metrics.(matlab.lang.makeValidName(cols{c})) = v;
            end
        end

        subjects_meta(end+1) = struct( ...
            'subject_id', id, ...
            'csv_path',   csv_path, ...
            'row_index',  height(t), ...
            'sub_dir',    sub_dir, ...
            'img_dir',    img_dir, ...
            'report_path', rpt_path, ...
            'metrics',    metrics); %#ok<AGROW>
    end
end

function T = outerjoin_tables(A, B)
% Append B's row(s) to A, padding missing columns with NaN/missing.
    colsA = A.Properties.VariableNames;
    colsB = B.Properties.VariableNames;
    all_cols = unique([colsA, colsB], 'stable');

    for c = 1:numel(all_cols)
        col = all_cols{c};
        if ~ismember(col, colsA)
            % Add column to A as missing
            A = addvars(A, repmat(missing_like(B.(col)), height(A), 1), 'NewVariableNames', col);
        end
        if ~ismember(col, colsB)
            B = addvars(B, repmat(missing_like(A.(col)), height(B), 1), 'NewVariableNames', col);
        end
    end
    % Reorder B to match A
    B = B(:, A.Properties.VariableNames);
    T = [A; B];
end

function m = missing_like(v)
    if isnumeric(v)
        m = NaN;
    elseif iscell(v)
        m = {''};
    elseif isstring(v)
        m = string(missing);
    else
        m = missing;
    end
end


%% =========================================================================
%  GROUP BOX-PLOT RENDERING
%  =========================================================================

function plot_path = render_group_boxplot(group_table, cols, title_str, out_path)
% Render a multi-panel box plot of the requested columns and save as PNG.
% Returns the saved path on success or '' on failure / no usable columns.

    plot_path = '';
    try
        avail = intersect(cols, group_table.Properties.VariableNames, 'stable');
        keep = false(1, numel(avail));
        for i = 1:numel(avail)
            v = group_table.(avail{i});
            keep(i) = isnumeric(v) && ~all(isnan(v));
        end
        avail = avail(keep);
        if isempty(avail), return; end

        fig = figure('Visible', 'off', 'Color', 'w', ...
                     'Position', [100 100 max(600, 140 * numel(avail)) 420]);
        ax = axes(fig); %#ok<LAXES>

        all_vals = [];
        all_grp  = [];
        for i = 1:numel(avail)
            v = group_table.(avail{i});
            v = v(~isnan(v));
            all_vals = [all_vals; v(:)]; %#ok<AGROW>
            all_grp  = [all_grp;  repmat(i, numel(v), 1)]; %#ok<AGROW>
        end

        if exist('boxchart', 'file') == 2 || exist('boxchart', 'builtin') == 5
            boxchart(ax, all_grp, all_vals);
        else
            boxplot(ax, all_vals, all_grp);
        end
        set(ax, 'XTick', 1:numel(avail), 'XTickLabel', avail, ...
                'XTickLabelRotation', 30, 'TickLabelInterpreter', 'none', ...
                'FontSize', 9);
        title(ax, title_str, 'Interpreter', 'none', 'FontSize', 11);
        grid(ax, 'on');

        if exist('exportgraphics', 'file') == 2 || exist('exportgraphics', 'builtin') == 5
            exportgraphics(ax, out_path, 'Resolution', 150);
        else
            print(fig, out_path, '-dpng', '-r150');
        end
        close(fig);
        plot_path = out_path;
    catch ME
        warning('generate_group_report:plotFailed', ...
            'Box plot rendering failed: %s', ME.message);
    end
end

function cols = get_intensity_plot_cols(is_layered)
    cols = {'Isppa', 'Ipa_target'};
    if is_layered
        cols = [cols, {'Isppa_brain', 'Isppa_skull', 'Isppa_skin'}];
    end
end

function cols = get_mi_plot_cols()
    cols = {'MI_tc', 'MI_brain', 'MI_skull', 'MI_skin'};
end

function cols = get_maxT_plot_cols()
    cols = {'maxT', 'maxT_brain', 'maxT_skull', 'maxT_skin'};
end

function cols = get_riseT_plot_cols()
    cols = {'riseT_brain', 'riseT_skull', 'riseT_skin'};
end

function cols = get_cem43_plot_cols()
    cols = {'CEM43_brain', 'CEM43_skull', 'CEM43_skin', ...
            'CEM43iso_brain', 'CEM43iso_skull', 'CEM43iso_skin'};
end


%% =========================================================================
%  SECTION BUILDERS
%  =========================================================================

function html = build_toc()
    html = ['<nav class="toc">' ...
        '<strong>Group Report:</strong>' ...
        '<a href="#header">Header</a>' ...
        '<a href="#filter">Filter</a>' ...
        '<a href="#dashboard">Dashboard</a>' ...
        '<a href="#roster">Roster</a>' ...
        '<a href="#group-acoustic">Acoustic</a>' ...
        '<a href="#group-thermal">Thermal</a>' ...
        '<a href="#subjects">Subjects</a>' ...
        '<a href="#config">Config</a>' ...
        '</nav>'];
end

function html = build_group_header(parameters, subject_list, medium, affix, is_layered)
    n = numel(subject_list);
    html = '<section class="report-section" id="header">';
    html = [html '<h1>PRESTUS Group Report</h1>'];
    html = [html '<table class="info-table">'];
    html = [html sprintf('<tr><th>Subjects</th><td><span id="n-active">%d</span> of %d active</td></tr>', n, n)];
    html = [html sprintf('<tr><th>Subject IDs</th><td>%s</td></tr>', ...
        html_utils.escape(strtrim(sprintf('%d ', subject_list))))];
    if is_layered
        html = [html sprintf('<tr><th>Medium</th><td>%s <span class="badge badge-amber">Layered</span></td></tr>', ...
            html_utils.escape(medium))];
    else
        html = [html sprintf('<tr><th>Medium</th><td>%s <span class="badge badge-blue">Water / Free-field</span></td></tr>', ...
            html_utils.escape(medium))];
    end
    if ~isempty(affix)
        html = [html sprintf('<tr><th>Affix</th><td>%s</td></tr>', html_utils.escape(affix))];
    end
    html = [html sprintf('<tr><th>Sim path</th><td>%s</td></tr>', html_utils.escape(parameters.path.sim))];
    html = [html sprintf('<tr><th>Generated</th><td>%s</td></tr>', datestr(now, 'yyyy-mm-dd HH:MM:SS'))];
    html = [html '</table>'];
    html = [html '<p class="note">Use the filter bar below to toggle subjects in/out. ' ...
        'Dashboard statistics and per-subject cards update live; group box plots reflect the full set.</p>'];
    html = [html '</section>'];
end

function html = build_filter_bar(subjects_meta)
    html = '<section class="report-section filter-bar" id="filter">';
    html = [html '<h2>Subject Filter</h2>'];
    html = [html '<div class="filter-controls">'];
    html = [html '<button type="button" class="filter-btn" onclick="setAllSubjects(true)">All</button>'];
    html = [html '<button type="button" class="filter-btn" onclick="setAllSubjects(false)">None</button>'];
    html = [html '<span class="filter-hint">Toggle subjects to update dashboard and per-subject cards.</span>'];
    html = [html '</div>'];
    html = [html '<div class="filter-checkboxes">'];
    for k = 1:numel(subjects_meta)
        id = subjects_meta(k).subject_id;
        html = [html sprintf( ...
            ['<label class="filter-chk"><input type="checkbox" class="subj-toggle" ' ...
             'value="%d" checked> sub-%03d</label>'], id, id)];
    end
    html = [html '</div>'];
    html = [html '</section>'];
end

function html = build_group_dashboard(subjects_meta, limits, is_layered)
% Cross-subject statistics per safety metric, with data-* attributes so the
% JS filter can recompute on toggle.
    html = '<section class="report-section" id="dashboard">';
    html = [html '<h2>Group Exposure Dashboard</h2>'];

    if ~is_layered
        html = [html '<div class="medium-banner medium-water">' ...
            '<strong>Water / Free-field</strong> &mdash; tissue-specific NSR limits are not applicable.</div>'];
    end

    metric_names = fieldnames(limits);

    html = [html '<div class="safety-grid">'];
    for i = 1:numel(metric_names)
        name = metric_names{i};
        info = limits.(name);
        js_key = matlab.lang.makeValidName(name);
        vals = collect_metric(subjects_meta, js_key);
        s = stat_summary(vals);

        % Colour driven by max (worst case across subjects); matches the
        % uncertainty report's "conservative variant" colour rule.
        color = risk_color(s.max, info.limit);

        % Axis baseline: 37 for absolute temperatures, 0 elsewhere
        if strncmp(name, 'maxT', 4) || strncmp(name, 'endT', 4)
            scale_min = 37;
        else
            scale_min = 0;
        end
        % Scale max: large enough to fit both worst-case sim and the NSR limit.
        if isinf(info.limit)
            scale_max = max(scale_min + 1, s.max);
        else
            scale_max = max([scale_min + 1, s.max, info.limit]);
        end
        if isnan(scale_max), scale_max = scale_min + 1; end
        scale_range = scale_max - scale_min;

        html = [html sprintf( ...
            '<div class="safety-card safety-%s" data-card-metric="%s" data-scale-min="%.6g" data-scale-max="%.6g" data-limit="%s">', ...
            color, html_utils.escape(js_key), scale_min, scale_max, ...
            iff(isinf(info.limit), 'null', sprintf('%.6g', info.limit)))];

        html = [html sprintf('<div class="safety-label">%s</div>', html_utils.escape(info.label))];

        % Main value (mean) + unit
        if isnan(s.mean)
            html = [html '<div class="safety-value" data-stat-mean>N/A</div>'];
        else
            html = [html sprintf('<div class="safety-value" data-stat-mean>%.3g</div>', s.mean)];
        end

        % ± SD
        if isnan(s.sd) || s.n < 1
            html = [html '<div class="safety-sd" data-stat-sd>&plusmn; &mdash;</div>'];
        else
            html = [html sprintf('<div class="safety-sd" data-stat-sd>&plusmn; %.3g %s</div>', s.sd, html_utils.escape(info.unit))];
        end

        % min – max (N=n)
        if s.n == 0
            range_str = '(N=0)';
        elseif s.n == 1
            range_str = sprintf('%.3g (N=1)', s.min);
        else
            range_str = sprintf('%.3g &ndash; %.3g (N=%d)', s.min, s.max, s.n);
        end
        html = [html sprintf('<div class="safety-range" data-stat-range>%s</div>', range_str)];

        % NSR limit annotation
        if isinf(info.limit)
            html = [html sprintf('<div class="risk-limit">%s (informational)</div>', html_utils.escape(info.unit))];
        else
            html = [html sprintf('<div class="risk-limit">NSR limit: %.3g %s</div>', info.limit, html_utils.escape(info.unit))];
        end

        % Axis range printed above the bar
        html = [html sprintf('<div class="safety-bar-range">%.3g &ndash; %.3g %s</div>', ...
            scale_min, scale_max, html_utils.escape(info.unit))];

        % Bar track: span (min→max), mean marker, NSR limit line
        bar_html = '<div class="safety-bar-track" data-bar-track>';
        if ~isnan(s.min) && ~isnan(s.max)
            bar_left  = clamp((s.min - scale_min) / scale_range * 100, 0, 100);
            bar_width = max(0, clamp((s.max - scale_min) / scale_range * 100, 0, 100) - bar_left);
        else
            bar_left = 0; bar_width = 0;
        end
        bar_html = [bar_html sprintf( ...
            '<div class="safety-bar" data-bar-span style="left:%.2f%%;width:%.2f%%"></div>', ...
            bar_left, bar_width)];
        if ~isnan(s.mean)
            mean_pct = clamp((s.mean - scale_min) / scale_range * 100, 0, 100);
            bar_html = [bar_html sprintf( ...
                '<div class="safety-bar-marker safety-bar-marker--default" data-bar-mean style="left:%.2f%%"></div>', mean_pct)];
        else
            bar_html = [bar_html '<div class="safety-bar-marker safety-bar-marker--default" data-bar-mean style="display:none"></div>'];
        end
        if ~isinf(info.limit)
            lim_pct = clamp((info.limit - scale_min) / scale_range * 100, 0, 100);
            bar_html = [bar_html sprintf( ...
                '<div class="safety-bar-limit" style="left:%.2f%%"></div>', lim_pct)];
        end
        bar_html = [bar_html '</div>'];
        html = [html bar_html];

        html = [html '</div>'];   % /safety-card
    end
    html = [html '</div>'];
    html = [html '<p class="safety-footnote">NSR limits per ITRUSST consensus ' ...
        '(Aubry et al., 2025). One tile per metric. ' ...
        '<strong>Number</strong> = mean across <em>currently selected</em> subjects; ' ...
        '<strong>&plusmn;</strong> = SD. ' ...
        '<strong>Range</strong> = min&ndash;max with N. ' ...
        '<strong>Bar</strong> spans min&ndash;max; <strong>&#x25A0;</strong> = mean; ' ...
        'dashed grey line = NSR limit. ' ...
        'Card colour reflects the <strong>worst-case subject</strong> (max) against the NSR limit: ' ...
        'green &lt;50% of limit, amber 50&ndash;100%, red exceeds.</p>'];
    html = [html '</section>'];
end

function y = clamp(x, lo, hi)
    if isnan(x), y = lo; return; end
    y = min(hi, max(lo, x));
end

function s = iff(cond, a, b)
    if cond, s = a; else, s = b; end
end

function vals = collect_metric(subjects_meta, key)
    vals = [];
    for k = 1:numel(subjects_meta)
        m = subjects_meta(k).metrics;
        if isfield(m, key) && isnumeric(m.(key)) && isscalar(m.(key))
            vals(end+1) = m.(key); %#ok<AGROW>
        end
    end
end

function s = stat_summary(v)
    v = v(~isnan(v));
    if isempty(v)
        s = struct('mean', NaN, 'median', NaN, 'min', NaN, 'max', NaN, 'sd', NaN, 'n', 0);
    else
        s = struct('mean', mean(v), 'median', median(v), 'min', min(v), ...
                   'max', max(v), 'sd', std(v), 'n', numel(v));
    end
end

function html = build_subject_roster(group_table, subjects_meta, parameters, medium, affix, limits, is_layered) %#ok<INUSL>
% One row per subject with key metrics + a link to the per-subject report
% (if present). Each row carries data-subject-id so the JS filter hides/shows.
    if isempty(group_table)
        html = '<p class="placeholder">No subject data.</p>';
        return
    end

    % Pick a small set of headline columns to keep the table scannable
    headline = {'subject_id', 'Isppa', 'Ipa_target', 'real_focal_distance_mm'};
    if is_layered
        headline = [headline, {'Isppa_brain', 'MI_tc', 'MI_brain'}];
        if any(ismember({'maxT', 'maxT_brain'}, group_table.Properties.VariableNames))
            headline = [headline, {'maxT_brain', 'CEM43_brain'}];
        end
    end
    avail = intersect(headline, group_table.Properties.VariableNames, 'stable');

    % Build HTML manually so each row gets data-subject-id (table2html can't do that)
    html = '<div class="table-wrapper"><table class="data-table"><thead><tr>';
    html = [html '<th>Subject</th>'];
    for c = 1:numel(avail)
        if strcmp(avail{c}, 'subject_id'), continue; end
        html = [html sprintf('<th>%s</th>', html_utils.escape(avail{c}))];
    end
    html = [html '<th>Report</th>'];
    html = [html '</tr></thead><tbody>'];

    for k = 1:numel(subjects_meta)
        id = subjects_meta(k).subject_id;
        % Find this subject's row in group_table
        row_idx = find(group_table.subject_id == id, 1, 'last');
        if isempty(row_idx)
            continue
        end
        html = [html sprintf('<tr data-subject-id="%d">', id)];
        html = [html sprintf('<td>sub-%03d</td>', id)];
        for c = 1:numel(avail)
            col = avail{c};
            if strcmp(col, 'subject_id'), continue; end
            val = group_table{row_idx, col};
            if iscell(val); val = val{1}; end
            val_str = html_utils.format_cell(val);

            cell_class = '';
            if isfield(limits, col) && isnumeric(val) && isscalar(val)
                color = risk_color(val, limits.(col).limit);
                if ~strcmp(color, 'info')
                    cell_class = sprintf(' class="cell-%s"', color);
                end
            end
            html = [html sprintf('<td%s>%s</td>', cell_class, val_str)];
        end

        rpt = subjects_meta(k).report_path;
        if isfile(rpt)
            html = [html sprintf('<td><a href="%s" target="_blank">open</a></td>', ...
                html_utils.escape(relative_path(parameters.path.sim, rpt)))];
        else
            html = [html '<td><span class="placeholder">&mdash;</span></td>'];
        end
        html = [html '</tr>'];
    end
    html = [html '</tbody></table></div>'];
end

function html = build_group_plot_section(plot_specs, fallback_msg)
% plot_specs: struct array with fields .path, .caption, .alt — one entry per box plot.
    valid = false(1, numel(plot_specs));
    for k = 1:numel(plot_specs)
        valid(k) = ~isempty(plot_specs(k).path) && isfile(plot_specs(k).path);
    end
    if ~any(valid)
        html = sprintf('<p class="placeholder">%s</p>', html_utils.escape(fallback_msg));
        return
    end

    html = '<div class="image-grid">';
    for k = 1:numel(plot_specs)
        if valid(k)
            html = [html html_utils.embed_image(plot_specs(k).path, ...
                plot_specs(k).alt, plot_specs(k).caption)];
        else
            html = [html sprintf( ...
                '<figure><div class="placeholder-img">(missing)</div><figcaption>%s</figcaption></figure>', ...
                html_utils.escape(plot_specs(k).caption))];
        end
    end
    html = [html '</div>'];
    html = [html '<p class="note">Group box plots reflect the full subject set; the JS filter only updates the dashboard and per-subject cards.</p>'];
end

function html = build_subject_cards(subjects_meta, parameters, medium, affix, is_layered) %#ok<INUSL>
% One collapsible card per subject, embedding key per-subject images. Each
% card has data-subject-id so the JS filter can hide/show them.
    html = '';
    for k = 1:numel(subjects_meta)
        id = subjects_meta(k).subject_id;
        img_dir = subjects_meta(k).img_dir;

        body = '<div class="image-grid">';
        body = [body try_embed(fullfile(img_dir, ...
            sprintf('sub-%03d_positioning%s.png', id, affix)), 'Positioning', 'Transducer positioning')];
        body = [body try_embed(fullfile(img_dir, ...
            sprintf('sub-%03d_%s_intensity%s.png', id, medium, affix)), 'Intensity', 'Intensity overlay (segmentation)')];
        body = [body try_embed(fullfile(img_dir, ...
            sprintf('sub-%03d_%s_intensity_t1%s.png', id, medium, affix)), 'Intensity T1', 'Intensity overlay (T1)')];
        if is_layered
            body = [body try_embed(fullfile(img_dir, ...
                sprintf('sub-%03d_%s_maxT%s.png', id, medium, affix)), 'maxT', 'Max temperature overlay')];
        end
        body = [body '</div>'];

        % Per-subject report link
        rpt = subjects_meta(k).report_path;
        if isfile(rpt)
            body = [body sprintf('<p class="note">Per-subject report: <a href="%s" target="_blank">%s</a></p>', ...
                html_utils.escape(relative_path(parameters.path.sim, rpt)), ...
                html_utils.escape(relative_path(parameters.path.sim, rpt)))];
        end

        html = [html sprintf('<details class="subject-card" data-subject-id="%d" open>', id)];
        html = [html sprintf('<summary><strong>sub-%03d</strong></summary>', id)];
        html = [html '<div class="section-content">'];
        html = [html body];
        html = [html '</div></details>'];
    end

    if isempty(html)
        html = '<p class="placeholder">No subjects.</p>';
    end
end

function s = try_embed(path, alt, cap)
    s = html_utils.embed_image(path, alt, cap);
    if isempty(s)
        s = sprintf('<figure><div class="placeholder-img">%s<br><small>(missing)</small></div><figcaption>%s</figcaption></figure>', ...
            html_utils.escape(alt), html_utils.escape(cap));
    end
end

function html = build_group_config_summary(parameters)
    rows = {};
    rows{end+1} = {'path.sim',           parameters.path.sim};
    rows{end+1} = {'simulation.medium',  parameters.simulation.medium};
    rows{end+1} = {'io.output_affix',    get_affix(parameters)};
    if isfield(parameters, 'grid') && isfield(parameters.grid, 'resolution_mm')
        rows{end+1} = {'grid.resolution_mm', num2str(parameters.grid.resolution_mm)};
    end
    if isfield(parameters, 'transducer') && ~isempty(parameters.transducer)
        td = parameters.transducer(1);
        if isfield(td, 'name')
            rows{end+1} = {'transducer(1).name', char(td.name)};
        end
        if isfield(td, 'freq_hz')
            rows{end+1} = {'transducer(1).freq_hz', num2str(td.freq_hz)};
        end
    end

    html = '<table class="info-table">';
    for i = 1:numel(rows)
        html = [html sprintf('<tr><th>%s</th><td>%s</td></tr>', ...
            html_utils.escape(rows{i}{1}), html_utils.escape(rows{i}{2}))]; %#ok<AGROW>
    end
    html = [html '</table>'];
    html = [html '<p class="note">Shown for the loaded parameters struct. Per-subject configs may differ; consult each subject''s own report for authoritative settings.</p>'];
end

function html = build_filter_script(subjects_meta, limits)
% Embed a JSON payload of per-subject metrics and a JS filter that
% recomputes dashboard stats / shows/hides subject cards and roster rows.

    % Build a slim JSON: one entry per subject with all numeric metrics
    payload = cell(1, numel(subjects_meta));
    for k = 1:numel(subjects_meta)
        payload{k} = subjects_meta(k).metrics;
    end
    json = jsonencode(payload);

    % Build a JS-friendly limits object
    lim_fields = fieldnames(limits);
    lim_js = cell(1, numel(lim_fields));
    for i = 1:numel(lim_fields)
        f = lim_fields{i};
        l = limits.(f).limit;
        if isinf(l)
            lim_js{i} = sprintf('"%s":null', matlab.lang.makeValidName(f));
        else
            lim_js{i} = sprintf('"%s":%.6g', matlab.lang.makeValidName(f), l);
        end
    end
    lim_json = ['{' strjoin(lim_js, ',') '}'];

    html = ['<script>' newline ...
        'const SUBJECTS = ' json ';' newline ...
        'const LIMITS   = ' lim_json ';' newline ...
        '' newline ...
        'function activeIds() {' newline ...
        '  return Array.from(document.querySelectorAll(".subj-toggle:checked"))' newline ...
        '              .map(cb => parseInt(cb.value, 10));' newline ...
        '}' newline ...
        '' newline ...
        'function setAllSubjects(on) {' newline ...
        '  document.querySelectorAll(".subj-toggle").forEach(cb => cb.checked = on);' newline ...
        '  recompute();' newline ...
        '}' newline ...
        '' newline ...
        'function mean(vs)   { return vs.reduce((a,b)=>a+b,0) / vs.length; }' newline ...
        'function sd(vs)     { if (vs.length < 2) return 0; const m = mean(vs); return Math.sqrt(vs.reduce((a,b)=>a+(b-m)*(b-m),0) / (vs.length-1)); }' newline ...
        'function vmin(vs)   { return Math.min.apply(null, vs); }' newline ...
        'function vmax(vs)   { return Math.max.apply(null, vs); }' newline ...
        'function fmt(v)     { return Number.isFinite(v) ? v.toPrecision(3) : "N/A"; }' newline ...
        'function clamp(x,lo,hi) { return Math.max(lo, Math.min(hi, x)); }' newline ...
        '' newline ...
        'function riskClass(value, limit) {' newline ...
        '  if (value === null || value === undefined || Number.isNaN(value)) return "safety-gray";' newline ...
        '  if (limit === null || limit === undefined)                        return "safety-info";' newline ...
        '  if (value > limit)         return "safety-red";' newline ...
        '  if (value >= 0.5 * limit)  return "safety-amber";' newline ...
        '  return "safety-green";' newline ...
        '}' newline ...
        '' newline ...
        'function recompute() {' newline ...
        '  const enabled = activeIds();' newline ...
        '  const enabledSet = new Set(enabled);' newline ...
        '' newline ...
        '  // Show/hide per-subject cards and roster rows' newline ...
        '  document.querySelectorAll("[data-subject-id]").forEach(el => {' newline ...
        '    const id = parseInt(el.dataset.subjectId, 10);' newline ...
        '    el.style.display = enabledSet.has(id) ? "" : "none";' newline ...
        '  });' newline ...
        '' newline ...
        '  const active = SUBJECTS.filter(s => enabledSet.has(s.subject_id));' newline ...
        '  const nActive = document.getElementById("n-active");' newline ...
        '  if (nActive) nActive.textContent = active.length;' newline ...
        '' newline ...
        '  // Recompute each dashboard tile' newline ...
        '  document.querySelectorAll(".safety-card[data-card-metric]").forEach(card => {' newline ...
        '    const metric    = card.dataset.cardMetric;' newline ...
        '    const scaleMin  = parseFloat(card.dataset.scaleMin);' newline ...
        '    const scaleMax  = parseFloat(card.dataset.scaleMax);' newline ...
        '    const scaleRng  = Math.max(1e-12, scaleMax - scaleMin);' newline ...
        '    const limitRaw  = card.dataset.limit;' newline ...
        '    const limit     = (limitRaw === "null" || limitRaw === undefined) ? null : parseFloat(limitRaw);' newline ...
        '' newline ...
        '    const vals = active.map(s => s[metric]).filter(v => v !== null && v !== undefined && !Number.isNaN(v));' newline ...
        '    const n    = vals.length;' newline ...
        '    const m    = n ? mean(vals) : NaN;' newline ...
        '    const sdv  = n ? sd(vals)   : NaN;' newline ...
        '    const lo   = n ? vmin(vals) : NaN;' newline ...
        '    const hi   = n ? vmax(vals) : NaN;' newline ...
        '' newline ...
        '    // Main value (mean)' newline ...
        '    const meanEl = card.querySelector("[data-stat-mean]");' newline ...
        '    if (meanEl) meanEl.textContent = n ? fmt(m) : "N/A";' newline ...
        '' newline ...
        '    // ± SD' newline ...
        '    const sdEl = card.querySelector("[data-stat-sd]");' newline ...
        '    if (sdEl) sdEl.innerHTML = n ? ("± " + fmt(sdv)) : "± —";' newline ...
        '' newline ...
        '    // Range and N' newline ...
        '    const rngEl = card.querySelector("[data-stat-range]");' newline ...
        '    if (rngEl) {' newline ...
        '      if (n === 0) rngEl.innerHTML = "(N=0)";' newline ...
        '      else if (n === 1) rngEl.innerHTML = fmt(lo) + " (N=1)";' newline ...
        '      else rngEl.innerHTML = fmt(lo) + " – " + fmt(hi) + " (N=" + n + ")";' newline ...
        '    }' newline ...
        '' newline ...
        '    // Bar geometry: span lo..hi' newline ...
        '    const span = card.querySelector("[data-bar-span]");' newline ...
        '    if (span) {' newline ...
        '      if (n) {' newline ...
        '        const left  = clamp((lo - scaleMin) / scaleRng * 100, 0, 100);' newline ...
        '        const right = clamp((hi - scaleMin) / scaleRng * 100, 0, 100);' newline ...
        '        span.style.left  = left.toFixed(2) + "%";' newline ...
        '        span.style.width = Math.max(0, right - left).toFixed(2) + "%";' newline ...
        '      } else {' newline ...
        '        span.style.width = "0%";' newline ...
        '      }' newline ...
        '    }' newline ...
        '' newline ...
        '    // Mean marker' newline ...
        '    const mk = card.querySelector("[data-bar-mean]");' newline ...
        '    if (mk) {' newline ...
        '      if (n) {' newline ...
        '        mk.style.display = "";' newline ...
        '        mk.style.left = clamp((m - scaleMin) / scaleRng * 100, 0, 100).toFixed(2) + "%";' newline ...
        '      } else {' newline ...
        '        mk.style.display = "none";' newline ...
        '      }' newline ...
        '    }' newline ...
        '' newline ...
        '    // Re-colour by max (worst case)' newline ...
        '    card.classList.remove("safety-green","safety-amber","safety-red","safety-gray","safety-info");' newline ...
        '    card.classList.add(riskClass(hi, limit));' newline ...
        '  });' newline ...
        '}' newline ...
        '' newline ...
        'document.querySelectorAll(".subj-toggle").forEach(cb => cb.addEventListener("change", recompute));' newline ...
        'recompute();  // initial render' newline ...
        '</script>'];
end


%% =========================================================================
%  CSS (extends css_styles_base)
%  =========================================================================

function css = group_css_extra()
    css = [...
        '/* Group filter bar */' newline ...
        '.filter-bar { position: sticky; top: 60px; z-index: 90; }' newline ...
        '.filter-controls { display: flex; align-items: center; gap: 10px; margin: 8px 0; flex-wrap: wrap; }' newline ...
        '.filter-btn { padding: 4px 12px; border-radius: 4px; border: 1px solid #cbd5e1; background: #f8fafc; cursor: pointer; font-size: 0.85em; }' newline ...
        '.filter-btn:hover { background: #e2e8f0; }' newline ...
        '.filter-hint { font-size: 0.8em; color: #94a3b8; }' newline ...
        '.filter-checkboxes { display: flex; flex-wrap: wrap; gap: 8px 14px; }' newline ...
        '.filter-chk { font-size: 0.85em; color: #334155; user-select: none; cursor: pointer; }' newline ...
        '.filter-chk input { margin-right: 4px; vertical-align: middle; }' newline ...
        '' newline ...
        '/* Subject cards */' newline ...
        'details.subject-card { background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 0; margin-bottom: 12px; }' newline ...
        'details.subject-card > summary { cursor: pointer; list-style: none; padding: 10px 16px; }' newline ...
        'details.subject-card > summary::-webkit-details-marker { display: none; }' newline ...
        'details.subject-card > summary::before { content: "\25B6"; font-size: 0.7em; margin-right: 8px; transition: transform 0.2s ease; }' newline ...
        'details[open].subject-card > summary::before { transform: rotate(90deg); }' newline ...
        'details.subject-card > .section-content { padding: 0 16px 16px 16px; }' newline ...
        '' newline ...
        '/* Group dashboard: compact tile additions */' newline ...
        '.safety-sd    { font-size: 0.78em; color: #475569; margin: 2px 0; }' newline ...
        '.safety-range { font-size: 0.78em; color: #64748b; margin-bottom: 4px; }' newline ...
        '' newline ...
        '/* Safety bar (lifted from uncertainty report) */' newline ...
        '.safety-bar-range { font-size: 0.72em; color: #64748b; margin-top: 6px; margin-bottom: 2px; }' newline ...
        '.safety-bar-track { position: relative; height: 4px; background: #e2e8f0; border-radius: 2px; overflow: visible; }' newline ...
        '.safety-bar { height: 100%; border-radius: 2px; position: absolute; z-index: 1; }' newline ...
        '.safety-bar-marker { position: absolute; top: -5px; width: 4px; height: 14px; border-radius: 1px; transform: translateX(-50%); z-index: 2; }' newline ...
        '.safety-bar-marker--default { background: #1e293b; }' newline ...
        '.safety-bar-limit { position: absolute; top: -3px; width: 2px; height: 10px; background: #94a3b8; border-radius: 0; transform: translateX(-50%); z-index: 1; border-left: 1px dashed #64748b; }' newline ...
        '.safety-green .safety-bar { background: var(--green); }' newline ...
        '.safety-amber .safety-bar { background: var(--amber); }' newline ...
        '.safety-red   .safety-bar { background: var(--red); }' newline ...
        '.safety-gray  .safety-bar { background: var(--gray); }' newline ...
        '.safety-info  .safety-bar { background: var(--info); }' newline ...
        '' newline ...
        '/* Cell colour classes used by table2html */' newline ...
        '.cell-green { background: #f0fdf4; }' newline ...
        '.cell-amber { background: #fffbeb; }' newline ...
        '.cell-red   { background: #fef2f2; font-weight: 600; }' newline ...
        '.cell-gray  { background: #f9fafb; color: #94a3b8; }' newline ...
        '' newline ...
        '/* Badges (compact) */' newline ...
        '.badge { display: inline-block; padding: 2px 8px; border-radius: 4px; font-size: 0.75em; margin-left: 6px; }' newline ...
        '.badge-amber { background: #fffbeb; color: #92400e; border: 1px solid #fde68a; }' newline ...
        '.badge-blue  { background: #eff6ff; color: #1e40af; border: 1px solid #bfdbfe; }' newline ...
        '.badge-gray  { background: #f9fafb; color: #475569; border: 1px solid #e2e8f0; }' newline ...
        ];
end


%% =========================================================================
%  MISC HELPERS
%  =========================================================================

function affix = get_affix(parameters)
    if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix') ...
            && ~isempty(parameters.io.output_affix)
        affix = parameters.io.output_affix;
    else
        affix = '';
    end
end

function rel = relative_path(base, target)
% Best-effort relative path so HTML links resolve when the report is opened
% under <sim>/. Falls back to the absolute path on any mismatch.
    base = char(base);
    target = char(target);
    if startsWith(target, base)
        rel = strrep(target(numel(base)+1:end), '\', '/');
        if startsWith(rel, '/'); rel = rel(2:end); end
    else
        rel = strrep(target, '\', '/');
    end
end
