function report_path = generate_group_report(parameters, runs, options)
% GENERATE_GROUP_REPORT  Generate a self-contained group-level HTML report
%
% Aggregates per-run PRESTUS outputs (CSV tables, targeting sidecars, PNG
% images, group MNI maps) and assembles a single portable HTML file with:
%   - a cross-run Exposure Dashboard (min / median / mean / max / SD,
%     colour-coded against ITRUSST consensus limits)
%   - study design: conditions and the per-run settings (placement, target,
%     transducer, depth, intensity, manifest overrides)      [study mode]
%   - condition comparison: per-condition statistics and box plots
%                                                        [> 1 condition]
%   - targeting accuracy: simulated focus vs planned focus / MNI target / ROI
%                                         [targeting sidecars present]
%   - group MNI maps: focal-volume overlap per condition     [maps given]
%   - outlier flags: robust z-score within each condition
%   - a run roster, group box plots, and one card per run
%   - an interactive filter bar (JavaScript) that toggles runs in/out and
%     recomputes dashboard stats client-side
%
% A run is one subject simulated with one set of settings (medium + output
% affix). Without a study, each subject is one run and the report covers all
% subjects sharing simulation.medium and io.output_affix. With a study
% (PRESTUS_GROUP_REPORT_START on a study config), runs come from the
% manifest and may differ in any setting; conditions group them.
%
% Mirrors the look-and-feel of GENERATE_SIMULATION_REPORT and reuses its
% HTML utilities (html_utils, table2html, css_styles_base, risk_color,
% get_risk_limits).
%
% Use as:
%   report_path = generate_group_report(parameters)
%   report_path = generate_group_report(parameters, subject_list)
%   report_path = generate_group_report(parameters, runs, options)
%
% Input:
%   parameters   - PRESTUS parameters struct (path.sim, simulation.medium,
%                  and io.output_affix are required)
%   runs         - (optional) numeric array of subject IDs, or a run table
%                  from STUDY_RUN_TABLE (included rows are reported). If
%                  omitted or empty, DISCOVER_GROUP_SUBJECTS is used to
%                  auto-scan path.sim for subjects sharing the medium / affix.
%   options      - (optional) struct, study mode:
%                    .name      study name (report file name)
%                    .out_dir   output folder (default path.sim)
%                    .maps      struct array from CREATE_GROUP_MAPS
%                    .manifest  manifest path (shown in the header)
%
% Output:
%   report_path - path to the generated HTML file, or '' if generation
%                 failed or no runs were found.
%
% Naming convention:
%   <path.sim>/group_<medium>_report<affix>.html        (subject list)
%   <out_dir>/group_<name>_report.html                  (study)
% Box plots go to <path.sim>/group_plots/ or <out_dir>/plots/; a tidy table
% with one row per run is written next to the report (…_runs.csv).
%
% See also: GENERATE_SIMULATION_REPORT, DISCOVER_GROUP_SUBJECTS,
%           PRESTUS_GROUP_REPORT_START, CREATE_GROUP_MAPS, STUDY_RUN_TABLE,
%           TARGETING_ANALYSIS, HTML_UTILS, TABLE2HTML

arguments
    parameters   (1,1) struct
    runs               = []
    options      (1,1) struct = struct()
end

report_path = '';

try
    %% Resolve runs (auto-discover subjects if not provided)
    if ~istable(runs)
        subject_list = runs;
        if isempty(subject_list)
            subject_list = discover_group_subjects(parameters);
        end
        subject_list = unique(subject_list(:).');  % dedup + row vector
        runs = legacy_run_table(parameters, subject_list);
    elseif ismember('include', runs.Properties.VariableNames)
        runs = runs(runs.include, :);
    end

    if isempty(runs)
        fprintf('generate_group_report: no subjects found for medium=%s, affix=%s in %s — skipping.\n', ...
            parameters.simulation.medium, get_affix(parameters), parameters.path.sim);
        return
    end

    %% Determine output path
    medium = parameters.simulation.medium;
    is_layered = contains(medium, {'layered', 'phantom'});
    affix = get_affix(parameters);
    study_mode = isfield(options, 'name') && ~isempty(options.name);

    if study_mode
        name = char(options.name);
        out_dir = char(parameters.path.sim);
        if isfield(options, 'out_dir') && ~isempty(options.out_dir); out_dir = char(options.out_dir); end
        report_path = fullfile(out_dir, sprintf('group_%s_report.html', name));
        group_plots_dir = fullfile(out_dir, 'plots');
        plot_suffix = ['_' name];
        export_path = fullfile(out_dir, sprintf('group_%s_runs.csv', name));
    else
        out_dir = char(parameters.path.sim);
        report_path = fullfile(out_dir, sprintf('group_%s_report%s.html', medium, affix));
        group_plots_dir = fullfile(out_dir, 'group_plots');
        plot_suffix = affix;
        export_path = fullfile(out_dir, sprintf('group_%s_runs%s.csv', medium, affix));
    end
    if ~exist(out_dir, 'dir'); mkdir(out_dir); end
    maps = struct([]);
    if isfield(options, 'maps'); maps = options.maps; end

    %% Aggregate per-run data
    [group_table, subjects_meta, missing] = aggregate_runs(runs);

    if isempty(group_table)
        fprintf('generate_group_report: no per-subject CSVs could be read — skipping.\n');
        report_path = '';
        return
    end
    conds = unique(string(group_table.condition), 'stable');
    multi_cond = numel(conds) > 1;
    [outliers, subjects_meta] = find_outliers(group_table, subjects_meta);

    %% Tidy export: one row per run
    try
        writetable(group_table, export_path);
    catch ME
        warning('generate_group_report:export', 'Could not write %s: %s', export_path, ME.message);
        export_path = '';
    end

    %% Render group-level box plots (best-effort)
    if ~exist(group_plots_dir, 'dir'); mkdir(group_plots_dir); end

    plot_path = @(tag) fullfile(group_plots_dir, ...
        sprintf('group_%s_%s%s.png', medium, tag, plot_suffix));

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

    has_thermal = is_layered && ((isfield(parameters.modules, 'run_heating_sims') && ...
        parameters.modules.run_heating_sims) || any(ismember({'maxT_brain', 'maxT_skull'}, ...
        group_table.Properties.VariableNames)));
    thermal_specs = struct('path', {}, 'caption', {}, 'alt', {});
    if has_thermal
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

    cond_plot = '';
    if multi_cond
        cond_plot = render_condition_boxplot(group_table, comparison_metrics(group_table), ...
            plot_path('conditions'));
    end
    has_targeting = any(ismember(targeting_cols(), group_table.Properties.VariableNames));
    targeting_plot = '';
    if has_targeting
        targeting_plot = render_targeting_scatter(group_table, plot_path('targeting'));
    end

    %% Build HTML
    limits = get_risk_limits(is_layered);
    n_runs = numel(subjects_meta);
    report_dir = fileparts(report_path);

    sections = {'header', 'Header'; 'filter', 'Filter'; 'dashboard', 'Safety'; 'other-group', 'Other'};
    if study_mode || multi_cond; sections(end+1, :) = {'design', 'Design'}; end
    if multi_cond;               sections(end+1, :) = {'compare', 'Compare'}; end
    if has_targeting;            sections(end+1, :) = {'targeting', 'Targeting'}; end
    if ~isempty(maps);           sections(end+1, :) = {'maps', 'Maps'}; end
    sections = [sections; {'outliers', 'Outliers'; 'analysis', 'Analysis'; 'roster', 'Roster'; ...
        'group-acoustic', 'Acoustic'}];
    if has_thermal; sections(end+1, :) = {'group-thermal', 'Thermal'}; end
    sections = [sections; {'subjects', 'Runs'; 'config', 'Config'}];

    if study_mode
        title_str = sprintf('PRESTUS Group Report — %s (N=%d)', name, n_runs);
    else
        title_str = sprintf('PRESTUS Group Report — %s%s (N=%d)', medium, affix, n_runs);
    end

    html_parts = {};
    html_parts{end+1} = '<!DOCTYPE html>';
    html_parts{end+1} = '<html lang="en" data-theme="light">';
    html_parts{end+1} = '<head>';
    html_parts{end+1} = '<meta charset="UTF-8">';
    html_parts{end+1} = '<meta name="viewport" content="width=device-width, initial-scale=1.0">';
    html_parts{end+1} = ['<script>try{var t=localStorage.getItem("prestus-theme")||' ...
        '((window.matchMedia&&matchMedia("(prefers-color-scheme:dark)").matches)?"dark":"light");' ...
        'document.documentElement.setAttribute("data-theme",t);}catch(e){}</script>'];
    html_parts{end+1} = sprintf('<title>%s</title>', html_utils.escape(title_str));
    html_parts{end+1} = '<style>';
    html_parts{end+1} = css_styles_base();
    html_parts{end+1} = group_css_extra();
    html_parts{end+1} = '</style>';
    html_parts{end+1} = '</head>';
    html_parts{end+1} = '<body>';
    html_parts{end+1} = '<a id="top"></a>';

    % Table of contents (sticky)
    try
        html_parts{end+1} = build_toc(sections);
    catch ME
        html_parts{end+1} = html_utils.section_error('Table of Contents', ME);
    end

    % Section 1: Header
    try
        html_parts{end+1} = build_group_header(parameters, subjects_meta, missing, medium, affix, ...
            is_layered, options, conds, report_dir, export_path);
    catch ME
        html_parts{end+1} = html_utils.section_error('Header', ME);
    end

    % Group verdict banner + runs exceeding NSR limits
    try
        html_parts{end+1} = build_group_safety_banner(subjects_meta, group_table, limits);
    catch ME
        html_parts{end+1} = html_utils.section_error('Group verdict', ME);
    end

    % Section 2: Filter bar (run checkboxes)
    try
        html_parts{end+1} = build_filter_bar(subjects_meta, conds);
    catch ME
        html_parts{end+1} = html_utils.section_error('Filter', ME);
    end

    % Section 3: Group Safety dashboard (ITRUSST-limited metrics only)
    try
        html_parts{end+1} = build_group_dashboard(subjects_meta, limits, is_layered);
    catch ME
        html_parts{end+1} = html_utils.section_error('Group Safety dashboard', ME);
    end

    % Section 3a: Other (informational) group metrics
    try
        html_parts{end+1} = build_group_other_metrics(subjects_meta, limits);
    catch ME
        html_parts{end+1} = html_utils.section_error('Other exposure metrics', ME);
    end

    % Section 3b: Study design (conditions, per-run settings)
    if study_mode || multi_cond
        try
            html_parts{end+1} = html_utils.collapsible('Study Design', ...
                build_study_design(runs, subjects_meta, missing, conds), true, 'design');
        catch ME
            html_parts{end+1} = html_utils.section_error('Study Design', ME);
        end
    end

    % Section 3c: Condition comparison
    if multi_cond
        try
            html_parts{end+1} = html_utils.collapsible('Condition Comparison', ...
                build_condition_comparison(group_table, conds, cond_plot, limits), true, 'compare');
        catch ME
            html_parts{end+1} = html_utils.section_error('Condition Comparison', ME);
        end
    end

    % Section 3d: Targeting accuracy
    if has_targeting
        try
            html_parts{end+1} = html_utils.collapsible('Targeting Accuracy', ...
                build_targeting_section(group_table, subjects_meta, conds, targeting_plot), true, 'targeting');
        catch ME
            html_parts{end+1} = html_utils.section_error('Targeting Accuracy', ME);
        end
    end

    % Section 3e: Group MNI maps
    if ~isempty(maps)
        try
            html_parts{end+1} = html_utils.collapsible('Group MNI Maps', ...
                build_maps_section(maps, report_dir), true, 'maps');
        catch ME
            html_parts{end+1} = html_utils.section_error('Group MNI Maps', ME);
        end
    end

    % Section 3f: Outliers
    try
        html_parts{end+1} = html_utils.collapsible('Outliers', ...
            build_outlier_section(outliers, subjects_meta), true, 'outliers');
    catch ME
        html_parts{end+1} = html_utils.section_error('Outliers', ME);
    end

    % Section 3g: Data analysis (client-side distribution / correlation / by-group)
    try
        html_parts{end+1} = build_analysis_section(subjects_meta, limits);
    catch ME
        html_parts{end+1} = html_utils.section_error('Data Analysis', ME);
    end

    % Section 4: Run roster table
    try
        html_parts{end+1} = html_utils.collapsible('Subject Roster', ...
            build_subject_roster(group_table, subjects_meta, report_dir, limits, is_layered, multi_cond), ...
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
    if has_thermal
        try
            html_parts{end+1} = html_utils.collapsible('Group Thermal Summary', ...
                build_group_plot_section(thermal_specs, 'No thermal columns available for plotting.'), ...
                true, 'group-thermal');
        catch ME
            html_parts{end+1} = html_utils.section_error('Group Thermal Summary', ME);
        end
    end

    % Section 7: Per-run cards
    try
        html_parts{end+1} = html_utils.collapsible('Per-Subject Details', ...
            build_subject_cards(subjects_meta, report_dir, is_layered), ...
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

    % Embedded JSON payload + filter JS (checkbox-driven dashboard recompute)
    try
        html_parts{end+1} = build_filter_script(subjects_meta, limits);
    catch ME
        html_parts{end+1} = html_utils.section_error('Filter Script', ME);
    end

    % Shared interactive scripts: theme/search/lightbox + sortable roster + analysis charts
    html_parts{end+1} = report_scripts.common();
    html_parts{end+1} = report_scripts.group_runtime();

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

function runs = legacy_run_table(parameters, subject_list)
% One run per subject, all with the medium and affix of parameters.
    n = numel(subject_list);
    subject_id = subject_list(:);
    runs = table((1:n)', subject_id, strings(n, 1), true(n, 1), ...
        repmat(string(get_affix(parameters)), n, 1), ...
        repmat(string(parameters.simulation.medium), n, 1), ...
        repmat(string(parameters.path.sim), n, 1), ...
        compose("sub-%03d", subject_id), ...
        'VariableNames', {'row', 'subject_id', 'condition', 'include', 'affix', 'medium', ...
                          'sim_dir', 'label'});
end

function [group_table, subjects_meta, missing] = aggregate_runs(runs)
% Read each run's CSV (+ targeting sidecar), vertically concatenate, and
% build a per-run struct array with the metrics needed for the dashboard
% and JS filter. Runs without a CSV are returned in `missing`.

    group_table = table();
    subjects_meta = struct('run_id', {}, 'subject_id', {}, 'label', {}, 'condition', {}, ...
                           'affix', {}, 'medium', {}, 'csv_path', {}, 'row_index', {}, ...
                           'sub_dir', {}, 'img_dir', {}, 'report_path', {}, ...
                           'descriptors', {}, 'overrides', {}, 'outlier', {}, 'metrics', {});
    missing = runs([], :);
    text_cols = {'placement_mode', 'target_label', 'transducer_name', 'intensity_quantity', ...
                 'target_mni_source'};

    for k = 1:height(runs)
        id = runs.subject_id(k);
        f = study_run_files(id, runs.medium(k), runs.affix(k), runs.sim_dir(k));

        t = [];
        if isfile(f.csv)
            try
                t = read_metrics_table(f.csv);
            catch
                t = [];
            end
        end
        if isempty(t)
            missing(end+1, :) = runs(k, :); %#ok<AGROW>
            continue
        end

        % Always pick the last row of each per-subject CSV (matches the
        % per-subject report's csv_value() convention).
        row = t(end, :);
        if ~ismember('subject_id', row.Properties.VariableNames)
            row = addvars(row, id, 'Before', 1, 'NewVariableNames', {'subject_id'});
        end
        row.subject_id = id;

        % Targeting sidecar: numeric columns joined, descriptors kept as text
        desc = struct('placement_mode', "", 'target_label', "", 'transducer_name', "", ...
                      'intensity_quantity', "", 'target_mni_source', "");
        if isfile(f.targeting_csv)
            try
                opts = detectImportOptions(f.targeting_csv, 'VariableNamingRule', 'preserve');
                tv = intersect([text_cols, {'condition', 'affix'}], opts.VariableNames);
                opts = setvartype(opts, tv, 'string');
                tt = readtable(f.targeting_csv, opts);
                tt = tt(end, :);
                for c = 1:numel(text_cols)
                    if ismember(text_cols{c}, tt.Properties.VariableNames)
                        v = tt.(text_cols{c});
                        if ~ismissing(v); desc.(text_cols{c}) = v; end
                    end
                end
                for c = tt.Properties.VariableNames
                    v = tt.(c{1});
                    if isnumeric(v) && isscalar(v) && ~ismember(c{1}, row.Properties.VariableNames)
                        row.(c{1}) = double(v);
                    end
                end
            catch ME
                warning('generate_group_report:targetingRead', ...
                    'sub-%03d: could not read %s (%s).', id, f.targeting_csv, ME.message);
            end
        end

        run_id = numel(subjects_meta) + 1;
        row = addvars(row, run_id, string(runs.label(k)), string(runs.condition(k)), ...
            string(runs.affix(k)), desc.placement_mode, desc.target_label, desc.transducer_name, ...
            'After', 'subject_id', 'NewVariableNames', {'run_id', 'run_label', 'condition', ...
            'affix', 'placement_mode', 'target_label', 'transducer_name'});

        % Append to combined table (outer-join semantics — union columns,
        % missing values become NaN/empty). If a per-subject row has an
        % incompatible column type (rare — schemas should match), skip the
        % table append for this run but keep its struct so it still
        % appears in subjects_meta / the JSON payload.
        try
            if isempty(group_table)
                group_table = row;
            else
                group_table = outerjoin_tables(group_table, row);
            end
        catch ME
            warning('generate_group_report:tableAppend', ...
                ['%s: CSV row could not be merged into group table (%s); ' ...
                 'run still included in JSON / dashboard.'], runs.label(k), ME.message);
        end

        % Build per-run metrics struct for JSON payload
        metrics = struct('subject_id', id, 'run_id', run_id);
        cols = row.Properties.VariableNames;
        for c = 1:numel(cols)
            v = row{1, c};
            if isnumeric(v) && isscalar(v)
                metrics.(matlab.lang.makeValidName(cols{c})) = v;
            end
        end
        metrics.subject_id = id;
        metrics.run_id = run_id;
        metrics.id = char(runs.label(k));
        metrics.subject = sprintf('sub-%03d', id);
        metrics.condition = char(runs.condition(k));
        metrics.placement_mode = char(desc.placement_mode);
        metrics.target_label = char(desc.target_label);
        metrics.transducer_name = char(desc.transducer_name);

        ov = {};
        if ismember('overrides', runs.Properties.VariableNames); ov = runs.overrides{k}; end
        subjects_meta(end+1) = struct( ...
            'run_id',     run_id, ...
            'subject_id', id, ...
            'label',      char(runs.label(k)), ...
            'condition',  char(runs.condition(k)), ...
            'affix',      char(runs.affix(k)), ...
            'medium',     char(runs.medium(k)), ...
            'csv_path',   f.csv, ...
            'row_index',  height(t), ...
            'sub_dir',    f.sub_dir, ...
            'img_dir',    f.img_dir, ...
            'report_path', f.report, ...
            'descriptors', desc, ...
            'overrides',  {ov}, ...
            'outlier',    {{}}, ...
            'metrics',    metrics); %#ok<AGROW>
    end
end

function T = outerjoin_tables(A, B)
% Append B's row(s) to A, padding missing columns with NaN/missing.
    colsA = A.Properties.VariableNames;
    colsB = B.Properties.VariableNames;
    all_cols = unique([colsA, colsB], 'stable');

    % Columns present in both tables are not padded, so a type mismatch here
    % (e.g. numeric for one subject, string for another under the same CSV
    % column name — possible if per-subject CSVs in path.sim were written by
    % different PRESTUS versions/schemas) would otherwise be silently
    % coerced or dropped by the concatenation below with no diagnostic.
    shared_cols = intersect(colsA, colsB, 'stable');
    for c = 1:numel(shared_cols)
        col = shared_cols{c};
        if ~strcmp(class(A.(col)), class(B.(col)))
            subj_str = '';
            if ismember('subject_id', colsB)
                subj_str = sprintf(' (sub-%03d)', B.subject_id(1));
            end
            warning('generate_group_report:columnTypeMismatch', ...
                ['Column "%s"%s has type "%s" here vs "%s" in the accumulated group table; ' ...
                 'values may be silently coerced or dropped.'], ...
                col, subj_str, class(B.(col)), class(A.(col)));
        end
    end

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
    cols = {'riseT_brain', 'riseT_skull', 'riseT_skin', ...
            'riseT37_brain', 'riseT37_skull', 'riseT37_skin'};
end
function cols = get_cem43_plot_cols()
    cols = {'CEM43_brain', 'CEM43_skull', 'CEM43_skin', ...
            'CEM43iso_brain', 'CEM43iso_skull', 'CEM43iso_skin'};
end


%% =========================================================================
%  SECTION BUILDERS
%  =========================================================================
function html = build_toc(sections)
    logo = html_utils.report_logo(26);
    if ~isempty(logo)
        brand = ['<a class="brand" href="#top">' logo '</a>'];
    else
        brand = '<a class="brand" href="#top"><span class="wm">PRE<b>STUS</b></span></a>';
    end
    links = '';
    for i = 1:size(sections, 1)
        links = [links sprintf('<a href="#%s">%s</a>', sections{i, 1}, sections{i, 2})]; %#ok<AGROW>
    end
    html = ['<nav class="toc">' brand links ...
        '<span class="toc-spacer"></span>' ...
        '<label class="report-search hide-sm">' ...
        '<svg width="14" height="14" viewBox="0 0 24 24" fill="none"><circle cx="11" cy="11" r="7" stroke="currentColor" stroke-width="2"/><path d="m20 20-3.5-3.5" stroke="currentColor" stroke-width="2" stroke-linecap="round"/></svg>' ...
        '<input id="report-search" type="search" placeholder="Search&hellip;" autocomplete="off"></label>' ...
        '<button class="theme-btn" id="theme-btn" type="button" title="Toggle dark mode" aria-label="Toggle theme">' ...
        '<svg id="theme-sun" width="17" height="17" viewBox="0 0 24 24" fill="none"><circle cx="12" cy="12" r="4.5" stroke="currentColor" stroke-width="2"/><path d="M12 2v2M12 20v2M2 12h2M20 12h2M5 5l1.5 1.5M17.5 17.5 19 19M19 5l-1.5 1.5M6.5 17.5 5 19" stroke="currentColor" stroke-width="2" stroke-linecap="round"/></svg>' ...
        '<svg id="theme-moon" width="17" height="17" viewBox="0 0 24 24" fill="none" style="display:none"><path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8Z" stroke="currentColor" stroke-width="2" stroke-linejoin="round"/></svg>' ...
        '</button>' ...
        '</nav>'];
end

function html = build_group_safety_banner(subjects_meta, group_table, limits) %#ok<INUSL>
% Group verdict: how many runs exceed any finite NSR limit, with a
% clickable list of the offenders (the participant-level flagging).
    names = fieldnames(limits);
    lim_names = {};
    for i = 1:numel(names)
        if ~isinf(limits.(names{i}).limit); lim_names{end+1} = names{i}; end %#ok<AGROW>
    end
    n = numel(subjects_meta);
    exceed_ids = []; exceed_info = {};
    for k = 1:n
        m = subjects_meta(k).metrics;
        over = {};
        for i = 1:numel(lim_names)
            key = matlab.lang.makeValidName(lim_names{i});
            if isfield(m, key) && isnumeric(m.(key)) && isscalar(m.(key)) && ~isnan(m.(key)) ...
                    && m.(key) > limits.(lim_names{i}).limit
                over{end+1} = limits.(lim_names{i}).label; %#ok<AGROW>
            end
        end
        if ~isempty(over)
            exceed_ids(end+1) = subjects_meta(k).run_id; %#ok<AGROW>
            exceed_info{end+1} = sprintf('%s &middot; %s', ...
                html_utils.escape(subjects_meta(k).label), html_utils.escape(strjoin(over, ', '))); %#ok<AGROW>
        end
    end
    nOver = numel(exceed_ids);
    if nOver == 0
        cls = 'pass';
        glyph = '<svg viewBox="0 0 24 24" fill="none"><path d="m5 13 4 4L19 7" stroke="currentColor" stroke-width="2.6" stroke-linecap="round" stroke-linejoin="round"/></svg>';
        title = sprintf('All %d %s within NSR limits', n, run_noun(subjects_meta));
        sub = 'None exceeds any ITRUSST non-significant-risk limit.';
        badges = '';
    else
        cls = 'review';
        glyph = '<svg viewBox="0 0 24 24" fill="none"><circle cx="10.5" cy="10.5" r="6.5" stroke="currentColor" stroke-width="2.2"/><path d="m15.5 15.5 5 5" stroke="currentColor" stroke-width="2.4" stroke-linecap="round"/><path d="M10.5 7.5v3.5" stroke="currentColor" stroke-width="2" stroke-linecap="round"/><circle cx="10.5" cy="13.4" r="1.1" fill="currentColor"/></svg>';
        title = sprintf('%d of %d %s above an ITRUSST NSR reference level: inspect NSR adherence', ...
            nOver, n, run_noun(subjects_meta));
        sub = ['This is not by itself a failed safety check. Some ITRUSST criteria are alternatives ' ...
            '(e.g. temperature rise or thermal dose), and study-specific NSR criteria may apply. ' ...
            'Click an entry to jump to its row and inspect it against the criteria that govern your protocol.'];
        b = '';
        for i = 1:numel(exceed_info)
            b = [b sprintf('<a class="exbadge" href="#row-%d">%s</a>', exceed_ids(i), exceed_info{i})]; %#ok<AGROW>
        end
        badges = ['<div class="exsub">' b '</div>'];
    end
    html = sprintf('<section class="verdict %s" id="verdict">', cls);
    html = [html sprintf('<div class="glyph">%s</div>', glyph)];
    html = [html sprintf('<div><h2>%s</h2><p>%s</p>%s</div>', html_utils.escape(title), sub, badges)];
    html = [html '<div class="tally">'];
    html = [html sprintf('<div class="t g"><div class="n">%d</div><div class="l">Within</div></div>', n - nOver)];
    html = [html sprintf('<div class="t r"><div class="n">%d</div><div class="l">Review</div></div>', nOver)];
    html = [html '</div></section>'];
end

function s = run_noun(subjects_meta)
% 'subjects' when every subject appears once, else 'runs'.
    ids = [subjects_meta.subject_id];
    if numel(unique(ids)) == numel(ids); s = 'subjects'; else, s = 'runs'; end
end

function html = build_analysis_section(subjects_meta, limits)
% Client-side analysis panel: distribution / correlation / by-group,
% rendered from an embedded JSON payload by report_scripts.group_runtime().
    subs = cell(1, numel(subjects_meta));
    for k = 1:numel(subjects_meta)
        subs{k} = subjects_meta(k).metrics;
    end
    if isempty(subs)
        subjects_json = '[]';
    else
        subjects_json = html_utils.script_json(subs);
    end

    % Metric registry (emitted as a JS literal so informational limits can use Infinity)
    specs = { ...
        'MI_tc','MI transcranial','','ac'; ...
        'TIC','Cranial TI','','ac'; ...
        'Isppa','ISPPA','W/cm2','ac'; ...
        'Ipa_target','Ipa at target','W/cm2','ac'; ...
        'real_focal_distance_bowl_mm','Focal distance (peak, apex-ref.)','mm','ac'; ...
        'real_focal_distance_ep_mm','Focal distance (peak, exit-plane-ref.)','mm','ac'; ...
        'maxT_skull','Max temp skull','C','th'; ...
        'riseT_skull','dT skull','C','th'; ...
        'riseT37_skull','dT from 37 skull','C','th'; ...
        'CEM43_skull','CEM43 skull','min','th'; ...
        'maxT_brain','Max temp brain','C','th'; ...
        'CEM43_brain','CEM43 brain','min','th'; ...
        'focus_error_mm','Focus error (vs planned)','mm','tg'; ...
        'focus_error_axial_mm','Focus error, axial','mm','tg'; ...
        'focus_error_lateral_mm','Focus error, lateral','mm','tg'; ...
        'target_error_mni_mm','MNI target error','mm','tg'; ...
        'roi_coverage','ROI coverage','','tg'; ...
        'roi_focal_fraction','Focal volume in ROI','','tg'};
    present = false(size(specs, 1), 1);
    for i = 1:size(specs, 1)
        present(i) = ~isempty(collect_metric(subjects_meta, specs{i, 1}));
    end
    specs = specs(present | ~strcmp(specs(:, 4), 'tg'), :);
    met = {};
    for i = 1:size(specs, 1)
        key = specs{i, 1};
        if isfield(limits, key) && ~isinf(limits.(key).limit)
            lim_s = sprintf('%.6g', limits.(key).limit);
        else
            lim_s = 'Infinity';
        end
        met{end+1} = sprintf('{"k":"%s","label":"%s","unit":"%s","limit":%s,"dom":"%s"}', ...
            key, specs{i,2}, specs{i,3}, lim_s, specs{i,4}); %#ok<AGROW>
    end
    metrics_json = ['[' strjoin(met, ',') ']'];

    opts = '';
    for i = 1:size(specs, 1)
        opts = [opts sprintf('<option value="%s">%s</option>', specs{i,1}, specs{i,2})]; %#ok<AGROW>
    end

    % Group-by choices: always subject; condition / placement / target /
    % transducer when they differ between runs.
    by_opts = '<option value="subject">subject</option>';
    by_keys = {'condition', 'placement_mode', 'target_label', 'transducer_name'};
    by_lbls = {'condition', 'placement mode', 'target', 'transducer'};
    for b = 1:numel(by_keys)
        vals = arrayfun(@(s) string(s.metrics.(by_keys{b})), subjects_meta);
        if numel(unique(vals)) > 1
            sel = '';
            if strcmp(by_keys{b}, 'condition'); sel = ' selected'; end
            by_opts = [by_opts sprintf('<option value="%s"%s>%s</option>', by_keys{b}, sel, by_lbls{b})]; %#ok<AGROW>
        end
    end

    html = '<section class="report-section" id="analysis">';
    html = [html '<h2>Data Analysis</h2>'];
    html = [html '<p class="note">Explore the group distribution of each metric and correlations between metrics; the dashed red line marks the NSR limit. Points are individual runs.</p>'];
    html = [html '<div class="an-tabs">' ...
        '<div class="an-tab on" data-view="dist">Distribution</div>' ...
        '<div class="an-tab" data-view="corr">Correlation</div>' ...
        '<div class="an-tab" data-view="cond">By group</div></div>'];
    html = [html sprintf('<div class="an-ctrls" id="an-c-dist"><label>Metric <select id="an-dist">%s</select></label></div>', opts)];
    html = [html sprintf('<div class="an-ctrls" id="an-c-corr" style="display:none"><label>X <select id="an-x">%s</select></label><label>Y <select id="an-y">%s</select></label></div>', opts, opts)];
    html = [html sprintf('<div class="an-ctrls" id="an-c-cond" style="display:none"><label>Metric <select id="an-cond">%s</select></label><label>Group by <select id="an-by">%s</select></label></div>', opts, by_opts)];
    html = [html '<div class="an-plot"><svg id="an-plot" viewBox="0 0 760 320" preserveAspectRatio="xMidYMid meet"></svg></div>'];
    html = [html '<div class="an-legend" id="an-legend"></div>'];
    html = [html '<script>window.PRESTUS_SUBJECTS=' subjects_json ';window.PRESTUS_METRICS=' metrics_json ';</script>'];
    html = [html '</section>'];
end

function html = build_group_header(parameters, subjects_meta, missing, medium, affix, is_layered, ...
        options, conds, report_dir, export_path)
    n = numel(subjects_meta);
    ids = unique([subjects_meta.subject_id]);
    html = '<section class="report-section" id="header">';
    html = [html '<h1>PRESTUS Group Report</h1>'];
    html = [html '<table class="info-table">'];
    if isfield(options, 'name') && ~isempty(options.name)
        html = [html sprintf('<tr><th>Study</th><td>%s</td></tr>', html_utils.escape(options.name))];
        if isfield(options, 'manifest') && ~isempty(options.manifest)
            html = [html sprintf('<tr><th>Manifest</th><td>%s</td></tr>', html_utils.escape(options.manifest))];
        end
    end
    html = [html sprintf('<tr><th>Runs</th><td><span id="n-active">%d</span> of %d active</td></tr>', n, n)];
    html = [html sprintf('<tr><th>Subjects</th><td>%d: %s</td></tr>', numel(ids), ...
        html_utils.escape(strtrim(sprintf('%d ', ids))))];
    if numel(conds) > 1 || (numel(conds) == 1 && conds ~= "")
        html = [html sprintf('<tr><th>Conditions</th><td>%s</td></tr>', ...
            html_utils.escape(strjoin(cellstr(conds), ', ')))];
    end
    if ~isempty(missing)
        html = [html sprintf('<tr><th>Missing</th><td>%d run(s) without outputs: %s</td></tr>', ...
            height(missing), html_utils.escape(strjoin(cellstr(missing.label), ', ')))];
    end
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
    if ~isempty(export_path) && isfile(export_path)
        html = [html sprintf('<tr><th>Data table</th><td><a href="%s">%s</a> (one row per run)</td></tr>', ...
            html_utils.escape(rel_link(report_dir, export_path)), html_utils.escape(rel_link(report_dir, export_path)))];
    end
    html = [html sprintf('<tr><th>Generated</th><td>%s</td></tr>', datestr(now, 'yyyy-mm-dd HH:MM:SS'))];
    html = [html '</table>'];
    html = [html '<p class="note">Use the filter bar below to toggle runs in/out. ' ...
        'Dashboard statistics and per-run cards update live; group box plots reflect the full set.</p>'];
    html = [html '</section>'];
end

function html = build_filter_bar(subjects_meta, conds)
    html = '<section class="report-section filter-bar" id="filter">';
    html = [html '<h2>Filter</h2>'];
    html = [html '<div class="filter-controls">'];
    html = [html '<button type="button" class="filter-btn" onclick="setAllSubjects(true)">All</button>'];
    html = [html '<button type="button" class="filter-btn" onclick="setAllSubjects(false)">None</button>'];
    if numel(conds) > 1
        for c = 1:numel(conds)
            lbl = char(conds(c)); if isempty(lbl); lbl = '(none)'; end
            html = [html sprintf('<button type="button" class="filter-btn" onclick="onlyCondition(''%s'')">%s</button>', ...
                html_utils.escape(strrep(char(conds(c)), '''', '')), html_utils.escape(lbl))]; %#ok<AGROW>
        end
    end
    html = [html '<span class="filter-hint">Toggle runs to update dashboard and per-run cards.</span>'];
    html = [html '</div>'];
    html = [html '<div class="filter-checkboxes">'];
    for k = 1:numel(subjects_meta)
        s = subjects_meta(k);
        flag = '';
        if ~isempty(s.outlier); flag = ' <span class="roster-flag">outlier</span>'; end
        html = [html sprintf( ...
            ['<label class="filter-chk"><input type="checkbox" class="subj-toggle" ' ...
             'value="%d" data-cond="%s" checked> %s%s</label>'], s.run_id, ...
             html_utils.escape(strrep(s.condition, '''', '')), html_utils.escape(s.label), flag)]; %#ok<AGROW>
    end
    html = [html '</div>'];
    html = [html '</section>'];
end

function html = build_group_dashboard(subjects_meta, limits, is_layered)
% Group ITRUSST safety section: only finite-limit metrics, split into a
% Mechanical (MI) and a Thermal sub-grid. Cards carry data-* attributes so the
% JS filter recomputes on toggle. Informational metrics live in
% build_group_other_metrics.
    html = '<section class="report-section" id="dashboard">';
    html = [html '<h2>Group Safety &mdash; ITRUSST exposure</h2>'];

    if ~is_layered
        html = [html '<div class="medium-banner medium-water">' ...
            '<strong>Water / Free-field</strong> &mdash; tissue-specific NSR limits are not applicable. ' ...
            'See <a href="#other-group">Other exposure metrics</a> below.</div>'];
    end

    metric_names = fieldnames(limits);
    mech_keys = {}; therm_keys = {};
    for i = 1:numel(metric_names)
        name = metric_names{i};
        switch metric_group(name, limits.(name))
            case 'mechanical'
                % Group Safety board reports one top-line MI tile per medium
                % (MI free-water or MI_tc transcranial) —
                % per-tissue MI_brain/MI_skull/MI_skin are intentionally
                % excluded here (still available in the CSV / box plots).
                if strcmp(name, 'MI_tc') || strcmp(name, 'MI')
                    mech_keys{end+1} = name; %#ok<AGROW>
                end
            case 'thermal'
                therm_keys{end+1} = name; %#ok<AGROW>
        end
    end
    html = [html render_group_grid(mech_keys,  'mech', 'Mechanical', 'mechanical index (MI)', limits, subjects_meta)];
    html = [html render_group_grid(therm_keys, 'th',   'Thermal',    'temperature &amp; dose', limits, subjects_meta)];
    if isempty(mech_keys) && isempty(therm_keys)
        html = [html '<p class="placeholder">No tissue-specific ITRUSST-limited metrics for this medium &mdash; see Other exposure metrics below.</p>'];
    end
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
function html = build_group_other_metrics(subjects_meta, limits)
% Informational group metrics (no ITRUSST limit): TIC, ISPPA.
    metric_names = fieldnames(limits);
    other_keys = {};
    for i = 1:numel(metric_names)
        if strcmp(metric_group(metric_names{i}, limits.(metric_names{i})), 'other')
            other_keys{end+1} = metric_names{i}; %#ok<AGROW>
        end
    end
    html = '<section class="report-section" id="other-group">';
    html = [html '<h2>Other exposure metrics</h2>'];
    if isempty(other_keys)
        html = [html '<p class="placeholder">No additional metrics.</p></section>'];
        return
    end
    html = [html render_group_grid(other_keys, 'other', '', '', limits, subjects_meta)];
    html = [html '<p class="safety-footnote">These intensity and output-display indices have ' ...
        'no ITRUSST non-significant-risk limit and are shown for reference only.</p>'];
    html = [html '</section>'];
end
function g = metric_group(name, info)
% Classify a metric: 'other' (informational / Inf limit), 'mechanical' (MI*),
% else 'thermal'. Mirrors the per-subject report split.
    if isinf(info.limit)
        g = 'other';
    elseif startsWith(name, 'MI')
        g = 'mechanical';
    else
        g = 'thermal';
    end
end
function html = render_group_grid(keys, dom_class, title, subtitle, limits, subjects_meta)
% A labelled sub-grid of group cards (skips the label when title is empty).
    html = '';
    if isempty(keys); return; end
    if ~isempty(title)
        html = sprintf(['<div class="dash-group-label %s"><span class="dot"></span>' ...
            '<span class="t">%s</span><span class="s">%s</span></div>'], dom_class, title, subtitle);
    end
    html = [html '<div class="safety-grid">'];
    for i = 1:numel(keys)
        html = [html render_group_card(keys{i}, limits.(keys{i}), subjects_meta)]; %#ok<AGROW>
    end
    html = [html '</div>'];
end
function html = render_group_card(name, info, subjects_meta)
% One group dashboard tile: mean / SD / range with a min-max bar, mean marker
% and NSR limit line. Carries data-* attributes for the JS filter recompute.
    js_key = matlab.lang.makeValidName(name);
    vals = collect_metric(subjects_meta, js_key);
    s = stat_summary(vals);

    % Colour driven by max (worst case across subjects).
    color = risk_color(s.max, info.limit);

    % Axis baseline: 37 for absolute temperatures, 0 elsewhere
    if strncmp(name, 'maxT', 4) || strncmp(name, 'endT', 4)
        scale_min = 37;
    else
        scale_min = 0;
    end
    if isinf(info.limit)
        scale_max = max(scale_min + 1, s.max);
    else
        scale_max = max([scale_min + 1, s.max, info.limit]);
    end
    if isnan(scale_max), scale_max = scale_min + 1; end
    scale_range = scale_max - scale_min;

    html = sprintf( ...
        '<div class="safety-card safety-%s" data-card-metric="%s" data-scale-min="%.6g" data-scale-max="%.6g" data-limit="%s">', ...
        color, html_utils.escape(js_key), scale_min, scale_max, ...
        iff(isinf(info.limit), 'null', sprintf('%.6g', info.limit)));

    html = [html sprintf('<div class="safety-label">%s</div>', html_utils.escape(info.label))];

    if isnan(s.mean)
        html = [html '<div class="safety-value" data-stat-mean>N/A</div>'];
    else
        html = [html sprintf('<div class="safety-value" data-stat-mean>%.3g</div>', s.mean)];
    end

    if isnan(s.sd) || s.n < 1
        html = [html '<div class="safety-sd" data-stat-sd>&plusmn; &mdash;</div>'];
    else
        html = [html sprintf('<div class="safety-sd" data-stat-sd>&plusmn; %.3g %s</div>', s.sd, html_utils.escape(info.unit))];
    end

    if s.n == 0
        range_str = '(N=0)';
    elseif s.n == 1
        range_str = sprintf('%.3g (N=1)', s.min);
    else
        range_str = sprintf('%.3g &ndash; %.3g (N=%d)', s.min, s.max, s.n);
    end
    html = [html sprintf('<div class="safety-range" data-stat-range>%s</div>', range_str)];

    if isinf(info.limit)
        html = [html sprintf('<div class="risk-limit">%s (informational)</div>', html_utils.escape(info.unit))];
    else
        html = [html sprintf('<div class="risk-limit">NSR limit: %.3g %s</div>', info.limit, html_utils.escape(info.unit))];
    end

    html = [html sprintf('<div class="safety-bar-range">%.3g &ndash; %.3g %s</div>', ...
        scale_min, scale_max, html_utils.escape(info.unit))];

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

    html = [html '</div>'];
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
function html = build_subject_roster(group_table, subjects_meta, report_dir, limits, is_layered, multi_cond)
% One row per run with key metrics + a link to the per-subject report
% (if present). Each row carries data-subject-id (the run id) so the JS
% filter hides/shows it.
    if isempty(group_table)
        html = '<p class="placeholder">No subject data.</p>';
        return
    end

    % Headline columns: ITRUSST-relevant metrics first (MI, then thermal NSR
    % metrics), then the informational intensity / pressure / distance and
    % targeting columns.
    headline = {'MI_tc'};
    if is_layered && any(ismember({'maxT_skull'}, group_table.Properties.VariableNames))
        headline = [headline, {'maxT_skull', 'riseT_skull', 'riseT37_skull', 'CEM43_skull'}];
    end
    headline = [headline, {'Isppa', 'Ipa_target', 'real_focal_distance_bowl_mm', 'real_focal_distance_ep_mm', ...
        'TIC', 'focus_error_mm', 'target_error_mni_mm'}];
    avail = intersect(headline, group_table.Properties.VariableNames, 'stable');
    keep = true(size(avail));
    for c = 1:numel(avail)
        v = group_table.(avail{c});
        keep(c) = ~(isnumeric(v) && all(isnan(v)));
    end
    avail = avail(keep);
    lim_names = fieldnames(limits);

    % Sortable header; each numeric <td> carries data-col + data-value so the
    % shared roster script can sort, and rows exceeding any NSR limit are flagged.
    html = ['<p class="note">Click a column header to sort. Rows exceeding any NSR limit are flagged ' ...
        '<em>over</em>; statistical outliers within their condition <em>outlier</em> (see Outliers).</p>'];
    html = [html '<div class="table-wrapper"><table class="data-table roster" id="roster-table"><thead><tr>'];
    html = [html '<th>Run</th>'];
    if multi_cond; html = [html '<th>Condition</th>']; end
    for c = 1:numel(avail)
        html = [html sprintf('<th data-col="%s">%s</th>', avail{c}, html_utils.escape(avail{c}))]; %#ok<AGROW>
    end
    html = [html '<th>Report</th>'];
    html = [html '</tr></thead><tbody>'];

    for k = 1:numel(subjects_meta)
        rid = subjects_meta(k).run_id;
        row_idx = find(group_table.run_id == rid, 1, 'last');
        if isempty(row_idx)
            continue
        end

        % Flag: does this run exceed any finite-limit metric in its row?
        is_exceed = false;
        for i = 1:numel(lim_names)
            ln = lim_names{i};
            if isinf(limits.(ln).limit); continue; end
            if ismember(ln, group_table.Properties.VariableNames)
                vv = group_table{row_idx, ln};
                if iscell(vv); vv = vv{1}; end
                if isnumeric(vv) && isscalar(vv) && ~isnan(vv) && vv > limits.(ln).limit
                    is_exceed = true; break
                end
            end
        end

        flag = '';
        if is_exceed
            html = [html sprintf('<tr data-subject-id="%d" class="exceed" id="row-%d">', rid, rid)]; %#ok<AGROW>
            flag = '<span class="roster-flag">over</span>';
        else
            html = [html sprintf('<tr data-subject-id="%d" id="row-%d">', rid, rid)]; %#ok<AGROW>
        end
        if ~isempty(subjects_meta(k).outlier)
            flag = [flag '<span class="roster-flag" title="' ...
                html_utils.escape(strjoin(subjects_meta(k).outlier, ', ')) '">outlier</span>']; %#ok<AGROW>
        end
        html = [html sprintf('<td>%s%s</td>', flag, html_utils.escape(sprintf('sub-%03d', subjects_meta(k).subject_id)))]; %#ok<AGROW>
        if multi_cond
            html = [html sprintf('<td>%s</td>', html_utils.escape(subjects_meta(k).condition))]; %#ok<AGROW>
        end

        for c = 1:numel(avail)
            col = avail{c};
            val = group_table{row_idx, col};
            if iscell(val); val = val{1}; end
            val_str = html_utils.format_cell(val);

            dval = ''; cell_class = '';
            if isnumeric(val) && isscalar(val) && ~isnan(val)
                dval = sprintf(' data-value="%.6g"', val);
                if isfield(limits, col)
                    color = risk_color(val, limits.(col).limit);
                    if strcmp(color, 'red'); cell_class = ' over';
                    elseif strcmp(color, 'amber'); cell_class = ' near'; end
                end
            end
            html = [html sprintf('<td data-col="%s"%s class="num%s">%s</td>', col, dval, cell_class, val_str)]; %#ok<AGROW>
        end

        rpt = subjects_meta(k).report_path;
        if isfile(rpt)
            html = [html sprintf('<td><a href="%s" target="_blank">open</a></td>', ...
                html_utils.escape(rel_link(report_dir, rpt)))]; %#ok<AGROW>
        else
            html = [html '<td><span class="placeholder">&mdash;</span></td>']; %#ok<AGROW>
        end
        html = [html '</tr>']; %#ok<AGROW>
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
function html = build_subject_cards(subjects_meta, report_dir, is_layered)
% One collapsible card per run, embedding key per-subject images. Each
% card has data-subject-id (the run id) so the JS filter can hide/show them.
    html = '';
    for k = 1:numel(subjects_meta)
        s = subjects_meta(k);
        id = s.subject_id; medium = s.medium; affix = s.affix;
        img_dir = s.img_dir;

        body = '';
        d = s.descriptors;
        facts = {};
        if d.placement_mode ~= "";  facts{end+1} = sprintf('placement %s', d.placement_mode); end %#ok<AGROW>
        if d.target_label ~= "";    facts{end+1} = sprintf('target %s', d.target_label); end %#ok<AGROW>
        if d.transducer_name ~= ""; facts{end+1} = sprintf('transducer %s', d.transducer_name); end %#ok<AGROW>
        if ~isempty(s.affix);       facts{end+1} = sprintf('affix %s', s.affix); end %#ok<AGROW>
        if ~isempty(facts)
            body = sprintf('<p class="note">%s</p>', html_utils.escape(strjoin(facts, ' · ')));
        end
        body = [body '<div class="image-grid">']; %#ok<AGROW>
        body = [body try_embed(fullfile(img_dir, ...
            sprintf('sub-%03d_positioning%s.png', id, affix)), 'Positioning', 'Transducer positioning')]; %#ok<AGROW>
        body = [body try_embed(fullfile(img_dir, ...
            sprintf('sub-%03d_%s_intensity%s.png', id, medium, affix)), 'Intensity', 'Intensity overlay (segmentation)')]; %#ok<AGROW>
        body = [body try_embed(fullfile(img_dir, ...
            sprintf('sub-%03d_%s_intensity_t1%s.png', id, medium, affix)), 'Intensity T1', 'Intensity overlay (T1)')]; %#ok<AGROW>
        if is_layered
            body = [body try_embed(fullfile(img_dir, ...
                sprintf('sub-%03d_%s_maxT%s.png', id, medium, affix)), 'maxT', 'Max temperature overlay')]; %#ok<AGROW>
        end
        body = [body '</div>']; %#ok<AGROW>

        % Per-subject report link
        rpt = s.report_path;
        if isfile(rpt)
            body = [body sprintf('<p class="note">Per-subject report: <a href="%s" target="_blank">%s</a></p>', ...
                html_utils.escape(rel_link(report_dir, rpt)), ...
                html_utils.escape(rel_link(report_dir, rpt)))]; %#ok<AGROW>
        end

        html = [html sprintf('<details class="subject-card" data-subject-id="%d" open>', s.run_id)]; %#ok<AGROW>
        html = [html sprintf('<summary><strong>%s</strong></summary>', html_utils.escape(s.label))]; %#ok<AGROW>
        html = [html '<div class="section-content">']; %#ok<AGROW>
        html = [html body]; %#ok<AGROW>
        html = [html '</div></details>']; %#ok<AGROW>
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
    json = html_utils.script_json(payload);

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
        'function onlyCondition(c) {' newline ...
        '  document.querySelectorAll(".subj-toggle").forEach(cb => cb.checked = (cb.dataset.cond === c));' newline ...
        '  recompute();' newline ...
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
        '  const active = SUBJECTS.filter(s => enabledSet.has(s.run_id));' newline ...
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
    % Consolidated into css_styles_base(); kept unreachable for reference.
    css = '';
    return; %#ok<UNRCH>
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
function rel = rel_link(from_dir, target)
% Relative link from the report folder to target, so links keep working
% when the whole simulation folder is moved. Falls back to the absolute
% path when the two share no common root.
    a = strsplit(strrep(char(from_dir), '\', '/'), '/');
    b = strsplit(strrep(char(target), '\', '/'), '/');
    a = a(~cellfun(@isempty, a)); b = b(~cellfun(@isempty, b));
    k = 0;
    while k < min(numel(a), numel(b)) && strcmp(a{k+1}, b{k+1})
        k = k + 1;
    end
    if k == 0
        rel = strrep(char(target), '\', '/');
        return
    end
    rel = strjoin([repmat({'..'}, 1, numel(a) - k), b(k+1:end)], '/');
end


%% =========================================================================
%  STUDY / CONDITION / TARGETING SECTIONS
%  =========================================================================

function cols = comparison_metrics(group_table)
% Key metrics compared between conditions (those present with data).
    cols = {'Isppa_brain', 'Ipa_target', 'MI_tc', 'real_focal_distance_ep_mm', 'focal_vol_3dB_mm3', ...
            'maxT_brain', 'maxT_skull', 'riseT_skull', 'CEM43_skull', ...
            'focus_error_mm', 'target_error_mni_mm', 'roi_coverage'};
    cols = cols(ismember(cols, group_table.Properties.VariableNames));
    keep = cellfun(@(c) isnumeric(group_table.(c)) && any(isfinite(group_table.(c))), cols);
    cols = cols(keep);
end

function cols = targeting_cols()
    cols = {'focus_error_mm', 'focus_error_axial_mm', 'focus_error_lateral_mm', 'peak_error_mm', ...
            'target_error_mni_mm', 'roi_focal_fraction', 'roi_coverage', 'roi_mean_intensity'};
end

function [outliers, subjects_meta] = find_outliers(group_table, subjects_meta)
% Robust outliers within each condition: modified z-score
% 0.6745 (x - median) / MAD, |z| > 3.5 (Iglewicz & Hoaglin, 1993), for
% conditions with at least 5 runs.
    outliers = table('Size', [0 7], ...
        'VariableTypes', {'double', 'string', 'string', 'string', 'double', 'double', 'double'}, ...
        'VariableNames', {'run_id', 'label', 'condition', 'metric', 'value', 'median', 'z'});
    metrics = comparison_metrics(group_table);
    conds = unique(string(group_table.condition), 'stable');
    for c = 1:numel(conds)
        in = string(group_table.condition) == conds(c);
        if nnz(in) < 5; continue; end
        for m = 1:numel(metrics)
            x = group_table.(metrics{m})(in);
            ok = isfinite(x);
            if nnz(ok) < 5; continue; end
            med = median(x(ok));
            mad = median(abs(x(ok) - med));
            if mad <= 0; continue; end
            z = 0.6745 * (x - med) / mad;
            rid = group_table.run_id(in);
            lbl = string(group_table.run_label(in));
            for i = find(abs(z) > 3.5 & ok)'
                outliers(end+1, :) = {rid(i), lbl(i), conds(c), string(metrics{m}), x(i), med, z(i)}; %#ok<AGROW>
                k = find([subjects_meta.run_id] == rid(i), 1);
                if ~isempty(k)
                    subjects_meta(k).outlier{end+1} = metrics{m};
                end
            end
        end
    end
end

function html = build_study_design(runs, subjects_meta, missing, conds)
% Conditions, the settings that vary, and every run's own settings.
    html = '<h3>Conditions</h3><div class="table-wrapper"><table class="data-table"><thead><tr>';
    html = [html '<th>Condition</th><th>Runs with outputs</th><th>Subjects</th><th>Missing</th></tr></thead><tbody>'];
    all_conds = unique([string(runs.condition); conds], 'stable');
    for c = 1:numel(all_conds)
        sm = subjects_meta(arrayfun(@(s) string(s.condition) == all_conds(c), subjects_meta));
        nm = 0;
        if ~isempty(missing); nm = nnz(string(missing.condition) == all_conds(c)); end
        lbl = char(all_conds(c)); if isempty(lbl); lbl = '(none)'; end
        html = [html sprintf('<tr><td>%s</td><td class="num">%d</td><td>%s</td><td class="num">%d</td></tr>', ...
            html_utils.escape(lbl), numel(sm), ...
            html_utils.escape(strtrim(sprintf('%d ', unique([sm.subject_id])))), nm)]; %#ok<AGROW>
    end
    html = [html '</tbody></table></div>'];

    % Manifest keys and how many distinct values each takes
    if ismember('overrides', runs.Properties.VariableNames)
        keys = {};
        for r = 1:height(runs)
            ov = runs.overrides{r};
            if ~isempty(ov); keys = [keys; ov(:, 1)]; end %#ok<AGROW>
        end
        keys = unique(keys, 'stable');
        if ~isempty(keys)
            parts = {};
            for k = 1:numel(keys)
                vals = strings(0, 1);
                for r = 1:height(runs)
                    ov = runs.overrides{r};
                    if isempty(ov); continue; end
                    hit = find(strcmp(ov(:, 1), keys{k}), 1, 'last');
                    if ~isempty(hit); vals(end+1) = string(value_text(ov{hit, 2})); end %#ok<AGROW>
                end
                parts{end+1} = sprintf('<code>%s</code> (%d value%s)', html_utils.escape(keys{k}), ...
                    numel(unique(vals)), iff(numel(unique(vals)) == 1, '', 's')); %#ok<AGROW>
            end
            html = [html '<p class="note"><strong>Per-run settings from the manifest:</strong> ' ...
                strjoin(parts, ', ') '</p>'];
        end
    end

    html = [html '<h3>Runs</h3><div class="table-wrapper"><table class="data-table"><thead><tr>'];
    html = [html '<th>Run</th><th>Affix</th><th>Placement</th><th>Target</th><th>Transducer</th>' ...
        '<th>Planned depth (mm)</th><th>Target Isppa (W/cm²)</th><th>Manifest settings</th><th>Outputs</th></tr></thead><tbody>'];
    for k = 1:numel(subjects_meta)
        s = subjects_meta(k); d = s.descriptors; m = s.metrics;
        html = [html sprintf('<tr data-subject-id="%d"><td>%s</td><td><code>%s</code></td><td>%s</td><td>%s</td><td>%s</td><td class="num">%s</td><td class="num">%s</td><td>%s</td><td>%s</td></tr>', ...
            s.run_id, html_utils.escape(s.label), html_utils.escape(s.affix), ...
            html_utils.escape(char(d.placement_mode)), html_utils.escape(char(d.target_label)), ...
            html_utils.escape(char(d.transducer_name)), num_field(m, 'focal_distance_ep_planned_mm'), ...
            num_field(m, 'target_isppa_wcm2'), overrides_text(s.overrides), ...
            iff(isfield(m, 'focus_error_mm'), 'results + targeting', 'results'))]; %#ok<AGROW>
    end
    for k = 1:height(missing)
        ov = {};
        if ismember('overrides', missing.Properties.VariableNames); ov = missing.overrides{k}; end
        html = [html sprintf(['<tr class="exceed"><td>%s</td><td><code>%s</code></td><td></td><td></td><td></td>' ...
            '<td></td><td></td><td>%s</td><td><span class="roster-flag">missing</span></td></tr>'], ...
            html_utils.escape(char(missing.label(k))), html_utils.escape(char(missing.affix(k))), ...
            overrides_text(ov))]; %#ok<AGROW>
    end
    html = [html '</tbody></table></div>'];
    html = [html '<p class="note">Placement, target, transducer, planned depth and target intensity come from ' ...
        'each run''s targeting sidecar (written by the pipeline); runs simulated before it existed show blanks.</p>'];
end

function html = build_condition_comparison(group_table, conds, plot_file, limits)
    metrics = comparison_metrics(group_table);
    html = '<div class="table-wrapper"><table class="data-table"><thead><tr><th>Metric</th>';
    for c = 1:numel(conds)
        lbl = char(conds(c)); if isempty(lbl); lbl = '(none)'; end
        html = [html sprintf('<th>%s</th>', html_utils.escape(lbl))]; %#ok<AGROW>
    end
    html = [html '</tr></thead><tbody>'];
    for m = 1:numel(metrics)
        label = metrics{m};
        if isfield(limits, label); label = sprintf('%s (%s)', limits.(label).label, limits.(label).unit); end
        html = [html sprintf('<tr><td>%s</td>', html_utils.escape(label))]; %#ok<AGROW>
        for c = 1:numel(conds)
            v = group_table.(metrics{m})(string(group_table.condition) == conds(c));
            html = [html sprintf('<td class="num">%s</td>', msd(v))]; %#ok<AGROW>
        end
        html = [html '</tr>']; %#ok<AGROW>
    end
    html = [html '</tbody></table></div>'];
    html = [html '<p class="note">Mean &plusmn; SD (n runs) per condition. Descriptive only: runs of one ' ...
        'subject under several conditions are not independent.</p>'];
    if ~isempty(plot_file) && isfile(plot_file)
        html = [html '<div class="image-grid--wide">' html_utils.embed_image(plot_file, ...
            'Metrics by condition', 'Key metrics by condition (box = IQR, points = runs).') '</div>'];
    end
end

function html = build_targeting_section(group_table, subjects_meta, conds, plot_file)
    cols = targeting_cols();
    cols = cols(ismember(cols, group_table.Properties.VariableNames));
    names = struct('focus_error_mm', 'Focus error (centroid vs planned focus)', ...
        'focus_error_axial_mm', 'Axial component (+ = deeper)', ...
        'focus_error_lateral_mm', 'Lateral component', ...
        'peak_error_mm', 'Peak vs planned focus', ...
        'target_error_mni_mm', 'MNI centroid vs MNI target', ...
        'roi_focal_fraction', 'Focal volume inside ROI (fraction)', ...
        'roi_coverage', 'ROI covered by focal volume (fraction)', ...
        'roi_mean_intensity', 'Mean intensity in ROI (W/cm²)');

    html = '<div class="table-wrapper"><table class="data-table"><thead><tr><th>Measure</th>';
    for c = 1:numel(conds)
        lbl = char(conds(c)); if isempty(lbl); lbl = 'all runs'; end
        html = [html sprintf('<th>%s</th>', html_utils.escape(lbl))]; %#ok<AGROW>
    end
    html = [html '</tr></thead><tbody>'];
    for m = 1:numel(cols)
        v_all = group_table.(cols{m});
        if ~any(isfinite(v_all)); continue; end
        html = [html sprintf('<tr><td>%s</td>', html_utils.escape(names.(cols{m})))]; %#ok<AGROW>
        for c = 1:numel(conds)
            v = v_all(string(group_table.condition) == conds(c));
            html = [html sprintf('<td class="num">%s</td>', msd(v))]; %#ok<AGROW>
        end
        html = [html '</tr>']; %#ok<AGROW>
    end
    html = [html '</tbody></table></div>'];
    html = [html '<p class="note">Mean &plusmn; SD (n runs); distances in mm. The simulated focus is the ' ...
        'centroid of the &minus;3 dB focal volume (brain-masked); the planned focus is the focus position ' ...
        'resolved by placement; axial/lateral split the error along/across the planned beam axis of ' ...
        'transducer 1. MNI targets come from <code>analysis.target_mni_mm</code> or the placement mode; ' ...
        'ROI measures need <code>analysis.target_roi_mni</code>.</p>'];

    % Per-run table
    per = {'focus_error_mm', 'focus_error_axial_mm', 'focus_error_lateral_mm', 'target_error_mni_mm', ...
           'roi_focal_fraction', 'roi_coverage'};
    per = per(ismember(per, group_table.Properties.VariableNames));
    html = [html '<div class="table-wrapper"><table class="data-table roster"><thead><tr><th>Run</th><th>Target</th>' ...
        '<th>MNI centroid (mm)</th><th>MNI target (mm)</th>'];
    for p = 1:numel(per)
        html = [html sprintf('<th>%s</th>', html_utils.escape(per{p}))]; %#ok<AGROW>
    end
    html = [html '</tr></thead><tbody>'];
    for k = 1:numel(subjects_meta)
        s = subjects_meta(k);
        i = find(group_table.run_id == s.run_id, 1);
        if isempty(i); continue; end
        html = [html sprintf('<tr data-subject-id="%d"><td>%s</td><td>%s</td><td>%s</td><td>%s</td>', ...
            s.run_id, html_utils.escape(s.label), html_utils.escape(char(s.descriptors.target_label)), ...
            xyz_text(group_table, i, 'focal_centroid_mni'), xyz_text(group_table, i, 'target_mni'))]; %#ok<AGROW>
        for p = 1:numel(per)
            html = [html sprintf('<td class="num">%s</td>', html_utils.format_cell(group_table.(per{p})(i)))]; %#ok<AGROW>
        end
        html = [html '</tr>']; %#ok<AGROW>
    end
    html = [html '</tbody></table></div>'];

    if ~isempty(plot_file) && isfile(plot_file)
        html = [html '<div class="image-grid--wide">' html_utils.embed_image(plot_file, 'Targeting', ...
            'MNI focal centroids per run (colour = condition, + = MNI targets) and the axial/lateral focus error.') '</div>'];
    end
end

function html = build_maps_section(maps, report_dir)
    html = ['<p class="note">Voxelwise group maps in MNI space from the per-run MNI NIfTIs ' ...
        '(<code>io.save_MNI</code>). The overlap map gives, per voxel, the fraction of runs whose ' ...
        '&minus;3 dB focal volume covers it.</p>'];
    for k = 1:numel(maps)
        m = maps(k);
        lbl = m.condition; if isempty(lbl); lbl = 'all runs'; end
        html = [html sprintf('<h3>%s <span class="badge badge-gray">N = %d</span></h3>', ...
            html_utils.escape(lbl), m.n)]; %#ok<AGROW>
        if m.n_missing > 0
            html = [html sprintf('<p class="note">%d run(s) without MNI focal masks are not included.</p>', m.n_missing)]; %#ok<AGROW>
        end
        if ~isempty(m.figure) && isfile(m.figure)
            html = [html '<div class="image-grid--wide">' html_utils.embed_image(m.figure, ...
                sprintf('Focal overlap %s', lbl), sprintf('Focal-volume overlap, %s.', lbl)) '</div>']; %#ok<AGROW>
        end
        fn = fieldnames(m.files);
        links = cell(1, numel(fn));
        for f = 1:numel(fn)
            p = m.files.(fn{f});
            links{f} = sprintf('<a href="%s">%s</a>', html_utils.escape(rel_link(report_dir, p)), ...
                html_utils.escape(regexprep(p, '^.*[\\/]', '')));
        end
        html = [html '<p class="note">NIfTI: ' strjoin(links, ' &middot; ') '</p>']; %#ok<AGROW>
    end
end

function html = build_outlier_section(outliers, subjects_meta)
    if isempty(outliers)
        conds = unique(string({subjects_meta.condition}));
        n_min = min(arrayfun(@(c) nnz(string({subjects_meta.condition}) == c), conds));
        if n_min < 5
            html = '<p class="placeholder">Outlier screening needs at least 5 runs per condition.</p>';
        else
            html = '<p class="placeholder">No run deviates strongly from its condition (|robust z| &le; 3.5).</p>';
        end
        return
    end
    html = ['<p class="note">Runs whose value lies far from their condition''s median: modified z-score ' ...
        '0.6745&middot;(x &minus; median)/MAD with |z| &gt; 3.5 (Iglewicz &amp; Hoaglin, 1993). Check the ' ...
        'run''s placement and inputs before interpreting it.</p>'];
    html = [html '<div class="table-wrapper"><table class="data-table"><thead><tr><th>Run</th><th>Condition</th>' ...
        '<th>Metric</th><th>Value</th><th>Condition median</th><th>z</th></tr></thead><tbody>'];
    for i = 1:height(outliers)
        html = [html sprintf(['<tr data-subject-id="%d"><td><a href="#row-%d">%s</a></td><td>%s</td>' ...
            '<td>%s</td><td class="num">%.4g</td><td class="num">%.4g</td><td class="num">%+.1f</td></tr>'], ...
            outliers.run_id(i), outliers.run_id(i), html_utils.escape(char(outliers.label(i))), ...
            html_utils.escape(char(outliers.condition(i))), html_utils.escape(char(outliers.metric(i))), ...
            outliers.value(i), outliers.median(i), outliers.z(i))]; %#ok<AGROW>
    end
    html = [html '</tbody></table></div>'];
end

function plot_path = render_condition_boxplot(group_table, metrics, out_path)
% One panel per metric, one box per condition.
    plot_path = '';
    if isempty(metrics); return; end
    try
        conds = unique(string(group_table.condition), 'stable');
        lbls = cellstr(conds); lbls(conds == "") = {'(none)'};
        nm = numel(metrics);
        ncol = min(4, nm); nrow = ceil(nm / ncol);
        fig = figure('Visible', 'off', 'Color', 'w', ...
            'Position', [100 100 max(700, 260 * ncol) 260 * nrow]);
        tl = tiledlayout(fig, nrow, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');
        for m = 1:nm
            ax = nexttile(tl);
            v = group_table.(metrics{m});
            g = zeros(size(v));
            for c = 1:numel(conds); g(string(group_table.condition) == conds(c)) = c; end
            ok = isfinite(v);
            if exist('boxchart', 'file') == 2 || exist('boxchart', 'builtin') == 5
                boxchart(ax, g(ok), v(ok));
            else
                boxplot(ax, v(ok), g(ok));
            end
            hold(ax, 'on');
            jitter = (mod((1:nnz(ok))', 7) - 3) * 0.05;
            plot(ax, g(ok) + jitter, v(ok), '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 10);
            set(ax, 'XTick', 1:numel(conds), 'XTickLabel', lbls, 'XLim', [0.4 numel(conds) + 0.6], ...
                'TickLabelInterpreter', 'none', 'FontSize', 8);
            title(ax, metrics{m}, 'Interpreter', 'none', 'FontSize', 9, 'FontWeight', 'normal');
            grid(ax, 'on');
        end
        export_figure(fig, out_path);
        plot_path = out_path;
    catch ME
        warning('generate_group_report:plotFailed', 'Condition plot failed: %s', ME.message);
    end
end

function plot_path = render_targeting_scatter(group_table, out_path)
% MNI centroid projections (colour = condition, + = targets) and the
% axial vs lateral focus error.
    plot_path = '';
    try
        cols = group_table.Properties.VariableNames;
        has_mni = all(ismember({'focal_centroid_mni_x', 'focal_centroid_mni_y', 'focal_centroid_mni_z'}, cols)) && ...
            any(isfinite(group_table.focal_centroid_mni_x));
        has_err = all(ismember({'focus_error_axial_mm', 'focus_error_lateral_mm'}, cols)) && ...
            any(isfinite(group_table.focus_error_axial_mm));
        if ~has_mni && ~has_err; return; end
        conds = unique(string(group_table.condition), 'stable');
        colors = lines(max(1, numel(conds)));
        nt = 3 * has_mni + has_err;
        fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 300 * nt 320]);
        tl = tiledlayout(fig, 1, nt, 'TileSpacing', 'compact', 'Padding', 'compact');
        if has_mni
            C = [group_table.focal_centroid_mni_x, group_table.focal_centroid_mni_y, group_table.focal_centroid_mni_z];
            G = nan(height(group_table), 3);
            if all(ismember({'target_mni_x', 'target_mni_y', 'target_mni_z'}, cols))
                G = [group_table.target_mni_x, group_table.target_mni_y, group_table.target_mni_z];
            end
            planes = {[2 3], 'Sagittal', 'y (mm)', 'z (mm)'; [1 3], 'Coronal', 'x (mm)', 'z (mm)'; ...
                      [1 2], 'Axial', 'x (mm)', 'y (mm)'};
            for p = 1:3
                ax = nexttile(tl); hold(ax, 'on');
                d = planes{p, 1};
                for c = 1:numel(conds)
                    in = string(group_table.condition) == conds(c);
                    lbl = char(conds(c)); if isempty(lbl); lbl = 'runs'; end
                    scatter(ax, C(in, d(1)), C(in, d(2)), 28, colors(c, :), 'filled', ...
                        'MarkerFaceAlpha', 0.75, 'DisplayName', lbl);
                end
                ok = all(isfinite(G), 2);
                if any(ok)
                    Gu = unique(round(G(ok, :), 2), 'rows');
                    plot(ax, Gu(:, d(1)), Gu(:, d(2)), 'k+', 'MarkerSize', 12, 'LineWidth', 1.8, ...
                        'DisplayName', 'target');
                end
                axis(ax, 'equal'); grid(ax, 'on');
                title(ax, planes{p, 2}, 'FontWeight', 'normal', 'FontSize', 9);
                xlabel(ax, planes{p, 3}); ylabel(ax, planes{p, 4}); set(ax, 'FontSize', 8);
                if p == 1 && numel(conds) > 1; legend(ax, 'Location', 'best', 'Interpreter', 'none', 'FontSize', 7); end
            end
        end
        if has_err
            ax = nexttile(tl); hold(ax, 'on');
            for c = 1:numel(conds)
                in = string(group_table.condition) == conds(c);
                scatter(ax, group_table.focus_error_lateral_mm(in), group_table.focus_error_axial_mm(in), ...
                    28, colors(c, :), 'filled', 'MarkerFaceAlpha', 0.75);
            end
            yline(ax, 0, 'k:');
            grid(ax, 'on'); set(ax, 'FontSize', 8);
            xlabel(ax, 'lateral error (mm)'); ylabel(ax, 'axial error (mm, + = deeper)');
            title(ax, 'Focus vs planned focus', 'FontWeight', 'normal', 'FontSize', 9);
        end
        export_figure(fig, out_path);
        plot_path = out_path;
    catch ME
        warning('generate_group_report:plotFailed', 'Targeting plot failed: %s', ME.message);
    end
end

function export_figure(fig, out_path)
    if exist('exportgraphics', 'file') == 2 || exist('exportgraphics', 'builtin') == 5
        exportgraphics(fig, out_path, 'Resolution', 150);
    else
        print(fig, out_path, '-dpng', '-r150');
    end
    close(fig);
end

function s = msd(v)
    v = v(isfinite(v));
    if isempty(v)
        s = '&mdash;';
    elseif numel(v) == 1
        s = sprintf('%.3g (n=1)', v);
    else
        s = sprintf('%.3g &plusmn; %.2g (n=%d)', mean(v), std(v), numel(v));
    end
end

function s = num_field(m, f)
    if isfield(m, f) && isnumeric(m.(f)) && isfinite(m.(f))
        s = sprintf('%.4g', m.(f));
    else
        s = '';
    end
end

function s = xyz_text(t, i, name)
    cols = {[name '_x'], [name '_y'], [name '_z']};
    if ~all(ismember(cols, t.Properties.VariableNames)); s = ''; return; end
    v = [t.(cols{1})(i), t.(cols{2})(i), t.(cols{3})(i)];
    if all(isfinite(v)); s = sprintf('[%.1f, %.1f, %.1f]', v); else, s = ''; end
end

function s = overrides_text(ov)
    if isempty(ov); s = ''; return; end
    parts = cell(1, size(ov, 1));
    for i = 1:size(ov, 1)
        parts{i} = sprintf('<code>%s</code> = %s', html_utils.escape(ov{i, 1}), ...
            html_utils.escape(value_text(ov{i, 2})));
    end
    s = strjoin(parts, '<br>');
end

function s = value_text(v)
    if ischar(v) || (isstring(v) && isscalar(v))
        s = char(v);
    elseif isnumeric(v) || islogical(v)
        if isempty(v); s = '[]';
        elseif isscalar(v); s = sprintf('%g', v);
        else, s = ['[' strjoin(compose('%g', double(v(:)')), ', ') ']'];
        end
    else
        s = class(v);
    end
end
