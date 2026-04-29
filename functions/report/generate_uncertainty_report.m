function report_path = generate_uncertainty_report(parameters, affixes)
% GENERATE_UNCERTAINTY_REPORT  Generate an uncertainty-range HTML simulation report
%
% Loads acoustic and thermal CSV tables and heating .mat files from three
% simulation variants (default, liberal, conservative), computes ranges
% across them, and writes a self-contained HTML report.
%
% Use as:
%   report_path = generate_uncertainty_report(parameters)
%   report_path = generate_uncertainty_report(parameters, affixes)
%
% Input:
%   parameters - (1,1) simulation parameters struct; parameters.io.dir_output
%                must point to the shared output folder
%   affixes    - struct with fields .default, .liberal, .conservative
%                (defaults: '', '_liberal', '_conservative')
%
% Output:
%   report_path - path to the generated HTML report file
%
% See also: GENERATE_SIMULATION_REPORT, GET_RISK_LIMITS

arguments
    parameters  (1,1) struct
    affixes     (1,1) struct = struct('default', '', 'liberal', '_desc-liberal', 'conservative', '_desc-conservative')
end

report_path = '';

try
    subject_id = parameters.subject_id;
    medium     = parameters.simulation.medium;
    output_dir = parameters.io.dir_output;

    % Determine whether this is a layered (tissue) simulation based on the
    % medium name. The modules.run_heating_sims flag reflects whether heating
    % was run in this particular pipeline call, NOT whether the medium is
    % tissue — so it must not be used here (the report stage always sets it
    % to 0, which would incorrectly trigger the free-water branch).
    is_layered = contains(medium, {'layered', 'phantom'});

    %% Load the three CSV tables
    variant_affixes = {affixes.liberal, affixes.default, affixes.conservative};
    variant_labels  = {'Liberal', 'Default', 'Conservative'};

    tables = cell(1, 3);
    for v = 1:3
        csv_path = fullfile(output_dir, ...
            sprintf('sub-%03d_%s%s.csv', subject_id, medium, variant_affixes{v}));
        if isfile(csv_path)
            try
                tables{v} = readtable(csv_path, 'VariableNamingRule', 'preserve');
            catch
                tables{v} = [];
            end
        end
    end

    if all(cellfun(@isempty, tables))
        warn('generate_uncertainty_report:noData', ...
            'No output tables found in %s — report will be empty.', output_dir);
    end

    %% Load heating timeseries from .mat files (if available)
    % Only load the fields needed for the SVG timeseries plot to limit memory use.
    heating_fields = {'timeseries', 't_array_s', 'focal_plume_T', 'T_time'};
    heating = cell(1, 3);
    for v = 1:3
        mat_path = fullfile(output_dir, 'cache', ...
            sprintf('sub-%03d_%s%s_heating_res.mat', subject_id, medium, variant_affixes{v}));
        if isfile(mat_path)
            try
                s = load(mat_path, 'results_heating');
                if isfield(s, 'results_heating')
                    % Copy only timeseries-relevant fields to avoid loading
                    % large temperature-volume arrays into memory.
                    h = struct();
                    for fi = 1:numel(heating_fields)
                        fn = heating_fields{fi};
                        if isfield(s.results_heating, fn)
                            h.(fn) = s.results_heating.(fn);
                        end
                    end
                    heating{v} = h;
                end
                clear s;
            catch
                heating{v} = [];
            end
        end
    end

    %% Build HTML
    html_parts = {};
    html_parts{end+1} = '<!DOCTYPE html>';
    html_parts{end+1} = '<html lang="en">';
    html_parts{end+1} = '<head>';
    html_parts{end+1} = '<meta charset="UTF-8">';
    html_parts{end+1} = '<meta name="viewport" content="width=device-width, initial-scale=1.0">';
    html_parts{end+1} = sprintf('<title>PRESTUS Uncertainty Report — sub-%03d %s</title>', subject_id, medium);
    html_parts{end+1} = '<style>';
    html_parts{end+1} = css_styles();
    html_parts{end+1} = '</style>';
    html_parts{end+1} = '</head>';
    html_parts{end+1} = '<body>';

    % TOC
    html_parts{end+1} = build_toc(is_layered, any(~cellfun(@isempty, heating)));

    % Header
    html_parts{end+1} = build_header(subject_id, medium, parameters, affixes);

    % Safety dashboard (ranges)
    try
        html_parts{end+1} = build_safety_dashboard(tables, variant_labels, is_layered);
    catch ME
        html_parts{end+1} = html_utils.section_error('Safety Dashboard', ME);
    end

    % Acoustic summary
    try
        html_parts{end+1} = html_utils.collapsible('Acoustic Results', ...
            build_acoustic_section(tables, variant_labels, variant_affixes, parameters, subject_id, medium, is_layered), ...
            true, 'acoustic');
    catch ME
        html_parts{end+1} = html_utils.section_error('Acoustic Results', ME);
    end

    % Thermal summary — show whenever heating data or thermal table columns exist,
    % regardless of medium name.
    has_thermal = any(~cellfun(@isempty, heating)) || ...
        (any(~cellfun(@isempty, tables)) && is_layered);
    if has_thermal
        try
            html_parts{end+1} = html_utils.collapsible('Thermal Results', ...
                build_thermal_section(tables, heating, variant_labels, variant_affixes, parameters, subject_id, medium), ...
                true, 'thermal');
        catch ME
            html_parts{end+1} = html_utils.section_error('Thermal Results', ME);
        end
    end

    % Side-by-side intensity images
    try
        html_parts{end+1} = html_utils.collapsible('Additional Maps', ...
            build_images_section(variant_affixes, variant_labels, parameters, subject_id, medium, is_layered), ...
            true, 'images');
    catch ME
        html_parts{end+1} = html_utils.section_error('Additional Maps', ME);
    end

    % Timing summary
    try
        timing_html = build_timing_section(parameters, variant_affixes, variant_labels);
        if ~isempty(timing_html)
            html_parts{end+1} = html_utils.collapsible('Pipeline Timing', timing_html, false, 'timing');
        end
    catch ME
        html_parts{end+1} = html_utils.section_error('Pipeline Timing', ME);
    end

    % Raw data tables
    try
        html_parts{end+1} = html_utils.collapsible('Raw Data Tables', ...
            build_tables_section(tables, variant_labels, is_layered), ...
            false, 'tables');
    catch ME
        html_parts{end+1} = html_utils.section_error('Raw Data Tables', ME);
    end

    % Footer
    html_parts{end+1} = '<footer>';
    html_parts{end+1} = sprintf('<p>Generated by PRESTUS uncertainty report — %s</p>', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    html_parts{end+1} = '</footer>';
    html_parts{end+1} = html_utils.lightbox();
    html_parts{end+1} = '</body>';
    html_parts{end+1} = '</html>';

    %% Write file
    report_filename = sprintf('sub-%03d_%s_desc-uncertainty_report.html', subject_id, medium);
    report_path = fullfile(output_dir, report_filename);
    fid = fopen(report_path, 'w', 'n', 'UTF-8');
    if fid == -1
        warn('generate_uncertainty_report:fileOpen', 'Cannot open %s for writing.', report_path);
        return
    end
    fprintf(fid, '%s\n', html_parts{:});
    fclose(fid);
    fprintf('Uncertainty report saved to: %s\n', report_path);

catch ME
    warn('generate_uncertainty_report:failed', ...
        'Report generation failed: %s\n%s', ME.message, getReport(ME, 'extended'));
    report_path = '';
end

end

%% ========================================================================
%  SECTION BUILDERS
%% ========================================================================

function html = build_header(subject_id, medium, parameters, affixes)
    html = '<section class="report-section" id="header">';
    html = [html '<h1>PRESTUS Uncertainty Report</h1>'];
    html = [html '<table class="info-table">'];
    html = [html sprintf('<tr><th>Subject</th><td>sub-%03d</td></tr>', subject_id)];
    html = [html sprintf('<tr><th>Medium</th><td>%s</td></tr>', html_utils.escape(medium))];
    html = [html sprintf('<tr><th>Default affix</th><td><code>%s</code></td></tr>', html_utils.escape(affixes.default))];
    html = [html sprintf('<tr><th>Liberal affix</th><td><code>%s</code></td></tr>', html_utils.escape(affixes.liberal))];
    html = [html sprintf('<tr><th>Conservative affix</th><td><code>%s</code></td></tr>', html_utils.escape(affixes.conservative))];
    html = [html sprintf('<tr><th>Output dir</th><td>%s</td></tr>', html_utils.escape(parameters.io.dir_output))];
    html = [html sprintf('<tr><th>Generated</th><td>%s</td></tr>', datestr(now, 'yyyy-mm-dd HH:MM:SS'))];
    html = [html '</table>'];
    html = [html '</section>'];
end

function html = build_safety_dashboard(tables, variant_labels, is_layered)
    limits = struct();
    if is_layered
        limits = get_risk_limits(true);
    end
    metric_names = fieldnames(limits);

    html = '<section class="report-section" id="safety">';
    html = [html '<h2>Safety Dashboard</h2>'];

    if ~is_layered
        html = [html '<div class="medium-banner medium-water">' ...
            '<strong>Water / Free-field</strong> &mdash; tissue-specific non-significant risk limits not applicable.</div>'];
    else
        html = [html '<div class="medium-banner medium-layered">' ...
            'Values shown for the <strong>default</strong> variant. ' ...
            'Range annotations span liberal to conservative.</div>'];
    end

    html = [html '<div class="safety-grid">'];

    % variant order: liberal=1, default=2, conservative=3
    for i = 1:length(metric_names)
        name = metric_names{i};
        info = limits.(name);

        % Extract values from all three variants
        vals = NaN(1, 3);
        for v = 1:3
            if ~isempty(tables{v}) && ismember(name, tables{v}.Properties.VariableNames)
                col = tables{v}.(name);
                if isnumeric(col) && ~isempty(col)
                    vals(v) = col(end);
                end
            end
        end

        % Display default (index 2); bar length and color driven by conservative (index 3)
        display_val = vals(2);
        liberal_val = vals(1);
        cons_val    = vals(3);

        display_unit = info.unit;
        color_val = cons_val;
        if isnan(color_val); color_val = display_val; end
        color = risk_color(color_val, info.limit);

        html = [html sprintf('<div class="safety-card safety-%s">', color)];
        html = [html sprintf('<div class="safety-label">%s</div>', html_utils.escape(info.label))];

        if isnan(display_val)
            html = [html '<div class="safety-value">N/A</div>'];
            html = [html '<div class="risk-limit">No data</div>'];
        else
            html = [html sprintf('<div class="safety-value">%.3g</div>', display_val)];

            % Range annotation: liberal to conservative
            if ~isnan(liberal_val) && ~isnan(cons_val)
                html = [html sprintf('<div class="risk-limit">%.3g &ndash; %.3g %s (lib&ndash;cons)</div>', ...
                    liberal_val, cons_val, display_unit)];
            elseif ~isinf(info.limit)
                html = [html sprintf('<div class="risk-limit">NSR limit: %.3g %s</div>', info.limit, display_unit)];
            else
                html = [html sprintf('<div class="risk-limit">%s (informational)</div>', display_unit)];
            end

            if ~isinf(info.limit)
                % Axis lower bound: hardcoded baseline
                if strncmp(name, 'maxT', 4) || strncmp(name, 'endT', 4)
                    scale_min = 37;
                else
                    scale_min = 0;
                end

                all_vals_finite = [liberal_val, display_val, cons_val];
                all_vals_finite = all_vals_finite(~isnan(all_vals_finite));
                min_sim = min(all_vals_finite);
                max_sim = max(all_vals_finite);

                % Scale runs from scale_min to max(max_sim, limit)
                scale_max = max([max_sim, info.limit]);
                if scale_max <= scale_min; scale_max = scale_min + 1; end
                scale_range = scale_max - scale_min;

                % Bar spans min to max of simulated values across variants
                bar_left  = max((min_sim - scale_min) / scale_range * 100, 0);
                bar_right = min((max_sim  - scale_min) / scale_range * 100, 100);
                bar_width = max(bar_right - bar_left, 0);

                lim_pct = (info.limit - scale_min) / scale_range * 100;
                limit_line = sprintf('<div class="safety-bar-limit" style="left:%.1f%%"></div>', lim_pct);

                markers = '';
                if ~isnan(liberal_val)
                    pct = (liberal_val - scale_min) / scale_range * 100;
                    markers = [markers sprintf('<div class="safety-bar-marker safety-bar-marker--liberal" style="left:%.1f%%"></div>', pct)];
                end
                if ~isnan(display_val)
                    pct = (display_val - scale_min) / scale_range * 100;
                    markers = [markers sprintf('<div class="safety-bar-marker safety-bar-marker--default" style="left:%.1f%%"></div>', pct)];
                end
                if ~isnan(cons_val)
                    pct = (cons_val - scale_min) / scale_range * 100;
                    markers = [markers sprintf('<div class="safety-bar-marker safety-bar-marker--conservative" style="left:%.1f%%"></div>', pct)];
                end

                % Scale range printed directly above the bar
                range_str = sprintf('%.3g &ndash; %.3g %s', scale_min, scale_max, display_unit);
                html = [html sprintf('<div class="safety-bar-range">%s</div>', range_str)];
                html = [html sprintf('<div class="safety-bar-track"><div class="safety-bar" style="left:%.1f%%;width:%.1f%%"></div>%s%s</div>', bar_left, bar_width, limit_line, markers)];
            end
        end

        html = [html '</div>'];
    end

    html = [html '</div>'];
    html = [html '<p class="safety-footnote">Non-significant risk limits: ITRUSST consensus (Aubry et al., 2025). ' ...
        '<strong>Card colour</strong> reflects the <strong>conservative</strong> variant: ' ...
        'green &lt;50% of limit, amber 50&ndash;100%, red exceeds limit. ' ...
        '<strong>Value</strong> shown for the <em>default</em> variant; ' ...
        '<strong>lib&ndash;cons range</strong> printed below the value. ' ...
        '<strong>Scale</strong> (above bar): fixed axis range — baseline (0 or 37&thinsp;&deg;C for temperature) to max(simulated, NSR limit). ' ...
        '<strong>Bar</strong>: simulated range from minimum (liberal) to maximum (conservative) variant. ' ...
        'Grey dashed line = NSR limit. ' ...
        'Ticks: &#x25A0;&nbsp;black&nbsp;=&nbsp;default, ' ...
        '<span style="color:#2563eb">&#x25A0;</span>&nbsp;blue&nbsp;=&nbsp;liberal, ' ...
        '<span style="color:#92400e">&#x25A0;</span>&nbsp;brown&nbsp;=&nbsp;conservative.</p>'];
    html = [html '</section>'];
end

function html = build_acoustic_section(tables, variant_labels, variant_affixes, parameters, subject_id, medium, is_layered)
    html = '';

    % Summary range table
    if is_layered
        key_metrics = {'Isppa', 'Isppa_brain', 'Isppa_skull', 'Isppa_skin', ...
                       'real_focal_distance_mm', 'MI_brain', 'MI_skull', 'MI_skin', 'MI_tc', ...
                       'Ipa_target', 'halfmax_ISPPA_volume_brain_mm3'};
    else
        key_metrics = {'Isppa', 'Isppa_after_exitplane', 'Psptp', 'Ptp_target', ...
                       'real_focal_distance_mm', 'Ipa_target'};
    end

    html = [html '<h3>Key Acoustic Metrics</h3>'];
    html = [html build_range_table(tables, variant_labels, key_metrics)];

    % Intensity images — y-slice only (3 columns: liberal | default | conservative)
    html = [html '<h3 style="margin-top:24px;">Intensity Maps (Liberal | Default | Conservative)</h3>'];
    html = [html '<div class="image-grid image-grid--3col">'];
    any_intensity = false;
    for v = 1:3
        aff = variant_affixes{v};
        lbl = variant_labels{v};
        img_path = fullfile(parameters.io.dir_img, ...
            sprintf('sub-%03d_%s%s_intensity_y.png', subject_id, medium, aff));
        img_html = html_utils.embed_image(img_path, lbl, lbl);
        if ~isempty(img_html)
            html = [html img_html];
            any_intensity = true;
        else
            html = [html sprintf('<figure><div class="placeholder-img">%s<br>Image not found</div></figure>', html_utils.escape(lbl))];
        end
    end
    if ~any_intensity
        % strip the heading if no images at all
        html = regexprep(html, '<h3[^>]*>Intensity Maps[^<]*</h3>', '');
    end
    html = [html '</div>'];
end

function html = build_thermal_section(tables, heating, variant_labels, variant_affixes, parameters, subject_id, medium)
    html = '';

    key_metrics = {'maxT', 'maxT_brain', 'maxT_skull', 'maxT_skin', ...
                   'riseT_brain', 'riseT_skull', 'riseT_skin', ...
                   'CEM43_brain', 'CEM43_skull', 'CEM43_skin'};

    html = [html '<h3>Thermal Metrics</h3>'];
    html = [html build_range_table(tables, variant_labels, key_metrics)];

    % Temperature vs. Time images (thermal_max, Liberal | Default | Conservative)
    html = [html '<h3 style="margin-top:24px;">Temperature vs. Time (Liberal | Default | Conservative)</h3>'];
    html = [html '<div class="image-grid image-grid--3col">'];
    any_thermal_max = false;
    for v = 1:3
        aff = variant_affixes{v};
        lbl = variant_labels{v};
        img_path = fullfile(parameters.io.dir_img, ...
            sprintf('sub-%03d_%s%s_thermal_max.png', subject_id, medium, aff));
        img_html = html_utils.embed_image(img_path, lbl, lbl);
        if ~isempty(img_html)
            html = [html img_html];
            any_thermal_max = true;
        else
            html = [html sprintf('<figure><div class="placeholder-img">%s<br>Image not found</div></figure>', html_utils.escape(lbl))];
        end
    end
    if ~any_thermal_max
        html = regexprep(html, '<h3[^>]*>Temperature vs\. Time[^<]*</h3>', '');
    end
    html = [html '</div>'];
end

function html = build_images_section(variant_affixes, variant_labels, parameters, subject_id, medium, is_layered)
    html = '';

    image_specs = {};
    if is_layered
        image_specs{end+1} = {fullfile(parameters.io.dir_img,  sprintf('sub-%03d_%s%%s_maxT_y.png',   subject_id, medium)), 'Max Temperature (y-slice)'};
        image_specs{end+1} = {fullfile(parameters.io.dir_img,  sprintf('sub-%03d_%s%%s_thermal.png', subject_id, medium)), 'Temperature vs. Time'};
        image_specs{end+1} = {fullfile(parameters.io.dir_img,  sprintf('sub-%03d_%s%%s_CEM.png',     subject_id, medium)), 'CEM43 vs. Time'};
    end

    for k = 1:numel(image_specs)
        pattern = image_specs{k}{1};
        title   = image_specs{k}{2};
        html = [html sprintf('<h3>%s</h3>', html_utils.escape(title))];
        html = [html '<div class="image-grid image-grid--3col">'];
        any_found = false;
        for v = 1:3
            aff = variant_affixes{v};
            lbl = variant_labels{v};
            img_path = sprintf(pattern, aff);
            img_html = html_utils.embed_image(img_path, lbl, lbl);
            if ~isempty(img_html)
                html = [html img_html];
                any_found = true;
            else
                html = [html sprintf('<figure><div class="placeholder-img">%s<br>Not found</div></figure>', html_utils.escape(lbl))];
            end
        end
        if ~any_found
            html = strrep(html, sprintf('<h3>%s</h3>', html_utils.escape(title)), '');
        end
        html = [html '</div>'];
    end
end

function html = build_tables_section(tables, variant_labels, is_layered)
    html = '';
    for v = 1:3
        if isempty(tables{v})
            html = [html sprintf('<h3>%s</h3><p class="placeholder">Table not found.</p>', ...
                html_utils.escape(variant_labels{v}))];
            continue;
        end
        html = [html sprintf('<h3>%s</h3>', html_utils.escape(variant_labels{v}))];
        html = [html table2html(tables{v})];
    end
end

%% ========================================================================
%  RANGE TABLE
%% ========================================================================

function html = build_range_table(tables, variant_labels, metric_names)
% Build a table with one row per metric, columns: metric | conservative | default | liberal | range
    limits = struct();
    if is_layered_limits(metric_names)
        limits = get_risk_limits(true);
    end

    html = '<div class="table-wrapper"><table class="data-table"><thead><tr>';
    html = [html '<th>Metric</th>'];
    for v = 1:3
        html = [html sprintf('<th>%s</th>', html_utils.escape(variant_labels{v}))];
    end
    html = [html '<th>Range</th><th>Limit</th></tr></thead><tbody>'];

    for m = 1:numel(metric_names)
        name = metric_names{m};

        vals = NaN(1, 3);
        for v = 1:3
            if ~isempty(tables{v}) && ismember(name, tables{v}.Properties.VariableNames)
                col = tables{v}.(name);
                if isnumeric(col) && ~isempty(col)
                    vals(v) = col(end);
                end
            end
        end

        % Determine limit & color for worst-case (liberal = index 3)
        limit_val = Inf;
        limit_str = '—';
        if isfield(limits, name)
            limit_val = limits.(name).limit;
            if ~isinf(limit_val)
                limit_str = sprintf('%.3g %s', limit_val, limits.(name).unit);
            else
                limit_str = sprintf('inf %s', limits.(name).unit);
            end
        end
        worst = vals(3);
        if isnan(worst); worst = nanmax(vals); end
        row_color = risk_color(worst, limit_val);

        html = [html sprintf('<tr class="row-%s">', row_color)];
        html = [html sprintf('<td><strong>%s</strong></td>', html_utils.escape(name))];
        for v = 1:3
            if isnan(vals(v))
                html = [html '<td>N/A</td>'];
            else
                html = [html sprintf('<td>%.4g</td>', vals(v))];
            end
        end

        % Range column
        valid = vals(~isnan(vals));
        if numel(valid) >= 2
            html = [html sprintf('<td>%.3g &ndash; %.3g</td>', min(valid), max(valid))];
        elseif numel(valid) == 1
            html = [html sprintf('<td>%.3g</td>', valid)];
        else
            html = [html '<td>N/A</td>'];
        end

        html = [html sprintf('<td>%s</td>', html_utils.escape(limit_str))];
        html = [html '</tr>'];
    end

    html = [html '</tbody></table></div>'];
end

function tf = is_layered_limits(metric_names)
% Heuristic: if any metric is tissue-specific, use layered limits.
    layered_keys = {'MI_brain', 'MI_skull', 'MI_skin', 'MI_tc', ...
                    'riseT_brain', 'CEM43_brain', 'maxT_brain', ...
                    'Isppa_brain', 'Isppa_skull', 'Isppa_skin'};
    tf = any(ismember(metric_names, layered_keys));
end

%% ========================================================================
%  TIMESERIES SVG PLOT
%% ========================================================================

function html = build_timeseries_svg(heating, variant_labels)
% Build an inline SVG showing max temperature timeseries with uncertainty band.
% heating{v}.timeseries is a struct with fields T, Tdiff, CEM43, each containing
% per-layer vectors (e.g. timeseries.T.brain). We take the max across all layers
% to show the global maximum temperature over time.

    html = '';

    % Extract max-temperature timeseries per variant
    ts_data = cell(1, 3);
    n_pts = 0;
    max_temp_all = 0;

    for v = 1:3
        ts_line = [];
        if ~isempty(heating{v}) && isfield(heating{v}, 'timeseries')
            ts = heating{v}.timeseries;
            if isstruct(ts) && isfield(ts, 'T')
                % timeseries.T is a struct of per-layer vectors; take pointwise max
                layer_names = fieldnames(ts.T);
                for li = 1:numel(layer_names)
                    ldata = ts.T.(layer_names{li})(:);
                    if isempty(ts_line)
                        ts_line = ldata;
                    else
                        n_common = min(numel(ts_line), numel(ldata));
                        ts_line(1:n_common) = max(ts_line(1:n_common), ldata(1:n_common));
                    end
                end
            elseif isnumeric(ts)
                % fallback: numeric matrix, take max across columns
                if size(ts, 2) > 1
                    ts_line = max(ts, [], 2);
                else
                    ts_line = ts(:);
                end
            end
        end
        if ~isempty(ts_line)
            ts_data{v} = ts_line;
            n_pts = max(n_pts, numel(ts_line));
            max_temp_all = max(max_temp_all, max(ts_line(~isnan(ts_line))));
        end
    end

    if n_pts == 0
        html = '<p class="placeholder">No timeseries data available.</p>';
        return;
    end

    % Pad shorter timeseries with NaN
    for v = 1:3
        if isempty(ts_data{v})
            ts_data{v} = NaN(n_pts, 1);
        elseif numel(ts_data{v}) < n_pts
            ts_data{v}(end+1:n_pts) = NaN;
        end
    end

    % SVG dimensions
    W = 700; H = 300;
    pad_l = 60; pad_r = 20; pad_t = 20; pad_b = 50;
    plot_w = W - pad_l - pad_r;
    plot_h = H - pad_t - pad_b;

    y_min = 37;   % body temperature baseline
    y_max = max(max_temp_all * 1.02, y_min + 0.1);
    y_range = y_max - y_min;
    x_scale = plot_w / max(n_pts - 1, 1);
    y_scale = plot_h / y_range;

    % Build band polygon: liberal (upper, index 1) forward, conservative (lower, index 3) backward
    band_pts = '';
    for i = 1:n_pts
        x = pad_l + (i-1) * x_scale;
        val = ts_data{1}(i);   % liberal
        if isnan(val); val = 0; end
        y = pad_t + plot_h - (val - y_min) * y_scale;
        band_pts = [band_pts sprintf('%.1f,%.1f ', x, y)];
    end
    for i = n_pts:-1:1
        x = pad_l + (i-1) * x_scale;
        val = ts_data{3}(i);   % conservative
        if isnan(val); val = 0; end
        y = pad_t + plot_h - (val - y_min) * y_scale;
        band_pts = [band_pts sprintf('%.1f,%.1f ', x, y)];
    end

    % Build default line
    default_pts = '';
    for i = 1:n_pts
        x = pad_l + (i-1) * x_scale;
        val = ts_data{2}(i);
        if isnan(val); continue; end
        y = pad_t + plot_h - (val - y_min) * y_scale;
        default_pts = [default_pts sprintf('%.1f,%.1f ', x, y)];
    end

    % Y-axis ticks
    n_ticks = 5;
    tick_html = '';
    for k = 0:n_ticks
        val = y_min + y_range * k / n_ticks;
        y = pad_t + plot_h - (val - y_min) * y_scale;
        tick_html = [tick_html ...
            sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e2e8f0" stroke-width="1"/>', ...
                pad_l, y, W-pad_r, y) ...
            sprintf('<text x="%d" y="%.1f" text-anchor="end" font-size="10" fill="#94a3b8">%.3g</text>', ...
                pad_l-4, y+4, val)];
    end

    % X-axis ticks (5 evenly spaced)
    x_tick_html = '';
    for k = 0:4
        idx = round(k * (n_pts-1) / 4) + 1;
        x = pad_l + (idx-1) * x_scale;
        x_tick_html = [x_tick_html ...
            sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="10" fill="#94a3b8">%d</text>', ...
                x, H-pad_b+14, idx)];
    end

    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' num2str(W) ' ' num2str(H) '" ' ...
           'style="width:100%;max-width:' num2str(W) 'px;height:auto;display:block;margin:0 auto 16px;">'];
    % Background
    svg = [svg sprintf('<rect width="%d" height="%d" fill="white" rx="4" ry="4"/>', W, H)];
    % Axes
    svg = [svg tick_html x_tick_html];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#cbd5e1" stroke-width="1.5"/>', ...
        pad_l, pad_t, pad_l, pad_t+plot_h)];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#cbd5e1" stroke-width="1.5"/>', ...
        pad_l, pad_t+plot_h, pad_l+plot_w, pad_t+plot_h)];
    % Axis labels
    svg = [svg sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="11" fill="#64748b">Time step</text>', ...
        pad_l + plot_w/2, H-4)];
    svg = [svg sprintf('<text x="%d" y="%.1f" text-anchor="middle" font-size="11" fill="#64748b" ', ...
        12, pad_t + plot_h/2) ...
        'transform="rotate(-90,' num2str(12) ',' num2str(pad_t+plot_h/2) ')">Temp (&#176;C)</text>'];
    % Uncertainty band
    if ~isempty(strtrim(band_pts))
        svg = [svg sprintf('<polygon points="%s" fill="#3b82f6" fill-opacity="0.15" stroke="none"/>', band_pts)];
    end
    % Default line
    if ~isempty(strtrim(default_pts))
        svg = [svg sprintf('<polyline points="%s" fill="none" stroke="#1d4ed8" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>', ...
            default_pts)];
    end
    % Liberal line (index 1)
    lib_pts = '';
    for i = 1:n_pts
        x = pad_l + (i-1) * x_scale;
        val = ts_data{1}(i);
        if isnan(val); continue; end
        y = pad_t + plot_h - (val - y_min) * y_scale;
        lib_pts = [lib_pts sprintf('%.1f,%.1f ', x, y)];
    end
    if ~isempty(strtrim(lib_pts))
        svg = [svg sprintf('<polyline points="%s" fill="none" stroke="#ef4444" stroke-width="1.5" stroke-dasharray="4,3"/>', lib_pts)];
    end
    % Conservative line (index 3)
    cons_pts = '';
    for i = 1:n_pts
        x = pad_l + (i-1) * x_scale;
        val = ts_data{3}(i);
        if isnan(val); continue; end
        y = pad_t + plot_h - (val - y_min) * y_scale;
        cons_pts = [cons_pts sprintf('%.1f,%.1f ', x, y)];
    end
    if ~isempty(strtrim(cons_pts))
        svg = [svg sprintf('<polyline points="%s" fill="none" stroke="#94a3b8" stroke-width="1.5" stroke-dasharray="4,3"/>', cons_pts)];
    end
    % Legend
    lx = pad_l + plot_w - 150; ly = pad_t + 10;
    svg = [svg sprintf('<rect x="%d" y="%d" width="150" height="64" rx="3" fill="white" fill-opacity="0.85" stroke="#e2e8f0"/>', lx, ly)];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#1d4ed8" stroke-width="2"/><text x="%d" y="%d" font-size="10" fill="#475569">Default</text>', ...
        lx+6, ly+16, lx+26, ly+16, lx+30, ly+20)];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#94a3b8" stroke-width="1.5" stroke-dasharray="4,3"/><text x="%d" y="%d" font-size="10" fill="#475569">Conservative</text>', ...
        lx+6, ly+32, lx+26, ly+32, lx+30, ly+36)];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#ef4444" stroke-width="1.5" stroke-dasharray="4,3"/><text x="%d" y="%d" font-size="10" fill="#475569">Liberal</text>', ...
        lx+6, ly+48, lx+26, ly+48, lx+30, ly+52)];
    svg = [svg '</svg>'];

    html = svg;
end

%% ========================================================================
%  PIPELINE TIMING
%% ========================================================================

function timing = parse_timing_from_log(log_path)
% Parse a diary log file and return a struct array with fields:
%   label (char), seconds (double).
% Returns empty struct array if the file is missing or parsing fails.
    timing = struct('label', {}, 'seconds', {});
    if isempty(log_path) || ~isfile(log_path); return; end
    try
        fid = fopen(log_path, 'r', 'n', 'UTF-8');
        if fid == -1; return; end
        lines = {};
        while ~feof(fid)
            lines{end+1} = fgetl(fid); %#ok<AGROW>
        end
        fclose(fid);
    catch
        return;
    end
    % Match timer-stop lines: ⏱ <label>  <seconds>s | ...
    % The ⏱ character is U+23F1 (UTF-8: E2 8F B1).
    timer_char = native2unicode([0xE2, 0x8F, 0xB1], 'UTF-8');
    pat = ['^' timer_char '\s+(\S+)\s+([\d.]+)s\b'];
    for li = 1:numel(lines)
        ln = lines{li};
        if ~ischar(ln) || isempty(ln); continue; end
        tok = regexp(ln, pat, 'tokens');
        if ~isempty(tok)
            t.label   = tok{1}{1};
            t.seconds = str2double(tok{1}{2});
            timing(end+1) = t; %#ok<AGROW>
        end
    end
end

function html = build_timing_section(parameters, variant_affixes, variant_labels)
% Build an HTML table showing per-step wall times across all five pipeline
% stages, using the pre-assigned log file paths in
% parameters.uncertainty.log_files.
% Returns empty string if no timing data is found.

    if ~isfield(parameters, 'uncertainty') || ~isfield(parameters.uncertainty, 'log_files')
        html = '';
        return;
    end
    lf = parameters.uncertainty.log_files;  % struct: stage1, default, liberal, conservative, report

    % Stage 1: preprocessing steps (preproc, medium, source)
    t_stage1 = parse_timing_from_log(char(lf.stage1));

    % Stages 2–4: simulation variants (acoustic, thermal, total pipeline)
    % variant_affixes = {liberal, default, conservative}
    t_variants = cell(1, 3);
    variant_log_fields = {'liberal', 'default', 'conservative'};
    for v = 1:3
        t_variants{v} = parse_timing_from_log(char(lf.(variant_log_fields{v})));
    end

    % Stage 5: report
    t_report = parse_timing_from_log(char(lf.report));

    any_data = ~isempty(t_stage1) || any(~cellfun(@isempty, t_variants)) || ~isempty(t_report);
    if ~any_data
        html = '';
        return;
    end

    % Step definitions: {log_label, display_label, source}
    % source: 'stage1' | 'variant' | 'report' | 'total'
    % 'total' spans all variant columns and the report column.
    steps = { ...
        'preproc',           'Head preprocessing',  'stage1';  ...
        'medium',            'Medium setup',         'stage1';  ...
        'source',            'Source setup',         'stage1';  ...
        'acoustic',          'Acoustic simulation',  'variant'; ...
        'acoustic_analysis', 'Acoustic analysis',    'variant'; ...
        'thermal',           'Thermal simulation',   'variant'; ...
        'thermal_analysis',  'Thermal analysis',     'variant'; ...
        'nifti',             'NIfTI export',         'variant'; ...
        'prestus_pipeline',  'Total (pipeline)',     'total';   ...
    };

    % Column headers: Stage 1 | Liberal | Default | Conservative | Report
    col_labels = [{'Stage 1'}, variant_labels, {'Report'}];
    col_colors = {'#475569', '#1d4ed8', '#1e293b', '#92400e', '#475569'};

    html = '<p style="font-size:0.82em;color:#64748b;margin:0 0 12px;">Wall times from pre-assigned log files for each pipeline stage. Steps skipped due to caching are shown as —.</p>';
    html = [html '<div style="overflow-x:auto"><table class="range-table" style="min-width:600px">'];
    html = [html '<thead><tr><th>Step</th>'];
    for c = 1:numel(col_labels)
        html = [html sprintf('<th style="color:%s">%s</th>', col_colors{c}, col_labels{c})]; %#ok<AGROW>
    end
    html = [html '</tr></thead><tbody>'];

    function s = lookup(t, key)
        s = NaN;
        for ti = 1:numel(t)
            if strcmp(t(ti).label, key); s = t(ti).seconds; return; end
        end
    end

    function str = fmt(s)
        if isnan(s); str = ''; return; end
        if s >= 3600;     str = sprintf('%.1f h',   s / 3600);
        elseif s >= 60;   str = sprintf('%.1f min', s / 60);
        else;             str = sprintf('%.0f s',   s);
        end
    end

    function str = cell_html(s)
        if isnan(s)
            str = '<td style="color:#94a3b8">—</td>';
        else
            str = sprintf('<td>%s</td>', fmt(s));
        end
    end

    for si = 1:size(steps, 1)
        key    = steps{si, 1};
        lbl    = steps{si, 2};
        src    = steps{si, 3};
        is_tot = strcmp(key, 'prestus_pipeline');
        row_style = '';
        if is_tot; row_style = ' style="font-weight:600;border-top:2px solid #e2e8f0"'; end

        html = [html sprintf('<tr%s><td>%s</td>', row_style, lbl)]; %#ok<AGROW>

        % Stage 1 column
        if strcmp(src, 'stage1')
            html = [html cell_html(lookup(t_stage1, key))]; %#ok<AGROW>
        elseif strcmp(src, 'total')
            html = [html '<td style="color:#94a3b8">—</td>']; %#ok<AGROW>
        else
            html = [html '<td style="color:#e2e8f0">·</td>']; %#ok<AGROW>
        end

        % Variant columns (liberal, default, conservative)
        for v = 1:3
            if strcmp(src, 'variant') || strcmp(src, 'total')
                html = [html cell_html(lookup(t_variants{v}, key))]; %#ok<AGROW>
            else
                html = [html '<td style="color:#e2e8f0">·</td>']; %#ok<AGROW>
            end
        end

        % Report column
        if strcmp(src, 'report') || strcmp(src, 'total')
            html = [html cell_html(lookup(t_report, key))]; %#ok<AGROW>
        else
            html = [html '<td style="color:#e2e8f0">·</td>']; %#ok<AGROW>
        end

        html = [html '</tr>']; %#ok<AGROW>
    end

    html = [html '</tbody></table></div>'];
end

%% ========================================================================
%  UTILITIES
%% ========================================================================

function html = build_toc(is_layered, has_thermal)
    html = '<nav class="toc" id="toc">';
    html = [html '<strong style="font-size:0.85em;color:#475569;">Uncertainty</strong>'];
    html = [html '<a href="#header">Header</a>'];
    html = [html '<a href="#safety">Safety</a>'];
    html = [html '<a href="#acoustic">Acoustic</a>'];
    if has_thermal && is_layered
        html = [html '<a href="#thermal">Thermal</a>'];
    end
    html = [html '<a href="#images">Images</a>'];
    html = [html '<a href="#timing">Timing</a>'];
    html = [html '<a href="#tables">Tables</a>'];
    html = [html '</nav>'];
end

%% ========================================================================
%  CSS STYLES
%% ========================================================================

function css = css_styles()
% Uncertainty-report CSS: base rules plus uncertainty-specific additions.
css = [css_styles_base() ...
newline ...
'/* Uncertainty-specific: positioned safety bars with markers */' newline ...
'.safety-bar-range { font-size: 0.72em; color: #64748b; margin-top: 6px; margin-bottom: 2px; }' newline ...
'.safety-bar-track { position: relative; height: 4px; background: #e2e8f0; border-radius: 2px; overflow: visible; }' newline ...
'.safety-bar { height: 100%; border-radius: 2px; position: absolute; z-index: 1; }' newline ...
'.safety-bar-marker { position: absolute; top: -5px; width: 4px; height: 14px; border-radius: 1px; transform: translateX(-50%); z-index: 2; }' newline ...
'.safety-bar-marker--default      { background: #1e293b; }' newline ...
'.safety-bar-marker--liberal      { background: #2563eb; }' newline ...
'.safety-bar-marker--conservative { background: #92400e; }' newline ...
'.safety-bar-limit { position: absolute; top: -3px; width: 2px; height: 10px; background: #94a3b8; border-radius: 0; transform: translateX(-50%); z-index: 1; border-left: 1px dashed #64748b; }' newline ...
'.safety-green .safety-bar { background: var(--green); }' newline ...
'.safety-amber .safety-bar { background: var(--amber); }' newline ...
'.safety-red   .safety-bar { background: var(--red); }' newline ...
newline ...
'/* Uncertainty-specific: range table */' newline ...
'.range-table { width: 100%; border-collapse: collapse; font-size: 0.9em; }' newline ...
'.range-table th, .range-table td { padding: 8px 12px; text-align: left; border-bottom: 1px solid #e2e8f0; }' newline ...
'.range-table thead { background: #f1f5f9; }' newline ...
];
end
