function report_path = generate_uncertainty_report(parameters, affixes)
% GENERATE_UNCERTAINTY_REPORT  Generate an uncertainty-range HTML simulation report
%
% Loads acoustic and thermal CSV tables and heating .mat files from the
% simulation variants (UNCERTAINTY_VARIANTS: default, Temp-H, Temp-L, Trans-H,
% Trans-L) and writes a self-contained HTML report. Acoustic metrics are
% ranged over the Trans pair, thermal metrics over the Temp pair; colours
% follow the highest value of the default and that pair.
%
% Use as:
%   report_path = generate_uncertainty_report(parameters)
%   report_path = generate_uncertainty_report(parameters, affixes)
%
% Input:
%   parameters - (1,1) simulation parameters struct; parameters.io.dir_output
%                must point to the shared output folder; the optional
%                parameters.uncertainty.report_suffix (the default variant's
%                affix, e.g. '_desc-isppa030') is appended to the report file name
%   affixes    - struct with one field per variant name (default: each
%                variant's descriptor, e.g. .temp_high = '_desc-temphigh')
%
% Output:
%   report_path - path to the generated HTML report file
%
% See also: GENERATE_SIMULATION_REPORT, GET_RISK_LIMITS, UNCERTAINTY_VARIANTS

arguments
    parameters  (1,1) struct
    affixes     (1,1) struct = default_affixes()
end

report_path = '';

try
    parameters = resolve_io_dirs(parameters);

    subject_id = parameters.subject_id;
    medium     = parameters.simulation.medium;
    output_dir = parameters.io.dir_output;

    % Determine whether this is a layered (tissue) simulation based on the
    % medium name. The modules.run_heating_sims flag reflects whether heating
    % was run in this particular pipeline call, NOT whether the medium is
    % tissue — so it must not be used here (the report stage always sets it
    % to 0, which would incorrectly trigger the free-water branch).
    is_layered = contains(medium, {'layered', 'phantom'});

    %% Load the variants' CSV tables
    variants        = uncertainty_variants();
    n_var           = numel(variants);
    variant_affixes = cellfun(@(n) affixes.(n), {variants.name}, 'UniformOutput', false);
    variant_labels  = {variants.label};

    tables = cell(1, n_var);
    for v = 1:n_var
        csv_path = fullfile(output_dir, ...
            sprintf('sub-%03d_%s%s.csv', subject_id, medium, variant_affixes{v}));
        if isfile(csv_path)
            try
                tables{v} = read_metrics_table(csv_path);
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
    heating = cell(1, n_var);
    for v = 1:n_var
        mat_path = fullfile(output_dir, 'cache', ...
            sprintf('sub-%03d_%s_heating_res%s.mat', subject_id, medium, variant_affixes{v}));
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
    html_parts{end+1} = build_header(subject_id, medium, parameters, variants, variant_affixes);

    % Safety dashboard (ranges)
    try
        html_parts{end+1} = build_safety_dashboard(tables, variants, is_layered);
    catch ME
        html_parts{end+1} = html_utils.section_error('Safety Dashboard', ME);
    end

    % Pseudo-CT skull of the default variant against the bound variants' homogeneous skull
    if pct_skull_mapped(parameters)
        try
            html_parts{end+1} = html_utils.collapsible('Pseudo-CT Skull vs. Bounds', ...
                build_pct_bounds_section(parameters, variants, affixes.default), true, 'pct');
        catch ME
            html_parts{end+1} = html_utils.section_error('Pseudo-CT Skull vs. Bounds', ME);
        end
    end

    % Acoustic summary
    try
        html_parts{end+1} = html_utils.collapsible('Acoustic Results', ...
            build_acoustic_section(tables, variants, variant_affixes, parameters, subject_id, medium, is_layered), ...
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
                build_thermal_section(tables, heating, variants, variant_affixes, parameters, subject_id, medium), ...
                true, 'thermal');
        catch ME
            html_parts{end+1} = html_utils.section_error('Thermal Results', ME);
        end
    end

    % Side-by-side intensity images
    try
        html_parts{end+1} = html_utils.collapsible('Additional Maps', ...
            build_images_section(variants, variant_affixes, parameters, subject_id, medium, is_layered), ...
            true, 'images');
    catch ME
        html_parts{end+1} = html_utils.section_error('Additional Maps', ME);
    end

    % Timing summary
    try
        timing_html = build_timing_section(parameters, variants);
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
    % Named after the default variant's affix: one report per base affix and,
    % for nested sweeps, per sweep point (UNCERTAINTY_PIPELINE).
    report_suffix = '';
    if isfield(parameters, 'uncertainty') && isfield(parameters.uncertainty, 'report_suffix')
        report_suffix = parameters.uncertainty.report_suffix;
    end
    report_filename = sprintf('sub-%03d_%s_desc-uncertainty_report%s.html', ...
        subject_id, medium, report_suffix);
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

function affixes = default_affixes()
    affixes = struct();
    for v = uncertainty_variants()
        affixes.(v.name) = v.desc;
    end
end

function trio = pair_trio(variants, pair)
% Indices [high, default, low] of a pair ('temp' | 'trans') in variants.
    trio = [find(strcmp({variants.pair}, pair) & strcmp({variants.bound}, 'high')), ...
            find(strcmp({variants.name}, 'default')), ...
            find(strcmp({variants.pair}, pair) & strcmp({variants.bound}, 'low'))];
end

function pair = metric_pair(name)
% Thermal metrics are ranged over the Temp pair, all others over the Trans pair.
    if ~isempty(regexp(name, '^(maxT|riseT|endT|CEM43)', 'once'))
        pair = 'temp';
    else
        pair = 'trans';
    end
end

function label = pair_label(pair)
    label = [upper(pair(1)) pair(2:end)];
end

function vals = metric_values(tables, name)
% Last row of column name in each variant's table (NaN when absent).
    vals = NaN(1, numel(tables));
    for v = 1:numel(tables)
        if ~isempty(tables{v}) && ismember(name, tables{v}.Properties.VariableNames)
            col = tables{v}.(name);
            if isnumeric(col) && ~isempty(col)
                vals(v) = col(end);
            end
        end
    end
end

function html = build_header(subject_id, medium, parameters, variants, variant_affixes)
    html = '<section class="report-section" id="header">';
    html = [html '<h1>PRESTUS Uncertainty Report</h1>'];
    html = [html '<table class="info-table">'];
    html = [html sprintf('<tr><th>Subject</th><td>sub-%03d</td></tr>', subject_id)];
    html = [html sprintf('<tr><th>Medium</th><td>%s</td></tr>', html_utils.escape(medium))];
    for v = 1:numel(variants)
        html = [html sprintf('<tr><th>%s affix</th><td><code>%s</code></td></tr>', ...
            html_utils.escape(variants(v).label), html_utils.escape(variant_affixes{v}))]; %#ok<AGROW>
    end
    edits = edited_bounds(parameters, variants);
    if isempty(edits)
        html = [html '<tr><th>Variant bounds</th><td>ITRUSST bounds (config/uncertainty/), unchanged</td></tr>'];
    else
        html = [html sprintf(['<tr><th>Variant bounds</th><td style="color:#b91c1c;font-weight:600;">' ...
            'Edited (simulation.uncertainty_overrides): %s. These variants no longer bracket the ' ...
            'ITRUSST bounds.</td></tr>'], html_utils.escape(strjoin(edits, '; ')))];
    end
    html = [html sprintf('<tr><th>Output dir</th><td>%s</td></tr>', html_utils.escape(parameters.io.dir_output))];
    html = [html sprintf('<tr><th>Generated</th><td>%s</td></tr>', datestr(now, 'yyyy-mm-dd HH:MM:SS'))];
    html = [html '</table>'];
    html = [html '</section>'];
end

function edits = edited_bounds(parameters, variants)
% 'Temp-H brain.sound_speed = 1560' per value in simulation.uncertainty_overrides.
    edits = {};
    if ~isfield(parameters, 'simulation') || ~isfield(parameters.simulation, 'uncertainty_overrides') || ...
            ~isstruct(parameters.simulation.uncertainty_overrides)
        return
    end
    ov = parameters.simulation.uncertainty_overrides;
    for v = variants(:)'
        if ~isfield(ov, v.name) || ~isstruct(ov.(v.name))
            continue
        end
        for t = fieldnames(ov.(v.name))'
            props = ov.(v.name).(t{1});
            if ~isstruct(props); continue; end
            for q = fieldnames(props)'
                edits{end+1} = sprintf('%s %s.%s = %g', v.label, t{1}, q{1}, props.(q{1})); %#ok<AGROW>
            end
        end
    end
end

function html = build_safety_dashboard(tables, variants, is_layered)
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
            'Ranges span the Trans pair (acoustic metrics) or the Temp pair (thermal metrics).</div>'];
    end

    html = [html '<div class="safety-grid">'];

    for i = 1:length(metric_names)
        name = metric_names{i};
        info = limits.(name);

        % Default, high and low variant of the pair that bounds this metric
        pair = metric_pair(name);
        trio = pair_trio(variants, pair);
        vals = metric_values(tables, name);
        high_val    = vals(trio(1));
        display_val = vals(trio(2));
        low_val     = vals(trio(3));

        % Colour follows the highest of default, high and low
        display_unit = info.unit;
        color_val = max([high_val, display_val, low_val], [], 'omitnan');
        color = risk_color(color_val, info.limit);

        html = [html sprintf('<div class="safety-card safety-%s">', color)];
        html = [html sprintf('<div class="safety-label">%s</div>', html_utils.escape(info.label))];

        if isnan(display_val)
            html = [html '<div class="safety-value">N/A</div>'];
            html = [html '<div class="risk-limit">No data</div>'];
        else
            html = [html sprintf('<div class="safety-value">%.3g</div>', display_val)];

            % Range annotation: low to high variant of the pair
            if ~isnan(high_val) && ~isnan(low_val)
                html = [html sprintf('<div class="risk-limit">%.3g &ndash; %.3g %s (%s)</div>', ...
                    min(low_val, high_val), max(low_val, high_val), display_unit, pair_label(pair))];
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

                all_vals_finite = [high_val, display_val, low_val];
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
                if ~isnan(low_val)
                    pct = (low_val - scale_min) / scale_range * 100;
                    markers = [markers sprintf('<div class="safety-bar-marker safety-bar-marker--low" style="left:%.1f%%"></div>', pct)];
                end
                if ~isnan(display_val)
                    pct = (display_val - scale_min) / scale_range * 100;
                    markers = [markers sprintf('<div class="safety-bar-marker safety-bar-marker--default" style="left:%.1f%%"></div>', pct)];
                end
                if ~isnan(high_val)
                    pct = (high_val - scale_min) / scale_range * 100;
                    markers = [markers sprintf('<div class="safety-bar-marker safety-bar-marker--high" style="left:%.1f%%"></div>', pct)];
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
        '<strong>Card colour</strong> reflects the highest value of the default and the metric''s pair ' ...
        '(Trans-H/L for acoustic metrics, Temp-H/L for thermal metrics): ' ...
        'green &lt;50% of limit, amber 50&ndash;100%, red exceeds limit. ' ...
        '<strong>Value</strong> shown for the <em>default</em> variant; ' ...
        'the pair''s <strong>range</strong> is printed below the value. ' ...
        '<strong>Scale</strong> (above bar): fixed axis range from baseline (0 or 37&thinsp;&deg;C for temperature) to max(simulated, NSR limit). ' ...
        '<strong>Bar</strong>: simulated range across the default and the pair. ' ...
        'Grey dashed line = NSR limit. ' ...
        'Ticks: &#x25A0;&nbsp;black&nbsp;=&nbsp;default, ' ...
        '<span style="color:#b91c1c">&#x25A0;</span>&nbsp;red&nbsp;=&nbsp;high (Temp-H / Trans-H), ' ...
        '<span style="color:#2563eb">&#x25A0;</span>&nbsp;blue&nbsp;=&nbsp;low (Temp-L / Trans-L).</p>'];
    html = [html '</section>'];
end

function html = build_acoustic_section(tables, variants, variant_affixes, parameters, subject_id, medium, is_layered)
    html = '';

    % Summary range table
    if is_layered
        key_metrics = {'Isppa', 'Isppa_brain', 'Isppa_skull', 'Isppa_skin', ...
                       'real_focal_distance_bowl_mm', 'MI_brain', 'MI_skull', 'MI_skin', 'MI_tc', ...
                       'Ipa_target', 'halfmax_ISPPA_volume_brain_mm3'};
    else
        key_metrics = {'Isppa', 'Isppa_after_exitplane', 'Psptp', 'Ptp_target', ...
                       'real_focal_distance_bowl_mm', 'Ipa_target'};
    end

    html = [html '<h3>Key Acoustic Metrics</h3>'];
    html = [html build_range_table(tables, variants, key_metrics)];

    % Intensity images — y-slice only (3 columns: Trans-H | Default | Trans-L)
    trio = pair_trio(variants, 'trans');
    html = [html sprintf('<h3 style="margin-top:24px;">Intensity Maps (%s | %s | %s)</h3>', variants(trio).label)];
    html = [html '<div class="image-grid image-grid--3col">'];
    any_intensity = false;
    for v = trio
        aff = variant_affixes{v};
        lbl = variants(v).label;
        img_path = fullfile(parameters.io.dir_img, ...
            sprintf('sub-%03d_%s_intensity_y%s.png', subject_id, medium, aff));
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

function html = build_thermal_section(tables, heating, variants, variant_affixes, parameters, subject_id, medium)
    html = '';

    key_metrics = {'maxT', 'maxT_brain', 'maxT_skull', 'maxT_skin', ...
                   'riseT_brain', 'riseT_skull', 'riseT_skin', ...
                   'CEM43_brain', 'CEM43_skull', 'CEM43_skin'};

    html = [html '<h3>Thermal Metrics</h3>'];
    html = [html build_range_table(tables, variants, key_metrics)];

    % Temperature vs. Time images (thermal_max, Temp-H | Default | Temp-L)
    trio = pair_trio(variants, 'temp');
    html = [html sprintf('<h3 style="margin-top:24px;">Temperature vs. Time (%s | %s | %s)</h3>', variants(trio).label)];
    html = [html '<div class="image-grid image-grid--3col">'];
    any_thermal_max = false;
    for v = trio
        aff = variant_affixes{v};
        lbl = variants(v).label;
        img_path = fullfile(parameters.io.dir_img, ...
            sprintf('sub-%03d_%s_thermal_max%s.png', subject_id, medium, aff));
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

    % Uncertainty band timeseries plots (separate plots, generated from .mat data)
    % over the Temp pair, in [high, default, low] order.
    heating = heating(trio);
    has_ts = any(~cellfun(@isempty, heating));
    if has_ts
        caption = ['<p style="font-size:0.85em;color:#64748b;margin:4px 0 12px;">' ...
            'Shaded area: Temp-H&#8211;Temp-L range. ' ...
            'Solid line: default variant. Dashed lines: Temp-H (red) / Temp-L (gray).</p>'];

        html = [html '<h3 style="margin-top:28px;">Temperature Timeseries — Global Maximum (Uncertainty)</h3>'];
        html = [html caption];
        html = [html build_metric_uncertainty_svg(heating, 'T', 37, 'Temp (&#176;C)', false)];

        html = [html '<h3 style="margin-top:20px;">CEM43 Timeseries — Global Maximum (Uncertainty)</h3>'];
        html = [html caption];
        html = [html build_metric_uncertainty_svg(heating, 'CEM43', 0, 'CEM43 (min)', false)];

        % Per-layer plots when layered data is present
        layer_names = get_ts_layer_names(heating, 'T');
        if ~isempty(layer_names)
            html = [html '<h3 style="margin-top:20px;">Per-Layer Temperature Timeseries (Uncertainty)</h3>'];
            html = [html caption];
            for li = 1:numel(layer_names)
                lname = layer_names{li};
                html = [html sprintf('<h4 style="margin:12px 0 4px;color:#475569;">%s</h4>', ...
                    html_utils.escape(strrep(lname, '_', ' ')))];
                html = [html build_metric_uncertainty_svg(heating, 'T', 37, 'Temp (&#176;C)', lname)];
            end

            layer_names_cem = get_ts_layer_names(heating, 'CEM43');
            if ~isempty(layer_names_cem)
                html = [html '<h3 style="margin-top:20px;">Per-Layer CEM43 Timeseries (Uncertainty)</h3>'];
                html = [html caption];
                for li = 1:numel(layer_names_cem)
                    lname = layer_names_cem{li};
                    html = [html sprintf('<h4 style="margin:12px 0 4px;color:#475569;">%s</h4>', ...
                        html_utils.escape(strrep(lname, '_', ' ')))];
                    html = [html build_metric_uncertainty_svg(heating, 'CEM43', 0, 'CEM43 (min)', lname)];
                end
            end
        end
    end
end

function layer_names = get_ts_layer_names(heating, metric)
% Return layer names found in heating{v}.timeseries.<metric> for any variant v.
    layer_names = {};
    for v = 1:numel(heating)
        if ~isempty(heating{v}) && isfield(heating{v}, 'timeseries') && ...
                isstruct(heating{v}.timeseries) && isfield(heating{v}.timeseries, metric)
            layer_names = fieldnames(heating{v}.timeseries.(metric));
            return
        end
    end
end

function html = build_images_section(variants, variant_affixes, parameters, subject_id, medium, is_layered)
% Thermal maps of the Temp pair and the default (the Trans variants run no heating).
    html = '';
    trio = pair_trio(variants, 'temp');

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
        for v = trio
            aff = variant_affixes{v};
            lbl = variants(v).label;
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
    for v = 1:numel(tables)
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

function html = build_range_table(tables, variants, metric_names)
% Build a table with one row per metric, columns: metric | one per variant |
% range over the metric's pair (Trans for acoustic, Temp for thermal) | limit
    limits = struct();
    if is_layered_limits(metric_names)
        limits = get_risk_limits(true);
    end

    html = '<div class="table-wrapper"><table class="data-table"><thead><tr>';
    html = [html '<th>Metric</th>'];
    for v = 1:numel(variants)
        html = [html sprintf('<th>%s</th>', html_utils.escape(variants(v).label))];
    end
    html = [html '<th>Range</th><th>Limit</th></tr></thead><tbody>'];

    for m = 1:numel(metric_names)
        name = metric_names{m};

        vals = metric_values(tables, name);
        trio = pair_trio(variants, metric_pair(name));

        % Limit, and colour from the highest value of default and pair
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
        worst = max(vals(trio), [], 'omitnan');
        row_color = risk_color(worst, limit_val);

        html = [html sprintf('<tr class="row-%s">', row_color)];
        html = [html sprintf('<td><strong>%s</strong></td>', html_utils.escape(name))];
        for v = 1:numel(vals)
            if isnan(vals(v))
                html = [html '<td>N/A</td>'];
            else
                html = [html sprintf('<td>%.4g</td>', vals(v))];
            end
        end

        % Range column: default and the metric's pair
        valid = vals(trio);
        valid = valid(~isnan(valid));
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
%  PSEUDO-CT SKULL VS. BOUNDS
%% ========================================================================

function html = build_pct_bounds_section(parameters, variants, default_affix)
% Histograms of the default variant's pseudo-CT skull properties (the
% medium_stats cache of MEDIUM_SETUP, else PCT_PREDICT_SKULL_PROPERTIES) with
% the band of the bound variants' homogeneous skull values.
    stats_file = fullfile(parameters.io.dir_output, 'cache', sprintf('sub-%03d_%s_medium_stats%s.json', ...
        parameters.subject_id, parameters.simulation.medium, default_affix));
    if isfile(stats_file)
        stats = jsondecode(fileread(stats_file));
        origin = sprintf('mapped skull of the default run (%s)', html_utils.escape(stats_file));
    else
        stats = pct_predict_skull_properties(parameters, write=false);
        origin = sprintf('predicted from the pseudo-CT (%s voxels), no medium_stats cache found', stats.source);
    end

    bounds = struct();
    for v = variants(~strcmp({variants.name}, 'default'))
        if isfield(parameters, 'uncertainty') && isfield(parameters.uncertainty, 'bound_skull') && ...
                isfield(parameters.uncertainty.bound_skull, v.name)
            bounds.(v.name) = parameters.uncertainty.bound_skull.(v.name);
            continue
        end
        try
            bounds.(v.name) = uncertainty_medium_layer(parameters, v.name).medium_properties.skull;
        catch
        end
    end

    props = {'density', 'Density [kg/m&#179;]'; 'sound_speed', 'Sound speed [m/s]'; ...
             'alpha_coeff', 'Attenuation [dB/(cm&#183;MHz^y)]'};
    html = ['<p style="font-size:0.85em;color:#64748b;margin:0 0 12px;">The default variant maps the skull ' ...
        'from the pseudo-CT; the bound variants give the same skull voxels the homogeneous values of their ' ...
        'bounds. Histogram: ' origin '. Shaded band: range of the bound variants'' skull values. ' ...
        sprintf('Mapping: density %s, sound speed %s, attenuation %s.</p>', ...
        html_utils.escape(char(stats.mapping.density)), html_utils.escape(char(stats.mapping.sound_speed)), ...
        html_utils.escape(char(stats.mapping.attenuation)))];
    html = [html '<div class="image-grid image-grid--3col">'];
    for k = 1:size(props, 1)
        key = props{k, 1};
        vals = cellfun(@(n) bounds.(n).(key), fieldnames(bounds));
        if isfield(stats.properties, key)
            html = [html '<figure>' pct_histogram_svg(stats.properties.(key), vals, props{k, 2}) '</figure>']; %#ok<AGROW>
        end
    end
    html = [html '</div>'];
end

function svg = pct_histogram_svg(dist, band_vals, x_label)
% Inline SVG histogram with a shaded band [min, max] of band_vals and the
% share of voxels in bins whose centre lies outside it.
    W = 320; H = 200; pl = 10; pr = 10; pt = 22; pb = 34;
    counts = double(dist.histogram.counts(:)');
    edges  = double(dist.histogram.edges(:)');
    lo = min(band_vals); hi = max(band_vals);
    x0 = min([edges(1), lo]); x1 = max([edges(end), hi]);
    if x1 <= x0; x1 = x0 + 1; end
    pad = 0.04 * (x1 - x0); x0 = x0 - pad; x1 = x1 + pad;
    sx = @(x) pl + (x - x0) / (x1 - x0) * (W - pl - pr);
    cmax = max([counts, 1]);
    sy = @(c) H - pb - c / cmax * (H - pt - pb);
    svg = sprintf(['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" ' ...
        'style="width:100%%;height:auto;background:white;border-radius:4px;">'], W, H);
    if ~isempty(band_vals)
        svg = [svg sprintf('<rect x="%.1f" y="%d" width="%.1f" height="%d" fill="#d98a1f" fill-opacity="0.22"/>', ...
            sx(lo), pt, max(sx(hi) - sx(lo), 1.5), H - pt - pb)];
    end
    for i = 1:numel(counts)
        svg = [svg sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#8d9db6"/>', ...
            sx(edges(i)), sy(counts(i)), max(sx(edges(i+1)) - sx(edges(i)), 0.5), H - pb - sy(counts(i)))]; %#ok<AGROW>
    end
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#cbd5e1"/>', pl, H - pb, W - pr, H - pb)];
    for t = [x0 + pad, (x0 + x1) / 2, x1 - pad]
        svg = [svg sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="9" fill="#64748b">%.3g</text>', ...
            sx(t), H - pb + 12, t)]; %#ok<AGROW>
    end
    svg = [svg sprintf('<text x="%d" y="%d" text-anchor="middle" font-size="10" fill="#475569">%s</text>', ...
        W / 2, H - 4, x_label)];
    if sum(counts) > 0 && ~isempty(band_vals)
        centres = (edges(1:end-1) + edges(2:end)) / 2;
        outside = sum(counts(centres < lo | centres > hi)) / sum(counts);
        svg = [svg sprintf('<text x="%d" y="%d" text-anchor="end" font-size="10" fill="#1e293b">%.0f%% outside</text>', ...
            W - pr, 14, 100 * outside)];
    end
    svg = [svg '</svg>'];
end

%% ========================================================================
%  UNCERTAINTY TIMESERIES SVG PLOT
%% ========================================================================

function html = build_metric_uncertainty_svg(heating, metric, y_floor, y_label, layer)
% Build an inline SVG showing a timeseries metric with uncertainty band across
% three variants (high=1, default=2, low=3; the Temp pair).
%
% metric  - field name in heating{v}.timeseries, e.g. 'T' or 'CEM43'
% y_floor - lower y-axis bound (37 for temperature, 0 for CEM43)
% y_label - y-axis label string (may include HTML entities)
% layer   - char layer name to plot (e.g. 'brain'), or false for global max

    html = '';
    ts_data = cell(1, 3);
    n_pts = 0;
    y_max_all = y_floor;

    for v = 1:3
        ts_line = [];
        if ~isempty(heating{v}) && isfield(heating{v}, 'timeseries')
            ts = heating{v}.timeseries;
            if isstruct(ts) && isfield(ts, metric)
                m = ts.(metric);
                if ischar(layer) || isstring(layer)
                    % Single named layer
                    if isfield(m, layer)
                        ts_line = m.(layer)(:);
                    end
                else
                    % Global maximum across all layers (layer == false)
                    if isstruct(m)
                        fns = fieldnames(m);
                        for li = 1:numel(fns)
                            ldata = m.(fns{li})(:);
                            if isempty(ts_line)
                                ts_line = ldata;
                            else
                                n_c = min(numel(ts_line), numel(ldata));
                                ts_line(1:n_c) = max(ts_line(1:n_c), ldata(1:n_c));
                            end
                        end
                    elseif isnumeric(m)
                        ts_line = max(m, [], 2);
                    end
                end
            end
        end
        if ~isempty(ts_line)
            finite_vals = ts_line(isfinite(ts_line));
            if ~isempty(finite_vals)
                y_max_all = max(y_max_all, max(finite_vals));
            end
            ts_data{v} = ts_line;
            n_pts = max(n_pts, numel(ts_line));
        end
    end

    if n_pts == 0
        html = '<p class="placeholder">No timeseries data available.</p>';
        return;
    end

    for v = 1:3
        if isempty(ts_data{v})
            ts_data{v} = NaN(n_pts, 1);
        elseif numel(ts_data{v}) < n_pts
            ts_data{v}(end+1:n_pts) = NaN;
        end
    end

    W = 700; H = 300;
    pad_l = 62; pad_r = 20; pad_t = 20; pad_b = 50;
    plot_w = W - pad_l - pad_r;
    plot_h = H - pad_t - pad_b;

    y_min   = y_floor;
    y_range = max(y_max_all * 1.05, y_min + 0.01) - y_min;
    x_scale = plot_w / max(n_pts - 1, 1);
    y_scale = plot_h / y_range;

    % Band polygon: high (index 1) forward then low (index 3) backward
    band_pts = '';
    for i = 1:n_pts
        x   = pad_l + (i-1) * x_scale;
        val = ts_data{1}(i);
        if isnan(val); val = y_min; end
        y   = pad_t + plot_h - (val - y_min) * y_scale;
        band_pts = [band_pts sprintf('%.1f,%.1f ', x, y)];
    end
    for i = n_pts:-1:1
        x   = pad_l + (i-1) * x_scale;
        val = ts_data{3}(i);
        if isnan(val); val = y_min; end
        y   = pad_t + plot_h - (val - y_min) * y_scale;
        band_pts = [band_pts sprintf('%.1f,%.1f ', x, y)];
    end

    % Helper: build polyline point string for one variant
    function pts = variant_pts(v_idx)
        pts = '';
        for ii = 1:n_pts
            xx  = pad_l + (ii-1) * x_scale;
            vv  = ts_data{v_idx}(ii);
            if isnan(vv); continue; end
            yy  = pad_t + plot_h - (vv - y_min) * y_scale;
            pts = [pts sprintf('%.1f,%.1f ', xx, yy)];
        end
    end

    % Y-axis ticks
    n_ticks  = 5;
    tick_svg = '';
    for k = 0:n_ticks
        val = y_min + y_range * k / n_ticks;
        yp  = pad_t + plot_h - (val - y_min) * y_scale;
        tick_svg = [tick_svg ...
            sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e2e8f0" stroke-width="1"/>', ...
                pad_l, yp, W-pad_r, yp) ...
            sprintf('<text x="%d" y="%.1f" text-anchor="end" font-size="10" fill="#94a3b8">%.3g</text>', ...
                pad_l-4, yp+4, val)];
    end

    % X-axis ticks
    xtick_svg = '';
    for k = 0:4
        idx = round(k * (n_pts-1) / 4) + 1;
        xp  = pad_l + (idx-1) * x_scale;
        xtick_svg = [xtick_svg ...
            sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="10" fill="#94a3b8">%d</text>', ...
                xp, H-pad_b+14, idx)];
    end

    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' num2str(W) ' ' num2str(H) '" ' ...
           'style="width:100%;max-width:' num2str(W) 'px;height:auto;display:block;margin:0 auto 12px;">'];
    svg = [svg sprintf('<rect width="%d" height="%d" fill="white" rx="4" ry="4"/>', W, H)];
    svg = [svg tick_svg xtick_svg];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#cbd5e1" stroke-width="1.5"/>', ...
        pad_l, pad_t, pad_l, pad_t+plot_h)];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#cbd5e1" stroke-width="1.5"/>', ...
        pad_l, pad_t+plot_h, pad_l+plot_w, pad_t+plot_h)];
    svg = [svg sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="11" fill="#64748b">Time step</text>', ...
        pad_l + plot_w/2, H-4)];
    svg = [svg sprintf('<text x="%d" y="%.1f" text-anchor="middle" font-size="11" fill="#64748b" ', ...
        12, pad_t + plot_h/2) ...
        'transform="rotate(-90,' num2str(12) ',' num2str(pad_t+plot_h/2) ')">' y_label '</text>'];

    % Uncertainty band
    if ~isempty(strtrim(band_pts))
        svg = [svg sprintf('<polygon points="%s" fill="#3b82f6" fill-opacity="0.15" stroke="none"/>', band_pts)];
    end
    % High (dashed red)
    lib_pts = variant_pts(1);
    if ~isempty(strtrim(lib_pts))
        svg = [svg sprintf('<polyline points="%s" fill="none" stroke="#ef4444" stroke-width="1.5" stroke-dasharray="4,3"/>', lib_pts)];
    end
    % Low (dashed gray)
    cons_pts = variant_pts(3);
    if ~isempty(strtrim(cons_pts))
        svg = [svg sprintf('<polyline points="%s" fill="none" stroke="#94a3b8" stroke-width="1.5" stroke-dasharray="4,3"/>', cons_pts)];
    end
    % Default (solid blue, drawn last so it sits on top)
    def_pts = variant_pts(2);
    if ~isempty(strtrim(def_pts))
        svg = [svg sprintf('<polyline points="%s" fill="none" stroke="#1d4ed8" stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>', ...
            def_pts)];
    end

    % Legend
    lx = pad_l + plot_w - 150; ly = pad_t + 10;
    svg = [svg sprintf('<rect x="%d" y="%d" width="150" height="64" rx="3" fill="white" fill-opacity="0.85" stroke="#e2e8f0"/>', lx, ly)];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#1d4ed8" stroke-width="2"/><text x="%d" y="%d" font-size="10" fill="#475569">Default</text>', ...
        lx+6, ly+16, lx+26, ly+16, lx+30, ly+20)];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#ef4444" stroke-width="1.5" stroke-dasharray="4,3"/><text x="%d" y="%d" font-size="10" fill="#475569">Temp-H</text>', ...
        lx+6, ly+32, lx+26, ly+32, lx+30, ly+36)];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#94a3b8" stroke-width="1.5" stroke-dasharray="4,3"/><text x="%d" y="%d" font-size="10" fill="#475569">Temp-L</text>', ...
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

function html = build_timing_section(parameters, variants)
% Build an HTML table showing per-step wall times across all pipeline
% stages, using the pre-assigned log file paths in
% parameters.uncertainty.log_files.
% Returns empty string if no timing data is found.

    if ~isfield(parameters, 'uncertainty') || ~isfield(parameters.uncertainty, 'log_files')
        html = '';
        return;
    end
    lf = parameters.uncertainty.log_files;  % struct: stage1, report and one field per variant name

    % Stage 1: preprocessing steps (preproc, medium, source)
    t_stage1 = parse_timing_from_log(char(lf.stage1));

    % Stages 2–6: simulation variants (acoustic, thermal, total pipeline)
    t_variants = cell(1, numel(variants));
    for v = 1:numel(variants)
        if isfield(lf, variants(v).name)
            t_variants{v} = parse_timing_from_log(char(lf.(variants(v).name)));
        end
    end

    % Stage 7: report
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

    % Column headers: Stage 1 | one per variant | Report
    col_labels = [{'Stage 1'}, {variants.label}, {'Report'}];
    pair_colors = struct('temp', '#b91c1c', 'trans', '#1d4ed8');
    col_colors = cell(1, numel(variants));
    for v = 1:numel(variants)
        if isempty(variants(v).pair)
            col_colors{v} = '#1e293b';
        else
            col_colors{v} = pair_colors.(variants(v).pair);
        end
    end
    col_colors = [{'#475569'}, col_colors, {'#475569'}];

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

        % Variant columns
        for v = 1:numel(variants)
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
'.safety-bar-marker--high         { background: #b91c1c; }' newline ...
'.safety-bar-marker--low          { background: #2563eb; }' newline ...
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
