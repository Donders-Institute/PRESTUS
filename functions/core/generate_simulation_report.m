function report_path = generate_simulation_report(parameters)
% GENERATE_SIMULATION_REPORT creates a self-contained HTML report consolidating
% all PRESTUS simulation outputs (CSV data, images, logs) into a single portable file.
%
% Features a safety dashboard at the top with color-coded metrics against
% ITRUSST consensus limits (Aubry et al., 2025).
%
% Use as:
%   report_path = generate_simulation_report(parameters)
%
% Inputs:
%   parameters       - PRESTUS parameters struct (required)
%
% Output:
%   report_path      - Path to the generated HTML file

arguments
    parameters      struct
end

report_path = '';

try
    %% Determine output path
    subject_id = parameters.subject_id;
    medium = parameters.simulation.medium;
    is_layered = contains(medium, {'layered', 'phantom'});
    affix = '';
    if isfield(parameters.io, 'output_affix')
        affix = parameters.io.output_affix;
    end

    report_filename = sprintf('sub-%03d_%s_report%s.html', subject_id, medium, affix);
    report_path = fullfile(parameters.io.output_dir, report_filename);

    %% Load CSV data if available
    csv_table = [];
    if isfield(parameters.io, 'filename_output_table') && isfile(parameters.io.filename_output_table)
        try
            csv_table = readtable(parameters.io.filename_output_table, 'VariableNamingRule', 'preserve');
        catch
            csv_table = [];
        end
    end

    %% Build HTML sections
    html_parts = {};

    html_parts{end+1} = '<!DOCTYPE html>';
    html_parts{end+1} = '<html lang="en">';
    html_parts{end+1} = '<head>';
    html_parts{end+1} = '<meta charset="UTF-8">';
    html_parts{end+1} = '<meta name="viewport" content="width=device-width, initial-scale=1.0">';
    html_parts{end+1} = sprintf('<title>PRESTUS Report — sub-%03d %s%s</title>', subject_id, medium, affix);
    html_parts{end+1} = '<style>';
    html_parts{end+1} = css_styles();
    html_parts{end+1} = '</style>';
    html_parts{end+1} = '</head>';
    html_parts{end+1} = '<body>';

    % Table of contents (sticky nav bar)
    try
        html_parts{end+1} = build_toc(parameters, is_layered);
    catch ME
        html_parts{end+1} = section_error('Table of Contents', ME);
    end

    % Section 1: Safety Dashboard (always open, no toggle)
    try
        html_parts{end+1} = build_safety_dashboard(csv_table, parameters, is_layered);
    catch ME
        html_parts{end+1} = section_error('Safety Dashboard', ME);
    end

    % Section 2: Header (always open, no toggle)
    try
        html_parts{end+1} = build_header(subject_id, medium, affix, parameters, is_layered);
    catch ME
        html_parts{end+1} = section_error('Header', ME);
    end

    % Section 3: Simulation Summary (always open, no toggle)
    try
        html_parts{end+1} = build_simulation_summary(csv_table, parameters, is_layered);
    catch ME
        html_parts{end+1} = section_error('Simulation Summary', ME);
    end

    % Section 4: Configuration Summary (collapsed by default)
    try
        html_parts{end+1} = collapsible_section('Configuration Summary', ...
            build_config_summary(parameters), false, 'config');
    catch ME
        html_parts{end+1} = section_error('Configuration Summary', ME);
    end

    % Section 5: Medium Properties (layered only, collapsed by default)
    if is_layered
        try
            html_parts{end+1} = collapsible_section('Medium Properties', ...
                build_medium_properties_section(parameters), false, 'medium');
        catch ME
            html_parts{end+1} = section_error('Medium Properties', ME);
        end
    end

    % Section 6: pseudoCT [Optional] (collapsed by default)
    try
        if isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && parameters.pct.enabled
            html_parts{end+1} = collapsible_section('pseudo-CT', ...
                build_pseudoCT_section(parameters), ...
                false, 'pseudoCT');
        end
    catch ME
        html_parts{end+1} = section_error('pseudo-CT', ME);
    end

    % Section 7: Positioning (open by default)
    try
        html_parts{end+1} = collapsible_section('Transducer Positioning', ...
            build_positioning_section(parameters, subject_id, affix), true, 'positioning');
    catch ME
        html_parts{end+1} = section_error('Positioning', ME);
    end

    % Section 8: Acoustic Results (open by default)
    try
        html_parts{end+1} = collapsible_section('Acoustic Results', ...
            build_acoustic_section(csv_table, parameters, subject_id, medium, affix, is_layered), true, 'acoustic');
    catch ME
        html_parts{end+1} = section_error('Acoustic Results', ME);
    end

    % Section 9: Thermal Results (open by default, conditional)
    try
        if isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims
            html_parts{end+1} = collapsible_section('Thermal Results', ...
                build_thermal_section(csv_table, parameters, subject_id, medium, affix, is_layered), true, 'thermal');
        end
    catch ME
        html_parts{end+1} = section_error('Thermal Results', ME);
    end

    % Section 10: Debug Information (collapsed by default, layered only)
    try
        if isfield(parameters.simulation, 'debug') && parameters.simulation.debug && is_layered
            html_parts{end+1} = collapsible_section('Debug Information', ...
                build_debug_section(parameters, subject_id, medium, affix), false, 'debug');
        end
    catch ME
        html_parts{end+1} = section_error('Debug Information', ME);
    end

    % Section 11: Post-hoc Water Simulation (collapsed by default, layered only)
    if is_layered
        try
            html_parts{end+1} = collapsible_section('Post-hoc Water Simulation', ...
                build_posthoc_section(parameters), false, 'posthoc');
        catch ME
            html_parts{end+1} = section_error('Post-hoc Water Simulation', ME);
        end
    end

    % Section 12: Performance (collapsed by default)
    try
        html_parts{end+1} = collapsible_section('Performance', ...
            build_performance_section(parameters, subject_id, medium, affix), false, 'performance');
    catch ME
        html_parts{end+1} = section_error('Performance', ME);
    end

    % Section 13: Log (collapsed by default)
    try
        html_parts{end+1} = collapsible_section('Log', ...
            build_log_section(parameters, subject_id, medium, affix), false, 'log');
    catch ME
        html_parts{end+1} = section_error('Log', ME);
    end

    % Footer
    html_parts{end+1} = '<footer>';
    html_parts{end+1} = sprintf('<p>Generated by PRESTUS %s</p>', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    html_parts{end+1} = '</footer>';
    % Lightbox overlay and script
    html_parts{end+1} = lightbox_html();

    html_parts{end+1} = '</body>';
    html_parts{end+1} = '</html>';

    %% Write HTML file
    fid = fopen(report_path, 'w', 'n', 'UTF-8');
    if fid == -1
        warning('generate_simulation_report:fileOpen', 'Cannot open %s for writing.', report_path);
        return
    end
    fprintf(fid, '%s\n', html_parts{:});
    fclose(fid);

    fprintf('HTML report saved to: %s\n', report_path);

catch ME
    warning('generate_simulation_report:failed', ...
        'Report generation failed: %s\n%s', ME.message, getReport(ME, 'extended'));
    report_path = '';
end

end

%% ========================================================================
%  SAFETY LIMITS & DASHBOARD
%  ========================================================================

function limits = get_safety_limits()
% Returns struct with ITRUSST consensus thresholds (Aubry et al., 2025).
% Green: < 50% of limit. Amber: 50-100% of limit. Red: exceeds limit.

    limits = struct();

    % MI limits (ITRUSST: MI <= 1.9)
    limits.MI_tc    = struct('label', 'MI (transcranial)', 'limit', 1.9, 'unit', '');
    limits.MI_brain = struct('label', 'MI (brain)',  'limit', 1.9, 'unit', '');
    limits.MI_skull = struct('label', 'MI (skull)',  'limit', 1.9, 'unit', '');
    limits.MI_skin  = struct('label', 'MI (skin)',   'limit', 1.9, 'unit', '');

    % Temperature rise limits (ITRUSST: <= 2 C rise)
    limits.riseT_brain = struct('label', 'Temp rise (brain)', 'limit', 2.0, 'unit', [char(176) 'C']);
    limits.riseT_skull = struct('label', 'Temp rise (skull)', 'limit', 2.0, 'unit', [char(176) 'C']);
    limits.riseT_skin  = struct('label', 'Temp rise (skin)',  'limit', 2.0, 'unit', [char(176) 'C']);

    % CEM43 limits (ITRUSST: brain <= 2, skull <= 16, skin <= 21)
    limits.CEM43_brain = struct('label', 'CEM43 (brain)', 'limit', 2.0,  'unit', 'min');
    limits.CEM43_skull = struct('label', 'CEM43 (skull)', 'limit', 16.0, 'unit', 'min');
    limits.CEM43_skin  = struct('label', 'CEM43 (skin)',  'limit', 21.0, 'unit', 'min');

    % Absolute temperature limits (ITRUSST: <= 39 C)
    limits.maxT_brain = struct('label', 'Max temp (brain)', 'limit', 39.0, 'unit', [char(176) 'C']);
    limits.maxT_skull = struct('label', 'Max temp (skull)', 'limit', 39.0, 'unit', [char(176) 'C']);
    limits.maxT_skin  = struct('label', 'Max temp (skin)',  'limit', 39.0, 'unit', [char(176) 'C']);

    % ISPPA (informational, no ITRUSST limit)
    limits.Isppa_brain = struct('label', 'ISPPA (brain)', 'limit', Inf, 'unit', 'W/cm²');
    limits.Isppa_skull = struct('label', 'ISPPA (skull)', 'limit', Inf, 'unit', 'W/cm²');
    limits.Isppa_skin  = struct('label', 'ISPPA (skin)',  'limit', Inf, 'unit', 'W/cm²');
end

function limits = get_safety_limits_water()
% Returns global-only safety limits for water/free-field simulations.
% Tissue-specific limits are not applicable.
    limits = struct();
    limits.Isppa  = struct('label', 'ISPPA (global)', 'limit', Inf, 'unit', 'W/cm²');
    limits.Psptp = struct('label', 'Psptp',         'limit', Inf, 'unit', 'Pa');
end

function color = safety_color(value, limit)
% Returns 'green', 'amber', or 'red' based on value vs ITRUSST limit.
% Green: < 50% of limit. Amber: 50-100%. Red: exceeds limit.
    if isnan(value) || isempty(value)
        color = 'gray';
    elseif isinf(limit)
        color = 'info'; % informational, no limit
    elseif value > limit
        color = 'red';
    elseif value >= 0.5 * limit
        color = 'amber';
    else
        color = 'green';
    end
end

function html = build_safety_dashboard(csv_table, parameters, is_layered)
    if is_layered
        limits = get_safety_limits();
    else
        limits = get_safety_limits_water();
    end
    metric_names = fieldnames(limits);

    html = '<section class="report-section" id="safety">';
    html = [html '<h2>Safety Dashboard</h2>'];

    if ~is_layered
        html = [html '<div class="medium-banner medium-water">' ...
            '<strong>Water / Free-field</strong> &mdash; tissue-specific safety limits not applicable. ' ...
            'Only global metrics are shown.</div>'];
    end

    html = [html '<div class="safety-grid">'];

    for i = 1:length(metric_names)
        name = metric_names{i};
        info = limits.(name);

        % Extract value from CSV
        value = NaN;
        display_unit = info.unit;
        if ~isempty(csv_table) && ismember(name, csv_table.Properties.VariableNames)
            val = csv_table.(name);
            if isnumeric(val) && ~isempty(val)
                value = val(end); % last row if multiple
            end
        end

        % Special handling for pressure: use dynamic unit scaling
        if strcmp(name, 'Psptp') && ~isnan(value)
            [value, display_unit] = scale_pressure(value);
        end

        % Determine color
        color = safety_color(value, info.limit);

        % Build card
        html = [html sprintf('<div class="safety-card safety-%s">', color)];
        html = [html sprintf('<div class="safety-label">%s</div>', html_escape(info.label))];

        if isnan(value)
            html = [html '<div class="safety-value">N/A</div>'];
            html = [html '<div class="safety-limit">No data</div>'];
        else
            html = [html sprintf('<div class="safety-value">%.3g</div>', value)];
            if isinf(info.limit)
                html = [html sprintf('<div class="safety-limit">%s (informational)</div>', display_unit)];
            else
                html = [html sprintf('<div class="safety-limit">Limit: %.3g %s</div>', info.limit, display_unit)];
                % Progress bar showing value as % of limit
                pct = min(100, max(0, (value / info.limit) * 100));
                html = [html sprintf('<div class="safety-bar-track"><div class="safety-bar" style="width:%.0f%%"></div></div>', pct)];
            end
        end

        html = [html '</div>'];
    end

    html = [html '</div>'];
    html = [html '<p class="safety-footnote">Limits based on ITRUSST consensus ' ...
            '(Aubry et al., 2025; ' ...
            '<a href="https://doi.org/10.1016/j.brs.2025.10.007">DOI: 10.1016/j.brs.2025.10.007</a>). ' ...
            'Green: &lt; 50% of limit. ' ...
            'Amber: 50&ndash;100% of limit. Red: exceeds limit.</p>'];
    html = [html '</section>'];
end

%% ========================================================================
%  SECTION BUILDERS
%  ========================================================================

function html = build_header(subject_id, medium, affix, parameters, is_layered)
    html = '<section class="report-section" id="header">';
    if is_layered
        html = [html '<h1>PRESTUS Simulation Report</h1>'];
    else
        html = [html '<h1>PRESTUS Free-water Report</h1>'];
    end
    html = [html '<table class="info-table">'];
    html = [html sprintf('<tr><th>Subject</th><td>sub-%03d</td></tr>', subject_id)];
    if is_layered && isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && parameters.pct.enabled
        html = [html sprintf('<tr><th>Medium</th><td>%s <span class="badge badge-amber">Layered (Realistic Head)</span> <span class="badge badge-gray">pseudo-CT (Continuous Skull)</span></td></tr>', html_escape(medium))];
    elseif is_layered
        html = [html sprintf('<tr><th>Medium</th><td>%s <span class="badge badge-amber">Layered (Realistic Head)</span></td></tr>', html_escape(medium))];
    else
        html = [html sprintf('<tr><th>Medium</th><td>%s <span class="badge badge-blue">Water / Free-field</span></td></tr>', html_escape(medium))];
    end
    if ~isempty(affix)
        html = [html sprintf('<tr><th>Affix</th><td>%s</td></tr>', html_escape(affix))];
    end
    html = [html sprintf('<tr><th>Generated</th><td>%s</td></tr>', datestr(now, 'yyyy-mm-dd HH:MM:SS'))];
    if isfield(parameters, 'output_dir')
        html = [html sprintf('<tr><th>Output dir</th><td>%s</td></tr>', html_escape(parameters.io.output_dir))];
    end
    html = [html '</table>'];
    html = [html '</section>'];
end

function html = build_simulation_summary(csv_table, parameters, is_layered)
% Quick-glance summary cards showing key simulation outcomes.
    html = '<section class="report-section" id="summary">';
    html = [html '<h2>Simulation Summary</h2>'];
    html = [html '<div class="summary-grid">'];

    % Always show: IPA at target, focal distance, max ISPPA
    html = [html summary_card('Ipa at target', csv_value(csv_table, 'Ipa_target'), 'W/cm²')];
    html = [html summary_card('Focal distance', csv_value(csv_table, 'real_focal_distance_mm'), 'mm')];
    html = [html summary_card('ISPPA', csv_value(csv_table, 'Isppa'), 'W/cm²')];

    if is_layered
        html = [html summary_card('ISPPA brain', csv_value(csv_table, 'Isppa_brain'), 'W/cm²')];
        html = [html summary_card('-6dB vol. brain', csv_value(csv_table, 'halfmax_ISPPA_volume_brain_mm3'), 'mm³')];
    else
        [p_val, p_unit] = scale_pressure(csv_value(csv_table, 'Psptp'));
        html = [html summary_card('Psptp', p_val, p_unit)];
    end

    % Thermal summary if available
    if is_layered && isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims
        html = [html summary_card('Max temp.', csv_value(csv_table, 'maxT'), [char(176) 'C'])];
    end

    html = [html '</div>'];
    html = [html '</section>'];
end

function html = summary_card(label, value, unit)
% Single summary metric card.
    if isnan(value)
        val_str = 'N/A';
    else
        val_str = sprintf('%.4g', value);
    end
    html = '<div class="summary-card">';
    html = [html sprintf('<div class="summary-label">%s</div>', html_escape(label))];
    html = [html sprintf('<div class="summary-value">%s</div>', val_str)];
    html = [html sprintf('<div class="summary-unit">%s</div>', html_escape(unit))];
    html = [html '</div>'];
end

function [val, unit] = scale_pressure(val_Pa)
% Auto-scale pressure value to Pa, kPa, or MPa for display.
    if isnan(val_Pa)
        val  = NaN;
        unit = 'Pa';
    elseif abs(val_Pa) >= 1e6
        val  = val_Pa / 1e6;
        unit = 'MPa';
    elseif abs(val_Pa) >= 1e3
        val  = val_Pa / 1e3;
        unit = 'kPa';
    else
        val  = val_Pa;
        unit = 'Pa';
    end
end

function val = csv_value(csv_table, col_name)
% Extract last-row numeric value from CSV table, or NaN if missing.
    val = NaN;
    if ~isempty(csv_table) && ismember(col_name, csv_table.Properties.VariableNames)
        v = csv_table.(col_name);
        if isnumeric(v) && ~isempty(v)
            val = v(end);
        end
    end
end

function html = build_medium_properties_section(parameters)
    html = '';
    if ~isfield(parameters, 'medium_properties') || ~isfield(parameters, 'layers')
        html = '<p class="placeholder">Missing medium or layers data.</p>';
        return
    end

    med = parameters.medium_properties;
    layer_names = fieldnames(parameters.layers);
    
    % === CORRECT FILTER: layers in BOTH parameters.layers AND parameters.simulation.medium ===
    valid_tissues = {};
    tissue_labels = {};
    for i = 1:length(layer_names)
        layer_name = layer_names{i};
        if isfield(med, layer_name)  % Only include if medium has properties for this layer
            valid_tissues{end+1} = layer_name;
            
            % Human-readable labels
            switch layer_name
                case 'skull_cortical'
                    tissue_labels{end+1} = 'Skull (cortical)';
                case 'skull_trabecular' 
                    tissue_labels{end+1} = 'Skull (trabecular)';
                otherwise
                    tissue_labels{end+1} = strrep(layer_name, '_', ' ');
            end
        end
    end
    
    if isempty(valid_tissues)
        html = '<p class="placeholder">No matching layers found in medium properties.</p>';
        return
    end

    % Rest unchanged...
    props = {'sound_speed', 'Sound speed', 'm/s'; ...
             'density', 'Density', 'kg/m³'; ...
             'alpha_coeff', 'Attenuation coeff.', 'dB/(MHz^y·cm)'; ...
             'alpha_power', 'Attenuation power', ''; ...
             'thermal_conductivity', 'Thermal cond.', 'W/(m·K)'; ...
             'specific_heat_capacity', 'Specific heat', 'J/(kg·K)'; ...
             'perfusion', 'Perfusion rate', 'kg/(m³·s)'; ...
             'absorption_fraction', 'Absorption fraction', 'Fraction'};

    html = '<div class="table-wrapper"><table class="data-table">';
    html = [html '<thead><tr><th>Property</th><th>Unit</th>'];
    for t = 1:length(valid_tissues)
        html = [html sprintf('<th>%s</th>', html_escape(tissue_labels{t}))];
    end
    html = [html '</tr></thead><tbody>'];

    for p = 1:size(props, 1)
        field = props{p, 1}; label = props{p, 2}; unit = props{p, 3};
        html = [html sprintf('<tr><td><strong>%s</strong></td><td>%s</td>', ...
            html_escape(label), html_escape(unit))];
        for t = 1:length(valid_tissues)
            tissue = valid_tissues{t};
            if isfield(med.(tissue), field) && isnumeric(med.(tissue).(field)) && isscalar(med.(tissue).(field))
                html = [html sprintf('<td>%.4g</td>', med.(tissue).(field))];
            else
                html = [html '<td>&mdash;</td>'];
            end
        end
        html = [html '</tr>'];
    end
    html = [html '</tbody></table></div>'];
    
    html = [html sprintf('<p class="note"><strong>Active layers:</strong> %s</p>', ...
        strjoin(tissue_labels, ', '))];
end

function html = build_config_summary(parameters)
    html = '<table class="config-table">';

    % Simulation type
    html = config_row(html, 'Simulation medium', safe_field(parameters.simulation, 'medium', 'N/A'));
    html = config_row(html, 'Dimensions', sprintf('%dD', safe_field(parameters.grid, 'n_dims', 3)));

    % Grid
    grid_step = safe_field(parameters.grid, 'resolution_mm', NaN);
    if ~isnan(grid_step)
        html = config_row(html, 'Grid step', sprintf('%.2f mm', grid_step));
    end
    grid_dims = safe_field(parameters.grid, 'default_dims', []);
    if ~isempty(grid_dims)
        html = config_row(html, 'Grid dims', sprintf('[%s]', strtrim(sprintf('%g ', grid_dims))));
    end
    pml = safe_field(parameters.grid, 'pml_size', NaN);
    if ~isnan(pml)
        html = config_row(html, 'PML size', sprintf('%d', pml));
    end

    % Transducer (first transducer)
    if isfield(parameters, 'transducer') && ~isempty(parameters.transducer)
        td = parameters.transducer(1);
        html = config_row(html, 'Frequency', sprintf('%.0f Hz', safe_field(td, 'freq_hz', NaN)));
        html = config_row(html, 'Elements', sprintf('%d', safe_field(td.(td.type), 'elem_n', NaN)));
        html = config_row(html, 'Curvature radius', sprintf('%.0f mm', safe_field(td.(td.type), 'curv_radius_mm', NaN)));
        html = config_row(html, 'Source amplitude', sprintf('%.1f Pa', safe_field(td.(td.type), 'elem_amp', NaN)));
    end
    focal = safe_field(parameters.transducer(1), 'focal_distance_bowl', NaN);
    if ~isnan(focal)
        html = config_row(html, 'Expected focal distance', sprintf('%.1f mm', focal));
    end

    % Modules
    modules = {'run_source_setup', 'Source setup'; ...
               'run_acoustic_sims', 'Acoustic sims'; ...
               'run_heating_sims', 'Thermal sims'; ...
               'run_posthoc_water_sims', 'Post-hoc water'};
    for i = 1:size(modules, 1)
        val = safe_field(parameters.modules, modules{i,1}, 0);
        if val
            status = 'Enabled';
        else
            status = 'Disabled';
        end
        html = config_row(html, modules{i,2}, status);
    end

    % Thermal protocol (if enabled)
    if isfield(parameters, 'thermal') && isfield(parameters, 'modules') && safe_field(parameters.modules, 'run_heating_sims', 0)
        th = parameters.thermal;
        thermal_fields = {'pd', 'Pulse duration'; 'pri', 'Pulse repetition interval'; ...
                          'ptd', 'Pulse train duration'; 'ptri', 'Pulse train rep. interval'; ...
                          'ptrd', 'Pulse train rep. duration'; 'post_ptri_dur', 'Steady-state duration'};
        for i = 1:size(thermal_fields, 1)
            val = safe_field(th, thermal_fields{i,1}, NaN);
            if ~isnan(val)
                html = config_row(html, thermal_fields{i,2}, sprintf('%.3f s', val));
            end
        end
    end

    html = [html '</table>'];
end

function html = build_positioning_section(parameters, subject_id, affix)
    html = '';

    n_trans = 1;
    if isfield(parameters, 'transducer')
        n_trans = numel(parameters.transducer);
    end

    found_any = false;
    html = [html '<div class="image-grid image-grid--wide">'];
    for t = 1:n_trans
        img_path = fullfile(parameters.io.output_dir, ...
            sprintf('sub-%03d_positioning_T%02d%s.png', subject_id, t, affix));
        img_html = embed_image(img_path, sprintf('Positioning T%02d', t), ...
            sprintf('Transducer %d positioning', t));
        if ~isempty(img_html)
            html = [html img_html];
            found_any = true;
        end
    end
    html = [html '</div>'];

    if ~found_any
        html = '<p class="placeholder">No positioning images found.</p>';
    end
end

function html = build_acoustic_section(csv_table, parameters, subject_id, medium, affix, is_layered)
    html = '';

    % CSV table with color coding
    if ~isempty(csv_table)
        if is_layered
            limits = get_safety_limits();
            acoustic_cols = get_acoustic_columns();
            html = [html table2html(csv_table, limits, acoustic_cols)];
        else
            % Water: filter to water-relevant columns only, no safety color coding
            water_cols = get_acoustic_columns_water();
            avail_cols = intersect(water_cols, csv_table.Properties.VariableNames, 'stable');
            if ~isempty(avail_cols)
                sub_table = csv_table(:, avail_cols);
                % Scale pressure values dynamically (Pa -> kPa/MPa based on magnitude)
                if ismember('Psptp', sub_table.Properties.VariableNames)
                    pressure_vals = sub_table{:, 'Psptp'};
                    if ~all(isnan(pressure_vals))
                        [scaled_vals, display_unit] = scale_pressure(pressure_vals);
                        sub_table{:, 'Psptp'} = scaled_vals;
                        % Rename column header to reflect actual unit
                        sub_table.Properties.VariableNames{'Psptp'} = ['Psptp_' display_unit];
                    end
                end
                html = [html table2html(sub_table, struct(), {})];
            else
                html = [html '<p class="placeholder">No acoustic columns found in CSV.</p>'];
            end
        end
    else
        html = [html '<p class="placeholder">No acoustic CSV data found.</p>'];
    end

    % Intensity images
    n_trans = 1;
    if isfield(parameters, 'transducer')
        n_trans = numel(parameters.transducer);
    end

    html = [html '<div class="image-grid">'];
    for t = 1:n_trans
        trans_suffix = '';
        if n_trans > 1; trans_suffix = sprintf('_T%02d', t); end

        % Intensity on segmentation
        img_path = fullfile(parameters.io.output_dir, ...
            sprintf('sub-%03d_%s_intensity%s%s.png', subject_id, medium, trans_suffix, affix));
        img_html = embed_image(img_path, sprintf('Intensity on segmentation%s', trans_suffix), ...
            sprintf('Intensity overlay (segmentation)%s', trans_suffix));
        if ~isempty(img_html)
            html = [html img_html];
        end

        % Intensity on T1
        img_path = fullfile(parameters.io.output_dir, ...
            sprintf('sub-%03d_%s_intensity_t1%s%s.png', subject_id, medium, trans_suffix, affix));
        img_html = embed_image(img_path, sprintf('Intensity on T1%s', trans_suffix), ...
            sprintf('Intensity overlay (T1)%s', trans_suffix));
        if ~isempty(img_html)
            html = [html img_html];
        end
    end
    html = [html '</div>'];
end

function html = build_thermal_section(csv_table, parameters, subject_id, medium, affix, is_layered)
    html = '';

    % CSV thermal columns with color coding
    if ~isempty(csv_table)
        if is_layered
            limits = get_safety_limits();
            thermal_cols = get_thermal_columns();
        else
            limits = get_safety_limits_water();
            thermal_cols = get_thermal_columns_water();
        end
        avail_cols = intersect(thermal_cols, csv_table.Properties.VariableNames, 'stable');
        if ~isempty(avail_cols)
            sub_table = csv_table(:, avail_cols);
            if is_layered
                html = [html table2html(sub_table, limits, avail_cols)];
            else
                html = [html table2html(sub_table, limits, fieldnames(limits)')];
            end
        else
            html = [html '<p class="placeholder">No thermal columns found in CSV.</p>'];
        end
    end

    % Thermal images
    thermal_images = {
        'sub-%03d_%s_maxT%s.png',            'Max temperature overlay';
        'sub-%03d_%s_thermal%s.png',         'Temperature vs time (focal)';
        'sub-%03d_%s_thermalrise%s.png',     'Temperature rise vs time (focal)';
        'sub-%03d_%s_CEM%s.png',             'CEM43 vs time (focal)';
        'sub-%03d_%s_thermal_max%s.png',     'Temperature vs time (layer max)';
        'sub-%03d_%s_thermalrise_max%s.png', 'Temperature rise vs time (layer max)';
        'sub-%03d_%s_CEM_max%s.png',         'CEM43 vs time (layer max)';
        'sub-%03d_%s_thermal_protocol%s.png','Thermal protocol diagram';
    };

    html = [html '<div class="image-grid">'];
    for i = 1:size(thermal_images, 1)
        img_path = fullfile(parameters.io.output_dir, ...
            sprintf(thermal_images{i,1}, subject_id, medium, affix));
        img_html = embed_image(img_path, thermal_images{i,2}, thermal_images{i,2});
        if ~isempty(img_html)
            html = [html img_html];
        end
    end
    html = [html '</div>'];

    % Note about heating animation
    avi_path = fullfile(parameters.io.output_dir, ...
        sprintf('sub-%03d_%s_heating_animation%s.avi', subject_id, medium, affix));
    if isfile(avi_path)
        html = [html sprintf('<p class="note">Heating animation available: <code>%s</code></p>', ...
            html_escape(avi_path))];
    end

end

function html = build_debug_section(parameters, subject_id, medium, affix)
    html = '';

    debug_dir = '';
    if isfield(parameters.io, 'debug_dir')
        debug_dir = parameters.io.debug_dir;
    end
    if isempty(debug_dir) || ~isfolder(debug_dir)
        html = '<p class="placeholder">Debug directory not found.</p>';
        return
    end

    n_trans = 1;
    if isfield(parameters, 'transducer')
        n_trans = numel(parameters.transducer);
    end

    % Debug images (non-transducer-specific)
    debug_images = {
        'sub-%03d_t1_after_rotating_and_scaling%s.png',          'T1 rotation & scaling';
        'sub-%03d_segmented_after_rotating_and_scaling%s.png',   'Segmentation rotation & scaling';
        'sub-%03d_after_rotating_and_scaling_orig%s.png',        'Bone mask rotation & scaling';
        'sub-%03d_t1_skin_skull%s.png',                          'T1 skin/skull overlay';
    };

    html = [html '<div class="image-grid image-grid--narrow">'];
    for i = 1:size(debug_images, 1)
        img_path = fullfile(debug_dir, sprintf(debug_images{i,1}, subject_id, affix));
        img_html = embed_image(img_path, debug_images{i,2}, debug_images{i,2});
        if ~isempty(img_html)
            html = [html img_html];
        end
    end

    % Per-transducer debug images
    for t = 1:n_trans
        img_path = fullfile(debug_dir, ...
            sprintf('sub-%03d_%s_segmented_brain_final_T%02d%s.png', subject_id, medium, t, affix));
        img_html = embed_image(img_path, sprintf('Final segmentation T%02d', t), ...
            sprintf('Final segmentation — T%02d', t));
        if ~isempty(img_html)
            html = [html img_html];
        end
    end

    % Attenuation fit
    img_path = fullfile(debug_dir, sprintf('attenuation_fit%s.png', affix));
    img_html = embed_image(img_path, 'Attenuation fit', 'Attenuation fit');
    if ~isempty(img_html)
        html = [html img_html];
    end

    html = [html '</div>'];
end

function html = build_posthoc_section(parameters)
    if isfield(parameters.modules, 'run_posthoc_water_sims') && parameters.modules.run_posthoc_water_sims
        html = '<span class="badge badge-green">Requested</span>';
    else
        html = '<span class="badge badge-gray">Not requested</span>';
    end
end

function html = build_performance_section(parameters, subject_id, medium, affix)
    html = '';

    log_path = find_log_file(parameters, subject_id, medium, affix);
    if isempty(log_path)
        html = '<p class="placeholder">No log file found.</p>';
        return
    end

    timer_data = parse_timer_lines(log_path);
    if isempty(timer_data)
        html = '<p class="placeholder">No timer data found in log.</p>';
        return
    end

    html = [html '<div class="table-wrapper"><table class="data-table">'];
    html = [html '<thead><tr><th>Stage</th><th>Duration (s)</th><th>Peak RAM (GB)</th>' ...
            '<th>&Delta;RAM (GB)</th><th>&Delta;HDD (GB)</th></tr></thead>'];
    html = [html '<tbody>'];
    for i = 1:length(timer_data)
        d = timer_data(i);
        html = [html sprintf('<tr><td>%s</td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%.2f</td></tr>', ...
            html_escape(d.label), d.duration_s, d.peak_ram_gb, d.delta_ram_gb, d.delta_hdd_gb)];
    end
    html = [html '</tbody></table></div>'];
end

function html = build_log_section(parameters, subject_id, medium, affix)
    html = '';

    log_path = find_log_file(parameters, subject_id, medium, affix);
    if isempty(log_path)
        html = '<p class="placeholder">No log file found.</p>';
        return
    end

    % Read all lines in the log (beginning crucial for parameters)
    fid = fopen(log_path, 'r');
    if fid == -1
        html = '<p class="placeholder">Cannot read log file.</p>';
        return
    end
    lines = {};
    line = fgetl(fid);
    while ischar(line)
        lines{end+1} = line; %#ok<AGROW>
        line = fgetl(fid);
    end
    fclose(fid);

    n = length(lines);
    tail_lines = lines(end-n+1:end);
    log_text = strjoin(tail_lines, newline);

    html = [html sprintf('<p class="note">Log: <code>%s</code></p>', ...
        html_escape(log_path))];
    html = [html '<pre class="log-block">' html_escape(log_text) '</pre>'];
end

function html = build_pseudoCT_section(parameters)
    html = '';
    
    % Extract from parameters
    debug_dir = parameters.io.output_dir;
    affix = '';
    if isfield(parameters.io, 'output_affix')
        affix = parameters.io.output_affix;
    end
    debug_path = fullfile(debug_dir, 'debug');
    
    % 1. MAPPING ALGORITHMS TABLE FIRST
    html = [html, '<dl class="mapping-list">'];
    seg_fields = {'density', 'soundspeed', 'attenuation'};
    labels = {'Density', 'Sound speed', 'Attenuation'};

    for i = 1:3
        if isfield(parameters, 'pct') && isfield(parameters.pct, ['mapping_' seg_fields{i}])
            val = html_escape(char(parameters.pct.(['mapping_' seg_fields{i}])));
            html = [html, sprintf('<dt class="inline-dt">%s mapping: <code>%s</code></dt>', labels{i}, val)];
        end
    end
    html = [html, '</dl>'];
    
    % NEW LINE BEFORE IMAGES
    html = [html, '<div style="height: 20px;"></div>'];
    
    % 2. pCT histograms
    pct_hist_path = fullfile(debug_path, sprintf('pCT_histograms%s.png', affix));
    if isfile(pct_hist_path)
        rel_path = sprintf('debug/pCT_histograms%s.png', affix);
        html = [html, sprintf('<div class="image-row"><figure class="image-container">')];
        html = [html, sprintf('<img src="%s" alt="pCT histograms" ', rel_path)];
        html = [html, sprintf('onclick="openLightbox(%d)" style="max-width: 400px; height: auto;">', 1)];
        html = [html, sprintf('<figcaption>pCT histograms [%s]</figcaption></figure></div>', affix)];
        
        % NEW LINE AFTER histogram image
        html = [html, '<div style="height: 15px;"></div>'];
    end
    
    % 3. KPLAN mapping - ONLY when pct_mapping_density == 'k-plan'
    lightbox_idx = 2;
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'mapping_density') && strcmp(parameters.pct.mapping_density, 'k-plan')
        kplan_path = fullfile(debug_path, 'pCT_hounsfield-density_kplan.png');
        if isfile(kplan_path)
            html = [html, sprintf('<div class="image-row"><figure class="image-container">')];
            html = [html, sprintf('<img src="debug/pCT_hounsfield-density_kplan.png" alt="KPLAN mapping" ')];
            html = [html, sprintf('onclick="openLightbox(%d)" style="max-width: 400px; height: auto;">', lightbox_idx)];
            html = [html, sprintf('<figcaption>Hounsfield → Density (KPLAN)</figcaption></figure></div>')];
            lightbox_idx = lightbox_idx + 1;
            
            % NEW LINE AFTER KPLAN image
            html = [html, '<div style="height: 15px;"></div>'];
        end
    end
    
    % 4. k-Wave mapping - ONLY when pct_mapping_density == 'k-wave'
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'mapping_density') && strcmp(parameters.pct.mapping_density, 'k-wave')
        kwave_path = fullfile(debug_path, 'pCT_hounsfield-density_kwave.png');
        if isfile(kwave_path)
            html = [html, sprintf('<div class="image-row"><figure class="image-container">')];
            html = [html, sprintf('<img src="debug/pCT_hounsfield-density_kwave.png" alt="k-Wave mapping" ')];
            html = [html, sprintf('onclick="openLightbox(%d)" style="max-width: 400px; height: auto;">', lightbox_idx)];
            html = [html, sprintf('<figcaption>Hounsfield → Density (k-Wave)</figcaption></figure></div>')];
            
            % NEW LINE AFTER k-Wave image
            html = [html, '<div style="height: 15px;"></div>'];
        end
    end
    
    % Fallback if no images
    if isempty(strtrim(html)) && ~any(ismember(fields, fieldnames(parameters)))
        html = '<p class="placeholder">No pCT configuration or images found</p>';
    end
    
    html = [html, '<p class="note"><strong>pCT:</strong> Hounsfield → acoustic properties via mappings</p>'];
end

%% ========================================================================
%  IMAGE ENCODING
%  ========================================================================

function b64 = encode_png_base64(filepath)
% Read image file and return base64 string.
    b64 = '';
    if ~isfile(filepath), return; end
    fid = fopen(filepath, 'r');
    if fid == -1, return; end
    raw = fread(fid, '*uint8');
    fclose(fid);
    b64 = matlab.net.base64encode(raw);
end

function html = embed_image(filepath, alt_text, caption)
% Returns <figure> with base64-encoded <img>, or empty string if file missing.
    html = '';
    if ~isfile(filepath), return; end

    b64 = encode_png_base64(filepath);
    if isempty(b64), return; end

    % Detect MIME type
    [~, ~, ext] = fileparts(filepath);
    switch lower(ext)
        case '.png', mime = 'image/png';
        case {'.jpg', '.jpeg'}, mime = 'image/jpeg';
        otherwise, mime = 'image/png';
    end

    html = '<figure>';
    html = [html sprintf('<img src="data:%s;base64,%s" alt="%s">', ...
        mime, b64, html_escape(alt_text))];
    html = [html sprintf('<figcaption>%s</figcaption>', html_escape(caption))];
    html = [html '</figure>'];
end

%% ========================================================================
%  TABLE RENDERING
%  ========================================================================

function html = table2html(tbl, limits, color_columns)
% Convert MATLAB table to HTML <table> with optional color-coded cells.
    if nargin < 2, limits = struct(); end
    if nargin < 3, color_columns = {}; end

    cols = tbl.Properties.VariableNames;

    html = '<div class="table-wrapper"><table class="data-table"><thead><tr>';
    for c = 1:length(cols)
        html = [html sprintf('<th>%s</th>', html_escape(cols{c}))];
    end
    html = [html '</tr></thead><tbody>'];

    for r = 1:height(tbl)
        html = [html '<tr>'];
        for c = 1:length(cols)
            col_name = cols{c};
            val = tbl{r, c};
            val_str = format_cell_value(val);

            % Color coding for safety-relevant columns
            cell_class = '';
            if ismember(col_name, color_columns) && isfield(limits, col_name) && isnumeric(val) && isscalar(val)
                color = safety_color(val, limits.(col_name).limit);
                if ~strcmp(color, 'info')
                    cell_class = sprintf(' class="cell-%s"', color);
                end
            end

            html = [html sprintf('<td%s>%s</td>', cell_class, val_str)];
        end
        html = [html '</tr>'];
    end
    html = [html '</tbody></table></div>'];
end

function str = format_cell_value(val)
% Format a table cell value for HTML display.
    if isnumeric(val)
        if isscalar(val)
            if isnan(val)
                str = 'N/A';
            elseif val == round(val) && abs(val) < 1e6
                str = sprintf('%d', val);
            else
                str = sprintf('%.4g', val);
            end
        else
            str = ['[' strtrim(sprintf('%.4g ', val)) ']'];
        end
    elseif ischar(val) || isstring(val)
        str = html_escape(char(val));
    elseif iscell(val)
        str = html_escape(char(val{1}));
    else
        str = '—';
    end
end

%% ========================================================================
%  COLUMN DEFINITIONS
%  ========================================================================

function cols = get_acoustic_columns()
    cols = {'MI_tc', 'MI_brain', 'MI_skull', 'MI_skin', ...
            'Isppa_brain', 'Isppa_skull', 'Isppa_skin', ...
            'Psptp_brain', 'Psptp_skull', 'Psptp_skin', 'Ptp_target'};
end

function cols = get_thermal_columns()
    cols = {'maxT', 'maxT_brain', 'maxT_skull', 'maxT_skin', ...
            'riseT_brain', 'riseT_skull', 'riseT_skin', ...
            'CEM43_brain', 'CEM43_skull', 'CEM43_skin', ...
            'endT', 'endT_brain', 'endT_skull', 'endT_skin', ...
            'rise_endT_brain', 'rise_endT_skull', 'rise_endT_skin', ...
            'maxCEM43', 'maxCEM43end', ...
            'CEM43_end_brain', 'CEM43_end_skull', 'CEM43_end_skin'};
end

function cols = get_acoustic_columns_water()
    cols = {'subject_id', 'freq_Hz', 'Isppa', 'Isppa_after_exitplane', ...
            'Psptp', 'Ptp_target', 'real_focal_distance_mm', ...
            'Ipa_target', 'Ipa_target_radius'};
end

function cols = get_thermal_columns_water()
    cols = {'maxT', 'endT', 'maxCEM43', 'maxCEM43end'};
end

%% ========================================================================
%  LOG PARSING
%  ========================================================================

function log_path = find_log_file(parameters, subject_id, medium, affix)
% Find the most recent diary log file.
    log_path = '';
    if ~isfield(parameters, 'io') || ~isfield(parameters.io, 'output_dir') || ~isfolder(parameters.io.output_dir)
        return
    end

    pattern = sprintf('sub-%03d_%s%s_*.txt', subject_id, medium, affix);
    files = dir(fullfile(parameters.io.output_dir, pattern));
    if isempty(files), return; end

    [~, idx] = max([files.datenum]);
    log_path = fullfile(files(idx).folder, files(idx).name);
end

function data = parse_timer_lines(log_filepath)
% Parse timer lines from diary log. Returns struct array.
    data = struct('label', {}, 'duration_s', {}, 'delta_ram_gb', {}, ...
                  'peak_ram_gb', {}, 'delta_hdd_gb', {});

    fid = fopen(log_filepath, 'r');
    if fid == -1, return; end
    content = fread(fid, '*char')';
    fclose(fid);

    % Match timer stop lines
    % Format: ⏱ <label>  <time>s | ΔRAM <val>GB | PEAKRAM <val>GB | ΔHDD <val>GB (<val>MB) [#<n>]
    % char(9201) = U+23F1 = ⏱ (STOPWATCH), char(916) = U+0394 = Δ (DELTA)
    expr = [char(9201) '\s+(\S+)\s+([\d.]+)s\s+\|\s+' char(916) 'RAM\s+([-\d.]+)GB\s+\|\s+PEAKRAM\s+([\d.]+)GB\s+\|\s+' char(916) 'HDD\s+([-\d.]+)GB'];
    tokens = regexp(content, expr, 'tokens');

    for i = 1:length(tokens)
        t = tokens{i};
        entry.label = t{1};
        entry.duration_s = str2double(t{2});
        entry.delta_ram_gb = str2double(t{3});
        entry.peak_ram_gb = str2double(t{4});
        entry.delta_hdd_gb = str2double(t{5});
        data(end+1) = entry; %#ok<AGROW>
    end
end

%% ========================================================================
%  UTILITY FUNCTIONS
%  ========================================================================

function html = lightbox_html()
% Returns the lightbox overlay markup and JavaScript.
    html = ['<div id="lightbox" class="lb-overlay" style="display:none;" onclick="closeLightbox()">' ...
            '<span class="lb-close" onclick="closeLightbox()">&times;</span>' ...
            '<img id="lb-img" src="" alt="Enlarged image" onclick="event.stopPropagation()">' ...
            '</div>' ...
            '<script>' ...
            'function openLightbox(src){var lb=document.getElementById("lightbox");' ...
            'document.getElementById("lb-img").src=src;lb.style.display="flex";}' ...
            'function closeLightbox(){document.getElementById("lightbox").style.display="none";}' ...
            'document.addEventListener("keydown",function(e){if(e.key==="Escape")closeLightbox();});' ...
            'document.querySelectorAll(".image-grid figure img").forEach(function(img){' ...
            'img.addEventListener("click",function(){openLightbox(this.src);});});' ...
            '</script>'];
end

function html = collapsible_section(title, content_html, is_open, section_id)
% Wraps content in a <details>/<summary> collapsible section.
% is_open: true to default open, false to default collapsed.
    if is_open
        html = sprintf('<details class="report-section" id="%s" open>', section_id);
    else
        html = sprintf('<details class="report-section" id="%s">', section_id);
    end
    html = [html sprintf('<summary><h2>%s</h2></summary>', html_escape(title))];
    html = [html '<div class="section-content">'];
    html = [html content_html];
    html = [html '</div></details>'];
end

function html = build_toc(parameters, is_layered)
% Build a sticky table-of-contents navigation bar.
    html = '<nav class="toc" id="toc">';

    % Safety status dot: determine worst-case color
    worst = 'green';
    try
        if isfield(parameters.io, 'filename_output_table') && isfile(parameters.io.filename_output_table)
            csv_table = readtable(parameters.io.filename_output_table, 'VariableNamingRule', 'preserve');
            if is_layered
                limits = get_safety_limits();
            else
                limits = get_safety_limits_water();
            end
            metric_names = fieldnames(limits);
            for i = 1:length(metric_names)
                name = metric_names{i};
                info = limits.(name);
                if ~isempty(csv_table) && ismember(name, csv_table.Properties.VariableNames)
                    val = csv_table.(name);
                    if isnumeric(val) && ~isempty(val)
                        c = safety_color(val(end), info.limit);
                        if strcmp(c, 'red')
                            worst = 'red';
                        elseif strcmp(c, 'amber') && ~strcmp(worst, 'red')
                            worst = 'amber';
                        end
                    end
                end
            end
        end
    catch
    end
    html = [html sprintf('<span class="toc-status" style="background:var(--%s);" title="Overall safety: %s"></span>', worst, worst)];

    % Always-present links
    links = {
        'safety',      'Safety';
        'header',      'Header';
        'summary',     'Summary';
        'config',      'Config';
    };
    for i = 1:size(links, 1)
        html = [html sprintf('<a href="#%s">%s</a>', links{i,1}, links{i,2})];
    end

    % Medium properties (layered only)
    if is_layered
        html = [html '<a href="#medium">Medium</a>'];
    end

    % pCT (pCT only)
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && parameters.pct.enabled
        html = [html '<a href="#pseudoCT">pseudo-CT</a>'];
    end

    html = [html '<a href="#positioning">Positioning</a>'];
    html = [html '<a href="#acoustic">Acoustic</a>'];

    % Conditional links
    if isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims
        html = [html '<a href="#thermal">Thermal</a>'];
    end
    if isfield(parameters.simulation, 'debug') && parameters.simulation.debug && is_layered
        html = [html '<a href="#debug">Debug</a>'];
    end
    if is_layered
        html = [html '<a href="#posthoc">Post-hoc</a>'];
    end
    html = [html '<a href="#performance">Performance</a>'];
    html = [html '<a href="#log">Log</a>'];

    html = [html '</nav>'];
end

function val = safe_field(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

function html = config_row(html, label, value)
    html = [html sprintf('<tr><th>%s</th><td>%s</td></tr>', ...
        html_escape(label), html_escape(char(string(value))))];
end

function html = section_error(section_name, ME)
    html = sprintf(['<section class="report-section error-section">' ...
        '<h2>%s</h2><p class="error-notice">Section failed: %s</p></section>'], ...
        html_escape(section_name), html_escape(ME.message));
end

function str = html_escape(str)
% Escape <, >, & for safe HTML embedding.
    if ~ischar(str), str = char(string(str)); end
    str = strrep(str, '&', '&amp;');
    str = strrep(str, '<', '&lt;');
    str = strrep(str, '>', '&gt;');
    str = strrep(str, '"', '&quot;');
end

%% ========================================================================
%  CSS STYLES
%  ========================================================================

function css = css_styles()
css = [...
':root { --green: #22c55e; --amber: #f59e0b; --red: #ef4444; --gray: #9ca3af; --info: #3b82f6; }' newline ...
'html { scroll-behavior: smooth; scroll-padding-top: 60px; }' newline ...
'* { box-sizing: border-box; margin: 0; padding: 0; }' newline ...
'body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; ' ...
  'max-width: 1200px; margin: 0 auto; padding: 20px; background: #f8fafc; color: #1e293b; line-height: 1.6; }' newline ...
'h1 { font-size: 1.8em; margin-bottom: 0.3em; }' newline ...
'h2 { font-size: 1.3em; margin-bottom: 0.8em; color: #334155; border-bottom: 2px solid #e2e8f0; padding-bottom: 0.3em; }' newline ...
'.report-section { background: white; border-radius: 8px; padding: 24px; margin-bottom: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }' newline ...
'.error-section { border-left: 4px solid var(--red); }' newline ...
'.error-notice { color: var(--red); font-style: italic; }' newline ...
newline ...
'/* Safety dashboard */' newline ...
'.safety-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(160px, 1fr)); gap: 12px; margin-bottom: 12px; }' newline ...
'.safety-card { padding: 14px; border-radius: 6px; }' newline ...
'.safety-green  { border-left: 4px solid var(--green); background: #f0fdf4; }' newline ...
'.safety-amber  { border-left: 4px solid var(--amber); background: #fffbeb; }' newline ...
'.safety-red    { border-left: 4px solid var(--red);   background: #fef2f2; }' newline ...
'.safety-gray   { border-left: 4px solid var(--gray);  background: #f9fafb; }' newline ...
'.safety-info   { border-left: 4px solid var(--info);  background: #eff6ff; }' newline ...
'.safety-label  { font-size: 0.8em; color: #64748b; text-transform: uppercase; letter-spacing: 0.05em; }' newline ...
'.safety-value  { font-size: 1.5em; font-weight: 700; margin: 4px 0; }' newline ...
'.safety-limit  { font-size: 0.75em; color: #94a3b8; }' newline ...
'.safety-footnote { font-size: 0.8em; color: #94a3b8; margin-top: 8px; }' newline ...
'.safety-footnote a { color: #64748b; text-decoration: none; border-bottom: 1px dotted #94a3b8; }' newline ...
'.safety-footnote a:hover { color: #1e293b; border-bottom-style: solid; }' newline ...
newline ...
'/* Tables */' newline ...
'.info-table, .config-table { border-collapse: collapse; }' newline ...
'.info-table th, .config-table th { text-align: left; padding: 4px 16px 4px 0; color: #64748b; font-weight: 500; }' newline ...
'.info-table td, .config-table td { padding: 4px 0; }' newline ...
'.table-wrapper { overflow-x: auto; -webkit-overflow-scrolling: touch; margin: 12px 0; }' newline ...
'.data-table { width: 100%; border-collapse: collapse; font-size: 0.9em; }' newline ...
'.data-table thead { background: #f1f5f9; }' newline ...
'.data-table th, .data-table td { padding: 8px 12px; text-align: left; border-bottom: 1px solid #e2e8f0; white-space: nowrap; }' newline ...
'.data-table tbody tr:nth-child(even) { background: #f8fafc; }' newline ...
'.cell-green { background: #f0fdf4 !important; }' newline ...
'.cell-amber { background: #fffbeb !important; }' newline ...
'.cell-red   { background: #fef2f2 !important; font-weight: 600; }' newline ...
'.cell-gray  { background: #f9fafb !important; color: #94a3b8; }' newline ...
newline ...
'/* Images */' newline ...
'.image-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(min(400px, 100%), 1fr)); gap: 16px; margin: 16px 0; }' newline ...
'.image-grid--wide { grid-template-columns: repeat(auto-fill, minmax(min(600px, 100%), 1fr)); }' newline ...
'.image-grid--narrow { grid-template-columns: repeat(auto-fill, minmax(min(280px, 100%), 1fr)); }' newline ...
'figure { margin: 0; min-width: 0; overflow: hidden; }' newline ...
'figure img { display: block; width: 100%; height: auto; border-radius: 4px; border: 1px solid #e2e8f0; cursor: zoom-in; }' newline ...
'figcaption { font-size: 0.85em; color: #64748b; margin-top: 4px; text-align: center; }' newline ...
newline ...
'/* Badges */' newline ...
'.badge { display: inline-block; padding: 4px 12px; border-radius: 12px; font-size: 0.85em; font-weight: 500; }' newline ...
'.badge-green { background: #f0fdf4; color: #166534; border: 1px solid #bbf7d0; }' newline ...
'.badge-gray  { background: #f9fafb; color: #6b7280; border: 1px solid #e5e7eb; }' newline ...
'.badge-blue  { background: #eff6ff; color: #1e40af; border: 1px solid #bfdbfe; }' newline ...
'.badge-amber { background: #fffbeb; color: #92400e; border: 1px solid #fde68a; }' newline ...
newline ...
'/* Medium banners */' newline ...
'.medium-banner { padding: 12px 16px; border-radius: 6px; margin-bottom: 16px; font-size: 0.9em; }' newline ...
'.medium-water { background: #eff6ff; border-left: 4px solid var(--info); color: #1e40af; }' newline ...
'.medium-layered { background: #fffbeb; border-left: 4px solid var(--amber); color: #92400e; }' newline ...
newline ...
'/* Summary cards */' newline ...
'.summary-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr)); gap: 12px; }' newline ...
'.summary-card { background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 14px; text-align: center; }' newline ...
'.summary-label { font-size: 0.8em; color: #64748b; text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 4px; }' newline ...
'.summary-value { font-size: 1.4em; font-weight: 700; color: #1e293b; }' newline ...
'.summary-unit { font-size: 0.75em; color: #94a3b8; margin-top: 2px; }' newline ...
newline ...
'/* Log & misc */' newline ...
'.log-block { background: #1e293b; color: #e2e8f0; padding: 16px; border-radius: 6px; ' ...
  'overflow-x: auto; font-size: 0.8em; line-height: 1.5; max-height: 600px; overflow-y: auto; white-space: pre-wrap; word-wrap: break-word; }' newline ...
'.placeholder { color: #94a3b8; font-style: italic; }' newline ...
'.note { font-size: 0.85em; color: #64748b; margin: 8px 0; }' newline ...
'code { background: #f1f5f9; padding: 2px 6px; border-radius: 3px; font-size: 0.9em; }' newline ...
'footer { text-align: center; color: #94a3b8; font-size: 0.8em; padding: 20px 0; }' newline ...
newline ...
'/* Collapsible sections */' newline ...
'details.report-section { padding: 0; }' newline ...
'details.report-section > summary { cursor: pointer; list-style: none; padding: 16px 24px; display: flex; align-items: center; gap: 8px; }' newline ...
'details.report-section > summary::-webkit-details-marker { display: none; }' newline ...
'details.report-section > summary::before { content: "\25B6"; font-size: 0.7em; transition: transform 0.2s ease; flex-shrink: 0; }' newline ...
'details[open].report-section > summary::before { transform: rotate(90deg); }' newline ...
'details.report-section > summary > h2 { border-bottom: none; margin-bottom: 0; padding-bottom: 0; }' newline ...
'details.report-section > .section-content { padding: 0 24px 24px 24px; }' newline ...
newline ...
'/* Sticky TOC */' newline ...
'.toc { position: sticky; top: 0; z-index: 100; background: white; border-radius: 8px; ' ...
  'padding: 10px 20px; margin-bottom: 20px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); ' ...
  'display: flex; align-items: center; gap: 12px; flex-wrap: wrap; }' newline ...
'.toc a { text-decoration: none; color: #475569; font-size: 0.85em; padding: 4px 10px; ' ...
  'border-radius: 4px; white-space: nowrap; transition: background 0.15s; }' newline ...
'.toc a:hover { background: #f1f5f9; color: #1e293b; }' newline ...
'.toc-status { width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }' newline ...
newline ...
'/* Lightbox */' newline ...
'.lb-overlay { position: fixed; inset: 0; background: rgba(0,0,0,0.85); z-index: 1000; ' ...
  'display: flex; align-items: center; justify-content: center; cursor: zoom-out; }' newline ...
'.lb-overlay img { max-width: 90vw; max-height: 85vh; object-fit: contain; cursor: default; border-radius: 4px; }' newline ...
'.lb-close { position: fixed; top: 16px; right: 24px; color: white; font-size: 2em; cursor: pointer; z-index: 1001; line-height: 1; }' newline ...
newline ...
'/* Safety progress bars */' newline ...
'.safety-bar-track { height: 4px; background: #e2e8f0; border-radius: 2px; margin-top: 6px; overflow: hidden; }' newline ...
'.safety-bar { height: 100%; border-radius: 2px; transition: width 0.3s; }' newline ...
'.safety-green .safety-bar { background: var(--green); }' newline ...
'.safety-amber .safety-bar { background: var(--amber); }' newline ...
'.safety-red .safety-bar { background: var(--red); }' newline ...
newline ...
'/* Print */' newline ...
'@media print { html { scroll-behavior: auto; } body { background: white; } ' ...
  '.report-section { box-shadow: none; border: 1px solid #e2e8f0; break-inside: avoid; } ' ...
  'details { display: block; } details > summary { display: none; } details > .section-content { padding: 0 24px 24px 24px; } ' ...
  '.log-block { max-height: none; } .toc { display: none; } .lb-overlay { display: none; } ' ...
  'figure { break-inside: avoid; } .safety-grid { break-inside: avoid; } ' ...
  'figure img { cursor: default; } ' ...
  '.image-grid { grid-template-columns: repeat(2, 1fr); } .image-grid--wide { grid-template-columns: 1fr; } }' newline ...
];
end

