function report_path = generate_simulation_report(parameters)
% GENERATE_SIMULATION_REPORT  Generate a self-contained HTML simulation report
%
% Consolidates all PRESTUS simulation outputs (CSV data, images, logs)
% into a single portable HTML file with a safety dashboard showing
% colour-coded metrics against ITRUSST consensus limits.
%
% Use as:
%   report_path = generate_simulation_report(parameters)
%
% Input:
%   parameters - (1,1) simulation parameters struct
%
% Output:
%   report_path - path to the generated HTML file
%
% See also: GENERATE_UNCERTAINTY_REPORT, CSS_STYLES_BASE

arguments
    parameters (1,1) struct
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
        html_parts{end+1} = html_utils.section_error('Table of Contents', ME);
    end

    % Section 1: Safety Dashboard (always open, no toggle)
    try
        html_parts{end+1} = build_safety_dashboard(csv_table, parameters, is_layered);
    catch ME
        html_parts{end+1} = html_utils.section_error('Safety Dashboard', ME);
    end

    % Section 2: Header (always open, no toggle)
    try
        html_parts{end+1} = build_header(subject_id, medium, affix, parameters, is_layered);
    catch ME
        html_parts{end+1} = html_utils.section_error('Header', ME);
    end

    % Section 3: Methods Boilerplate (collapsed by default)
    try
        html_parts{end+1} = html_utils.collapsible('Methods', ...
            build_methods_boilerplate(parameters, is_layered), false, 'methods');
    catch ME
        html_parts{end+1} = html_utils.section_error('Methods', ME);
    end

    % Section 4 (formerly 3): Simulation Summary (always open, no toggle)
    try
        html_parts{end+1} = build_simulation_summary(csv_table, parameters, is_layered);
    catch ME
        html_parts{end+1} = html_utils.section_error('Simulation Summary', ME);
    end

    % Section 5: Configuration Summary (collapsed by default)
    try
        html_parts{end+1} = html_utils.collapsible('Configuration Summary', ...
            build_config_summary(parameters), false, 'config');
    catch ME
        html_parts{end+1} = html_utils.section_error('Configuration Summary', ME);
    end

    % Section 5: Medium Properties (layered only, collapsed by default)
    if is_layered
        try
            html_parts{end+1} = html_utils.collapsible('Medium Properties', ...
                build_medium_properties_section(parameters), false, 'medium');
        catch ME
            html_parts{end+1} = html_utils.section_error('Medium Properties', ME);
        end
    end

    % Section 6: pseudoCT [Optional] (collapsed by default)
    try
        if isfield(parameters, 'pct') && isfield(parameters.pct, 'enabled') && parameters.pct.enabled
            html_parts{end+1} = html_utils.collapsible('pseudo-CT', ...
                build_pseudoCT_section(parameters), ...
                false, 'pseudoCT');
        end
    catch ME
        html_parts{end+1} = html_utils.section_error('pseudo-CT', ME);
    end

    % Section 7: Positioning (open by default)
    try
        html_parts{end+1} = html_utils.collapsible('Transducer Positioning', ...
            build_positioning_section(parameters, subject_id, affix), true, 'positioning');
    catch ME
        html_parts{end+1} = html_utils.section_error('Positioning', ME);
    end

    % Section 8: Acoustic Results (open by default)
    try
        html_parts{end+1} = html_utils.collapsible('Acoustic Results', ...
            build_acoustic_section(csv_table, parameters, subject_id, medium, affix, is_layered), true, 'acoustic');
    catch ME
        html_parts{end+1} = html_utils.section_error('Acoustic Results', ME);
    end

    % Section 9: Thermal Results (open by default, conditional)
    try
        if isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims
            html_parts{end+1} = html_utils.collapsible('Thermal Results', ...
                build_thermal_section(csv_table, parameters, subject_id, medium, affix, is_layered), true, 'thermal');
        end
    catch ME
        html_parts{end+1} = html_utils.section_error('Thermal Results', ME);
    end

    % Section 10: Debug Information (collapsed by default, layered only)
    try
        if isfield(parameters.simulation, 'debug') && parameters.simulation.debug && is_layered
            html_parts{end+1} = html_utils.collapsible('Debug Information', ...
                build_debug_section(parameters, subject_id, medium, affix), false, 'debug');
        end
    catch ME
        html_parts{end+1} = html_utils.section_error('Debug Information', ME);
    end

    % Section 11: Post-hoc Water Simulation (collapsed by default, layered only)
    if is_layered
        try
            html_parts{end+1} = html_utils.collapsible('Post-hoc Water Simulation', ...
                build_posthoc_section(parameters), false, 'posthoc');
        catch ME
            html_parts{end+1} = html_utils.section_error('Post-hoc Water Simulation', ME);
        end
    end

    % Section 12: Performance (collapsed by default)
    try
        html_parts{end+1} = html_utils.collapsible('Performance', ...
            build_performance_section(parameters, subject_id, medium, affix), false, 'performance');
    catch ME
        html_parts{end+1} = html_utils.section_error('Performance', ME);
    end

    % Section 13: Log (collapsed by default)
    try
        html_parts{end+1} = html_utils.collapsible('Log', ...
            build_log_section(parameters, subject_id, medium, affix), false, 'log');
    catch ME
        html_parts{end+1} = html_utils.section_error('Log', ME);
    end

    % Footer
    html_parts{end+1} = '<footer>';
    html_parts{end+1} = sprintf('<p>Generated by PRESTUS %s</p>', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    html_parts{end+1} = '</footer>';
    % Lightbox overlay and script
    html_parts{end+1} = html_utils.lightbox();

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
%  SAFETY DASHBOARD
%  ========================================================================

function html = build_safety_dashboard(csv_table, parameters, is_layered)
    if is_layered
        limits = get_risk_limits(is_layered);
    else
    end
    metric_names = fieldnames(limits);

    % When ISO CEM43 was requested, remap CSV column lookups to iso variants
    use_cem43_iso = isfield(parameters, 'thermal') && ...
                   isfield(parameters.thermal, 'cem43_iso') && ...
                   parameters.thermal.cem43_iso == 1;
    cem43_iso_map = struct( ...
        'CEM43_brain', 'CEM43iso_brain', ...
        'CEM43_skull', 'CEM43iso_skull', ...
        'CEM43_skin',  'CEM43iso_skin');

    html = '<section class="report-section" id="safety">';
    html = [html '<h2>Safety Dashboard</h2>'];

    if ~is_layered
        html = [html '<div class="medium-banner medium-water">' ...
            '<strong>Water / Free-field</strong> &mdash; tissue-specific non-significant risk limits not applicable. ' ...
            'Only global metrics are shown.</div>'];
    end

    html = [html '<div class="safety-grid">'];

    for i = 1:length(metric_names)
        name = metric_names{i};
        info = limits.(name);

        % Extract value from CSV (remap CEM43 fields to ISO variant when requested)
        csv_name = name;
        if use_cem43_iso && isfield(cem43_iso_map, name)
            csv_name = cem43_iso_map.(name);
        end
        value = NaN;
        display_unit = info.unit;
        if ~isempty(csv_table) && ismember(csv_name, csv_table.Properties.VariableNames)
            val = csv_table.(csv_name);
            if isnumeric(val) && ~isempty(val)
                value = val(end); % last row if multiple
            end
        end

        % Special handling for pressure: use dynamic unit scaling
        if strcmp(name, 'Psptp') && ~isnan(value)
            [value, display_unit] = scale_pressure(value);
        end

        % Determine color
        color = risk_color(value, info.limit);

        % Build card
        html = [html sprintf('<div class="safety-card safety-%s">', color)];
        html = [html sprintf('<div class="safety-label">%s</div>', html_utils.escape(info.label))];

        if isnan(value)
            html = [html '<div class="safety-value">N/A</div>'];
            html = [html '<div class="risk-limit">No data</div>'];
        else
            html = [html sprintf('<div class="safety-value">%.3g</div>', value)];
            if isinf(info.limit)
                html = [html sprintf('<div class="risk-limit">%s (informational)</div>', display_unit)];
            else
                html = [html sprintf('<div class="risk-limit">NSR limit: %.3g %s</div>', info.limit, display_unit)];
                % Progress bar showing value as % of limit
                pct = min(100, max(0, (value / info.limit) * 100));
                html = [html sprintf('<div class="safety-bar-track"><div class="safety-bar" style="width:%.0f%%"></div></div>', pct)];
            end
        end

        html = [html '</div>'];
    end

    html = [html '</div>'];
    html = [html '<p class="safety-footnote">Non-significant risk limits: ITRUSST consensus ' ...
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
        html = [html sprintf('<tr><th>Medium</th><td>%s <span class="badge badge-amber">Layered (Realistic Head)</span> <span class="badge badge-gray">pseudo-CT (Continuous Skull)</span></td></tr>', html_utils.escape(medium))];
    elseif is_layered
        html = [html sprintf('<tr><th>Medium</th><td>%s <span class="badge badge-amber">Layered (Realistic Head)</span></td></tr>', html_utils.escape(medium))];
    else
        html = [html sprintf('<tr><th>Medium</th><td>%s <span class="badge badge-blue">Water / Free-field</span></td></tr>', html_utils.escape(medium))];
    end
    if ~isempty(affix)
        html = [html sprintf('<tr><th>Affix</th><td>%s</td></tr>', html_utils.escape(affix))];
    end
    html = [html sprintf('<tr><th>Generated</th><td>%s</td></tr>', datestr(now, 'yyyy-mm-dd HH:MM:SS'))];
    if isfield(parameters, 'output_dir')
        html = [html sprintf('<tr><th>Output dir</th><td>%s</td></tr>', html_utils.escape(parameters.io.output_dir))];
    end
    html = [html '</table>'];
    html = [html '</section>'];
end

function html = build_methods_boilerplate(parameters, is_layered)
% Build a CC0-licensed, fMRIPrep-style methods boilerplate aligned with
% ITRUSST standardised reporting guidelines (Martin et al., 2024,
% Brain Stimulation 17:607-615; DOI:10.1016/j.brs.2024.04.013).
%
% Terminology strictly follows ITRUSST conventions:
%   - "pulse duration (PD)", "pulse repetition interval (PRI)",
%     "pulse repetition frequency (PRF)", "pulse train"  (not "burst")
%   - MI_tc for transcranial mechanical index
%   - I_SPPA and I_SPTA for intensity parameters
%   - DC = PD/PRI × 100% (only for rectangular pulses)
%
% The generated text is released under CC0 1.0 and may be copied
% verbatim into manuscript Methods sections.

    % --- inline helpers ---
    function s = fmt(val, fmtstr, fallback)
        if isnumeric(val) && isscalar(val) && ~isnan(val)
            s = sprintf(fmtstr, val);
        else
            s = fallback;
        end
    end
    function s = fmt_pml(val)
        % Format PML size for report: scalar → '%d', vector → '[a b c]', non-numeric → as-is
        if isnumeric(val) && ~any(isnan(val))
            s = ['[' strtrim(num2str(val)) ']'];
        elseif ischar(val) || isstring(val)
            s = char(val);
        else
            s = '?';
        end
    end
    p  = @(text) sprintf('<p>%s</p>\n', text);
    b  = @(text) sprintf('<strong>%s</strong>', text);
    em = @(text) sprintf('<em>%s</em>', text);

    % ------------------------------------------------------------------ %
    %  Collect parameters
    % ------------------------------------------------------------------ %
    grid_res  = safe_field(parameters.grid, 'resolution_mm', NaN);
    n_dims    = safe_field(parameters.grid, 'n_dims', 3);
    pml       = safe_field(parameters.grid, 'pml_size_effective', ...
                  safe_field(parameters.grid, 'pml_size', NaN));
    code_type = safe_field(parameters.simulation, 'code_type', 'MATLAB');

    % Transducer (first element)
    td       = parameters.transducer(1);
    freq_hz  = safe_field(td, 'freq_hz', NaN);
    freq_khz = fmt(freq_hz / 1e3, '%.0f', '?');
    td_type  = safe_field(td, 'type', '');
    n_elem   = NaN;  curv_mm = NaN;  ap_mm = NaN;
    if ~isempty(td_type) && isfield(td, td_type)
        n_elem  = safe_field(td.(td_type), 'elem_n',           NaN);
        curv_mm = safe_field(td.(td_type), 'curv_radius_mm',   NaN);
        ap_mm   = safe_field(td.(td_type), 'aperture_diameter_mm', NaN);
    end
    focal_mm = safe_field(td, 'focal_distance_bowl', NaN);

    % Tissue layers
    if is_layered && isfield(parameters, 'layers')
        layer_names  = fieldnames(parameters.layers);
        tissue_list  = strrep(strjoin(layer_names, ', '), '_', ' ');
        n_layers     = numel(layer_names);
    else
        tissue_list  = '';
        n_layers     = 0;
    end

    % Thermal / pulse timing
    run_thermal = isfield(parameters, 'modules') && ...
                  isfield(parameters.modules, 'run_heating_sims') && ...
                  parameters.modules.run_heating_sims;

    pd_val = NaN; pri_val = NaN; ptd_val = NaN; ptri_val = NaN;
    if isfield(parameters, 'thermal')
        th      = parameters.thermal;
        pd_val  = safe_field(th, 'pd',   NaN);
        pri_val = safe_field(th, 'pri',  NaN);
        ptd_val = safe_field(th, 'ptd',  NaN);
        ptri_val= safe_field(th, 'ptri', NaN);
    end
    pd_ms   = fmt(pd_val  * 1e3, '%.3g', '?');
    pri_ms  = fmt(pri_val * 1e3, '%.3g', '?');
    prf_hz  = fmt(1 / pri_val,   '%.3g', '?');
    ptd_s   = fmt(ptd_val,       '%.3g', '?');
    ptri_s  = fmt(ptri_val,      '%.3g', '?');
    if ~isnan(pd_val) && ~isnan(pri_val) && pri_val > 0
        dc_pct = fmt(pd_val / pri_val * 100, '%.1f', '?');
    else
        dc_pct = '?';
    end

    % pseudoCT
    use_pct    = isfield(parameters, 'pct') && ...
                 isfield(parameters.pct, 'enabled') && parameters.pct.enabled;
    pct_method = '';
    if use_pct
        pct_method = safe_field(parameters.pct, 'skull_mapping', 'unknown');
    end

    % ------------------------------------------------------------------ %
    %  1. Transducer and drive system
    % ------------------------------------------------------------------ %
    geom_parts = {};
    if ~isnan(curv_mm); geom_parts{end+1} = sprintf('radius of curvature: %.0f mm', curv_mm); end
    if ~isnan(ap_mm);   geom_parts{end+1} = sprintf('aperture diameter: %.0f mm', ap_mm);     end
    if ~isnan(focal_mm);geom_parts{end+1} = sprintf('nominal focal distance: %.0f mm', focal_mm); end
    geom_str = iff(~isempty(geom_parts), [' (' strjoin(geom_parts, '; ') ')'], '');

    para1 = sprintf([...
        '%s The transducer was a %s-element %s array operating at a centre '   ...
        'frequency of %s kHz%s. '                                              ...
        'The transducer was modelled within the k-Wave Toolbox '               ...
        '(Treeby &amp; Cox, 2010) using the kWaveArray discrete source '       ...
        'approach (Bell et al., 2022).'],                                       ...
        b('Transducer and drive system.'),                                      ...
        fmt(n_elem, '%d', '?'), html_utils.escape(td_type),                    ...
        freq_khz, geom_str);

    % ------------------------------------------------------------------ %
    %  2. Simulation grid (in situ approach)
    % ------------------------------------------------------------------ %
    if is_layered
        medium_desc = sprintf(['a layered, subject-specific head model '       ...
            'comprising %d tissue compartments (%s)'], n_layers, tissue_list);
        seg_sentence = [' Tissue segmentation was derived from T1-weighted '   ...
            'MRI using the SimNIBS ' em('charm') ' pipeline '                  ...
            '(Thielscher et al., 2015; Puonti et al., 2020).'];
    else
        medium_desc  = 'a homogeneous water medium (free-field reference condition)';
        seg_sentence = '';
    end

    skull_sentence = '';
    if is_layered && use_pct
        skull_sentence = sprintf([' Skull bone acoustic properties were '      ...
            'mapped from pseudo-CT Hounsfield unit values using the '           ...
            '%s approach.'],                                                    ...
            html_utils.escape(pct_method));
    end

    % Tissue properties sentence — refer to the Medium Properties table in the report
    if is_layered
        props_sentence = [' Tissue acoustic properties (sound speed, density, ' ...
            'and power-law attenuation coefficient) are listed in the '         ...
            '<a href="#medium">Medium Properties</a> section.'];
    else
        props_sentence = '';
    end

    para2 = sprintf([...
        '%s '                                                                   ...
        'Acoustic pressure field propagation was modelled in %dD '             ...
        'using the pseudospectral time-domain method '                          ...
        'on a %s&#8201;mm isotropic Cartesian grid '                           ...
        '(perfectly matched layer: %s grid points). '                          ...
        'The medium was represented as %s.%s%s%s '                             ...
        'The method used to obtain ' em('in situ') ' exposure estimates '      ...
        'was full-wave acoustic simulation through the individual skull '       ...
        'and brain geometry.'],                                                 ...
        b(['Acoustic simulation and ' em('in situ') ' exposure estimates.']),   ...
        n_dims,                                                                 ...
        fmt(grid_res, '%.2g', '?'),                                            ...
        fmt_pml(pml),                                                          ...
        medium_desc, seg_sentence, skull_sentence, props_sentence);

    % ------------------------------------------------------------------ %
    %  3. Pulse timing parameters
    % ------------------------------------------------------------------ %
    has_pulse_train = ~isnan(ptd_val) && ~isnan(ptri_val);

    dc_str = '';
    if ~strcmp(dc_pct, '?')
        dc_str = sprintf(' (duty cycle: %s&#37;)', dc_pct);
    end

    if ~strcmp(pd_ms, '?') && ~strcmp(pri_ms, '?')
        timing_prose = sprintf(['Pulses of %s&#8201;ms duration were delivered '   ...
            'at a pulse repetition frequency of %s&#8201;Hz '                      ...
            '(pulse repetition interval: %s&#8201;ms%s).'],                        ...
            pd_ms, prf_hz, pri_ms, dc_str);
    else
        timing_prose = 'Pulse timing parameters were not fully specified.';
    end

    if has_pulse_train
        timing_prose = [timing_prose sprintf([...
            ' Pulses were organised into pulse trains of %s&#8201;s duration '     ...
            'repeated every %s&#8201;s.'], ptd_s, ptri_s)];
    end

    para3 = [b('Pulse timing parameters.') ' ' timing_prose];

    % ------------------------------------------------------------------ %
    %  4. In situ exposure and safety metrics
    % ------------------------------------------------------------------ %
    if run_thermal
        thermal_limits_str = [...
            'MI<sub>tc</sub> &#8804; 1.9; &#916;T &#8804; 2&#8201;&#176;C; '   ...
            'T<sub>peak</sub> &#8804; 39&#8201;&#176;C; '                       ...
            'CEM43 &#8804; 2&#8201;min (brain), &#8804; 16&#8201;min (skull), ' ...
            '&#8804; 21&#8201;min (skin).'];
    else
        thermal_limits_str = 'MI<sub>tc</sub> &#8804; 1.9.';
    end
    para4 = sprintf([...
        '%s '                                                                   ...
        'The ' em('in situ') ' spatial-peak pulse-average intensity '          ...
        '(I<sub>SPPA</sub>) and the transcranial mechanical index '             ...
        '(MI<sub>tc</sub>) were derived from the simulated pressure fields. '   ...
        'Non-significant risk limits followed the ITRUSST consensus on '        ...
        'biophysical safety (Aubry et al., 2025): %s'],                         ...
        b([em('In situ') ' exposure and safety metrics.']), thermal_limits_str);

    % ------------------------------------------------------------------ %
    %  5. Thermal simulations (conditional)
    % ------------------------------------------------------------------ %
    use_cem43_iso = isfield(parameters, 'thermal') && ...
                   isfield(parameters.thermal, 'cem43_iso') && ...
                   parameters.thermal.cem43_iso == 1;

    para5 = '';
    if run_thermal
        if use_cem43_iso
            cem43_formula = [...
                'CEM43 was computed using the ISO form: '                           ...
                'CEM43 = &#8747; R<sup>(43&#8722;T)</sup> dt, '                    ...
                'where R&#8201;=&#8201;0 for T&#8201;&lt;&#8201;39&#8201;&#176;C '  ...
                '(and T<sub>max</sub>&#8201;&lt;&#8201;43&#8201;&#176;C), '         ...
                'R&#8201;=&#8201;0.25 for 39&#8201;&#8804;&#8201;T&#8201;&lt;&#8201;43&#8201;&#176;C '  ...
                '(and T<sub>max</sub>&#8201;&lt;&#8201;43&#8201;&#176;C), '         ...
                'R&#8201;=&#8201;0.5 for T&#8201;&#8805;&#8201;43&#8201;&#176;C '  ...
                'or where T<sub>max</sub>&#8201;&#8805;&#8201;43&#8201;&#176;C, '   ...
                'and CEM43&#8201;=&#8201;&#8734; for T&#8201;&#8805;&#8201;57&#8201;&#176;C '...
                '(Martin et al., 2024, Eq.&#8201;5).'];
        else
            cem43_formula = [...
                'CEM43 was computed as CEM = &#8747; R<sup>(43&#8722;T)</sup> dt, ' ...
                'where R&#8201;=&#8201;0.5 for T&#8201;&#8805;&#8201;43&#8201;&#176;C '  ...
                'and R&#8201;=&#8201;0.25 for T&#8201;&lt;&#8201;43&#8201;&#176;C '        ...
                '(Sapareto &amp; Dewey, 1984; Martin et al., 2024, Eq.&#8201;5).'];
        end
        para5 = [b('Thermal simulations.') ...
            ' Temperature elevation and cumulative thermal dose (CEM43) '      ...
            'were estimated by solving the Pennes bioheat equation '           ...
            '(Pennes, 1948) over the full sonication protocol described above. '...
            cem43_formula];
    end

    % ------------------------------------------------------------------ %
    %  Assemble single text block
    % ------------------------------------------------------------------ %
    paras = {para1, para2, para3, para4};
    if run_thermal
        paras{end+1} = para5;
    end
    txt = strjoin(cellfun(@(s) p(s), paras, 'UniformOutput', false), '');

    % ------------------------------------------------------------------ %
    %  CC0 notice
    % ------------------------------------------------------------------ %
    cc0 = ['<div class="boilerplate-license">'                                 ...
        '<strong>License.</strong> '                                            ...
        'The above methods text was automatically generated by PRESTUS '        ...
        'and conforms to ITRUSST standardised reporting (Martin et al., 2024). '...
        'It has been released under the '                                       ...
        '<a href="https://creativecommons.org/publicdomain/zero/1.0/" '        ...
        'target="_blank" rel="noopener">CC0 1.0 Universal</a> license '        ...
        'and may be copied verbatim into manuscript Methods sections.'          ...
        '</div>'];

    html = [txt cc0];
end

function s = iff(cond, a, b)
% Inline conditional: return a if cond is true, else b.
    if cond; s = a; else; s = b; end
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
    html = [html sprintf('<div class="summary-label">%s</div>', html_utils.escape(label))];
    html = [html sprintf('<div class="summary-value">%s</div>', val_str)];
    html = [html sprintf('<div class="summary-unit">%s</div>', html_utils.escape(unit))];
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
        html = [html sprintf('<th>%s</th>', html_utils.escape(tissue_labels{t}))];
    end
    html = [html '</tr></thead><tbody>'];

    for p = 1:size(props, 1)
        field = props{p, 1}; label = props{p, 2}; unit = props{p, 3};
        html = [html sprintf('<tr><td><strong>%s</strong></td><td>%s</td>', ...
            html_utils.escape(label), html_utils.escape(unit))];
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
    pml = safe_field(parameters.grid, 'pml_size_effective', ...
              safe_field(parameters.grid, 'pml_size', NaN));
    if isnumeric(pml) && ~isnan(pml)
        html = config_row(html, 'PML size', sprintf('[%s]', strtrim(num2str(pml))));
    elseif ischar(pml) || isstring(pml)
        html = config_row(html, 'PML size', char(pml));
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
        img_html = html_utils.embed_image(img_path, sprintf('Positioning T%02d', t), ...
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
            limits = get_risk_limits(is_layered);
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
        img_html = html_utils.embed_image(img_path, sprintf('Intensity on segmentation%s', trans_suffix), ...
            sprintf('Intensity overlay (segmentation)%s', trans_suffix));
        if ~isempty(img_html)
            html = [html img_html];
        end

        % Intensity on T1
        img_path = fullfile(parameters.io.output_dir, ...
            sprintf('sub-%03d_%s_intensity_t1%s%s.png', subject_id, medium, trans_suffix, affix));
        img_html = html_utils.embed_image(img_path, sprintf('Intensity on T1%s', trans_suffix), ...
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
            limits = get_risk_limits(is_layered);
            thermal_cols = get_thermal_columns();
        else
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
        'sub-%03d_%s_CEM%s.png',             'CEM43 (kWave) vs time (focal)';
        'sub-%03d_%s_CEM_iso%s.png',         'CEM43 (ISO) vs time (focal)';
        'sub-%03d_%s_thermal_max%s.png',     'Temperature vs time (layer max)';
        'sub-%03d_%s_thermalrise_max%s.png', 'Temperature rise vs time (layer max)';
        'sub-%03d_%s_CEM_max%s.png',         'CEM43 (kWave) vs time (layer max)';
        'sub-%03d_%s_CEM_iso_max%s.png',     'CEM43 (ISO) vs time (layer max)';
        'sub-%03d_%s_thermal_protocol%s.png','Thermal protocol diagram';
    };

    html = [html '<div class="image-grid">'];
    for i = 1:size(thermal_images, 1)
        img_path = fullfile(parameters.io.output_dir, ...
            sprintf(thermal_images{i,1}, subject_id, medium, affix));
        img_html = html_utils.embed_image(img_path, thermal_images{i,2}, thermal_images{i,2});
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
            html_utils.escape(avi_path))];
    end

end

function html = build_debug_section(parameters, subject_id, medium, affix)
    html = '';

    debug_dir_preproc = '';
    debug_dir_medium  = '';
    if isfield(parameters.io, 'debug_dir_preproc')
        debug_dir_preproc = parameters.io.debug_dir_preproc;
    elseif isfield(parameters.io, 'debug_dir')
        debug_dir_preproc = fullfile(parameters.io.debug_dir, 'preproc');
    end
    if isfield(parameters.io, 'debug_dir_medium')
        debug_dir_medium = parameters.io.debug_dir_medium;
    elseif isfield(parameters.io, 'debug_dir')
        debug_dir_medium = fullfile(parameters.io.debug_dir, 'medium');
    end
    if isempty(debug_dir_preproc) || ~isfolder(debug_dir_preproc)
        html = '<p class="placeholder">Debug directory not found.</p>';
        return
    end

    n_trans = 1;
    if isfield(parameters, 'transducer')
        n_trans = numel(parameters.transducer);
    end

    % Preprocessing debug images
    preproc_images = {
        'sub-%03d_t1_after_rotating_and_scaling%s.png',          'T1 rotation & scaling';
        'sub-%03d_segmented_after_rotating_and_scaling%s.png',   'Segmentation rotation & scaling';
        'sub-%03d_after_rotating_and_scaling_orig%s.png',        'Bone mask rotation & scaling';
        'sub-%03d_t1_skin_skull%s.png',                          'T1 skin/skull overlay';
    };

    html = [html '<div class="image-grid image-grid--narrow">'];
    for i = 1:size(preproc_images, 1)
        img_path = fullfile(debug_dir_preproc, sprintf(preproc_images{i,1}, subject_id, affix));
        img_html = html_utils.embed_image(img_path, preproc_images{i,2}, preproc_images{i,2});
        if ~isempty(img_html)
            html = [html img_html];
        end
    end

    % Per-transducer debug images
    for t = 1:n_trans
        img_path = fullfile(debug_dir_preproc, ...
            sprintf('sub-%03d_%s_segmented_brain_final_T%02d%s.png', subject_id, medium, t, affix));
        img_html = html_utils.embed_image(img_path, sprintf('Final segmentation T%02d', t), ...
            sprintf('Final segmentation — T%02d', t));
        if ~isempty(img_html)
            html = [html img_html];
        end
    end

    % Attenuation fit (medium stage)
    if ~isempty(debug_dir_medium) && isfolder(debug_dir_medium)
        img_path = fullfile(debug_dir_medium, sprintf('attenuation_fit%s.png', affix));
        img_html = html_utils.embed_image(img_path, 'Attenuation fit', 'Attenuation fit');
        if ~isempty(img_html)
            html = [html img_html];
        end
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
            html_utils.escape(d.label), d.duration_s, d.peak_ram_gb, d.delta_ram_gb, d.delta_hdd_gb)];
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
        html_utils.escape(log_path))];
    html = [html '<pre class="log-block">' html_utils.escape(log_text) '</pre>'];
end

function html = build_pseudoCT_section(parameters)
    html = '';
    
    % Extract from parameters
    affix = '';
    if isfield(parameters.io, 'output_affix')
        affix = parameters.io.output_affix;
    end
    if isfield(parameters.io, 'debug_dir_medium')
        debug_path = parameters.io.debug_dir_medium;
    else
        debug_path = fullfile(parameters.io.output_dir, 'debug', 'medium');
    end
    
    % 1. MAPPING ALGORITHMS TABLE FIRST
    html = [html, '<dl class="mapping-list">'];
    seg_fields = {'density', 'soundspeed', 'attenuation'};
    labels = {'Density', 'Sound speed', 'Attenuation'};

    for i = 1:3
        if isfield(parameters, 'pct') && isfield(parameters.pct, ['mapping_' seg_fields{i}])
            val = html_utils.escape(char(parameters.pct.(['mapping_' seg_fields{i}])));
            html = [html, sprintf('<dt class="inline-dt">%s mapping: <code>%s</code></dt>', labels{i}, val)];
        end
    end
    html = [html, '</dl>'];
    
    % NEW LINE BEFORE IMAGES
    html = [html, '<div style="height: 20px;"></div>'];
    
    % 2. pCT histograms
    pct_hist_path = fullfile(debug_path, sprintf('pCT_histograms%s.png', affix));
    img_html = html_utils.embed_image(pct_hist_path, 'pCT histograms', sprintf('pCT histograms [%s]', affix));
    if ~isempty(img_html)
        html = [html, img_html, '<div style="height: 15px;"></div>'];
    end

    % 3. KPLAN mapping - ONLY when pct_mapping_density == 'k-plan'
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'mapping_density') && strcmp(parameters.pct.mapping_density, 'k-plan')
        kplan_path = fullfile(debug_path, 'pCT_hounsfield-density_kplan.png');
        img_html = html_utils.embed_image(kplan_path, 'KPLAN mapping', 'Hounsfield → Density (KPLAN)');
        if ~isempty(img_html)
            html = [html, img_html, '<div style="height: 15px;"></div>'];
        end
    end

    % 4. k-Wave mapping - ONLY when pct_mapping_density == 'k-wave'
    if isfield(parameters, 'pct') && isfield(parameters.pct, 'mapping_density') && strcmp(parameters.pct.mapping_density, 'k-wave')
        kwave_path = fullfile(debug_path, 'pCT_hounsfield-density_kwave.png');
        img_html = html_utils.embed_image(kwave_path, 'k-Wave mapping', 'Hounsfield → Density (k-Wave)');
        if ~isempty(img_html)
            html = [html, img_html, '<div style="height: 15px;"></div>'];
        end
    end
    
    % Fallback if no images
    if isempty(strtrim(html)) && ~any(ismember(fields, fieldnames(parameters)))
        html = '<p class="placeholder">No pCT configuration or images found</p>';
    end
    
    html = [html, '<p class="note"><strong>pCT:</strong> Hounsfield → acoustic properties via mappings</p>'];
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
            'CEM43iso_brain', 'CEM43iso_skull', 'CEM43iso_skin', ...
            'endT', 'endT_brain', 'endT_skull', 'endT_skin', ...
            'rise_endT_brain', 'rise_endT_skull', 'rise_endT_skin', ...
            'maxCEM43', 'maxCEM43end', 'maxCEM43iso', 'maxCEM43isoend', ...
            'CEM43_end_brain', 'CEM43_end_skull', 'CEM43_end_skin', ...
            'CEM43iso_end_brain', 'CEM43iso_end_skull', 'CEM43iso_end_skin'};
end

function cols = get_acoustic_columns_water()
    cols = {'subject_id', 'freq_Hz', 'Isppa', 'Isppa_after_exitplane', ...
            'Psptp', 'Ptp_target', 'real_focal_distance_mm', ...
            'Ipa_target', 'Ipa_target_radius'};
end

function cols = get_thermal_columns_water()
    cols = {'maxT', 'endT', 'maxCEM43', 'maxCEM43end', 'maxCEM43iso', 'maxCEM43isoend'};
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

function html = build_toc(parameters, is_layered)
% Build a sticky table-of-contents navigation bar.
    html = '<nav class="toc" id="toc">';

    % Safety status dot: determine worst-case color
    worst = 'green';
    try
        if isfield(parameters.io, 'filename_output_table') && isfile(parameters.io.filename_output_table)
            csv_table = readtable(parameters.io.filename_output_table, 'VariableNamingRule', 'preserve');
            if is_layered
                limits = get_risk_limits(is_layered);
            else
            end
            metric_names = fieldnames(limits);
            for i = 1:length(metric_names)
                name = metric_names{i};
                info = limits.(name);
                if ~isempty(csv_table) && ismember(name, csv_table.Properties.VariableNames)
                    val = csv_table.(name);
                    if isnumeric(val) && ~isempty(val)
                        c = risk_color(val(end), info.limit);
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
        'methods',     'Methods';
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
        html_utils.escape(label), html_utils.escape(char(string(value))))];
end

%% ========================================================================
%  CSS STYLES
%  ========================================================================

function css = css_styles()
% Simulation-report CSS: base rules plus simulation-specific additions.
css = [css_styles_base() ...
newline ...
'/* Simulation-specific: tables */' newline ...
'.info-table, .config-table { border-collapse: collapse; }' newline ...
'.info-table th, .config-table th { text-align: left; padding: 4px 16px 4px 0; color: #64748b; font-weight: 500; }' newline ...
'.info-table td, .config-table td { padding: 4px 0; }' newline ...
'.table-wrapper { -webkit-overflow-scrolling: touch; }' newline ...
'.data-table tbody tr:nth-child(even) { background: #f8fafc; }' newline ...
'.cell-green { background: #f0fdf4 !important; }' newline ...
'.cell-amber { background: #fffbeb !important; }' newline ...
'.cell-red   { background: #fef2f2 !important; font-weight: 600; }' newline ...
'.cell-gray  { background: #f9fafb !important; color: #94a3b8; }' newline ...
newline ...
'/* Simulation-specific: image grid variants */' newline ...
'.image-grid--wide { grid-template-columns: repeat(auto-fill, minmax(min(600px, 100%), 1fr)); }' newline ...
'.image-grid--narrow { grid-template-columns: repeat(auto-fill, minmax(min(280px, 100%), 1fr)); }' newline ...
newline ...
'/* Badges */' newline ...
'.badge { display: inline-block; padding: 4px 12px; border-radius: 12px; font-size: 0.85em; font-weight: 500; }' newline ...
'.badge-green { background: #f0fdf4; color: #166534; border: 1px solid #bbf7d0; }' newline ...
'.badge-gray  { background: #f9fafb; color: #6b7280; border: 1px solid #e5e7eb; }' newline ...
'.badge-blue  { background: #eff6ff; color: #1e40af; border: 1px solid #bfdbfe; }' newline ...
'.badge-amber { background: #fffbeb; color: #92400e; border: 1px solid #fde68a; }' newline ...
newline ...
'/* Summary cards */' newline ...
'.summary-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr)); gap: 12px; }' newline ...
'.summary-card { background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 14px; text-align: center; }' newline ...
'.summary-label { font-size: 0.8em; color: #64748b; text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 4px; }' newline ...
'.summary-value { font-size: 1.4em; font-weight: 700; color: #1e293b; }' newline ...
'.summary-unit { font-size: 0.75em; color: #94a3b8; margin-top: 2px; }' newline ...
newline ...
'/* Simulation-specific: safety footnote links */' newline ...
'.safety-footnote a { color: #64748b; text-decoration: none; border-bottom: 1px dotted #94a3b8; }' newline ...
'.safety-footnote a:hover { color: #1e293b; border-bottom-style: solid; }' newline ...
newline ...
'/* Log & misc */' newline ...
'.log-block { background: #1e293b; color: #e2e8f0; padding: 16px; border-radius: 6px; ' ...
  'overflow-x: auto; font-size: 0.8em; line-height: 1.5; max-height: 600px; overflow-y: auto; white-space: pre-wrap; word-wrap: break-word; }' newline ...
'.boilerplate-license { margin-top: 20px; padding: 12px 16px; background: #f8fafc; ' ...
  'border-left: 3px solid #94a3b8; border-radius: 4px; font-size: 0.82em; color: #475569; }' newline ...
'.timing-table { margin: 12px 0; font-size: 0.9em; }' newline ...
'.timing-table td:first-child { white-space: nowrap; }' newline ...
newline ...
'/* TOC status dot */' newline ...
'.toc a { transition: background 0.15s; }' newline ...
'.toc a:hover { color: #1e293b; }' newline ...
'.toc-status { width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }' newline ...
newline ...
'/* Simulation-specific: simple safety bars (no markers) */' newline ...
'.safety-bar-track { height: 4px; background: #e2e8f0; border-radius: 2px; margin-top: 6px; overflow: hidden; }' newline ...
'.safety-bar { height: 100%; border-radius: 2px; transition: width 0.3s; }' newline ...
'.safety-green .safety-bar { background: var(--green); }' newline ...
'.safety-amber .safety-bar { background: var(--amber); }' newline ...
'.safety-red .safety-bar { background: var(--red); }' newline ...
newline ...
'/* Simulation-specific: lightbox cursor override */' newline ...
'.lb-overlay img { cursor: default; }' newline ...
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

