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
    parameters = resolve_io_dirs(parameters);

    %% Determine output path
    subject_id = parameters.subject_id;
    medium = parameters.simulation.medium;
    is_layered = contains(medium, {'layered', 'phantom'});
    affix = '';
    if isfield(parameters.io, 'output_affix')
        affix = parameters.io.output_affix;
    end

    report_filename = sprintf('sub-%03d_%s_report%s.html', subject_id, medium, affix);
    report_path = fullfile(parameters.io.dir_reports, report_filename);

    %% Load CSV data if available
    csv_table = [];
    if isfield(parameters.io, 'filename_table') && isfile(parameters.io.filename_table)
        try
            csv_table = readtable(parameters.io.filename_table, 'VariableNamingRule', 'preserve');
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

    % Section 1: Exposure Dashboard (always open, no toggle)
    try
        html_parts{end+1} = build_safety_dashboard(csv_table, parameters, is_layered);
    catch ME
        html_parts{end+1} = html_utils.section_error('Exposure Dashboard', ME);
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
            build_methods_boilerplate(parameters, is_layered, csv_table), false, 'methods');
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
        warn('generate_simulation_report:fileOpen', 'Cannot open %s for writing.', report_path);
        return
    end
    fprintf(fid, '%s\n', html_parts{:});
    fclose(fid);

    fprintf('HTML report saved to: %s\n', report_path);

catch ME
    warn('generate_simulation_report:failed', ...
        'Report generation failed: %s\n%s', ME.message, getReport(ME, 'extended'));
    report_path = '';
end

end

%% ========================================================================
%  SAFETY DASHBOARD
%  ========================================================================

function html = build_safety_dashboard(csv_table, parameters, is_layered)
    limits = get_risk_limits(is_layered);
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
    html = [html '<h2>Exposure Dashboard</h2>'];

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
    html = [html '<p class="safety-footnote">NSR limits per ITRUSST consensus ' ...
            '(Aubry et al., 2025; ' ...
            '<a href="https://doi.org/10.1016/j.brs.2025.10.007">DOI: 10.1016/j.brs.2025.10.007</a>). ' ...
            'For thermal metrics, NSR may be indicated by either of multiple indices (temperature rise, absolute temperature, CEM43, or exposure duration). ' ...
            'A red tile for any individual metric does therefore not indicate that NSR cannot be assumed. ' ...
            'Always refer to current community guidelines when interpreting values and allow for uncertainty about simulated outcomes. ' ...
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
    prestus_v = get_prestus_version(parameters);
    if ~isempty(prestus_v)
        html = [html sprintf('<tr><th>PRESTUS version</th><td>%s</td></tr>', html_utils.escape(prestus_v))];
    end
    kwave_v = get_kwave_version(parameters);
    if ~isempty(kwave_v)
        html = [html sprintf('<tr><th>k-Wave version</th><td>%s</td></tr>', html_utils.escape(kwave_v))];
    end
    if isfield(parameters, 'output_dir')
        html = [html sprintf('<tr><th>Output dir</th><td>%s</td></tr>', html_utils.escape(parameters.io.dir_output))];
    end
    html = [html '</table>'];
    html = [html '</section>'];
end

function html = build_methods_boilerplate(parameters, is_layered, csv_table)
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
    function s = detect_mri_sequence(params)
        % Infer MRI sequence label from path field names in parameters.paths.
        % Returns a string like 'T1-weighted' or 'T2-weighted'.
        s = 'MRI';  % safe fallback
        if ~isfield(params, 'paths'); return; end
        flds = fieldnames(params.paths);
        has_t2 = any(~cellfun(@isempty, regexpi(flds, 'T2', 'once')));
        has_t1 = any(~cellfun(@isempty, regexpi(flds, 'T1', 'once')));
        if has_t2 && ~has_t1
            s = 'T2-weighted MRI';
        elseif has_t1
            s = 'T1-weighted MRI';
        end
    end
    function s = fmt(val, fmtstr, fallback)
        if isnumeric(val) && isscalar(val) && ~isnan(val)
            s = sprintf(fmtstr, val);
        else
            s = fallback;
        end
    end
    function s = build_timing_sentence(Nt, dt_us, dur_us, ppw, cfl, ref_sos, ...
                                        gx_mm, gy_mm, gz_mm, grid_res_mm)
        % Compose the temporal-discretisation sentence for the methods block.
        % Any unavailable (NaN) value is omitted gracefully.
        parts = {};

        % Grid domain extent (mm)
        if ~isnan(gx_mm) && ~isnan(gy_mm) && ~isnan(gz_mm)
            parts{end+1} = sprintf('total grid extent: %.0f&#8201;&#215;&#8201;%.0f&#8201;&#215;&#8201;%.0f&#8201;mm', ...
                gx_mm, gy_mm, gz_mm);
        end

        % Spatial sampling
        if ~isnan(ppw)
            if ~isnan(grid_res_mm)
                parts{end+1} = sprintf('%.1f points per wavelength at %.2g&#8201;mm resolution', ...
                    ppw, grid_res_mm);
            else
                parts{end+1} = sprintf('%.1f points per wavelength', ppw);
            end
        end

        % Temporal discretisation
        if ~isnan(Nt) && ~isnan(dt_us)
            parts{end+1} = sprintf('%d time steps of %.2f&#8201;&#956;s', Nt, dt_us);
        end
        if ~isnan(dur_us)
            parts{end+1} = sprintf('total simulation duration %.1f&#8201;&#956;s', dur_us);
        end

        % CFL and reference sound speed
        if ~isnan(cfl)
            if ~isnan(ref_sos)
                parts{end+1} = sprintf('CFL&#8201;=&#8201;%.2f (reference sound speed: %.0f&#8201;m/s)', ...
                    cfl, ref_sos);
            else
                parts{end+1} = sprintf('CFL&#8201;=&#8201;%.2f', cfl);
            end
        end

        if isempty(parts)
            s = '';
        else
            s = ['Simulation parameters: ' strjoin(parts, '; ') '. '];
        end
    end
    function v = csv_num(tbl, col)
        % Read a scalar numeric from csv_table; return NaN if unavailable.
        v = NaN;
        if ~isempty(tbl) && istable(tbl) && ismember(col, tbl.Properties.VariableNames)
            raw = tbl.(col);
            if isnumeric(raw) && ~isempty(raw)
                v = raw(1);
            end
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

    % Grid timing & discretisation (written by acoustic_analysis into the CSV)
    sim_Nt      = csv_num(csv_table, 'sim_Nt');
    sim_dt_us   = csv_num(csv_table, 'sim_dt_us');
    sim_dur_us  = csv_num(csv_table, 'sim_duration_us');
    sim_ppw     = csv_num(csv_table, 'sim_ppw');
    sim_cfl     = csv_num(csv_table, 'sim_cfl');
    sim_ref_sos = csv_num(csv_table, 'sim_ref_sound_speed_mps');
    sim_gx_mm   = csv_num(csv_table, 'sim_grid_x_mm');
    sim_gy_mm   = csv_num(csv_table, 'sim_grid_y_mm');
    sim_gz_mm   = csv_num(csv_table, 'sim_grid_z_mm');

    % Transducer (first element)
    td       = parameters.transducer(1);
    freq_hz  = safe_field(td, 'freq_hz', NaN);
    freq_khz = fmt(freq_hz / 1e3, '%.0f', '?');
    td_type  = safe_field(td, 'type', '');
    td_name  = safe_field(td, 'name', '');
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

    % Acoustic / thermal backend
    use_babelbrain_acoustic  = isfield(parameters, 'acoustic_backend') && ...
                               strcmpi(parameters.acoustic_backend, 'babelbrain');
    use_babelbrain_thermal   = isfield(parameters, 'thermal_backend') && ...
                               strcmpi(parameters.thermal_backend, 'babelbrain');
    use_jwave_acoustic       = isfield(parameters, 'acoustic_backend') && ...
                               strcmpi(parameters.acoustic_backend, 'jwave');

    jw_solver      = '';
    jw_precision   = '';
    jw_pml         = NaN;
    jw_fd_accuracy = NaN;
    jw_compute     = '';
    if use_jwave_acoustic && isfield(parameters, 'jwave')
        jw = parameters.jwave;
        if isfield(jw, 'solver') && ~isempty(jw.solver)
            jw_solver = jw.solver;
        end
        if isfield(jw, 'precision') && ~isempty(jw.precision)
            jw_precision = jw.precision;
        end
        if isfield(jw, 'pml_size') && ~isempty(jw.pml_size)
            jw_pml = jw.pml_size;
        end
        if isfield(jw, 'fd_accuracy') && ~isempty(jw.fd_accuracy)
            jw_fd_accuracy = jw.fd_accuracy;
        end
        % Resolve compute backend: explicit use_gpu flag or auto-detect hint
        if isfield(jw, 'use_gpu')
            if isequal(jw.use_gpu, true) || isequal(jw.use_gpu, 1)
                jw_compute = 'GPU (CUDA)';
            elseif isequal(jw.use_gpu, false) || isequal(jw.use_gpu, 0)
                jw_compute = 'CPU';
            end
        end
    end

    bb_ppw    = '';
    bb_preset = '';
    bb_gpu    = '';
    if use_babelbrain_acoustic && isfield(parameters, 'babelbrain')
        bb = parameters.babelbrain;
        if isfield(bb, 'ppw') && ~isempty(bb.ppw)
            bb_ppw = fmt(bb.ppw, '%d', '');
        end
        if isfield(bb, 'tx_system') && ~isempty(bb.tx_system)
            bb_preset = html_utils.escape(bb.tx_system);
        end
        if isfield(bb, 'computing_backend')
            switch bb.computing_backend
                case 1; bb_gpu = 'CUDA';
                case 2; bb_gpu = 'OpenCL';
                case 3; bb_gpu = 'Metal (Apple Silicon)';
                case 4; bb_gpu = 'CPU';
            end
        end
    end

    % Software versions
    prestus_ver  = get_prestus_version(parameters);
    kw_ver       = get_kwave_version(parameters);
    kw_ver_str   = iff(~isempty(kw_ver),      sprintf('v%s ', kw_ver), '');
    bb_ver_str   = '';   % BabelBrain version not yet captured at runtime
    jw_ver_str   = '';   % j-Wave version not yet captured at runtime

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

    if ischar(td_name) && ~isempty(td_name)
        td_name_str = sprintf(' (%s)', html_utils.escape(td_name));
    else
        td_name_str = '';
    end
    if use_babelbrain_acoustic
        preset_str = iff(~isempty(bb_preset), sprintf(', BabelBrain preset %s', bb_preset), '');
        para1 = sprintf([...
            '%s The transducer%s was a %s-element %s array operating at a centre ' ...
            'frequency of %s kHz%s. '                                              ...
            'The transducer geometry was specified via a named BabelBrain preset%s. ' ...
            'Element phases from hydrophone calibration were applied as '           ...
            'corrections relative to geometric (Rayleigh) focusing phases.'],       ...
            b('Transducer and drive system.'),                                      ...
            td_name_str,                                                            ...
            fmt(n_elem, '%d', '?'), html_utils.escape(td_type),                    ...
            freq_khz, geom_str, preset_str);
    elseif use_jwave_acoustic
        para1 = sprintf([...
            '%s The transducer%s was a %s-element %s array operating at a centre ' ...
            'frequency of %s kHz%s. '                                              ...
            'The transducer was modelled within the k-Wave Toolbox '               ...
            '(Treeby &amp; Cox, 2010) using direct grid-point sources '            ...
            '(kWaveArray not supported with the jwave backend).'],                  ...
            b('Transducer and drive system.'),                                      ...
            td_name_str,                                                            ...
            fmt(n_elem, '%d', '?'), html_utils.escape(td_type),                    ...
            freq_khz, geom_str);
    else
        para1 = sprintf([...
            '%s The transducer%s was a %s-element %s array operating at a centre ' ...
            'frequency of %s kHz%s. '                                              ...
            'The transducer was modelled within the k-Wave Toolbox '               ...
            '(Treeby &amp; Cox, 2010) using the kWaveArray discrete source '       ...
            'approach (Wise et al., 2019).'],                                        ...
            b('Transducer and drive system.'),                                      ...
            td_name_str,                                                            ...
            fmt(n_elem, '%d', '?'), html_utils.escape(td_type),                    ...
            freq_khz, geom_str);
    end

    % ------------------------------------------------------------------ %
    %  2. Simulation grid (in situ approach)
    % ------------------------------------------------------------------ %
    if is_layered
        medium_desc = sprintf(['a layered, subject-specific head model '       ...
            'comprising %d tissue compartments (%s)'], n_layers, tissue_list);
        seg_sentence = sprintf([' Tissue segmentation was derived from %s '    ...
            'using the SimNIBS ' em('charm') ' pipeline '                      ...
            '(Thielscher et al., 2015; Puonti et al., 2020).'],                 ...
            detect_mri_sequence(parameters));
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

    if use_babelbrain_acoustic
        gpu_str = iff(~isempty(bb_gpu), sprintf(' GPU acceleration used the %s backend.', bb_gpu), '');
        % BabelBrain has no PRESTUS-side kgrid, so Nt/dt/CFL are NaN.
        % Grid extent and PPW can be derived from parameters directly.
        bb_ppw_num = str2double(bb_ppw);   % bb_ppw is a fmt'd string; convert back
        bb_gx_mm = NaN; bb_gy_mm = NaN; bb_gz_mm = NaN;
        if isfield(parameters, 'grid') && isfield(parameters.grid, 'dims') && ...
                isfield(parameters.grid, 'resolution_mm') && numel(parameters.grid.dims) >= 3
            bb_gx_mm = parameters.grid.dims(1) * parameters.grid.resolution_mm;
            bb_gy_mm = parameters.grid.dims(2) * parameters.grid.resolution_mm;
            bb_gz_mm = parameters.grid.dims(3) * parameters.grid.resolution_mm;
        end
        bb_timing_sentence = build_timing_sentence(NaN, NaN, NaN, ...
            bb_ppw_num, NaN, NaN, bb_gx_mm, bb_gy_mm, bb_gz_mm, grid_res);
        % BabelBrain material properties at simulation frequency
        % All attenuations are linear in frequency (power-law exponent 1).
        % Coefficients from BabelIntegrationBASE.py:
        %   Cortical longitudinal: 203.25 * 0.8 = 162.6 Np/m/MHz
        %   Bone shear: (57/0.27 + 373/0.836)/2 = 328.65 Np/m/MHz
        %   Brain: 3.45/0.5 = 6.9 Np/m/MHz
        % 1 Np/m = 0.08686 dB/cm
        np_to_db_cm = 0.08686;
        freq_mhz    = freq_hz / 1e6;
        cort_long_att_np = 162.6 * freq_mhz;   % Np/m at simulation frequency
        cort_shear_att_np = 328.65 * freq_mhz; % Np/m at simulation frequency
        brain_att_np = 6.9 * freq_mhz;         % Np/m at simulation frequency
        cort_long_att_db = cort_long_att_np * np_to_db_cm;
        cort_shear_att_db = cort_shear_att_np * np_to_db_cm;
        brain_att_db = brain_att_np * np_to_db_cm;

        if ~isnan(freq_hz)
            mat_prop_str = sprintf([...
                'At %s&#8201;kHz, the applied material properties included: '           ...
                'cortical bone longitudinal attenuation %.1f&#8201;dB/cm '              ...
                '(%.0f&#8201;Np/m; coefficient: 162.6&#8201;Np/m/MHz), '                ...
                'cortical/trabecular bone shear attenuation %.1f&#8201;dB/cm '          ...
                '(%.0f&#8201;Np/m; coefficient: 328.7&#8201;Np/m/MHz), '                ...
                'and brain longitudinal attenuation %.2f&#8201;dB/cm '                  ...
                '(%.1f&#8201;Np/m; coefficient: 6.9&#8201;Np/m/MHz). '                  ...
                'All attenuation coefficients scale linearly with frequency.'],          ...
                freq_khz,                                                                ...
                cort_long_att_db,  cort_long_att_np,                                    ...
                cort_shear_att_db, cort_shear_att_np,                                   ...
                brain_att_db,      brain_att_np);
        else
            mat_prop_str = [...
                'Frequency-dependent material properties were assigned '                 ...
                'automatically by BabelBrain from published literature fits '            ...
                '(Pichardo et al., 2011; Pichardo et al., 2017; Goss et al., 1978). '   ...
                'All attenuation coefficients scale linearly with frequency.'];
        end

        para2 = sprintf([...
            '%s '                                                                       ...
            'Acoustic pressure field propagation was modelled in 3D '                  ...
            'using BabelBrain %s(Pichardo, 2023), which solves the full viscoelastic ' ...
            'wave equation via an explicit staggered-grid finite-difference '           ...
            'time-domain method (BabelViscoFDTD), including both longitudinal and '    ...
            'shear wave propagation in bone and mode conversion at the skull. '        ...
            'A Rayleigh integral pre-step was used to propagate the transducer '       ...
            'field through the water coupling layer prior to the skull. '              ...
            '%s'                                                                        ...
            '%s '                                                                       ...
            'The medium was represented as %s.%s%s%s%s '                               ...
            'The method used to obtain ' em('in situ') ' exposure estimates '          ...
            'was full-wave acoustic simulation through the individual skull '           ...
            'and brain geometry.'],                                                     ...
            b(['Acoustic simulation and ' em('in situ') ' exposure estimates.']),       ...
            bb_ver_str,                                                                 ...
            bb_timing_sentence,                                                         ...
            mat_prop_str,                                                               ...
            medium_desc, seg_sentence, skull_sentence, props_sentence, gpu_str);
    elseif use_jwave_acoustic
        solver_str    = iff(~isempty(jw_solver),    html_utils.escape(jw_solver),    'helmholtz');
        precision_str = iff(~isempty(jw_precision), html_utils.escape(jw_precision), 'single');
        % PML: prefer the jwave-specific field; fall back to the global pml variable
        jw_pml_val = iff(~isnan(jw_pml), jw_pml, pml);
        pml_str = fmt(jw_pml_val, '%d', '?');
        solver_desc = 'frequency-domain Helmholtz (GMRES)';
        compute_str = iff(~isempty(jw_compute), sprintf(' Computation was performed on %s.', jw_compute), '');
        % j-Wave Helmholtz is frequency-domain: no time steps to report
        jw_timing_sentence = '';
        jw_window_sentence = '';
        para2 = sprintf([...
            '%s '                                                                        ...
            'Acoustic pressure field propagation was modelled in %dD '                  ...
            'using j-Wave %s(Stanziola et al., 2023), a JAX-based open-source '         ...
            'differentiable wave simulator. '                                            ...
            'The %s solver was used with %s-precision (float%s) arithmetic, '           ...
            'on a %s&#8201;mm isotropic Cartesian grid '                                ...
            '(perfectly matched layer thickness: %s grid points). %s%s%s'               ...
            'The medium was represented as %s.%s%s%s '                                  ...
            'The method used to obtain ' em('in situ') ' exposure estimates '           ...
            'was full-wave acoustic simulation through the individual skull '            ...
            'and brain geometry.'],                                                      ...
            b(['Acoustic simulation and ' em('in situ') ' exposure estimates.']),        ...
            n_dims,                                                                      ...
            jw_ver_str,                                                                  ...
            solver_desc, precision_str,                                                  ...
            iff(strcmpi(precision_str,'single'),'32','64'),                              ...
            fmt(grid_res, '%.2g', '?'),                                                 ...
            pml_str, jw_timing_sentence, jw_window_sentence, compute_str,               ...
            medium_desc, seg_sentence, skull_sentence, props_sentence);
    else
        % k-Wave pseudo-spectral time-domain
        timing_sentence = build_timing_sentence(sim_Nt, sim_dt_us, sim_dur_us, ...
                                                sim_ppw, sim_cfl, sim_ref_sos, ...
                                                sim_gx_mm, sim_gy_mm, sim_gz_mm, ...
                                                grid_res);
        para2 = sprintf([...
            '%s '                                                                   ...
            'Acoustic pressure field propagation was modelled in %dD '             ...
            'using the pseudospectral time-domain method implemented in the '       ...
            'k-Wave Toolbox %s(Treeby &amp; Cox, 2010), '                          ...
            'on a %s&#8201;mm isotropic Cartesian grid '                           ...
            '(perfectly matched layer: %s grid points). '                          ...
            '%s'                                                                    ...
            'Steady-state pressure was recorded over the final 3 wave periods. '   ...
            'The medium was represented as %s.%s%s%s '                             ...
            'The method used to obtain ' em('in situ') ' exposure estimates '      ...
            'was full-wave acoustic simulation through the individual skull '       ...
            'and brain geometry.'],                                                 ...
            b(['Acoustic simulation and ' em('in situ') ' exposure estimates.']),   ...
            n_dims,                                                                 ...
            kw_ver_str,                                                             ...
            fmt(grid_res, '%.2g', '?'),                                            ...
            fmt_pml(pml),                                                          ...
            timing_sentence,                                                       ...
            medium_desc, seg_sentence, skull_sentence, props_sentence);
    end

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
        if use_babelbrain_thermal
            thermal_solver_str = [...
                ' Thermal simulations were performed using BabelBrain&#x2019;s '   ...
                'GPU-accelerated Pennes bioheat equation solver '                   ...
                '(Pichardo, 2023; Pennes, 1948), operating on the same simulation ' ...
                'domain as the BabelBrain acoustic step.'];
        else
            thermal_solver_str = [...
                ' Temperature elevation and cumulative thermal dose (CEM43) '      ...
                'were estimated by solving the Pennes bioheat equation '           ...
                '(Pennes, 1948) over the full sonication protocol described above.'];
        end
        para5 = [b('Thermal simulations.') thermal_solver_str ' ' cem43_formula];
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
    %  References block (expandable) + CC0 notice
    % ------------------------------------------------------------------ %
    ref_flags = struct( ...
        'use_babelbrain_acoustic', use_babelbrain_acoustic, ...
        'use_jwave_acoustic',      use_jwave_acoustic,      ...
        'use_kwave',               ~use_babelbrain_acoustic && ~use_jwave_acoustic, ...
        'is_layered',              is_layered,              ...
        'run_thermal',             run_thermal,             ...
        'use_babelbrain_thermal',  use_babelbrain_thermal,  ...
        'use_pct',                 use_pct);
    refs_html = build_references_block(ref_flags);

    if ~isempty(prestus_ver)
        prestus_ver_str = sprintf(' (PRESTUS %s)', prestus_ver);
    else
        prestus_ver_str = '';
    end
    cc0 = ['<div class="boilerplate-license">'                                 ...
        '<strong>License.</strong> '                                            ...
        sprintf('The above methods text was automatically generated by PRESTUS%s ', prestus_ver_str) ...
        'and conforms to ITRUSST standardised reporting (Martin et al., 2024). '...
        'It has been released under the '                                       ...
        '<a href="https://creativecommons.org/publicdomain/zero/1.0/" '        ...
        'target="_blank" rel="noopener">CC0 1.0 Universal</a> license '        ...
        'and may be copied verbatim into manuscript Methods sections.'          ...
        '</div>'];

    html = [txt refs_html cc0];
end

function html = build_references_block(f)
% Build an expandable <details> block listing all references cited in the
% methods boilerplate, showing only those relevant to the current run.

    % Each entry: {tag, full_text_html}
    % Always-present references
    refs = {};

    % --- Acoustic solver ---
    if f.use_babelbrain_acoustic
        refs{end+1} = ref( ...
            'Pichardo, 2023', ...
            ['Pichardo, S. (2023). BabelBrain: An Open-Source Application for Prospective Modeling of Transcranial Focused Ultrasound for Neuromodulation Applications. ' ...
            '<em>IEEE Trans Ultrason Ferroelectr Freq Control</em>, <strong>70</strong>, 587&#8211;599. ' ...
            '<a href="https://doi.org/10.1109/TUFFC.2023.3274046" target="_blank" rel="noopener">https://doi.org/10.1109/TUFFC.2023.3274046</a>']);
        refs{end+1} = ref( ...
            'Pichardo et al., 2011', ...
            ['Pichardo, S., Sin, V. W., &amp; Hynynen, K. (2011). Multi-frequency characterization of the speed of sound and attenuation coefficient for longitudinal transmission of freshly excised human skulls. ' ...
            '<em>Phys Med Biol</em>, <strong>56</strong>, 219&#8211;250. ' ...
            '<a href="https://doi.org/10.1088/0031-9155/56/1/014" target="_blank" rel="noopener">https://doi.org/10.1088/0031-9155/56/1/014</a>']);
        refs{end+1} = ref( ...
            'Pichardo et al., 2017', ...
            ['Pichardo, S., Bhatt, C., Huang, Y., &amp; Hynynen, K. (2017). Treatment of the skull in transcranial focused ultrasound. ' ...
            '<em>Phys Med Biol</em>, <strong>62</strong>, 7834&#8211;7857. ' ...
            '<a href="https://doi.org/10.1088/1361-6560/aa8d6c" target="_blank" rel="noopener">https://doi.org/10.1088/1361-6560/aa8d6c</a>']);
        refs{end+1} = ref( ...
            'Goss et al., 1978', ...
            ['Goss, S. A., Johnston, R. L., &amp; Dunn, F. (1978). Comprehensive compilation of empirical ultrasonic properties of mammalian tissues. ' ...
            '<em>J Acoust Soc Am</em>, <strong>64</strong>, 423&#8211;457. ' ...
            '<a href="https://doi.org/10.1121/1.382016" target="_blank" rel="noopener">https://doi.org/10.1121/1.382016</a>']);
    elseif f.use_jwave_acoustic
        refs{end+1} = ref( ...
            'Stanziola et al., 2023', ...
            ['Stanziola, A., Arridge, S. R., Cox, B. T., &amp; Treeby, B. E. (2023). j-Wave: An open-source differentiable wave simulator. ' ...
            '<em>SoftwareX</em>, <strong>22</strong>, 101338. ' ...
            '<a href="https://doi.org/10.1016/j.softx.2023.101338" target="_blank" rel="noopener">https://doi.org/10.1016/j.softx.2023.101338</a>']);
    end

    % k-Wave (always cited for transducer modelling; also the kwave backend)
    refs{end+1} = ref( ...
        'Treeby &amp; Cox, 2010', ...
        ['Treeby, B. E., &amp; Cox, B. T. (2010). k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave fields. ' ...
        '<em>J Biomed Opt</em>, <strong>15</strong>, 021314. ' ...
        '<a href="https://doi.org/10.1117/1.3360308" target="_blank" rel="noopener">https://doi.org/10.1117/1.3360308</a>']);

    if f.use_kwave
        refs{end+1} = ref( ...
            'Wise et al., 2019', ...
            ['Wise, E. S., Cox, B. T., Jaros, J., &amp; Treeby, B. E. (2019). Representing arbitrary acoustic source and sensor distributions in Fourier collocation methods. ' ...
            '<em>J Acoust Soc Am</em>, <strong>146</strong>, 278&#8211;288. ' ...
            '<a href="https://doi.org/10.1121/1.5116132" target="_blank" rel="noopener">https://doi.org/10.1121/1.5116132</a>']);
    end

    % Segmentation
    if f.is_layered
        refs{end+1} = ref( ...
            'Thielscher et al., 2015', ...
            ['Thielscher, A., Antunes, A., &amp; Saturnino, G. B. (2015). Field modeling for transcranial magnetic stimulation: A useful tool to understand the physiological effects of TMS? ' ...
            '<em>Proc IEEE EMBC</em>, 222&#8211;225. ' ...
            '<a href="https://doi.org/10.1109/EMBC.2015.7318340" target="_blank" rel="noopener">https://doi.org/10.1109/EMBC.2015.7318340</a>']);
        refs{end+1} = ref( ...
            'Puonti et al., 2020', ...
            ['Puonti, O., Van Leemput, K., Saturnino, G. B., Siebner, H. R., Madsen, K. H., &amp; Thielscher, A. (2020). Reproducible and automated T1-only segmentation of the brain and skull using the SimNIBS <em>charm</em> pipeline. ' ...
            '<em>NeuroImage</em>, <strong>204</strong>, 116237. ' ...
            '<a href="https://doi.org/10.1016/j.neuroimage.2019.116237" target="_blank" rel="noopener">https://doi.org/10.1016/j.neuroimage.2019.116237</a>']);
    end

    % Thermal
    if f.run_thermal
        refs{end+1} = ref( ...
            'Pennes, 1948', ...
            ['Pennes, H. H. (1948). Analysis of tissue and arterial blood temperatures in the resting human forearm. ' ...
            '<em>J Appl Physiol</em>, <strong>1</strong>, 93&#8211;122. ' ...
            '<a href="https://doi.org/10.1152/jappl.1948.1.2.93" target="_blank" rel="noopener">https://doi.org/10.1152/jappl.1948.1.2.93</a>']);
        refs{end+1} = ref( ...
            'Sapareto &amp; Dewey, 1984', ...
            ['Sapareto, S. A., &amp; Dewey, W. C. (1984). Thermal dose determination in cancer therapy. ' ...
            '<em>Int J Radiat Oncol Biol Phys</em>, <strong>10</strong>, 787&#8211;800. ' ...
            '<a href="https://doi.org/10.1016/0360-3016(84)90379-1" target="_blank" rel="noopener">https://doi.org/10.1016/0360-3016(84)90379-1</a>']);
    end

    % ITRUSST (always cited)
    refs{end+1} = ref( ...
        'Martin et al., 2024', ...
        ['Martin, E., Treeby, B. E., Bhatt, M., et al. (2024). Recommendations for standardised reporting of focused ultrasound neuromodulation studies. ' ...
        '<em>Brain Stimul</em>, <strong>17</strong>, 140&#8211;151. ' ...
        '<a href="https://doi.org/10.1016/j.brs.2024.01.007" target="_blank" rel="noopener">https://doi.org/10.1016/j.brs.2024.01.007</a>']);
    refs{end+1} = ref( ...
        'Aubry et al., 2025', ...
        ['Aubry, J.-F., et al. (2025). ITRUSST consensus on biophysical safety for low-intensity transcranial ultrasound stimulation. ' ...
        '<em>Brain Stimul</em>. ' ...
        '<a href="https://doi.org/10.1016/j.brs.2024.11.012" target="_blank" rel="noopener">https://doi.org/10.1016/j.brs.2024.11.012</a>']);

    items_html = strjoin(refs, '');
    html = ['<details class="boilerplate-refs">' ...
            '<summary class="boilerplate-refs-toggle">References</summary>' ...
            '<ol class="boilerplate-ref-list">' items_html '</ol>' ...
            '</details>'];
end

function s = ref(tag, text_html)
% Format a single reference list item.
    s = sprintf('<li>%s</li>', text_html);
end


function v = get_prestus_version(parameters)
% Read PRESTUS version from VERSION file next to this function's directory.
% Returns version string (e.g. 'v0.6.1') or empty string if not found.
    v = '';
    try
        this_dir = fileparts(mfilename('fullpath'));   % functions/report/
        repo_root = fullfile(this_dir, '..', '..');    % PRESTUS root
        ver_file  = fullfile(repo_root, 'VERSION');
        if isfile(ver_file)
            raw = strtrim(fileread(ver_file));
            if ~isempty(raw)
                v = raw;
            end
        end
    catch
    end
    % Also allow override via parameters.version for non-standard installs
    if isfield(parameters, 'version') && ischar(parameters.version) && ~isempty(parameters.version)
        v = parameters.version;
    end
end

function v = get_kwave_version(parameters)
% Read k-Wave version from its Contents.m file, resolved via the toolbox path.
% Returns version string (e.g. '1.4') or empty string if not found.
    v = '';
    try
        % Prefer path stored in parameters if set
        kw_root = '';
        if isfield(parameters, 'paths') && isfield(parameters.paths, 'kwave_toolbox')
            kw_root = parameters.paths.kwave_toolbox;
        end
        if isempty(kw_root)
            % Fall back: find kWaveGrid on the MATLAB path
            kw_grid = which('kWaveGrid');
            if ~isempty(kw_grid)
                kw_root = fileparts(kw_grid);
            end
        end
        if ~isempty(kw_root)
            contents = fullfile(kw_root, 'Contents.m');
            if isfile(contents)
                txt = fileread(contents);
                tok = regexp(txt, 'Version\s+([\d.]+)', 'tokens', 'once');
                if ~isempty(tok)
                    v = tok{1};
                end
            end
        end
    catch
    end
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

    use_bb = isfield(parameters, 'acoustic_backend') && ...
             strcmpi(parameters.acoustic_backend, 'babelbrain');

    if use_bb
        html = build_medium_properties_babelbrain(parameters);
        return
    end

    if ~isfield(parameters, 'medium_properties') || ~isfield(parameters, 'layers')
        html = '<p class="placeholder">Missing medium or layers data.</p>';
        return
    end

    med = parameters.medium_properties;
    layer_names = fieldnames(parameters.layers);

    valid_tissues = {};
    tissue_labels = {};
    for i = 1:length(layer_names)
        layer_name = layer_names{i};
        if isfield(med, layer_name)
            valid_tissues{end+1} = layer_name;
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

function html = build_medium_properties_babelbrain(parameters)
% Build the Medium Properties section for BabelBrain runs.
% Reads the Material matrix from the DataForSim.h5 output and converts
% units for display. BabelBrain Material rows: Water(0), Skin(1),
% Cortical(2), Trabecular(3), Brain(4), [WhiteMatter(5), GrayMatter(6), CSF(7)].
% Columns: [Density(kg/m³), LongSoS(m/s), ShearSoS(m/s), LongAtt(Np/m), ShearAtt(Np/m)].
    html = '';

    % Locate DataForSim.h5 — same naming logic as babelbrain_wrapper
    h5_path = '';
    try
        out_dir  = parameters.io.dir_cache;
        subj_id  = sprintf('sub-%03d', parameters.subject_id);
        target   = strrep(parameters.io.output_affix, '_', '');
        if isempty(target), target = 'target'; end
        bb       = safe_field(parameters, 'babelbrain', struct());
        freq_hz  = parameters.transducer(1).freq_hz;
        freq_kHz = round(freq_hz / 1e3);
        ppw      = safe_field(bb, 'ppw', 6);
        h5_name  = sprintf('%s_%dkHz_%dPPW_DataForSim.h5', target, freq_kHz, ppw);
        h5_path  = fullfile(out_dir, subj_id, h5_name);
        if ~isfile(h5_path)
            candidates = dir(fullfile(out_dir, subj_id, '*DataForSim.h5'));
            if ~isempty(candidates)
                h5_path = fullfile(candidates(end).folder, candidates(end).name);
            else
                h5_path = '';
            end
        end
    catch
        h5_path = '';
    end

    mat = [];
    if ~isempty(h5_path) && isfile(h5_path)
        try
            ensure_hdf5_plugin_path('', 'babelbrain');
            mat_raw = h5read(h5_path, '/Material');
            % h5read returns [5 x N] (column-major); transpose to [N x 5]
            if size(mat_raw, 1) == 5 && size(mat_raw, 2) > 1
                mat_raw = mat_raw';
            end
            mat = double(mat_raw);
        catch
        end
    end

    if isempty(mat)
        html = ['<p class="placeholder">BabelBrain material properties not available ' ...
                '(DataForSim.h5 not found or /Material dataset missing).</p>'];
        return
    end

    % BabelBrain material row labels (fixed ordering from BabelIntegrationBASE.py)
    bb_row_labels  = {'Water', 'Skin', 'Skull (cortical)', 'Skull (trabecular)', 'Brain', ...
                      'White matter', 'Gray matter', 'CSF'};
    n_rows = min(size(mat, 1), numel(bb_row_labels));

    % Frequency for attenuation unit conversion (Np/m → dB/cm/MHz)
    freq_mhz = NaN;
    try
        freq_mhz = parameters.transducer(1).freq_hz / 1e6;
    catch
    end

    % CT-derived skull note
    use_ct = isfield(parameters, 'babelbrain') && ...
             isfield(parameters.babelbrain, 'ct_mapping') && ...
             ~isempty(parameters.babelbrain.ct_mapping);
    skull_note = '';
    if use_ct
        skull_note = ' Skull properties are CT-derived (voxel-by-voxel); values shown are the simulation-wide mean.';
    end

    html = '<div class="table-wrapper"><table class="data-table">';
    html = [html '<thead><tr>' ...
        '<th>Tissue</th>' ...
        '<th>Density<br><small>kg/m³</small></th>' ...
        '<th>Long. SoS<br><small>m/s</small></th>' ...
        '<th>Shear SoS<br><small>m/s</small></th>' ...
        '<th>Long. att.<br><small>Np/m</small></th>' ...
        '<th>Long. att.<br><small>dB/cm</small></th>'];
    if ~isnan(freq_mhz)
        html = [html sprintf('<th>Long. att.<br><small>dB/cm/MHz @ %.3g MHz</small></th>', freq_mhz)];
    end
    html = [html '<th>Shear att.<br><small>Np/m</small></th></tr></thead><tbody>'];

    NP_TO_DB_PER_CM = 20 / log(10) / 100;  % Np/m → dB/cm

    for r = 1:n_rows
        row = mat(r, :);
        den   = row(1);  lsos  = row(2);  ssos  = row(3);
        latt  = row(4);  satt  = row(5);
        latt_db = latt * NP_TO_DB_PER_CM;
        html = [html sprintf('<tr><td><strong>%s</strong></td>', ...
            html_utils.escape(bb_row_labels{r}))];
        html = [html sprintf('<td>%.4g</td><td>%.4g</td><td>%.4g</td><td>%.4g</td><td>%.4g</td>', ...
            den, lsos, ssos, latt, latt_db)];
        if ~isnan(freq_mhz) && freq_mhz > 0
            html = [html sprintf('<td>%.4g</td>', latt_db / freq_mhz)];
        end
        html = [html sprintf('<td>%.4g</td></tr>', satt)];
    end
    html = [html '</tbody></table></div>'];

    html = [html sprintf(['<p class="note">Properties as used by BabelBrain (frequency-dependent; ' ...
        'evaluated at %.3g MHz). Longitudinal attenuation in both Np/m and dB/cm for reference.%s</p>'], ...
        freq_mhz, skull_note)];
end

function html = build_config_summary(parameters)
    html = '<table class="config-table">';

    % Simulation type
    html = config_row(html, 'Simulation medium', safe_field(parameters.simulation, 'medium', 'N/A'));
    html = config_row(html, 'Dimensions', sprintf('%dD', safe_field(parameters.grid, 'n_dims', 3)));

    % Grid
    grid_step = safe_field(parameters.grid, 'resolution_mm', NaN);
    if isnumeric(grid_step) && isscalar(grid_step) && ~isnan(grid_step)
        html = config_row(html, 'Grid step', sprintf('%.2f mm', grid_step));
    end
    grid_dims = safe_field(parameters.grid, 'default_dims', []);
    if ~isempty(grid_dims)
        html = config_row(html, 'Grid dims', sprintf('[%s]', strtrim(sprintf('%g ', grid_dims))));
    end
    pml = safe_field(parameters.grid, 'pml_size_effective', ...
              safe_field(parameters.grid, 'pml_size', NaN));
    if isnumeric(pml) && ~any(isnan(pml))
        html = config_row(html, 'PML size', sprintf('[%s]', strtrim(num2str(pml))));
    elseif ischar(pml) || isstring(pml)
        html = config_row(html, 'PML size', char(pml));
    end

    % Transducer (first transducer)
    if isfield(parameters, 'transducer') && ~isempty(parameters.transducer)
        td = parameters.transducer(1);
        td_name_cfg = safe_field(td, 'name', '');
        if ischar(td_name_cfg) && ~isempty(td_name_cfg)
            html = config_row(html, 'Transducer name', td_name_cfg);
        end
        html = config_row(html, 'Frequency', sprintf('%.0f Hz', safe_field(td, 'freq_hz', NaN)));
        if isfield(td, 'type') && isfield(td, td.type)
            html = config_row(html, 'Elements', sprintf('%d', safe_field(td.(td.type), 'elem_n', NaN)));
            html = config_row(html, 'Curvature radius', sprintf('%.0f mm', safe_field(td.(td.type), 'curv_radius_mm', NaN)));
            html = config_row(html, 'Source amplitude', sprintf('%.1f Pa', safe_field(td.(td.type), 'elem_amp', NaN)));
        end
        focal = safe_field(td, 'focal_distance_bowl', NaN);
        if isnumeric(focal) && isscalar(focal) && ~isnan(focal)
            html = config_row(html, 'Expected focal distance', sprintf('%.1f mm', focal));
        end
    end

    % Modules
    modules = {'run_source_setup', 'Source setup'; ...
               'run_acoustic_sims', 'Acoustic sims'; ...
               'run_heating_sims', 'Thermal sims'; ...
               'run_posthoc_water_sims', 'Post-hoc water'};
    for i = 1:size(modules, 1)
        val = safe_field(parameters.modules, modules{i,1}, 0);
        if any(val)
            status = 'Enabled';
        else
            status = 'Disabled';
        end
        html = config_row(html, modules{i,2}, status);
    end

    % Thermal protocol (if enabled)
    if isfield(parameters, 'thermal') && isfield(parameters, 'modules') && any(safe_field(parameters.modules, 'run_heating_sims', 0))
        th = parameters.thermal;
        thermal_fields = {'pd', 'Pulse duration'; 'pri', 'Pulse repetition interval'; ...
                          'ptd', 'Pulse train duration'; 'ptri', 'Pulse train rep. interval'; ...
                          'ptrd', 'Pulse train rep. duration'; 'post_ptri_dur', 'Steady-state duration'};
        for i = 1:size(thermal_fields, 1)
            val = safe_field(th, thermal_fields{i,1}, NaN);
            if isnumeric(val) && isscalar(val) && ~isnan(val)
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
        tnn_suffix = trans_suffix(t, n_trans);
        if n_trans > 1; cap = sprintf('Transducer %d positioning', t); else; cap = 'Transducer positioning'; end
        img_path = fullfile(parameters.io.dir_img, ...
            sprintf('sub-%03d_positioning%s%s.png', subject_id, affix, tnn_suffix));
        img_html = html_utils.embed_image(img_path, sprintf('Positioning%s', tnn_suffix), cap);
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
        t_suffix = trans_suffix(t, n_trans);

        % Intensity on segmentation
        img_path = fullfile(parameters.io.dir_img, ...
            sprintf('sub-%03d_%s_intensity%s%s.png', subject_id, medium, t_suffix, affix));
        img_html = html_utils.embed_image(img_path, sprintf('Intensity on segmentation%s', t_suffix), ...
            sprintf('Intensity overlay (segmentation)%s', t_suffix));
        if ~isempty(img_html)
            html = [html img_html];
        end

        % Intensity on T1w
        img_path = fullfile(parameters.io.dir_img, ...
            sprintf('sub-%03d_%s_intensity_t1w%s%s.png', subject_id, medium, t_suffix, affix));
        img_html = html_utils.embed_image(img_path, sprintf('Intensity on T1w%s', t_suffix), ...
            sprintf('Intensity overlay (T1w)%s', t_suffix));
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

    % Thermal images — {filename_pattern, label}
    % Patterns use new naming: sub-NNN_medium-MEDIUM{affix}_{metric}.png
    td = parameters.io.dir_img;
    thermal_images = {
        fullfile(td, sprintf('sub-%03d_%s_maxT%s.png',            subject_id, medium, affix)), 'Max temperature overlay';
        fullfile(td, sprintf('sub-%03d_%s_thermal_max%s.png',     subject_id, medium, affix)), 'Temperature vs time (layer max)';
        fullfile(td, sprintf('sub-%03d_%s_thermalrise_max%s.png', subject_id, medium, affix)), 'Temperature rise vs time (layer max)';
        fullfile(td, sprintf('sub-%03d_%s_CEM_max%s.png',         subject_id, medium, affix)), 'CEM43 (kWave) vs time (layer max)';
        fullfile(td, sprintf('sub-%03d_%s_CEM_iso_max%s.png',     subject_id, medium, affix)), 'CEM43 (ISO) vs time (layer max)';
        fullfile(td, sprintf('sub-%03d_%s_thermal_protocol%s.png',subject_id, medium, affix)), 'Thermal protocol diagram';
    };

    html = [html '<div class="image-grid">'];
    for i = 1:size(thermal_images, 1)
        img_path = thermal_images{i,1};
        img_html = html_utils.embed_image(img_path, thermal_images{i,2}, thermal_images{i,2});
        if ~isempty(img_html)
            html = [html img_html];
        end
    end
    html = [html '</div>'];

    % Note about heating animation
    avi_path = fullfile(parameters.io.dir_img, ...
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
        debug_dir_preproc = parameters.io.dir_debug_preproc;
    elseif isfield(parameters.io, 'debug_dir')
        debug_dir_preproc = fullfile(parameters.io.dir_debug, 'preproc');
    end
    if isfield(parameters.io, 'debug_dir_medium')
        debug_dir_medium = parameters.io.dir_debug_medium;
    elseif isfield(parameters.io, 'debug_dir')
        debug_dir_medium = fullfile(parameters.io.dir_debug, 'medium');
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
        'sub-%03d_t1w_rotated_scaled%s.png',          'T1w rotation & scaling';
        'sub-%03d_segmented_rotated_scaled%s.png',   'Segmentation rotation & scaling';
        'sub-%03d_rotated_scaled_orig%s.png',        'Bone mask rotation & scaling';
        'sub-%03d_t1w_skin_skull%s.png',             'T1w skin/skull overlay';
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
        seg_trans_suffix = trans_suffix(t, n_trans);
        img_path = fullfile(debug_dir_preproc, ...
            sprintf('sub-%03d_%s_segmented_brain_final%s%s.png', subject_id, medium, seg_trans_suffix, affix));
        img_html = html_utils.embed_image(img_path, sprintf('Final segmentation%s', seg_trans_suffix), ...
            sprintf('Final segmentation%s', seg_trans_suffix));
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
        debug_path = parameters.io.dir_debug_medium;
    else
        debug_path = fullfile(parameters.io.dir_output, 'debug', 'medium');
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
    if ~isfield(parameters, 'io') || ~isfield(parameters.io, 'dir_output') || ~isfolder(parameters.io.dir_output)
        return
    end

    pattern = sprintf('sub-%03d_%s_*%s.txt', subject_id, medium, affix);
    if isfield(parameters.io, 'dir_logs')
        logs_dir = parameters.io.dir_logs;
    else
        logs_dir = fullfile(parameters.io.dir_output, 'log');
    end
    files = dir(fullfile(logs_dir, pattern));
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
        if isfield(parameters.io, 'filename_table') && isfile(parameters.io.filename_table)
            csv_table = readtable(parameters.io.filename_table, 'VariableNamingRule', 'preserve');
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
    html = [html sprintf('<span class="toc-status" style="background:var(--%s);" title="Overall NSR status: %s"></span>', worst, worst)];

    % Always-present links
    links = {
        'safety',      'Exposure';
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
'.boilerplate-refs { margin-top: 16px; border: 1px solid #e2e8f0; border-radius: 6px; overflow: hidden; }' newline ...
'.boilerplate-refs-toggle { display: flex; align-items: center; gap: 6px; padding: 8px 14px; ' ...
  'cursor: pointer; font-size: 0.85em; font-weight: 600; color: #475569; background: #f8fafc; ' ...
  'user-select: none; list-style: none; }' newline ...
'.boilerplate-refs-toggle::before { content: "▶"; font-size: 0.7em; transition: transform 0.15s; }' newline ...
'.boilerplate-refs[open] .boilerplate-refs-toggle::before { transform: rotate(90deg); }' newline ...
'.boilerplate-ref-list { margin: 0; padding: 10px 14px 10px 32px; font-size: 0.82em; ' ...
  'color: #334155; line-height: 1.6; }' newline ...
'.boilerplate-ref-list li { margin-bottom: 6px; }' newline ...
'.boilerplate-ref-list li:last-child { margin-bottom: 0; }' newline ...
'.ref-tag { font-weight: 600; color: #1e293b; }' newline ...
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
  'details { display: block; } details > summary { display: none; } ' ...
  'details > *:not(summary) { display: block !important; } ' ...
  'details > .section-content { padding: 0 24px 24px 24px; } ' ...
  '.boilerplate-refs { border: none; } ' ...
  '.boilerplate-ref-list { display: block !important; padding-left: 1.5em; } ' ...
  '.log-block { max-height: none; } .toc { display: none; } .lb-overlay { display: none; } ' ...
  'figure { break-inside: avoid; } .safety-grid { break-inside: avoid; } ' ...
  'figure img { cursor: default; } ' ...
  '.image-grid { grid-template-columns: repeat(2, 1fr); } .image-grid--wide { grid-template-columns: 1fr; } }' newline ...
];
end

