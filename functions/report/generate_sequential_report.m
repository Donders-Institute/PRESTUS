function report_path = generate_sequential_report(run_params_list, run_labels)
% GENERATE_SEQUENTIAL_REPORT  HTML summary report for a chain of sequential simulations
%
% Produces a self-contained HTML file with:
%   (1) Acoustic comparison table and voxelwise-max intensity NIfTI across all runs
%   (2) Thermal timeseries (temperature and CEM43) concatenated across all runs,
%       drawn as per-layer SVG line plots with run-boundary markers
%
% The voxelwise-max intensity NIfTI is saved alongside the report and also
% referenced in the HTML.
%
% Use as:
%   report_path = generate_sequential_report(run_params_list)
%   report_path = generate_sequential_report(run_params_list, run_labels)
%
% Input:
%   run_params_list - (1×N) cell array of PRESTUS parameters structs,
%                    one per run in chronological order (base run first)
%   run_labels      - (optional) (1×N) cell array of short string labels for each run
%
% Output:
%   report_path - path to the generated HTML file
%
% See also: GENERATE_SIMULATION_REPORT, HTML_UTILS, CSS_STYLES_BASE

arguments
    run_params_list (1,:) cell
    run_labels      (1,:) cell = {}
end

report_path = '';

try
    n_runs = numel(run_params_list);
    if n_runs == 0
        warn('generate_sequential_report:empty', 'run_params_list is empty — nothing to report.');
        return
    end

    % Use first run's parameters for shared metadata
    base_p       = run_params_list{1};
    subject_id   = base_p.subject_id;
    medium       = base_p.simulation.medium;
    is_layered   = contains(medium, {'layered', 'phantom'});

    % Default run labels
    if numel(run_labels) < n_runs
        run_labels = arrayfun(@(i) sprintf('Run %d', i), 1:n_runs, 'UniformOutput', false);
    end

    % Output path: alongside base run reports
    report_filename = sprintf('sub-%03d_%s_sequential_report.html', subject_id, medium);
    report_path     = fullfile(base_p.io.dir_reports, report_filename);

    % ------------------------------------------------------------------ %
    %% Load per-run data
    % ------------------------------------------------------------------ %
    run_data = cell(1, n_runs);
    for ri = 1:n_runs
        p = run_params_list{ri};
        d = struct();

        % CSV table (acoustic + thermal scalars)
        if isfield(p.io, 'filename_table') && isfile(p.io.filename_table)
            try
                d.csv = readtable(p.io.filename_table, 'VariableNamingRule', 'preserve');
            catch
                d.csv = [];
            end
        else
            d.csv = [];
        end

        % Resolve dirs (may be absent if path_log_setup hasn't run yet)
        if isfield(p.io, 'dir_cache')
            dir_cache = p.io.dir_cache;
        elseif isfield(p.io, 'dir_output')
            dir_cache = fullfile(p.io.dir_output, 'cache');
        else
            dir_cache = '';
        end
        if ~isfield(p.io, 'dir_nii_T1w')
            if isfield(p.io, 'dir_output')
                p.io.dir_nii_T1w = fullfile(p.io.dir_output, 'nii');
            else
                p.io.dir_nii_T1w = '';
            end
        end

        % Thermal .mat file — use thermal_cache_affix when available
        if isfield(p.io, 'thermal_cache_affix')
            t_affix = p.io.thermal_cache_affix;
        else
            t_affix = p.io.output_affix;
        end
        heating_mat = fullfile(dir_cache, ...
            sprintf('sub-%03d_%s_heating_res%s.mat', p.subject_id, medium, t_affix));
        d.timeseries = [];
        if isfile(heating_mat)
            try
                tmp = load(heating_mat, 'results_heating');
                if isfield(tmp, 'results_heating') && ...
                   isfield(tmp.results_heating, 'timeseries')
                    d.timeseries = tmp.results_heating.timeseries;
                end
            catch
            end
        end

        d.params    = p;
        d.affix     = p.io.output_affix;
        d.t_affix   = t_affix;
        run_data{ri} = d;
    end

    % ------------------------------------------------------------------ %
    %% Build voxelwise-max intensity NIfTI
    % ------------------------------------------------------------------ %
    max_intensity_nii = build_max_intensity_nii(run_data, base_p, subject_id, medium);

    % ------------------------------------------------------------------ %
    %% Build HTML
    % ------------------------------------------------------------------ %
    html_parts = {};

    html_parts{end+1} = '<!DOCTYPE html>';
    html_parts{end+1} = '<html lang="en">';
    html_parts{end+1} = '<head>';
    html_parts{end+1} = '<meta charset="UTF-8">';
    html_parts{end+1} = '<meta name="viewport" content="width=device-width, initial-scale=1.0">';
    html_parts{end+1} = sprintf('<title>PRESTUS Sequential Report — sub-%03d %s</title>', subject_id, medium);
    html_parts{end+1} = '<style>';
    html_parts{end+1} = css_sequential();
    html_parts{end+1} = '</style>';
    html_parts{end+1} = '</head>';
    html_parts{end+1} = '<body>';

    % TOC
    html_parts{end+1} = build_toc(is_layered);

    % Header
    html_parts{end+1} = build_header_section(subject_id, medium, n_runs, run_labels, base_p);

    % Acoustic overview
    try
        html_parts{end+1} = html_utils.collapsible('Acoustic Overview', ...
            build_acoustic_overview(run_data, run_labels, is_layered, max_intensity_nii, subject_id, medium), ...
            true, 'acoustic');
    catch ME
        html_parts{end+1} = html_utils.section_error('Acoustic Overview', ME);
    end

    % Thermal timeseries
    try
        html_parts{end+1} = html_utils.collapsible('Thermal Timeseries', ...
            build_thermal_timeseries(run_data, run_labels, is_layered, subject_id, medium), ...
            true, 'thermal');
    catch ME
        html_parts{end+1} = html_utils.section_error('Thermal Timeseries', ME);
    end

    % Per-run details
    try
        html_parts{end+1} = html_utils.collapsible('Per-Run Details', ...
            build_per_run_details(run_data, run_labels, subject_id, medium), ...
            false, 'details');
    catch ME
        html_parts{end+1} = html_utils.section_error('Per-Run Details', ME);
    end

    html_parts{end+1} = '<footer>';
    html_parts{end+1} = sprintf('<p>Generated by PRESTUS %s &mdash; %d runs</p>', ...
        datestr(now, 'yyyy-mm-dd HH:MM:SS'), n_runs);
    html_parts{end+1} = '</footer>';
    html_parts{end+1} = html_utils.lightbox();
    html_parts{end+1} = '</body>';
    html_parts{end+1} = '</html>';

    % Write file
    fid = fopen(report_path, 'w', 'n', 'UTF-8');
    if fid == -1
        warn('generate_sequential_report:fileOpen', 'Cannot open %s for writing.', report_path);
        return
    end
    fprintf(fid, '%s\n', html_parts{:});
    fclose(fid);

    fprintf('Sequential HTML report saved to: %s\n', report_path);

catch ME
    warn('generate_sequential_report:failed', ...
        'Sequential report generation failed: %s\n%s', ME.message, getReport(ME, 'extended'));
    report_path = '';
end

end

%% ========================================================================
%  MAX INTENSITY NIFTI
%% ========================================================================

function out_path = build_max_intensity_nii(run_data, base_p, subject_id, medium)
% Load each run's T1w intensity NIfTI and write the voxelwise maximum.
    out_path = '';
    vol_max  = [];
    hdr_ref  = [];

    for ri = 1:numel(run_data)
        p   = run_data{ri}.params;
        nii = fullfile(p.io.dir_nii_T1w, ...
            sprintf('sub-%03d_%s_T1w%s_intensity.nii.gz', ...
                subject_id, medium, p.io.output_affix));
        if ~isfile(nii), continue; end
        try
            vol = single(niftiread(nii));
            if isempty(vol_max)
                vol_max  = vol;
                hdr_ref  = niftiinfo(nii);
            else
                vol_max = max(vol_max, vol);
            end
        catch
        end
    end

    if isempty(vol_max), return; end

    out_path = fullfile(base_p.io.dir_nii_T1w, ...
        sprintf('sub-%03d_%s_T1w_sequential_maxIntensity.nii.gz', subject_id, medium));
    try
        hdr_ref.Filename = out_path;
        hdr_ref.Datatype = 'single';
        niftiwrite(vol_max, strrep(out_path, '.nii.gz', '.nii'), hdr_ref, 'Compressed', true);
        fprintf('Max-intensity NIfTI saved to: %s\n', out_path);
    catch ME
        warn('generate_sequential_report:nifti', 'Could not save max-intensity NIfTI: %s', ME.message);
        out_path = '';
    end
end

%% ========================================================================
%  SECTION BUILDERS
%% ========================================================================

function html = build_header_section(subject_id, medium, n_runs, run_labels, base_p)
    html = '<section class="report-section" id="header">';
    html = [html '<h1>PRESTUS Sequential Simulation Report</h1>'];
    html = [html '<table class="info-table">'];
    html = [html sprintf('<tr><th>Subject</th><td>sub-%03d</td></tr>', subject_id)];
    html = [html sprintf('<tr><th>Medium</th><td>%s</td></tr>', html_utils.escape(medium))];
    html = [html sprintf('<tr><th>Runs</th><td>%d</td></tr>', n_runs)];
    html = [html sprintf('<tr><th>Run labels</th><td>%s</td></tr>', ...
        html_utils.escape(strjoin(run_labels, ' &rarr; ')))];
    html = [html sprintf('<tr><th>Generated</th><td>%s</td></tr>', datestr(now, 'yyyy-mm-dd HH:MM:SS'))];
    if isfield(base_p, 'io') && isfield(base_p.io, 'dir_output')
        html = [html sprintf('<tr><th>Output dir</th><td>%s</td></tr>', html_utils.escape(base_p.io.dir_output))];
    end
    html = [html '</table></section>'];
end

function html = build_acoustic_overview(run_data, run_labels, is_layered, max_nii_path, subject_id, medium)
    html = '';

    % --- Comparison table ---
    acoustic_cols = {'MI_tc', 'Isppa', 'Isppa_brain', 'Isppa_skull', 'Isppa_skin', ...
                     'MI_brain', 'Psptp_brain', 'real_focal_distance_mm', 'Ipa_target'};
    limits = get_risk_limits(is_layered);

    table_html = '<div class="table-wrapper"><table class="data-table"><thead><tr>';
    table_html = [table_html '<th>Run</th>'];
    for ci = 1:numel(acoustic_cols)
        table_html = [table_html sprintf('<th>%s</th>', html_utils.escape(acoustic_cols{ci}))];
    end
    table_html = [table_html '</tr></thead><tbody>'];

    for ri = 1:numel(run_data)
        table_html = [table_html sprintf('<tr><td><strong>%s</strong></td>', ...
            html_utils.escape(run_labels{ri}))];
        csv = run_data{ri}.csv;
        for ci = 1:numel(acoustic_cols)
            col = acoustic_cols{ci};
            val = NaN;
            if ~isempty(csv) && ismember(col, csv.Properties.VariableNames)
                v = csv.(col);
                if isnumeric(v) && ~isempty(v)
                    val = v(end);
                end
            end
            cell_class = '';
            if isfield(limits, col) && ~isnan(val)
                cell_class = sprintf(' class="cell-%s"', risk_color(val, limits.(col).limit));
            end
            if isnan(val)
                table_html = [table_html sprintf('<td%s>—</td>', cell_class)];
            else
                table_html = [table_html sprintf('<td%s>%.4g</td>', cell_class, val)];
            end
        end
        table_html = [table_html '</tr>'];
    end
    table_html = [table_html '</tbody></table></div>'];
    html = [html table_html];

    % --- Max-intensity NIfTI note ---
    if ~isempty(max_nii_path) && isfile(max_nii_path)
        html = [html sprintf(['<p class="note"><strong>Voxelwise-max intensity NIfTI</strong> ' ...
            '(maximum across all runs): <code>%s</code></p>'], html_utils.escape(max_nii_path))];
    else
        html = [html '<p class="note">Voxelwise-max intensity NIfTI not available ' ...
            '(individual run intensity NIfTIs may be missing).</p>'];
    end

    % --- Per-run intensity images ---
    html = [html '<h3>Intensity overlays</h3>'];
    html = [html '<div class="run-grid">'];
    for ri = 1:numel(run_data)
        p     = run_data{ri}.params;
        affix = run_data{ri}.affix;
        img   = fullfile(p.io.dir_img, ...
            sprintf('sub-%03d_%s_intensity_t1%s.png', subject_id, medium, affix));
        if ~isfile(img)
            img = fullfile(p.io.dir_img, ...
                sprintf('sub-%03d_%s_intensity%s.png', subject_id, medium, affix));
        end
        html = [html '<div class="run-col">'];
        html = [html sprintf('<h4>%s</h4>', html_utils.escape(run_labels{ri}))];
        img_html = html_utils.embed_image(img, run_labels{ri}, run_labels{ri});
        if ~isempty(img_html)
            html = [html img_html];
        else
            html = [html '<p class="placeholder">No image</p>'];
        end
        html = [html '</div>'];
    end
    html = [html '</div>'];
end

function html = build_thermal_timeseries(run_data, run_labels, is_layered, subject_id, medium)
    html = '';

    if ~is_layered
        html = '<p class="placeholder">Thermal timeseries only available for layered simulations.</p>';
        return
    end

    % Collect per-layer timeseries across runs
    % Detect layer names from first run that has timeseries data
    layer_names = {};
    for ri = 1:numel(run_data)
        ts = run_data{ri}.timeseries;
        if ~isempty(ts) && isstruct(ts) && isfield(ts, 'T')
            layer_names = fieldnames(ts.T);
            break
        end
    end

    if isempty(layer_names)
        html = '<p class="placeholder">No thermal timeseries data found in any run.</p>';
        return
    end

    % Concatenate per-layer vectors; record run boundary indices
    T_layers      = struct();
    CEM43_layers  = struct();
    CEM43_iso_layers = struct();
    for li = 1:numel(layer_names)
        T_layers.(layer_names{li})         = [];
        CEM43_layers.(layer_names{li})     = [];
        CEM43_iso_layers.(layer_names{li}) = [];
    end
    run_boundaries = [0];  % step index where each new run starts

    for ri = 1:numel(run_data)
        ts = run_data{ri}.timeseries;
        n_this = 0;
        for li = 1:numel(layer_names)
            lname = layer_names{li};
            if ~isempty(ts) && isstruct(ts) && isfield(ts, 'T') && isfield(ts.T, lname)
                vec = ts.T.(lname)(:)';
                T_layers.(lname) = [T_layers.(lname), vec];
                n_this = max(n_this, numel(vec));
            end
            if ~isempty(ts) && isstruct(ts) && isfield(ts, 'CEM43') && isfield(ts.CEM43, lname)
                CEM43_layers.(lname) = [CEM43_layers.(lname), ts.CEM43.(lname)(:)'];
            end
            if ~isempty(ts) && isstruct(ts) && isfield(ts, 'CEM43_iso') && isfield(ts.CEM43_iso, lname)
                CEM43_iso_layers.(lname) = [CEM43_iso_layers.(lname), ts.CEM43_iso.(lname)(:)'];
            end
        end
        run_boundaries(end+1) = run_boundaries(end) + n_this;
    end
    run_boundaries = run_boundaries(2:end-1);  % exclude 0 and final end

    html = [html '<h3>Temperature (max per layer)</h3>'];
    html = [html build_ts_svg(T_layers, layer_names, run_boundaries, run_labels, ...
        'Temp (&#176;C)', 37, @(x) x)];

    html = [html '<h3>CEM43 (k-Wave, max per layer)</h3>'];
    html = [html build_ts_svg(CEM43_layers, layer_names, run_boundaries, run_labels, ...
        'CEM43 (min)', 0, @(x) x)];

    % ISO CEM43 (only if any data)
    has_iso = false;
    for li = 1:numel(layer_names)
        if ~isempty(CEM43_iso_layers.(layer_names{li}))
            has_iso = true; break
        end
    end
    if has_iso
        html = [html '<h3>CEM43<sub>ISO</sub> (max per layer)</h3>'];
        html = [html build_ts_svg(CEM43_iso_layers, layer_names, run_boundaries, run_labels, ...
            'CEM43<sub>ISO</sub> (min)', 0, @(x) x)];
    end
end

function html = build_per_run_details(run_data, run_labels, subject_id, medium)
    html = '';
    thermal_image_types = {'maxT', 'thermal_max', 'CEM_max', 'CEM_iso_max'};

    for ri = 1:numel(run_data)
        p     = run_data{ri}.params;
        affix = run_data{ri}.affix;
        html  = [html sprintf('<h3>%s</h3>', html_utils.escape(run_labels{ri}))];
        html  = [html '<div class="image-grid">'];
        found = false;
        for ti = 1:numel(thermal_image_types)
            img = fullfile(p.io.dir_img, ...
                sprintf('sub-%03d_%s_%s%s.png', subject_id, medium, thermal_image_types{ti}, affix));
            img_html = html_utils.embed_image(img, thermal_image_types{ti}, thermal_image_types{ti});
            if ~isempty(img_html)
                html  = [html img_html];
                found = true;
            end
        end
        if ~found
            html = [html '<p class="placeholder">No thermal images found for this run.</p>'];
        end
        html = [html '</div>'];
    end
end

%% ========================================================================
%  SVG TIMESERIES
%% ========================================================================

function svg_html = build_ts_svg(layers, layer_names, run_boundaries, run_labels, y_label, y_floor, transform_fn)
% Build an inline SVG of per-layer timeseries with run-boundary markers.

    % Palette: brain=blue, skull=amber, skin=green, others=gray shades
    palette = {'#3b82f6','#f59e0b','#10b981','#8b5cf6','#ef4444','#64748b'};

    % Determine total length and y range
    n_pts   = 0;
    y_max_v = y_floor;
    for li = 1:numel(layer_names)
        vec = layers.(layer_names{li});
        if ~isempty(vec)
            n_pts   = max(n_pts, numel(vec));
            finite  = vec(isfinite(vec));
            if ~isempty(finite)
                y_max_v = max(y_max_v, max(finite));
            end
        end
    end

    if n_pts == 0
        svg_html = '<p class="placeholder">No data.</p>';
        return
    end

    W = 700; H = 280;
    pad_l = 62; pad_r = 20; pad_t = 20; pad_b = 50;
    plot_w = W - pad_l - pad_r;
    plot_h = H - pad_t - pad_b;

    y_min   = y_floor;
    y_range = max(y_max_v * 1.05, y_min + 0.01) - y_min;
    x_scale = plot_w / max(n_pts - 1, 1);
    y_scale = plot_h / y_range;

    svg = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ' num2str(W) ' ' num2str(H) '" ' ...
           'style="width:100%;max-width:' num2str(W) 'px;height:auto;display:block;margin:0 auto 16px;">'];
    svg = [svg sprintf('<rect width="%d" height="%d" fill="white" rx="4" ry="4"/>', W, H)];

    % Grid lines + Y ticks
    n_ticks = 5;
    for k = 0:n_ticks
        val = y_min + y_range * k / n_ticks;
        yp  = pad_t + plot_h - (val - y_min) * y_scale;
        svg = [svg sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#e2e8f0" stroke-width="1"/>', ...
            pad_l, yp, W-pad_r, yp)];
        svg = [svg sprintf('<text x="%d" y="%.1f" text-anchor="end" font-size="10" fill="#94a3b8">%.3g</text>', ...
            pad_l-4, yp+4, transform_fn(val))];
    end

    % Run-boundary vertical lines
    for bi = 1:numel(run_boundaries)
        bx = pad_l + (run_boundaries(bi) - 0.5) * x_scale;
        svg = [svg sprintf(['<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" ' ...
            'stroke="#94a3b8" stroke-width="1" stroke-dasharray="4,3"/>'], ...
            bx, pad_t, bx, pad_t+plot_h)];
        if bi <= numel(run_labels) + 1
            lbl = '';
            if bi < numel(run_labels)
                lbl = run_labels{bi+1};
            end
            if ~isempty(lbl)
                svg = [svg sprintf(['<text x="%.1f" y="%d" text-anchor="middle" ' ...
                    'font-size="9" fill="#94a3b8">%s</text>'], ...
                    bx, pad_t + plot_h + 24, html_utils.escape(lbl))];
            end
        end
    end

    % Axes
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#cbd5e1" stroke-width="1.5"/>', ...
        pad_l, pad_t, pad_l, pad_t+plot_h)];
    svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#cbd5e1" stroke-width="1.5"/>', ...
        pad_l, pad_t+plot_h, pad_l+plot_w, pad_t+plot_h)];

    % Axis labels
    svg = [svg sprintf('<text x="%.1f" y="%d" text-anchor="middle" font-size="11" fill="#64748b">Time step</text>', ...
        pad_l + plot_w/2, H-4)];
    svg = [svg sprintf('<text x="%d" y="%.1f" text-anchor="middle" font-size="11" fill="#64748b" ', ...
        12, pad_t + plot_h/2) ...
        'transform="rotate(-90,' num2str(12) ',' num2str(pad_t+plot_h/2) ')">' y_label '</text>'];

    % Per-layer lines + legend
    legend_y = pad_t + 10;
    for li = 1:numel(layer_names)
        lname = layer_names{li};
        vec   = layers.(lname);
        if isempty(vec), continue; end
        color = palette{mod(li-1, numel(palette))+1};

        pts = '';
        for i = 1:numel(vec)
            val = vec(i);
            if ~isfinite(val), continue; end
            xp = pad_l + (i-1) * x_scale;
            yp = pad_t + plot_h - (val - y_min) * y_scale;
            pts = [pts sprintf('%.1f,%.1f ', xp, yp)];
        end
        if ~isempty(strtrim(pts))
            svg = [svg sprintf('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.8" stroke-linejoin="round" stroke-linecap="round"/>', ...
                pts, color)];
        end

        % Legend entry
        lx = pad_l + plot_w - 120;
        svg = [svg sprintf('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2"/>', ...
            lx, legend_y + (li-1)*14 + 6, lx+18, legend_y + (li-1)*14 + 6, color)];
        svg = [svg sprintf('<text x="%d" y="%d" font-size="10" fill="#475569">%s</text>', ...
            lx+22, legend_y + (li-1)*14 + 10, html_utils.escape(strrep(lname, '_', ' ')))];
    end

    svg = [svg '</svg>'];
    svg_html = svg;
end

%% ========================================================================
%  TOC
%% ========================================================================

function html = build_toc(is_layered)
    html = '<nav class="toc" id="toc">';
    html = [html '<a href="#header">Overview</a>'];
    html = [html '<a href="#acoustic">Acoustic</a>'];
    if is_layered
        html = [html '<a href="#thermal">Thermal</a>'];
    end
    html = [html '<a href="#details">Per-run</a>'];
    html = [html '</nav>'];
end

%% ========================================================================
%  CSS
%% ========================================================================

function css = css_sequential()
css = [css_styles_base() newline ...
'.info-table th { text-align:left; padding:4px 16px 4px 0; color:#64748b; font-weight:500; }' newline ...
'.info-table td { padding:4px 0; }' newline ...
'.table-wrapper { overflow-x:auto; -webkit-overflow-scrolling:touch; }' newline ...
'.data-table tbody tr:nth-child(even) { background:#f8fafc; }' newline ...
'.cell-green { background:#f0fdf4 !important; }' newline ...
'.cell-amber { background:#fffbeb !important; }' newline ...
'.cell-red   { background:#fef2f2 !important; font-weight:600; }' newline ...
'.cell-gray  { background:#f9fafb !important; color:#94a3b8; }' newline ...
'.run-grid   { display:grid; grid-template-columns:repeat(auto-fill,minmax(min(280px,100%),1fr)); gap:16px; margin:12px 0; }' newline ...
'.run-col    { border:1px solid #e2e8f0; border-radius:6px; padding:12px; background:#fafafa; }' newline ...
'.run-col h4 { margin:0 0 8px; font-size:0.95em; color:#475569; }' newline ...
];
end
