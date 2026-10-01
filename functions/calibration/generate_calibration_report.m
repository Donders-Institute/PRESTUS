function html_path = generate_calibration_report(combo_name, library_path, output_path, equip_param)
% GENERATE_CALIBRATION_REPORT  Write an HTML summary of a transducer calibration
%
% Reads the library YAML for a given equipment combination and generates a
% self-contained HTML report: provenance, model and fit quality
% (meta.quality), the settings the entry resolves to at its fitted depths
% (LIBRARY_SETTINGS_TABLE: phases, drive, gains, model focus, as a study gets
% them), a static phase correction if any, and the amplitude lookup.
%
% Use as:
%   html_path = generate_calibration_report(combo_name)
%   html_path = generate_calibration_report(combo_name, library_path)
%   html_path = generate_calibration_report(combo_name, library_path, output_path)
%   html_path = generate_calibration_report(combo_name, library_path, output_path, equip_param)
%
% Input:
%   combo_name   - equipment combination key, e.g.
%                  'IS_PCD15473_01001_IGT_32_ch_comb_10_ch'
%   library_path - (optional) path to calibration library folder;
%                  defaults to config/calibration/ under PRESTUS root
%   output_path  - (optional) file path for the output HTML;
%                  defaults to <library_path>/<combo_name>_report.html
%   equip_param  - (optional) struct from load_equipment_config; loaded if
%                  missing (needed for the resolved settings table)
%
% Output:
%   html_path - absolute path to the written HTML file
%
% See also: LIBRARY_SETTINGS_TABLE, PLOT_TRANSDUCER_LIBRARY, LOAD_TRANSDUCER_FROM_LIBRARY

    if nargin < 2 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'calibration');
    end

    yaml_path = fullfile(library_path, [combo_name, '.yaml']);
    if ~isfile(yaml_path)
        error('generate_calibration_report: library file not found:\n  %s', yaml_path);
    end
    lib = yaml.loadFile(yaml_path, 'ConvertToArray', true);

    if nargin < 3 || isempty(output_path)
        output_path = fullfile(library_path, [combo_name, '_report.html']);
    end
    if nargin < 4 || isempty(equip_param)
        equip_param = load_equipment_config();
    end

    %% ── Gather provenance ────────────────────────────────────────────────
    meta = lib.meta;
    created_at      = get_field(meta, 'created_at',      '—');
    prestus_version = get_field(meta, 'prestus_version', '—');
    tran_serial     = get_field(meta, 'tran_serial',     combo_name);

    cal_type  = get_field(meta, 'calibration_type', '—');
    fwd_model = '—'; lambda_str = '—'; n_depths = '—';
    lib_mode  = get_field(meta, 'library_mode', '—');
    blend_str = '—';
    if isfield(lib, 'global_model'); fwd_model = get_field(lib.global_model, 'forward_model', '—'); end
    if isfield(lib, 'refinement')
        n_depths  = num2str(numel(lib.refinement.depths_ep_mm));
        blend_str = sprintf('%s mm', get_field(lib.refinement, 'blend_mm', '—'));
    end
    if isfield(meta, 'calibration')
        mc = meta.calibration;
        cal_type  = get_field(mc, 'calibration_type', get_field(mc, 'phase_mode', '—'));
        fwd_model = get_field(mc, 'forward_model', '—');
        if isfield(mc, 'lambda')
            lambda_str = num2str(mc.lambda);
        end
        if isfield(mc, 'n_depths_used')
            n_depths = num2str(mc.n_depths_used);
        end
    end

    %% ── Phase correction table ──────────────────────────────────────────
    correction_html = '<p><em>No per-element phase correction in library.</em></p>';
    if isfield(lib, 'global_model') && isfield(lib.global_model, 'elem_phase_correction_deg') && ...
            ~isempty(lib.global_model.elem_phase_correction_deg)
        corr = lib.global_model.elem_phase_correction_deg(:)';
        rows = '';
        for k = 1:numel(corr)
            rows = [rows, sprintf('<tr><td>%d</td><td>%.2f</td></tr>\n', k, corr(k))]; %#ok<AGROW>
        end
        correction_html = ['<table class="data">', ...
            '<tr><th>Element</th><th>Correction (deg)</th></tr>', rows, '</table>'];
    end

    %% ── Depth × intensity amplitude table ────────────────────────────────
    amp_html = '<p><em>No per-depth amplitude data found.</em></p>';
    if isfield(lib, 'global_model') && isfield(lib.global_model, 'amplitude_scaling') && ...
            ~isempty(lib.global_model.amplitude_scaling)
        amp_map   = lib.global_model.amplitude_scaling;
        int_keys  = fieldnames(amp_map);
        int_vals  = cellfun(@(k) str2double(strrep(strrep(k,'i',''),'p','.')), int_keys);
        [int_vals, si] = sort(int_vals);
        int_keys = int_keys(si);
        header = '<tr><th>Intensity (W/cm²)</th><th>elem_p_amp (Pa)</th></tr>';
        rows = '';
        for k = 1:numel(int_keys)
            rows = [rows, sprintf('<tr><td>%.1f</td><td>%d</td></tr>\n', ...
                int_vals(k), amp_map.(int_keys{k}))]; %#ok<AGROW>
        end
        amp_html = ['<table class="data">', header, rows, '</table>'];
    elseif isfield(lib, 'calibration') && isfield(lib.calibration, 'focal_depths')
        cal_depths = lib.calibration.focal_depths;
        dkeys = fieldnames(cal_depths);
        header = '<tr><th>Focal depth (mm)</th><th>Intensity (W/cm²)</th><th>elem_p_amp (Pa)</th><th>Refinement?</th></tr>';
        rows = '';
        for di = 1:numel(dkeys)
            entry = cal_depths.(dkeys{di});
            depth_val = str2double(strrep(strrep(dkeys{di},'f',''),'p','.'));
            is_ref = get_field(entry, 'is_refinement', false);
            if isfield(entry, 'amplitude_scaling')
                ak = fieldnames(entry.amplitude_scaling);
                for ai = 1:numel(ak)
                    ival = str2double(strrep(strrep(ak{ai},'i',''),'p','.'));
                    amp  = entry.amplitude_scaling.(ak{ai});
                    rows = [rows, sprintf('<tr><td>%.1f</td><td>%.1f</td><td>%d</td><td>%s</td></tr>\n', ...
                        depth_val, ival, amp, tf2str(is_ref))]; %#ok<AGROW>
                end
            end
        end
        if isempty(rows)
            rows = '<tr><td colspan="4"><em>No amplitude data.</em></td></tr>';
        end
        amp_html = ['<table class="data">', header, rows, '</table>'];
    end

    %% ── Fit quality (meta.quality) ──────────────────────────────────────
    quality_html = '<p><em>No fit quality recorded.</em></p>';
    if isfield(meta, 'quality') && isstruct(meta.quality)
        qn = fieldnames(meta.quality);
        rows = '';
        for k = 1:numel(qn)
            v = meta.quality.(qn{k});
            if isnumeric(v) || islogical(v)
                v = strjoin(arrayfun(@(x) sprintf('%.3g', x), double(v(:)'), 'UniformOutput', false), ', ');
            elseif iscell(v) || isstring(v)
                v = strjoin(cellstr(v), ', ');
            end
            rows = [rows, sprintf('<tr><td>%s</td><td>%s</td></tr>\n', qn{k}, v)]; %#ok<AGROW>
        end
        quality_html = ['<table class="data"><tr><th>Measure</th><th>Value</th></tr>', rows, '</table>'];
    end

    %% ── Settings at the fitted depths (as resolved at run time) ─────────
    settings_html = '<p><em>No fitted depths in library.</em></p>';
    depths = [];
    if isfield(lib, 'refinement') && isfield(lib.refinement, 'depths_ep_mm')
        depths = lib.refinement.depths_ep_mm;
    elseif isfield(lib, 'global_model') && isfield(lib.global_model, 'fit_depths_ep_mm')
        depths = lib.global_model.fit_depths_ep_mm;
    elseif isfield(lib, 'global_model') && isfield(lib.global_model, 'depths_ep_mm')
        depths = lib.global_model.depths_ep_mm;
    end
    tran_key = matlab.lang.makeValidName(tran_serial);
    if ~isempty(depths) && isfield(equip_param.trans, tran_key)
        try
            T = library_settings_table(equip_param.trans.(tran_key), lib, double(depths(:)'));
            n_elem = equip_param.trans.(tran_key).transducer.annular.elem_n;
            ph = table2array(T(:, compose('elem_phase_deg_%d', 1:n_elem)));
            gw = table2array(T(:, compose('elem_amp_weight_%d', 1:n_elem)));
            header = ['<tr><th>Depth (mm EP)</th><th>Phase model</th><th>Drive source</th>', ...
                '<th>Drive (kPa)</th><th>Model focus (mm EP)</th><th>Phases (deg)</th><th>Gains</th></tr>'];
            rows = '';
            for k = 1:height(T)
                rows = [rows, sprintf(['<tr><td>%.1f</td><td>%s</td><td>%s</td><td>%.1f</td>', ...
                    '<td>%.1f</td><td>%s</td><td>%s</td></tr>\n'], T.depth_ep_mm(k), T.phase_model(k), ...
                    T.amplitude_source(k), T.drive_pa(k) / 1e3, T.model_focus_ep_mm(k), ...
                    num2str(round(mod(ph(k, :), 360), 1)), num2str(round(gw(k, :), 3)))]; %#ok<AGROW>
            end
            settings_html = ['<table class="data">', header, rows, '</table>'];
        catch ME
            settings_html = sprintf('<p><em>Settings could not be resolved: %s</em></p>', ME.message);
        end
    elseif ~isempty(depths)
        settings_html = sprintf('<p><em>Transducer %s not in the equipment config.</em></p>', tran_serial);
    end

    %% ── Render HTML ─────────────────────────────────────────────────────
    now_str = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
    html = sprintf(['<!DOCTYPE html>\n<html lang="en">\n<head>\n' ...
        '<meta charset="UTF-8">\n' ...
        '<meta name="viewport" content="width=device-width, initial-scale=1">\n' ...
        '<title>Calibration report — %s</title>\n' ...
        '<style>\n' ...
        'body{font-family:system-ui,sans-serif;max-width:900px;margin:2rem auto;color:#222;line-height:1.5}\n' ...
        'h1{color:#1a3a5c;border-bottom:2px solid #1a3a5c;padding-bottom:.4rem}\n' ...
        'h2{color:#2c5f8a;margin-top:2rem}\n' ...
        'table.data{border-collapse:collapse;margin-top:.5rem;min-width:320px}\n' ...
        'table.data th,table.data td{border:1px solid #ccc;padding:.35rem .7rem;text-align:left}\n' ...
        'table.data th{background:#e8f0f8}\n' ...
        '.badge{display:inline-block;padding:.15rem .5rem;border-radius:.3rem;font-size:.85rem;' ...
        'background:#dff0d8;color:#2d6a2d;border:1px solid #b2d8b2}\n' ...
        '.meta-grid{display:grid;grid-template-columns:auto 1fr;gap:.25rem 1rem;max-width:600px}\n' ...
        '.meta-grid dt{font-weight:600;color:#555}\n' ...
        '.meta-grid dd{margin:0}\n' ...
        'footer{margin-top:3rem;font-size:.8rem;color:#888;border-top:1px solid #ddd;padding-top:.5rem}\n' ...
        '</style>\n</head>\n<body>\n'], combo_name);

    html = [html, sprintf(['<h1>Calibration Report</h1>\n' ...
        '<p><span class="badge">%s</span></p>\n' ...
        '<dl class="meta-grid">\n' ...
        '<dt>Combo</dt><dd>%s</dd>\n' ...
        '<dt>Transducer serial</dt><dd>%s</dd>\n' ...
        '<dt>Created</dt><dd>%s</dd>\n' ...
        '<dt>PRESTUS version</dt><dd>%s</dd>\n' ...
        '</dl>\n'], combo_name, combo_name, tran_serial, created_at, prestus_version)];

    html = [html, sprintf(['<h2>Calibration Parameters</h2>\n' ...
        '<dl class="meta-grid">\n' ...
        '<dt>Calibration type</dt><dd>%s</dd>\n' ...
        '<dt>Library mode</dt><dd>%s</dd>\n' ...
        '<dt>Forward model</dt><dd>%s</dd>\n' ...
        '<dt>Lambda (regularisation)</dt><dd>%s</dd>\n' ...
        '<dt>Refined depths</dt><dd>%s</dd>\n' ...
        '<dt>Refinement blend</dt><dd>%s</dd>\n' ...
        '</dl>\n'], cal_type, lib_mode, fwd_model, lambda_str, n_depths, blend_str)];

    html = [html, '<h2>Fit Quality</h2>', quality_html];
    html = [html, '<h2>Settings at the Fitted Depths</h2>', ...
        '<p>Resolved as a study gets them (<code>library_model_at_depth</code>); phases are absolute.</p>', settings_html];
    html = [html, '<h2>Static Per-Element Phase Correction</h2>', correction_html];
    html = [html, '<h2>Amplitude Lookup</h2>', amp_html];

    html = [html, sprintf(['<footer>Generated by GENERATE_CALIBRATION_REPORT &bull; ' ...
        '%s &bull; Library: %s</footer>\n</body>\n</html>\n'], now_str, yaml_path)];

    %% ── Write ────────────────────────────────────────────────────────────
    fid = fopen(output_path, 'w');
    if fid < 0
        error('generate_calibration_report: cannot write to %s', output_path);
    end
    fwrite(fid, html, 'char');
    fclose(fid);

    html_path = output_path;
    fprintf('Calibration report written: %s\n', html_path);
end

%% ── Helpers ──────────────────────────────────────────────────────────────
function v = get_field(s, fname, default)
    if isfield(s, fname) && ~isempty(s.(fname))
        val = s.(fname);
        if isnumeric(val); v = num2str(val);
        elseif islogical(val); v = tf2str(val);
        else; v = val;
        end
    else
        v = default;
    end
end

function s = tf2str(v)
    if v; s = 'yes'; else; s = 'no'; end
end
