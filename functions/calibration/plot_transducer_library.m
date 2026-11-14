function fig = plot_transducer_library(combo_name, library_path, equip_param)
% PLOT_TRANSDUCER_LIBRARY  Visualise the settings a calibration library entry resolves to
%
% Resolves the entry over a range of focal depths exactly as a study does
% (LIBRARY_SETTINGS_TABLE, i.e. LIBRARY_MODEL_AT_DEPTH, the runtime path of
% LOAD_TRANSDUCER_FROM_LIBRARY): per-depth refinements interpolated between
% the refined depths and blended back to the global model outside them.
% Works for every library schema the loader reads. Three panels:
%   1. Ring phase relative to ring 1 vs depth (markers: refined / fitted depths)
%   2. Drive (mean source pressure) vs depth, at the stored reference ISPPA
%   3. Per-ring amplitude gain vs depth
%
% Use as:
%   plot_transducer_library(combo_name)
%   plot_transducer_library(combo_name, library_path)
%   plot_transducer_library(combo_name, library_path, equip_param)
%   fig = plot_transducer_library(...)
%
% Input:
%   combo_name   - equipment combination key (library filename without .yaml)
%   library_path - path to transducer library folder
%                  (default: config/calibration/ under PRESTUS root)
%   equip_param  - struct from load_equipment_config (optional; loaded if missing)
%
% Output:
%   fig - figure handle
%
% See also: LIBRARY_SETTINGS_TABLE, GENERATE_CALIBRATION_REPORT, VALIDATE_TRANSDUCER_LIBRARY

    if nargin < 2 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'calibration');
    end
    if nargin < 3 || isempty(equip_param)
        equip_param = load_equipment_config();
    end

    [lib, tran, knots] = load_library_entry(combo_name, library_path, equip_param);
    depths = library_plot_depths(tran, knots);
    T  = library_settings_table(tran, lib, depths);
    Tk = [];
    if ~isempty(knots); Tk = library_settings_table(tran, lib, knots); end

    n_elem = tran.transducer.annular.elem_n;
    cmap   = lines(n_elem);

    fig = figure('Name', ['Library: ' combo_name], 'NumberTitle', 'off', ...
                 'Position', [100 100 1200 400]);

    %% Phase vs depth
    ax1 = subplot(1, 3, 1, 'Parent', fig);
    hold(ax1, 'on');
    ph = relative_phase(table2array(T(:, compose('elem_phase_deg_%d', 1:n_elem))));
    for el = 1:n_elem
        plot(ax1, T.depth_ep_mm, ph(:, el), '-', 'Color', cmap(el, :), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('Ring %d', el));
    end
    if ~isempty(Tk)
        phk = relative_phase(table2array(Tk(:, compose('elem_phase_deg_%d', 1:n_elem))));
        for el = 1:n_elem
            % markers on the unwrapped curve: shift each knot by whole turns
            ref = interp1(T.depth_ep_mm, ph(:, el), Tk.depth_ep_mm);
            plot(ax1, Tk.depth_ep_mm, phk(:, el) + 360 * round((ref - phk(:, el)) / 360), 'o', ...
                'Color', cmap(el, :), 'MarkerFaceColor', cmap(el, :), 'HandleVisibility', 'off');
        end
    end
    xlabel(ax1, 'Focal depth from exit plane (mm)');
    ylabel(ax1, 'Phase relative to ring 1 (°)');
    title(ax1, 'Ring phases');
    legend(ax1, 'show', 'Location', 'best', 'FontSize', 7);
    grid(ax1, 'on');

    %% Drive vs depth
    ax2 = subplot(1, 3, 2, 'Parent', fig);
    hold(ax2, 'on');
    plot(ax2, T.depth_ep_mm, T.drive_pa / 1e3, 'k-', 'LineWidth', 1.4);
    if ~isempty(Tk)
        plot(ax2, Tk.depth_ep_mm, Tk.drive_pa / 1e3, 'ko', 'MarkerFaceColor', 'k');
    end
    xlabel(ax2, 'Focal depth from exit plane (mm)');
    ylabel(ax2, 'Mean source pressure (kPa)');
    title(ax2, sprintf('Drive (%s)', strjoin(unique(cellstr(T.amplitude_source))', ', ')), ...
        'Interpreter', 'none');
    grid(ax2, 'on');
    xlim(ax2, [min(T.depth_ep_mm) max(T.depth_ep_mm)]);
    if all(isnan(T.drive_pa))
        text(0.5, 0.5, 'No amplitude data', 'Parent', ax2, ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
    end

    %% Ring gains vs depth
    ax3 = subplot(1, 3, 3, 'Parent', fig);
    hold(ax3, 'on');
    g = table2array(T(:, compose('elem_amp_weight_%d', 1:n_elem)));
    for el = 1:n_elem
        plot(ax3, T.depth_ep_mm, g(:, el), '-', 'Color', cmap(el, :), 'LineWidth', 1.2);
    end
    xlabel(ax3, 'Focal depth from exit plane (mm)');
    ylabel(ax3, 'Gain (mean 1)');
    title(ax3, 'Ring gains');
    grid(ax3, 'on');

    cal_type = '';
    if isfield(lib, 'meta') && isfield(lib.meta, 'calibration_type')
        cal_type = sprintf(' (%s)', lib.meta.calibration_type);
    end
    sgtitle(fig, [combo_name cal_type], 'Interpreter', 'none');
end

% -------------------------------------------------------------------------
function [lib, tran, knots] = load_library_entry(combo_name, library_path, equip_param)
% Library struct, transducer struct and the depths the entry was fitted at.
    lib_file = library_entry_file(library_path, combo_name);
    if ~isfile(lib_file)
        error('plot_transducer_library: library file not found:\n  %s', lib_file);
    end
    lib = yaml.loadFile(lib_file, 'ConvertToArray', true);
    tran_serial = combo_name;
    if isfield(lib, 'meta') && isfield(lib.meta, 'tran_serial')
        tran_serial = lib.meta.tran_serial;
    end
    if ~isfield(equip_param.trans, tran_serial)
        error('plot_transducer_library: transducer %s not in the equipment config.', tran_serial);
    end
    tran  = equip_param.trans.(tran_serial);
    knots = [];
    if isfield(lib, 'refinement') && isfield(lib.refinement, 'depths_ep_mm')
        knots = lib.refinement.depths_ep_mm;
    elseif isfield(lib, 'global_model') && isfield(lib.global_model, 'fit_depths_ep_mm')
        knots = lib.global_model.fit_depths_ep_mm;
    elseif isfield(lib, 'global_model') && isfield(lib.global_model, 'depths_ep_mm')
        knots = lib.global_model.depths_ep_mm;
    end
    knots = double(knots(:)');
end

function d = library_plot_depths(tran, knots)
% Depth grid: the fitted span plus 10 mm either side, else around the
% geometric focus.
    if isempty(knots)
        f = tran.transducer.annular.dist_geom_ep_mm;
        lo = 0.5 * f; hi = 1.3 * f;
    else
        lo = min(knots) - 10; hi = max(knots) + 10;
    end
    d = linspace(max(lo, 5), hi, 80);
end

function ph = relative_phase(ph)
% Phase of each ring relative to ring 1, wrapped to (-180, 180] and then
% unwrapped along depth so curves stay continuous.
    ph = angle(exp(1i * (ph - ph(:, 1)) * pi / 180));
    ph = unwrap(ph, [], 1) * 180 / pi;
end
