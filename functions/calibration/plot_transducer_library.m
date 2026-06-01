function plot_transducer_library(combo_name, library_path, equip_param)
% PLOT_TRANSDUCER_LIBRARY  Visualise calibration library entry for a combo
%
% Creates a figure with two subplots:
%   1. Phase vs depth for each transducer element (lines + optional markers)
%   2. Amplitude vs ISPPA from the calibrated amplitude_scaling entries
%
% Use as:
%   plot_transducer_library(combo_name)
%   plot_transducer_library(combo_name, library_path)
%   plot_transducer_library(combo_name, library_path, equip_param)
%
% Input:
%   combo_name   - equipment combination key (library filename without .yaml)
%   library_path - path to transducer library folder
%                  (default: config/transducer/ under PRESTUS root)
%   equip_param  - struct from load_equipment_config (optional; loaded if missing)
%
% Output:
%   (none) — produces a MATLAB figure
%
% See also: UPDATE_TRANSDUCER_LIBRARY, VALIDATE_TRANSDUCER_LIBRARY

    if nargin < 2 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'transducer');
    end
    if nargin < 3 || isempty(equip_param)
        equip_param = load_equipment_config();
    end

    %% Load library
    lib_file = fullfile(library_path, [combo_name '.yaml']);
    if ~isfile(lib_file)
        error('plot_transducer_library: library file not found:\n  %s', lib_file);
    end
    lib = yaml.loadFile(lib_file, 'ConvertToArray', true);

    %% Determine tran_serial and load geometry
    if isfield(lib, 'meta') && isfield(lib.meta, 'tran_serial')
        tran_serial = lib.meta.tran_serial;
    else
        tran_serial = combo_name;
    end

    n_elem = NaN;
    if isfield(equip_param.trans, tran_serial)
        n_elem = equip_param.trans.(tran_serial).transducer.annular.elem_n;
    end

    %% ── Figure ──────────────────────────────────────────────────────────────
    fig = figure('Name', ['Library: ' combo_name], ...
                 'NumberTitle', 'off', ...
                 'Position', [100 100 900 420]);

    %% Subplot 1: Phase vs depth
    ax1 = subplot(1, 2, 1, 'Parent', fig);
    hold(ax1, 'on');

    % Collect per-depth data from global_model (lines) and per-depth entries (markers)
    gm_depths = [];
    gm_phases = [];  % [n_depths × n_elem]

    if isfield(lib, 'global_model') && isfield(lib.global_model, 'depths_ep_mm')
        gm_depths  = lib.global_model.depths_ep_mm(:)';
        phase_cell = lib.global_model.elem_phase_deg;
        if ~iscell(phase_cell)
            phase_cell = num2cell(phase_cell, 2);
        end
        gm_phases = cell2mat(cellfun(@(r) r(:)', phase_cell, 'UniformOutput', false));

        % If n_elem still unknown, infer from data
        if isnan(n_elem)
            n_elem = size(gm_phases, 2);
        end

        cmap = lines(n_elem);
        for el = 1:size(gm_phases, 2)
            plot(ax1, gm_depths, gm_phases(:, el), '-o', ...
                'Color', cmap(el, :), 'LineWidth', 1.2, ...
                'DisplayName', sprintf('Elem %d', el));
        end
    end

    % Overlay per-depth parametric entries as square markers
    if isfield(lib, 'calibration') && isfield(lib.calibration, 'focal_depths')
        dk     = fieldnames(lib.calibration.focal_depths);
        depths_pd = cellfun(@(k) str2double(strrep(strrep(k,'f',''),'p','.')), dk);
        if isnan(n_elem), n_elem = 4; end
        cmap = lines(n_elem);
        for di = 1:numel(dk)
            entry = lib.calibration.focal_depths.(dk{di});
            if isfield(entry, 'elem_phase_deg')
                ph_deg = entry.elem_phase_deg(:)';
            elseif isfield(entry, 'phase_start_deg') && isfield(entry, 'phase_step_deg')
                ph_deg = entry.phase_start_deg + (0:n_elem-1) * entry.phase_step_deg;
            else
                continue;
            end
            for el = 1:min(numel(ph_deg), n_elem)
                plot(ax1, depths_pd(di), ph_deg(el), 's', ...
                    'Color', cmap(el, :), 'MarkerSize', 8, 'MarkerFaceColor', cmap(el, :), ...
                    'HandleVisibility', 'off');
            end
        end
    end

    xlabel(ax1, 'Depth (mm)');
    ylabel(ax1, 'Phase (°)');
    title(ax1, 'Phase vs Depth');
    if ~isnan(n_elem) && n_elem <= 12
        legend(ax1, 'show', 'Location', 'best', 'FontSize', 7);
    end
    grid(ax1, 'on');

    %% Subplot 2: Amplitude vs ISPPA
    ax2 = subplot(1, 2, 2, 'Parent', fig);
    hold(ax2, 'on');

    amp_map = [];
    if isfield(lib, 'global_model') && isfield(lib.global_model, 'amplitude_scaling')
        amp_map = lib.global_model.amplitude_scaling;
    elseif isfield(lib, 'calibration') && isfield(lib.calibration, 'focal_depths')
        dk = fieldnames(lib.calibration.focal_depths);
        if ~isempty(dk)
            first_entry = lib.calibration.focal_depths.(dk{1});
            if isfield(first_entry, 'amplitude_scaling')
                amp_map = first_entry.amplitude_scaling;
            end
        end
    end

    if ~isempty(amp_map) && isstruct(amp_map)
        int_keys   = fieldnames(amp_map);
        int_values = cellfun(@(k) str2double(strrep(strrep(k, 'i', ''), 'p', '.')), int_keys);
        amp_values = cellfun(@(k) double(amp_map.(k)), int_keys);
        [int_values, srt] = sort(int_values);
        amp_values = amp_values(srt);

        plot(ax2, int_values, amp_values, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
        xlabel(ax2, 'ISPPA (W/cm²)');
        ylabel(ax2, 'Amplitude (a.u.)');
        title(ax2, 'Amplitude vs Intensity');
        grid(ax2, 'on');
    else
        text(0.5, 0.5, 'No amplitude scaling data', ...
            'Parent', ax2, 'HorizontalAlignment', 'center', 'Units', 'normalized');
    end

    sgtitle(fig, combo_name, 'Interpreter', 'none');
end
