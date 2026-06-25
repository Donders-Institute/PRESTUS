function parameters = load_transducer_from_library(combo_name, focal_distance_ep, desired_intensity, equip_param, library_path)
% LOAD_TRANSDUCER_FROM_LIBRARY  Load calibrated transducer parameters from library
%
% Reads the calibration library YAML for the given equipment combination and
% the geometry from equip_param (load_equipment_config). Reconstructs per-
% element phases from the fitted parametric model and looks up (or scales)
% the amplitude for the desired intensity.
%
%
% When the library YAML contains a global_model section (phases calibrated
% across all depths simultaneously), phases are interpolated linearly at the
% requested depth. Per-depth entries with is_refinement:true override the
% global model when the requested depth is within 1 mm of that entry.
%
% Use as:
%   parameters = load_transducer_from_library(combo_name, focal_distance_ep, ...
%                    desired_intensity, equip_param)
%   parameters = load_transducer_from_library(combo_name, focal_distance_ep, ...
%                    desired_intensity, equip_param, library_path)
%
% Input:
%   combo_name         - equipment combination key, e.g.
%                        'IS_PCD15473_01001_IGT_32_ch_comb_10_ch'
%   focal_distance_ep  - desired focal distance from exit plane [mm]
%   desired_intensity  - target ISPPA [W/cm²]
%   equip_param        - struct from load_equipment_config()
%   library_path       - (optional) path to transducer library folder;
%                        defaults to config/calibration/ under PRESTUS root
%
% Output:
%   parameters - struct with:
%     .transducer              — geometry + calibrated phases and amplitude
%     .transducer.tx_calibration.isppa_wcm2        — ISPPA the elem_p_amp produces [W/cm²]
%     .transducer.tx_calibration.focal_depth_ep_mm — calibrated depth used [mm]
%     .transducer.tx_calibration.phase_source      — 'library'
%     .transducer.tx_calibration.isppa_source      — 'library'
%     .transducer.tx_calibration.library_key       — combo name used for lookup
%     .transducer.tx_calibration.library_meta      — provenance from the library YAML
%
% See also: UPDATE_TRANSDUCER_LIBRARY, CALIBRATION_TRANSDUCER

    if nargin < 5 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'calibration');
    end

    %% Load calibration library entry
    yaml_path = fullfile(library_path, [combo_name, '.yaml']);
    if ~isfile(yaml_path)
        error('load_transducer_from_library: library file not found:\n  %s', yaml_path);
    end
    lib = yaml.loadFile(yaml_path, 'ConvertToArray', true);
    if ~isstruct(lib) || ~isfield(lib, 'meta')
        error('load_transducer_from_library:missingMeta', ...
            'Library file is missing a ''meta'' section:\n  %s', yaml_path);
    end
    if ~isfield(lib.meta, 'tran_serial')
        error('load_transducer_from_library:missingSerial', ...
            'Library file meta section is missing ''tran_serial'':\n  %s', yaml_path);
    end

    %% Load geometry from equipment config
    tran_serial = lib.meta.tran_serial;
    if ~isfield(equip_param.trans, tran_serial)
        error('load_transducer_from_library: transducer ''%s'' not found in equip_param.', tran_serial);
    end
    tran = equip_param.trans.(tran_serial);
    parameters.transducer = tran.transducer;

    %% Reconstruct phases — global model takes precedence when present
    if isfield(lib, 'global_model') && ~isempty(lib.global_model) && ...
            isfield(lib.global_model, 'depths_ep_mm') && ...
            isfield(lib.global_model, 'elem_phase_deg')

        gm     = lib.global_model;
        depths = gm.depths_ep_mm(:)';

        % elem_phase_deg is stored as a cell array of per-depth row vectors
        phase_cell = gm.elem_phase_deg;
        if ~iscell(phase_cell)
            phase_cell = num2cell(phase_cell, 2);
        end
        phase_mat = cell2mat(cellfun(@(r) r(:)', phase_cell, 'UniformOutput', false));

        % Check if a per-depth refinement overrides the global model
        [dist, idx]    = min(abs(depths - focal_distance_ep));
        used_depth_gm  = depths(idx);
        use_refinement = false;

        if isfield(lib, 'calibration') && isfield(lib.calibration, 'focal_depths')
            fkey = sprintf('f%s', strrep(num2str(used_depth_gm), '.', 'p'));
            if dist < 1 && isfield(lib.calibration.focal_depths, fkey) && ...
                    isfield(lib.calibration.focal_depths.(fkey), 'is_refinement') && ...
                    lib.calibration.focal_depths.(fkey).is_refinement
                use_refinement = true;
            end
        end

        if ~use_refinement
            % Interpolate phases from global model
            if numel(depths) == 1
                interp_phases_deg = phase_mat(1, :);
            else
                interp_phases_deg = interp1(depths, phase_mat, focal_distance_ep, 'linear', 'extrap');
            end
            elem_phase_rad = interp_phases_deg / 180 * pi;

            if dist > 5
                warning(['load_transducer_from_library: closest global-model depth is %.1f mm ' ...
                    '(requested %.1f mm). Consider running calibration at this depth.'], ...
                    used_depth_gm, focal_distance_ep);
            end

            % Apply per-element phase correction.
            % Library YAML takes precedence (written by multi_depth calibration).
            % Falls back to equipment YAML for backward compatibility with older
            % calibrations written by save_elem_correction.
            n_elem = tran.transducer.annular.elem_n;
            if isfield(gm, 'elem_phase_correction_deg') && ~isempty(gm.elem_phase_correction_deg)
                correction_rad = gm.elem_phase_correction_deg(:)' * pi/180;
                if numel(correction_rad) ~= n_elem
                    warning(['load_transducer_from_library: library elem_phase_correction_deg has ' ...
                        '%d elements but transducer has %d — correction skipped.'], ...
                        numel(correction_rad), n_elem);
                else
                    elem_phase_rad = mod(elem_phase_rad + correction_rad, 2*pi);
                    fprintf('Applied phase correction from library (multi_depth).\n');
                end
            elseif isfield(tran, 'elem_phase_correction') && ~isempty(tran.elem_phase_correction)
                correction_rad = tran.elem_phase_correction.deg(:)' * pi/180;
                if numel(correction_rad) ~= n_elem
                    warning(['load_transducer_from_library: equipment elem_phase_correction has ' ...
                        '%d elements but transducer has %d — correction skipped.'], ...
                        numel(correction_rad), n_elem);
                else
                    elem_phase_rad = mod(elem_phase_rad + correction_rad, 2*pi);
                    fprintf('Applied phase correction from equipment YAML (legacy).\n');
                end
            end

            parameters.transducer.annular.elem_phase_rad = elem_phase_rad;
            parameters.transducer.annular.elem_phase_deg = elem_phase_rad / pi * 180;

            % Amplitude from global_model.amplitude_scaling
            if isfield(gm, 'amplitude_scaling')
                amp_map    = gm.amplitude_scaling;
                int_keys   = fieldnames(amp_map);
                int_values = cellfun(@(k) str2double(strrep(strrep(k, 'i', ''), 'p', '.')), int_keys);

                [int_dist, amp_idx] = min(abs(int_values - desired_intensity));
                int_used            = int_values(amp_idx);

                if int_dist / desired_intensity > 0.05
                    base_amp = amp_map.(int_keys{amp_idx});
                    elem_p_amp = round(base_amp * sqrt(desired_intensity / int_used));
                    warning(['load_transducer_from_library: scaling amplitude from %.1f to %.1f W/cm² ' ...
                        '(%.0f -> %d). Run calibration at this intensity for an exact value.'], ...
                        int_used, desired_intensity, base_amp, elem_p_amp);
                else
                    elem_p_amp = amp_map.(int_keys{amp_idx});
                end
                parameters.transducer.annular.elem_p_amp = repmat(elem_p_amp, 1, n_elem);
            end

            %% Focal distance bookkeeping
            parameters.transducer.focal_distance_ep   = focal_distance_ep;
            parameters.transducer.focal_distance_bowl = focal_distance_ep + ...
                (tran.transducer.annular.curv_radius_mm - ...
                 tran.transducer.annular.dist_geom_ep_mm);
            parameters.transducer.tx_calibration.isppa_wcm2        = desired_intensity;
            parameters.transducer.tx_calibration.focal_depth_ep_mm = used_depth_gm;
            parameters.transducer.tx_calibration.phase_source       = 'library';
            parameters.transducer.tx_calibration.isppa_source       = 'library';
            parameters.transducer.tx_calibration.library_key        = combo_name;
            parameters.transducer.tx_calibration.library_meta       = lib.meta;

            cal_type = '';
            if isfield(lib.meta, 'calibration_type'), cal_type = lib.meta.calibration_type; end
            fprintf('Library loaded (global model, %s): %s | focal %.1f mm (interp from %.1f mm) | amp %d @ %.1f W/cm²\n', ...
                cal_type, combo_name, focal_distance_ep, used_depth_gm, elem_p_amp(1), desired_intensity);
            return;
        end
        % fall through to per-depth if use_refinement == true
    end

    %% Find closest calibrated focal depth (refinement override or legacy per-depth)
    cal = lib.calibration.focal_depths;
    if isempty(fieldnames(cal))
        error(['load_transducer_from_library: no calibrated focal depths in library ' ...
            'for combo ''%s''. Run calibration_standalone first.'], combo_name);
    end

    depth_keys   = fieldnames(cal);
    depth_values = cellfun(@(k) str2double(strrep(strrep(k, 'f', ''), 'p', '.')), depth_keys);

    [dist, idx]    = min(abs(depth_values - focal_distance_ep));
    used_key       = depth_keys{idx};
    focal_used     = depth_values(idx);

    if dist > 5
        warning(['load_transducer_from_library: closest calibrated depth is %.1f mm ' ...
            '(requested %.1f mm). Consider running calibration at this depth.'], ...
            focal_used, focal_distance_ep);
    end

    entry  = cal.(used_key);
    n_elem = tran.transducer.annular.elem_n;

    %% Reconstruct per-element phases from parametric model
    phase_start_rad = entry.phase_start_deg / 180 * pi;

    switch entry.precession_mode
        case 'linear'
            phase_step_rad = entry.phase_step_deg / 180 * pi;
            elem_phase_rad = mod(phase_start_rad + (0:n_elem-1) * phase_step_rad, 2*pi);
        case 'monotonic'
            increments_rad = entry.phase_increments_deg / 180 * pi;
            elem_phase_rad = mod(phase_start_rad + [0, cumsum(increments_rad)], 2*pi);
        otherwise
            error('load_transducer_from_library: unknown precession_mode ''%s''.', ...
                entry.precession_mode);
    end

    % Per-depth refinement entries store phases from a single-depth global search
    % which already incorporates any correction — do not apply again.

    parameters.transducer.annular.elem_phase_rad = elem_phase_rad;
    parameters.transducer.annular.elem_phase_deg = elem_phase_rad / pi * 180;

    %% Look up amplitude for desired intensity
    amp_map     = entry.amplitude_scaling;
    int_keys    = fieldnames(amp_map);
    int_values  = cellfun(@(k) str2double(strrep(strrep(k, 'i', ''), 'p', '.')), int_keys);

    [int_dist, amp_idx] = min(abs(int_values - desired_intensity));
    int_used            = int_values(amp_idx);

    if int_dist / desired_intensity > 0.05
        % Scale analytically: amp ∝ sqrt(I) since I ∝ v² ∝ amp²
        base_amp = amp_map.(int_keys{amp_idx});
        elem_p_amp = round(base_amp * sqrt(desired_intensity / int_used));
        warning(['load_transducer_from_library: scaling amplitude from %.1f to %.1f W/cm² ' ...
            '(%.0f -> %d). Run calibration at this intensity for an exact value.'], ...
            int_used, desired_intensity, base_amp, elem_p_amp);
    else
        elem_p_amp = amp_map.(int_keys{amp_idx});
    end

    parameters.transducer.annular.elem_p_amp = repmat(elem_p_amp, 1, n_elem);

    %% Focal distance bookkeeping
    parameters.transducer.focal_distance_ep   = focal_distance_ep;
    parameters.transducer.focal_distance_bowl = focal_distance_ep + ...
        (tran.transducer.annular.curv_radius_mm - ...
         tran.transducer.annular.dist_geom_ep_mm);
    parameters.transducer.tx_calibration.isppa_wcm2        = desired_intensity;
    parameters.transducer.tx_calibration.focal_depth_ep_mm = focal_used;
    parameters.transducer.tx_calibration.phase_source       = 'library';
    parameters.transducer.tx_calibration.isppa_source       = 'library';
    parameters.transducer.tx_calibration.library_key        = combo_name;
    parameters.transducer.tx_calibration.library_meta       = lib.meta;

    tier = 'per-depth';
    if use_refinement, tier = 'refinement override'; end
    cal_type = '';
    if isfield(lib.meta, 'calibration_type'), cal_type = lib.meta.calibration_type; end
    fprintf('Library loaded (%s, %s): %s | focal %.1f mm (calibrated %.1f mm) | amp %d @ %.1f W/cm²\n', ...
        tier, cal_type, combo_name, focal_distance_ep, focal_used, elem_p_amp(1), desired_intensity);
end
