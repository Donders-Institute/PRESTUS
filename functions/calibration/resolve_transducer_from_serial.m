function tr = resolve_transducer_from_serial(tr, equip_param, library_path, t_i)
% RESOLVE_TRANSDUCER_FROM_SERIAL  Populate transducer struct from equipment and calibration library
%
% Given a transducer entry with a 'serial' field, loads geometry from the
% equipment config and—if elem_phase_deg or elem_amp are absent—resolves
% them from the calibration library.  Fields already set in the study config
% take precedence over library values.
%
% Errors with actionable instructions when:
%   - The serial is not found in the equipment config
%   - Calibration is needed (phases/amplitude absent) but no library file exists
%
% Use as:
%   tr = resolve_transducer_from_serial(tr, equip_param, library_path, t_i)
%
% Input:
%   tr           - single transducer struct (parameters.transducer(i))
%   equip_param  - struct from load_equipment_config()
%   library_path - path to calibration library folder (config/transducer/)
%   t_i          - transducer index (for error messages)
%
% Output:
%   tr - transducer struct populated with geometry and, if available,
%        calibrated phases and amplitude
%
% See also: LOAD_EQUIPMENT_CONFIG, LOAD_TRANSDUCER_FROM_LIBRARY,
%           LOAD_TRANSDUCER_PARAMETERS

    serial = char(tr.serial);

    %% ── Geometry from equipment config ───────────────────────────────────
    if ~isfield(equip_param.trans, serial)
        error(['resolve_transducer_from_serial: transducer serial ''%s'' ' ...
            '(transducer %d) not found in equipment library.\n' ...
            'Available serials: %s\n' ...
            'Add a YAML file to config/equipment/ — see doc/doc_transducer_library.md.'], ...
            serial, t_i, strjoin(fieldnames(equip_param.trans), ', '));
    end

    tran_entry = equip_param.trans.(serial);
    geom       = tran_entry.transducer;   % contains .annular (or .matrix) + .freq_hz

    % Merge geometry: equipment YAML provides defaults; study config overrides.
    tr = mergestruct_shallow(geom, tr);

    % Ensure type is set if only serial was provided.
    if ~isfield(tr, 'type') || isempty(tr.type)
        if isfield(geom, 'annular') && ~isempty(geom.annular)
            tr.type = 'annular';
        elseif isfield(geom, 'matrix') && ~isempty(geom.matrix)
            tr.type = 'matrix';
        end
    end

    %% ── Focal depth range check ─────────────────────────────────────────
    if isfield(tr, 'focal_distance_ep') && ~isempty(tr.focal_distance_ep)
        if isfield(tran_entry, 'min_foc') && ~isempty(tran_entry.min_foc) && ...
                tr.focal_distance_ep < tran_entry.min_foc
            warning('resolve_transducer_from_serial:outOfRange', ...
                'Transducer %d (%s): requested focal distance %.1f mm is below min_foc %.1f mm.', ...
                t_i, serial, tr.focal_distance_ep, tran_entry.min_foc);
        end
        if isfield(tran_entry, 'max_foc') && ~isempty(tran_entry.max_foc) && ...
                tr.focal_distance_ep > tran_entry.max_foc
            warning('resolve_transducer_from_serial:outOfRange', ...
                'Transducer %d (%s): requested focal distance %.1f mm exceeds max_foc %.1f mm.', ...
                t_i, serial, tr.focal_distance_ep, tran_entry.max_foc);
        end
    end

    %% ── Calibration library lookup ───────────────────────────────────────
    % Only annular transducers currently support library-based calibration.
    if ~strcmp(tr.type, 'annular')
        return;
    end

    phases_present = isfield(tr, 'annular') && isfield(tr.annular, 'elem_phase_deg') ...
        && ~isempty(tr.annular.elem_phase_deg) && ~all(isnan(tr.annular.elem_phase_deg(:)));
    amp_present    = isfield(tr, 'annular') && isfield(tr.annular, 'elem_amp') ...
        && ~isempty(tr.annular.elem_amp) && ~all(isnan(tr.annular.elem_amp(:)));

    if phases_present && amp_present
        return;  % fully specified inline — nothing to resolve
    end

    % Need focal distance and target intensity for library lookup.
    has_focal = isfield(tr, 'focal_distance_ep') && ~isempty(tr.focal_distance_ep) ...
        && ~isnan(tr.focal_distance_ep);
    has_isppa = isfield(tr, 'target_isppa_wcm2') && ~isempty(tr.target_isppa_wcm2) ...
        && ~isnan(tr.target_isppa_wcm2);

    if ~has_focal || ~has_isppa
        % Cannot resolve without these — let downstream validation catch it.
        return;
    end

    % Determine combo name (transducer + driving system).
    combo_name = resolve_combo_name(serial, tr, equip_param, t_i);

    % Check that a library YAML exists for this combo.
    yaml_path = fullfile(library_path, [combo_name, '.yaml']);
    if ~isfile(yaml_path)
        error(['resolve_transducer_from_serial: no calibration library found for ''%s'' ' ...
            '(transducer %d).\n' ...
            'Expected: %s\n\n' ...
            'Options:\n' ...
            '  1. Run calibration and deposit to config/transducer/ — see doc/doc_transducer_library.md\n' ...
            '  2. Provide elem_phase_deg and elem_amp inline in the study config to bypass the library.\n' ...
            '  3. Set combo.ds_serial to use a driving-system-specific calibration ' ...
            '(key: %s_{ds_serial}.yaml).'], ...
            combo_name, t_i, yaml_path, serial);
    end

    % Load from library and merge phases/amplitude (without overwriting inline values).
    lib_params = load_transducer_from_library(combo_name, tr.focal_distance_ep, ...
        tr.target_isppa_wcm2, equip_param, library_path);

    if ~phases_present
        tr.annular.elem_phase_deg = lib_params.transducer.annular.elem_phase_deg;
        tr.annular.elem_phase_rad = lib_params.transducer.annular.elem_phase_rad;
    end
    if ~amp_present
        tr.annular.elem_amp = lib_params.transducer.annular.elem_amp;
    end

    fprintf(['[transducer %d] Resolved from library: %s | focal %.1f mm ' ...
        '(calibrated %.1f mm) | %.1f W/cm²\n'], ...
        t_i, combo_name, tr.focal_distance_ep, lib_params.library_focal_used, ...
        tr.target_isppa_wcm2);
end

% ── Helpers ──────────────────────────────────────────────────────────────

function combo_name = resolve_combo_name(serial, tr, ~, ~)
% Return the library lookup key.
% With ds_serial: '{serial}_{ds_serial}'  (combo-specific calibration)
% Without:        '{serial}'              (generic / DS-agnostic calibration)
    if isfield(tr, 'combo') && isstruct(tr.combo) && ...
            isfield(tr.combo, 'ds_serial') && ~isempty(tr.combo.ds_serial)
        combo_name = [serial, '_', char(tr.combo.ds_serial)];
    else
        combo_name = serial;
    end
end

function dst = mergestruct_shallow(base, override)
% Merge two structs: fields from override take precedence over base.
% Only one level deep — nested structs are handled field-by-field.
    dst = base;
    fn  = fieldnames(override);
    for i = 1:numel(fn)
        f = fn{i};
        if isstruct(override.(f)) && isfield(base, f) && isstruct(base.(f))
            dst.(f) = mergestruct_shallow(base.(f), override.(f));
        else
            dst.(f) = override.(f);
        end
    end
end
