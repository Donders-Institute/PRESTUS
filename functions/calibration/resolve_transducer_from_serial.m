function tr = resolve_transducer_from_serial(tr, equip_param, library_path, t_i)
% RESOLVE_TRANSDUCER_FROM_SERIAL  Populate transducer struct from equipment and calibration library
%
% Given a transducer entry with a 'serial' field, loads geometry from the
% equipment config and—if elem_phase_deg or elem_p_amp are absent—resolves
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
%   library_path - path to calibration library folder (config/calibration/)
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
    amp_present    = isfield(tr, 'annular') && isfield(tr.annular, 'elem_p_amp') ...
        && ~isempty(tr.annular.elem_p_amp) && ~all(isnan(tr.annular.elem_p_amp(:)));

    if phases_present && amp_present
        return;  % fully specified inline — nothing to resolve
    end

    % Need focal distance and target intensity for library lookup.
    has_focal = isfield(tr, 'focal_distance_ep') && ~isempty(tr.focal_distance_ep) ...
        && ~isnan(tr.focal_distance_ep);
    has_isppa = isfield(tr, 'target_isppa_wcm2') && ~isempty(tr.target_isppa_wcm2) ...
        && all(~isnan(tr.target_isppa_wcm2(:)));
    if has_isppa
        % A vector requests a multi-ISPPA sweep: the acoustic run uses the
        % first value as drive level and multi_isppa_pipeline rescales the
        % pressure to every target (sqrt(target / free-water baseline)).
        drive_isppa = tr.target_isppa_wcm2(1);
    end

    % EXPERIMENTAL: a dephased library condition replaces the focal depth
    % (focal_distance_ep then defaults to the condition's nominal setpoint).
    condition = '';
    if isfield(tr, 'condition') && ~isempty(tr.condition) && ~isa(tr.condition, 'yaml.Null')
        condition = char(tr.condition);
    end
    if (~has_focal && isempty(condition)) || ~has_isppa
        % Cannot resolve without these — let downstream validation catch it.
        return;
    end
    focal_req = [];
    if has_focal; focal_req = tr.focal_distance_ep; end

    % Determine combo name (transducer + driving system).
    combo_name = resolve_combo_name(serial, tr, equip_param, t_i);

    % A user library (transducer.library_path, e.g. written by
    % refine_transducer_library) takes precedence over the shipped one.
    if isfield(tr, 'library_path') && ~isempty(tr.library_path)
        user_path = char(tr.library_path);
        if isfile(fullfile(user_path, [combo_name, '.yaml'])) || isfile(fullfile(user_path, [serial, '.yaml']))
            library_path = user_path;
        else
            warning('resolve_transducer_from_serial:userLibraryMiss', ...
                'Transducer %d: no entry for ''%s'' in transducer.library_path (%s); using the shipped library.', ...
                t_i, combo_name, user_path);
        end
    end

    % Resolve library key: prefer DS-specific, fall back to DS-agnostic ({serial}.yaml).
    yaml_path_combo   = fullfile(library_path, [combo_name, '.yaml']);
    yaml_path_generic = fullfile(library_path, [serial, '.yaml']);

    if isfile(yaml_path_combo)
        lib_key = combo_name;
    elseif isfile(yaml_path_generic)
        lib_key = serial;
        warning(['resolve_transducer_from_serial:dsAgnosticFallback', ...
            'No DS-specific library found for ''%s'' (transducer %d); ' ...
            'falling back to DS-agnostic entry ''%s.yaml''.\n' ...
            'For best accuracy, run build_transducer_library with your DS combo.'], ...
            combo_name, t_i, serial);
    else
        error(['resolve_transducer_from_serial: no calibration library found for ''%s'' ' ...
            '(transducer %d).\n' ...
            'Checked:\n  %s\n  %s\n\n' ...
            'Options:\n' ...
            '  1. Run examples/build_transducer_library.m to build or refresh the library.\n' ...
            '  2. Provide elem_phase_deg and elem_p_amp inline in the study config to bypass the library.'], ...
            combo_name, t_i, yaml_path_combo, yaml_path_generic);
    end

    % Load from library and merge phases/amplitude (without overwriting inline values).
    lib_params = load_transducer_from_library(lib_key, focal_req, ...
        drive_isppa, equip_param, library_path, condition);
    if ~has_focal
        tr.focal_distance_ep = lib_params.transducer.focal_distance_ep;
    end

    if ~phases_present
        tr.annular.elem_phase_deg = lib_params.transducer.annular.elem_phase_deg;
        tr.annular.elem_phase_rad = lib_params.transducer.annular.elem_phase_rad;
    end
    if ~amp_present
        tr.annular.elem_p_amp = lib_params.transducer.annular.elem_p_amp;
    end

    % Ring profile: the library drive refers to sub-rings, so the geometry
    % comes from the entry too (inline drives refer to the equipment rings)
    rp_lib = [];
    if isfield(lib_params.transducer.tx_calibration, 'ring_profile')
        rp_lib = lib_params.transducer.tx_calibration.ring_profile;
    end
    if ~isempty(rp_lib)
        if ~phases_present && ~amp_present
            for f = {'elem_n', 'elem_id_mm', 'elem_od_mm', 'steer_r_mid_mm'}
                tr.annular.(f{1}) = lib_params.transducer.annular.(f{1});
            end
            tr.tx_calibration.ring_profile = rp_lib;
        elseif phases_present && amp_present
            warning('resolve_transducer_from_serial:ringProfileInline', ...
                ['Transducer %d: the library entry ''%s'' has a ring profile, but elem_phase_deg and ' ...
                 'elem_p_amp are given inline for the equipment rings; the ring profile is not used.'], t_i, lib_key);
        else
            error('resolve_transducer_from_serial:ringProfileInline', ...
                ['Transducer %d: the library entry ''%s'' has a ring profile (drives per sub-ring). Give both ' ...
                 'elem_phase_deg and elem_p_amp inline (equipment rings), or neither.'], t_i, lib_key);
        end
    end

    fprintf(['[transducer %d] Resolved from library: %s | focal %.1f mm ' ...
        '(%s) | %.1f W/cm²\n'], ...
        t_i, lib_key, tr.focal_distance_ep, lib_params.transducer.tx_calibration.phase_model, ...
        drive_isppa);
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
% Empty override values count as unset: after load_parameters the study
% transducer carries every config_default field, mostly [] (elem_n,
% freq_hz, ...), which must not wipe the library geometry.
    dst = base;
    fn  = fieldnames(override);
    for i = 1:numel(fn)
        f = fn{i};
        if isstruct(override.(f)) && isfield(base, f) && isstruct(base.(f))
            dst.(f) = mergestruct_shallow(base.(f), override.(f));
        elseif (isempty(override.(f)) || isa(override.(f), 'yaml.Null')) && ...
                isfield(base, f) && ~isempty(base.(f))
            % keep the library value
        else
            dst.(f) = override.(f);
        end
    end
end
