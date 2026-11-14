function report = validate_transducer_library(library_path, equipment_path)
% VALIDATE_TRANSDUCER_LIBRARY  Validate calibration library against equipment configs
%
% For each YAML in library_path: checks consistency with the equipment config
% (static offset and refinement element counts against the current geometry,
% amplitude, recorded forward model, legacy per-depth forms). For each YAML in
% equipment_path: reports whether a library entry exists.
%
% Use as:
%   report = validate_transducer_library()
%   report = validate_transducer_library(library_path)
%   report = validate_transducer_library(library_path, equipment_path)
%
% Input:
%   library_path   - path to transducer library folder
%                    (default: config/calibration/ under PRESTUS root)
%   equipment_path - path to equipment config folder
%                    (default: config/equipment/ under PRESTUS root)
%
% Output:
%   report - struct array with fields:
%     .serial       - equipment serial
%     .combo        - combo_name (library filename without .yaml)
%     .lib_file     - full path to library YAML, or '' if missing
%     .n_depths     - number of refined depths (0 for a static-only entry)
%     .depth_range  - [min max] mm of refined depths, or []
%     .monotonic    - kept for compatibility (always true)
%     .status       - 'ok' | 'warning' | 'no_library' | 'error'
%     .issues       - cell array of issue strings
%
% See also: UPDATE_TRANSDUCER_LIBRARY, CHECK_EQUIPMENT_CONFIG

    if nargin < 1 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'calibration');
    end
    if nargin < 2 || isempty(equipment_path)
        equipment_path = fullfile(get_prestus_path(), 'config', 'equipment');
    end

    %% Load equipment config
    eq = load_equipment_config(equipment_path);

    %% Build list of all known combos from equipment
    combo_keys = fieldnames(eq.combos);
    tran_keys  = fieldnames(eq.trans);

    % Also include trans-only keys (no combo) as potential library entries
    [~, ~, lib_index] = library_entry_file(library_path);   % flat and nested entries
    lib_names = {lib_index.key};

    % Build report entries for all equipment combos
    all_serials = union(combo_keys, tran_keys);
    all_serials = union(all_serials, lib_names(:));

    report = struct('serial', {}, 'combo', {}, 'lib_file', {}, ...
                    'n_depths', {}, 'depth_range', {}, 'monotonic', {}, ...
                    'status', {}, 'issues', {});

    for ei = 1:numel(all_serials)
        entry_name = all_serials{ei};
        rec        = struct();
        rec.serial = entry_name;
        rec.combo  = entry_name;
        rec.issues = {};

        [lib_file, has_entry] = library_entry_file(library_path, entry_name);
        if ~has_entry
            rec.lib_file    = '';
            rec.n_depths    = 0;
            rec.depth_range = [];
            rec.monotonic   = false;
            rec.status      = 'no_library';
            rec.issues{end+1} = 'No library file found';
            report(end+1) = rec; %#ok<AGROW>
            continue;
        end

        rec.lib_file = lib_file;

        lib = yaml.loadFile(lib_file, 'ConvertToArray', true);

        % Determine tran_serial
        if isfield(lib, 'meta') && isfield(lib.meta, 'tran_serial')
            tran_serial = lib.meta.tran_serial;
        else
            tran_serial = entry_name;
            rec.issues{end+1} = 'meta.tran_serial missing';
        end
        rec.serial = tran_serial;

        % Check equipment entry exists
        n_elem_equip = NaN;
        if isfield(eq.trans, tran_serial)
            n_elem_equip = eq.trans.(tran_serial).transducer.annular.elem_n;
        elseif isfield(eq.combos, tran_serial)
            % combo — tran_serial is the combo itself; get n_elem from tran
            ts = eq.combos.(tran_serial).tran_serial;
            if isfield(eq.trans, ts)
                n_elem_equip = eq.trans.(ts).transducer.annular.elem_n;
            end
        else
            rec.issues{end+1} = sprintf('Equipment entry for ''%s'' not found', tran_serial);
        end

        % A ring profile splits every ring into n_sub sub-rings (RING_PROFILE_TRANSDUCER)
        if ~isnan(n_elem_equip) && isfield(lib, 'global_model') && isfield(lib.global_model, 'ring_profile') ...
                && isstruct(lib.global_model.ring_profile) && isfield(lib.global_model.ring_profile, 'n_sub')
            n_elem_equip = n_elem_equip * double(lib.global_model.ring_profile.n_sub);
        end

        %% Static offset --------------------------------------------------
        gm = struct();
        if isfield(lib, 'global_model') && isstruct(lib.global_model); gm = lib.global_model; end
        has_static = isfield(gm, 'elem_phase_correction_deg') && ~isempty(gm.elem_phase_correction_deg);
        is_error   = false;
        if has_static && ~isnan(n_elem_equip) && numel(gm.elem_phase_correction_deg) ~= n_elem_equip
            rec.issues{end+1} = sprintf(['global_model.elem_phase_correction_deg has %d elements but ' ...
                'equipment has %d — entry is stale for the current geometry, rebuild it'], ...
                numel(gm.elem_phase_correction_deg), n_elem_equip);
            is_error = true;
        end
        if ~isfield(gm, 'amplitude_scaling') || isempty(gm.amplitude_scaling)
            rec.issues{end+1} = 'No global_model.amplitude_scaling';
        end
        if ~isfield(gm, 'forward_model') && has_static
            rec.issues{end+1} = 'Forward model used for the fit not recorded (global_model.forward_model)';
        end

        %% Refinements (optional) -----------------------------------------
        depths_cal = [];
        if isfield(lib, 'refinement') && isstruct(lib.refinement) && isfield(lib.refinement, 'depths_ep_mm')
            depths_cal = double(lib.refinement.depths_ep_mm(:)');
            offs = lib.refinement.elem_phase_correction_deg;
            if ~iscell(offs)
                if numel(depths_cal) == 1; offs = {offs}; else; offs = num2cell(offs, 2); end
            end
            if numel(offs) ~= numel(depths_cal)
                rec.issues{end+1} = sprintf('refinement: %d depths but %d offset rows', numel(depths_cal), numel(offs));
                is_error = true;
            elseif ~isnan(n_elem_equip) && any(cellfun(@numel, offs) ~= n_elem_equip)
                rec.issues{end+1} = sprintf('refinement offsets do not have %d elements', n_elem_equip);
                is_error = true;
            end
        end
        if isfield(gm, 'depths_ep_mm') && ~isempty(gm.depths_ep_mm)
            rec.issues{end+1} = 'Legacy per-depth phase table in global_model (read as refinement)';
            depths_cal = [depths_cal, double(gm.depths_ep_mm(:)')];
        end
        if isfield(lib, 'calibration') && isfield(lib.calibration, 'focal_depths')
            rec.issues{end+1} = 'Legacy calibration.focal_depths entries (read as refinement)';
        end
        if ~has_static && isempty(depths_cal)
            rec.issues{end+1} = 'No offset stored — geometric steering only';
        end

        depths_cal      = unique(depths_cal);
        rec.n_depths    = numel(depths_cal);
        rec.depth_range = [];
        if ~isempty(depths_cal); rec.depth_range = [min(depths_cal) max(depths_cal)]; end
        rec.monotonic   = true;   % offsets are not expected to be monotonic in depth

        if is_error
            rec.status = 'error';
        elseif isempty(rec.issues)
            rec.status = 'ok';
        else
            rec.status = 'warning';
        end
        report(end+1) = rec; %#ok<AGROW>
    end

    %% Print formatted table
    fprintf('\n%-40s %-8s %-8s %-14s %-10s %s\n', ...
        'Serial/Combo', 'N_refine', 'Range(mm)', 'Monotonic', 'Status', 'Issues');
    fprintf('%s\n', repmat('-', 1, 100));
    for ri = 1:numel(report)
        r = report(ri);
        if isempty(r.depth_range)
            range_str = '—';
        else
            range_str = sprintf('%.0f-%.0f', r.depth_range(1), r.depth_range(2));
        end
        mono_str   = 'N/A';
        if r.n_depths > 1
            mono_str = ternary_str(r.monotonic, 'yes', 'NO');
        end
        issue_str = strjoin(r.issues, '; ');
        fprintf('%-40s %-8d %-8s %-14s %-10s %s\n', ...
            r.combo, r.n_depths, range_str, mono_str, r.status, issue_str);
    end
    fprintf('\n');
end

function s = ternary_str(cond, a, b)
    if cond, s = a; else, s = b; end
end
