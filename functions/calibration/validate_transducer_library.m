function report = validate_transducer_library(library_path, equipment_path)
% VALIDATE_TRANSDUCER_LIBRARY  Validate calibration library against equipment configs
%
% For each YAML in library_path: checks consistency with equipment config
% (element count, depth coverage, phase monotonicity). For each YAML in
% equipment_path: reports whether a library entry exists.
%
% Use as:
%   report = validate_transducer_library()
%   report = validate_transducer_library(library_path)
%   report = validate_transducer_library(library_path, equipment_path)
%
% Input:
%   library_path   - path to transducer library folder
%                    (default: config/transducer/ under PRESTUS root)
%   equipment_path - path to equipment config folder
%                    (default: config/equipment/ under PRESTUS root)
%
% Output:
%   report - struct array with fields:
%     .serial       - equipment serial
%     .combo        - combo_name (library filename without .yaml)
%     .lib_file     - full path to library YAML, or '' if missing
%     .n_depths     - number of calibrated depths (0 if no library)
%     .depth_range  - [min max] mm, or [] if no depths
%     .monotonic    - logical; true if all element phases monotone across depths
%     .status       - 'ok' | 'no_library' | 'error'
%     .issues       - cell array of issue strings
%
% See also: UPDATE_TRANSDUCER_LIBRARY, CHECK_EQUIPMENT_CONFIG

    if nargin < 1 || isempty(library_path)
        library_path = fullfile(get_prestus_path(), 'config', 'transducer');
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
    all_lib_files = dir(fullfile(library_path, '*.yaml'));
    lib_names     = cellfun(@(f) strrep(f, '.yaml', ''), {all_lib_files.name}, 'UniformOutput', false);

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

        lib_file = fullfile(library_path, [entry_name '.yaml']);
        if ~isfile(lib_file)
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

        %% Analyse calibrated depths
        depths_cal = [];
        if isfield(lib, 'global_model') && isfield(lib.global_model, 'depths_ep_mm')
            depths_cal = lib.global_model.depths_ep_mm(:)';
        elseif isfield(lib, 'calibration') && isfield(lib.calibration, 'focal_depths')
            dk = fieldnames(lib.calibration.focal_depths);
            depths_cal = sort(cellfun(@(k) str2double(strrep(strrep(k,'f',''),'p','.')), dk));
        end

        rec.n_depths    = numel(depths_cal);
        rec.depth_range = [];
        rec.monotonic   = true;

        if rec.n_depths > 0
            rec.depth_range = [min(depths_cal) max(depths_cal)];

            % Check gaps > 10 mm
            if rec.n_depths > 1
                gaps = diff(sort(depths_cal));
                big_gaps = gaps(gaps > 10);
                if ~isempty(big_gaps)
                    rec.issues{end+1} = sprintf('Depth gaps > 10 mm: [%s] mm', num2str(big_gaps, '%.1f '));
                end
            end

            % Check global_model element count and monotonicity
            if isfield(lib, 'global_model') && isfield(lib.global_model, 'elem_phase_deg')
                phase_cell = lib.global_model.elem_phase_deg;
                if ~iscell(phase_cell)
                    phase_cell = num2cell(phase_cell, 2);
                end
                phase_mat = cell2mat(cellfun(@(r) r(:)', phase_cell, 'UniformOutput', false));

                % Check depth vs phase_mat row count
                if size(phase_mat, 1) ~= rec.n_depths
                    rec.issues{end+1} = sprintf( ...
                        'global_model: depths_ep_mm has %d entries but elem_phase_deg has %d rows', ...
                        rec.n_depths, size(phase_mat, 1));
                end

                % Check n_elem matches equipment
                if ~isnan(n_elem_equip) && size(phase_mat, 2) ~= n_elem_equip
                    rec.issues{end+1} = sprintf( ...
                        'global_model: %d phase columns but equipment has %d elements', ...
                        size(phase_mat, 2), n_elem_equip);
                end

                % Monotonicity check per element
                if size(phase_mat, 1) > 1
                    for el = 1:size(phase_mat, 2)
                        diffs = diff(phase_mat(:, el));
                        if ~all(diffs >= 0) && ~all(diffs <= 0)
                            rec.monotonic = false;
                            break;
                        end
                    end
                    if ~rec.monotonic
                        rec.issues{end+1} = 'Non-monotonic phase progression in global_model';
                    end
                end
            end
        else
            rec.issues{end+1} = 'No calibrated depths found';
        end

        rec.status = 'ok';
        if ~isempty(rec.issues)
            rec.status = 'warning';
        end
        report(end+1) = rec; %#ok<AGROW>
    end

    %% Print formatted table
    fprintf('\n%-40s %-8s %-8s %-14s %-10s %s\n', ...
        'Serial/Combo', 'N_depths', 'Range(mm)', 'Monotonic', 'Status', 'Issues');
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
