function report = check_equipment_config(equipment_path)
% CHECK_EQUIPMENT_CONFIG  Validate equipment configuration YAML files
%
% Checks every transducer entry in equipment_path (as load_equipment_config
% reads it, so a unit includes what it takes from its model) for required fields,
% physical plausibility (ring overlap, focal range, frequency sanity), and
% dimension consistency. Prints a formatted report and returns a struct array.
%
% Use as:
%   report = check_equipment_config()
%   report = check_equipment_config(equipment_path)
%
% Input:
%   equipment_path - path to equipment config folder
%                    (default: config/equipment/ under PRESTUS root)
%
% Output:
%   report - struct array with fields:
%     .serial  - device serial string
%     .file    - full path to YAML file
%     .issues  - cell array of issue description strings
%     .valid   - logical; true if no issues found
%
% See also: LOAD_EQUIPMENT_CONFIG, VALIDATE_TRANSDUCER_LIBRARY

    if nargin < 1 || isempty(equipment_path)
        equipment_path = fullfile(get_prestus_path(), 'config', 'equipment');
    end

    % Check entries as PRESTUS uses them: a unit with the fields it takes from its model
    equip = load_equipment_config(equipment_path);
    serials = fieldnames(equip.trans);
    if isempty(serials)
        error('check_equipment_config: no transducer entries found in:\n  %s', equipment_path);
    end

    report = struct('serial', {}, 'file', {}, 'issues', {}, 'valid', {}, 'entry', {});

    for fi = 1:numel(serials)
        d          = equip.trans.(serials{fi});
        rec        = struct();
        rec.file   = equip.files.(serials{fi});
        rec.serial = serials{fi};
        rec.issues = {};
        rec.entry  = d;

        if isfield(d, 'serial')
            rec.serial = d.serial;
        else
            rec.issues{end+1} = 'Missing required field: serial';
        end

        %% Required top-level fields
        if ~isfield(d, 'type')
            rec.issues{end+1} = 'Missing required field: type';
        end

        %% Transducer sub-struct
        if ~isfield(d, 'transducer')
            rec.issues{end+1} = 'Missing required field: transducer';
            rec.valid = false;
            report(end+1) = rec; %#ok<AGROW>
            continue;
        end
        tr = d.transducer;

        if ~isfield(tr, 'freq_hz')
            rec.issues{end+1} = 'Missing required field: transducer.freq_hz';
        else
            freq = tr.freq_hz;
            if freq < 100e3 || freq > 5e6
                rec.issues{end+1} = sprintf('transducer.freq_hz = %.0f Hz is outside sane range [100 kHz, 5 MHz]', freq);
            end
        end

        %% Matrix-specific checks
        if ~isfield(tr, 'annular') && isfield(tr, 'matrix')
            rec.issues = [rec.issues, check_matrix(tr.matrix)];
            rec.valid = isempty(rec.issues);
            report(end+1) = rec; %#ok<AGROW>
            continue;
        end

        %% Annular-specific checks
        if ~isfield(tr, 'annular')
            rec.issues{end+1} = 'Missing required field: transducer.annular or transducer.matrix';
            rec.valid = isempty(rec.issues);
            report(end+1) = rec; %#ok<AGROW>
            continue;
        end
        ann = tr.annular;

        required_ann = {'elem_n', 'elem_id_mm', 'elem_od_mm', 'curv_radius_mm', 'dist_geom_ep_mm'};
        for ri = 1:numel(required_ann)
            if ~isfield(ann, required_ann{ri})
                rec.issues{end+1} = sprintf('Missing required field: transducer.annular.%s', required_ann{ri});
            end
        end

        % Dimension consistency
        if isfield(ann, 'elem_n') && isfield(ann, 'elem_id_mm') && isfield(ann, 'elem_od_mm')
            n = ann.elem_n;
            if numel(ann.elem_id_mm) ~= n
                rec.issues{end+1} = sprintf('elem_id_mm has %d entries but elem_n = %d', numel(ann.elem_id_mm), n);
            end
            if numel(ann.elem_od_mm) ~= n
                rec.issues{end+1} = sprintf('elem_od_mm has %d entries but elem_n = %d', numel(ann.elem_od_mm), n);
            end

            % Each ring: inner diameter below its outer diameter
            if numel(ann.elem_od_mm) == n && numel(ann.elem_id_mm) == n
                inverted = ann.elem_od_mm(:)' <= ann.elem_id_mm(:)';
                if any(inverted)
                    rec.issues{end+1} = sprintf('Ring(s) %s: elem_od_mm must exceed elem_id_mm', ...
                        num2str(find(inverted)));
                end
            end
            % No ring overlap: od(i) < id(i+1)
            if numel(ann.elem_od_mm) == n && numel(ann.elem_id_mm) == n && n > 1
                od = ann.elem_od_mm(:)';
                id = ann.elem_id_mm(:)';
                overlaps = od(1:end-1) > id(2:end);     % touching sub-rings (ring profile) are allowed
                if any(overlaps)
                    idx_str = num2str(find(overlaps));
                    rec.issues{end+1} = sprintf('Ring overlap detected at gap(s): %s', idx_str);
                end
            end
        end

        % Physical plausibility for curvature/geometry
        if isfield(ann, 'curv_radius_mm') && isfield(ann, 'elem_od_mm')
            if ann.curv_radius_mm <= max(ann.elem_od_mm(:)) / 2
                rec.issues{end+1} = sprintf( ...
                    'curv_radius_mm (%.1f) <= max(elem_od_mm)/2 (%.1f)', ...
                    ann.curv_radius_mm, max(ann.elem_od_mm(:)) / 2);
            end
        end

        if isfield(ann, 'curv_radius_mm') && isfield(ann, 'dist_geom_ep_mm')
            if ann.dist_geom_ep_mm > ann.curv_radius_mm
                rec.issues{end+1} = sprintf( ...
                    'dist_geom_ep_mm (%.1f) > curv_radius_mm (%.1f)', ...
                    ann.dist_geom_ep_mm, ann.curv_radius_mm);
            end
        end

        % Optional focal range
        if isfield(d, 'min_foc') && ~isempty(d.min_foc) && d.min_foc < 0
            rec.issues{end+1} = sprintf('min_foc (%.1f) < 0', d.min_foc);
        end
        % min_foc/max_foc are distances from the exit plane; steering beyond
        % the geometric focus is allowed, so they are not bounded by curv_radius_mm.
        if isfield(d, 'min_foc') && isfield(d, 'max_foc') && ~isempty(d.min_foc) && ...
                ~isempty(d.max_foc) && d.min_foc > d.max_foc
            rec.issues{end+1} = sprintf('min_foc (%.1f) > max_foc (%.1f)', d.min_foc, d.max_foc);
        end
        if isfield(d, 'status') && ~any(strcmp(char(d.status), {'active', 'deprecated', 'excluded'}))
            rec.issues{end+1} = sprintf('status ''%s'' is not active, deprecated or excluded', char(d.status));
        end

        rec.valid = isempty(rec.issues);
        report(end+1) = rec; %#ok<AGROW>
    end

    %% Print formatted report
    n_elem_str = @(r) '?';
    fprintf('\n%-35s %-8s %-8s %-10s %-10s %s\n', ...
        'Serial', 'Type', 'N_elem', 'Freq_kHz', 'Status', 'Issues');
    fprintf('%s\n', repmat('-', 1, 100));
    for ri = 1:numel(report)
        r   = report(ri);
        d   = r.entry;
        type_str = '';
        if isfield(d, 'type'), type_str = d.type; end
        n_el_str = n_elem_str(r);
        freq_str = '';
        if isfield(d, 'transducer') && isfield(d.transducer, 'freq_hz')
            freq_str = sprintf('%.0f', d.transducer.freq_hz / 1e3);
        end
        if isfield(d, 'transducer') && isfield(d.transducer, 'annular') && ...
                isfield(d.transducer.annular, 'elem_n')
            n_el_str = num2str(d.transducer.annular.elem_n);
        elseif isfield(d, 'transducer') && isfield(d.transducer, 'matrix') && isfield(d, 'n_elem')
            n_el_str = num2str(d.n_elem);
        end
        status_str = 'ok';
        if ~r.valid, status_str = 'ISSUES'; end
        issue_str = strjoin(r.issues, '; ');
        fprintf('%-35s %-8s %-8s %-10s %-10s %s\n', ...
            r.serial, type_str, n_el_str, freq_str, status_str, issue_str);
    end
    fprintf('\n');
end

function issues = check_matrix(m)
% Required fields and plausibility for a transducer.matrix block
    issues = {};
    req = {'elem_shape', 'elem_width_mm', 'elem_height_mm', 'outer_diameter_mm', 'is_curved', 'matrix_shape'};
    for k = 1:numel(req)
        if ~isfield(m, req{k}) || isempty(m.(req{k}))
            issues{end+1} = sprintf('Missing required field: transducer.matrix.%s', req{k}); %#ok<AGROW>
        end
    end
    if ~isempty(issues); return; end
    if ~any(strcmp(char(m.elem_shape), {'rect', 'disc', 'bowl'}))
        issues{end+1} = sprintf('matrix.elem_shape ''%s'' is not rect, disc or bowl', char(m.elem_shape));
    end
    if m.elem_width_mm <= 0 || m.elem_height_mm <= 0
        issues{end+1} = 'matrix element width and height must be positive';
    end
    if m.is_curved
        if ~isfield(m, 'curv_radius_mm') || isempty(m.curv_radius_mm)
            issues{end+1} = 'Missing required field: transducer.matrix.curv_radius_mm (is_curved: true)';
        elseif m.curv_radius_mm <= m.outer_diameter_mm / 2
            issues{end+1} = sprintf('matrix.curv_radius_mm (%.1f) <= outer_diameter_mm/2 (%.1f)', ...
                m.curv_radius_mm, m.outer_diameter_mm / 2);
        elseif isfield(m, 'dist_geom_ep_mm') && ~isempty(m.dist_geom_ep_mm) && m.dist_geom_ep_mm > m.curv_radius_mm
            issues{end+1} = sprintf('matrix.dist_geom_ep_mm (%.1f) > curv_radius_mm (%.1f)', ...
                m.dist_geom_ep_mm, m.curv_radius_mm);
        end
    end
    ms = m.matrix_shape;
    if ~isfield(ms, 'type') || ~any(strcmp(char(ms.type), {'define_here', 'extract_from_file'}))
        issues{end+1} = 'matrix.matrix_shape.type must be define_here or extract_from_file';
    elseif strcmp(char(ms.type), 'define_here')
        gs = ms.define_here.grid_shape;
        switch char(gs.type)
            case 'rect'
                f = {'elem_n_row', 'elem_n_col', 'elem_spacing_height_mm', 'elem_spacing_width_mm', 'sparsity_factor'};
            case {'fibonacci', 'fermat'}
                f = {'elem_n'};
            otherwise
                issues{end+1} = sprintf('matrix grid_shape.type ''%s'' is not rect, fibonacci or fermat', char(gs.type));
                return;
        end
        for k = 1:numel(f)
            if ~isfield(gs.(char(gs.type)), f{k})
                issues{end+1} = sprintf('Missing grid_shape.%s.%s', char(gs.type), f{k}); %#ok<AGROW>
            end
        end
    end
end
