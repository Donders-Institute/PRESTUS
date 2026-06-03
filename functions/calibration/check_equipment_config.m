function report = check_equipment_config(equipment_path)
% CHECK_EQUIPMENT_CONFIG  Validate equipment configuration YAML files
%
% Iterates over all *.yaml files in equipment_path and checks required fields,
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

    files = dir(fullfile(equipment_path, '*.yaml'));
    if isempty(files)
        error('check_equipment_config: no YAML files found in:\n  %s', equipment_path);
    end

    report = struct('serial', {}, 'file', {}, 'issues', {}, 'valid', {});

    for fi = 1:numel(files)
        fpath = fullfile(equipment_path, files(fi).name);

        rec        = struct();
        rec.file   = fpath;
        rec.serial = strrep(files(fi).name, '.yaml', '');
        rec.issues = {};

        d = yaml.loadFile(fpath, 'ConvertToArray', true);

        % Skip non-transducer entries (e.g. equipment_info.yaml)
        if ~isfield(d, 'type') || ~strcmp(d.type, 'transducer')
            continue;
        end

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

        %% Annular-specific checks
        if ~isfield(tr, 'annular')
            rec.issues{end+1} = 'Missing required field: transducer.annular (expected for type=transducer)';
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

            % No ring overlap: od(i) < id(i+1)
            if numel(ann.elem_od_mm) == n && numel(ann.elem_id_mm) == n && n > 1
                od = ann.elem_od_mm(:)';
                id = ann.elem_id_mm(:)';
                overlaps = od(1:end-1) >= id(2:end);
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
        if isfield(d, 'max_foc') && isfield(ann, 'curv_radius_mm') && ...
                ~isempty(d.max_foc) && d.max_foc > ann.curv_radius_mm
            rec.issues{end+1} = sprintf('max_foc (%.1f) > curv_radius_mm (%.1f)', ...
                d.max_foc, ann.curv_radius_mm);
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
        d   = yaml.loadFile(r.file, 'ConvertToArray', true);
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
        end
        status_str = 'ok';
        if ~r.valid, status_str = 'ISSUES'; end
        issue_str = strjoin(r.issues, '; ');
        fprintf('%-35s %-8s %-8s %-10s %-10s %s\n', ...
            r.serial, type_str, n_el_str, freq_str, status_str, issue_str);
    end
    fprintf('\n');
end
