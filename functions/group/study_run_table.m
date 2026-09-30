function runs = study_run_table(parameters, manifest)
% STUDY_RUN_TABLE  Resolve each manifest row to its output location (subject, affix, medium)
%
% Derives, without running anything, where each manifest row writes its
% outputs, so that the launcher and the group report agree on file names:
%
%   affix   - io.output_affix of the run. An explicit io.output_affix
%             column (or per-row config value) wins. Otherwise a non-empty
%             condition is appended to the study affix label with '-'
%             (study affix 'liberal' + condition 'M1' -> '_desc-liberal-M1';
%             no study affix -> '_desc-M1'). Sanitised like LOAD_PARAMETERS.
%   medium  - simulation.medium (row override, else study value)
%   sim_dir - path.sim (row override, else study value)
%
% Errors when two included rows would write to the same outputs (same
% subject, affix, medium and path.sim): give them distinct conditions.
%
% Use as:
%   runs = study_run_table(parameters, manifest)
%
% Input:
%   parameters - study config loaded with LOAD_PARAMETERS (no row applied)
%   manifest   - table from STUDY_READ_MANIFEST
%
% Output:
%   runs - table with columns row, subject_id, condition, include, affix,
%          medium, sim_dir, label (display name 'sub-003 · M1'), config,
%          overrides (copied from the manifest)
%
% See also: STUDY_READ_MANIFEST, STUDY_ROW_PARAMETERS, SANITIZE_OUTPUT_AFFIX

arguments
    parameters (1,1) struct
    manifest table
end

    n = height(manifest);
    affix   = strings(n, 1);
    medium  = strings(n, 1);
    sim_dir = strings(n, 1);
    label   = strings(n, 1);

    study_label = regexprep(char(parameters.io.output_affix), '^_desc-', '');

    for r = 1:n
        ov = manifest.overrides{r};
        row_cfg = struct();
        if strlength(manifest.config(r)) > 0
            if ~isfile(manifest.config(r))
                error('study_run_table:configMissing', ...
                    'Manifest row %d: config file not found: %s', manifest.row(r), manifest.config(r));
            end
            row_cfg = yaml.loadFile(char(manifest.config(r)), "ConvertToArray", true);
        end

        [has_affix, a] = row_value(ov, row_cfg, {'io', 'output_affix'});
        cond = char(manifest.condition(r));
        if has_affix
            lbl = char(string(a));
        elseif ~isempty(cond)
            if isempty(study_label)
                lbl = cond;
            else
                lbl = [study_label '-' cond];
            end
        else
            lbl = study_label;
        end
        affix(r) = string(sanitize_output_affix(lbl));

        [has_m, m] = row_value(ov, row_cfg, {'simulation', 'medium'});
        if has_m, medium(r) = string(m); else, medium(r) = string(parameters.simulation.medium); end
        [has_s, s] = row_value(ov, row_cfg, {'path', 'sim'});
        if has_s, sim_dir(r) = string(s); else, sim_dir(r) = string(parameters.path.sim); end

        label(r) = string(sprintf('sub-%03d', manifest.subject_id(r)));
        if ~isempty(cond)
            label(r) = label(r) + " · " + string(cond);
        end
    end

    runs = table(manifest.row, manifest.subject_id, manifest.condition, manifest.include, ...
        affix, medium, sim_dir, label, manifest.config, manifest.overrides, ...
        'VariableNames', {'row', 'subject_id', 'condition', 'include', 'affix', 'medium', ...
                          'sim_dir', 'label', 'config', 'overrides'});

    % Two included rows must never write to the same files.
    inc = runs(runs.include, :);
    keys = compose('%d|%s|%s|%s', inc.subject_id, inc.affix, inc.medium, inc.sim_dir);
    [u, ~, j] = unique(keys, 'stable');
    if numel(u) < numel(keys)
        dup = find(accumarray(j, 1) > 1, 1);
        rows = inc.row(j == dup);
        error('study_run_table:duplicateRun', ...
            ['Manifest rows %s write to the same outputs (sub-%03d, affix "%s"). ' ...
             'Give them distinct condition labels.'], ...
            strjoin(compose('%d', rows), ', '), inc.subject_id(find(j == dup, 1)), ...
            inc.affix(find(j == dup, 1)));
    end
end

% -------------------------------------------------------------------------
function [found, v] = row_value(ov, row_cfg, path)
% Value of a config key for this row: manifest column first, then row config.
    key = strjoin(path, '.');
    found = false; v = [];
    if ~isempty(ov)
        hit = find(strcmp(ov(:, 1), key), 1, 'last');
        if ~isempty(hit)
            found = true; v = ov{hit, 2};
            return
        end
    end
    s = row_cfg;
    for p = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{p}); return; end
        s = s.(path{p});
    end
    found = true; v = s;
end
