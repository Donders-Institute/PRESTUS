function report = prestus_config_upgrade(target_dir, options)
% PRESTUS_CONFIG_UPGRADE  Refresh a project's config_default.yaml from the toolbox
%
% prestus_config_init copies config_default.yaml into a project folder to pin
% its defaults. When the toolbox later adds keys, LOAD_PARAMETERS fills them
% from the toolbox default and warns (prestus:staleProjectDefault). This
% function compares the project copy with the current toolbox default and,
% unless DryRun is set, replaces the copy with the toolbox version after
% keeping a dated backup.
%
% Values that the project copy pins differently from the toolbox default are
% listed, so they can be moved into the study config (override YAML) if they
% were intentional.
%
% Use as:
%   report = prestus_config_upgrade(target_dir)
%   report = prestus_config_upgrade(target_dir, DryRun=true)
%
% Input:
%   target_dir     - folder holding the project config_default.yaml
%   options.DryRun - only report, do not change any file (default false)
%
% Output:
%   report - struct with fields
%              added    cellstr of keys missing from the project copy
%              removed  cellstr of keys no longer in the toolbox default
%              changed  table (key, project, toolbox) of differing values
%              backup   path of the backup ('' for a dry run / no change)
%
% See also: PRESTUS_CONFIG_INIT, LOAD_PARAMETERS

    arguments
        target_dir     (1,:) char
        options.DryRun (1,1) logical = false
    end

    this_dir     = fileparts(mfilename('fullpath'));
    prestus_root = fileparts(fileparts(this_dir));
    src_default  = fullfile(prestus_root, 'config', 'config_default.yaml');
    dst_default  = fullfile(target_dir, 'config_default.yaml');
    assert(isfile(dst_default), 'prestus_config_upgrade:notFound', ...
        'No config_default.yaml in %s (use prestus_config_init to create one).', target_dir);

    if ~exist('yaml.loadFile', 'file')
        addpath(genpath(fullfile(prestus_root, 'external')));
    end
    toolbox = yaml.loadFile(src_default, 'ConvertToArray', true);
    project = yaml.loadFile(dst_default, 'ConvertToArray', true);

    report = struct();
    report.added   = config_unknown_keys(toolbox, project);
    report.removed = config_unknown_keys(project, toolbox);
    report.changed = changed_values(project, toolbox);
    report.backup  = '';

    fprintf('Project defaults: %s\n', dst_default);
    fprintf('  %d key(s) missing from the project copy, %d no longer in the toolbox, %d pinned differently.\n', ...
        numel(report.added), numel(report.removed), height(report.changed));
    if height(report.changed) > 0
        fprintf('  Values the project copy pins differently from the toolbox default:\n');
        for i = 1:height(report.changed)
            fprintf('    %-45s project: %-20s toolbox: %s\n', report.changed.key{i}, ...
                report.changed.project{i}, report.changed.toolbox{i});
        end
        fprintf(['  To keep any of these, add them to your study config before upgrading; ' ...
                 'the refreshed copy holds the toolbox values.\n']);
    end

    if options.DryRun
        return
    end
    if isempty(report.added) && isempty(report.removed) && height(report.changed) == 0
        fprintf('  Already up to date.\n');
        return
    end
    report.backup = sprintf('%s.bak-%s', dst_default, char(datetime('now', 'Format', 'yyyyMMdd-HHmmss')));
    copyfile(dst_default, report.backup);
    copyfile(src_default, dst_default);
    fprintf('  Replaced with the toolbox default; backup: %s\n', report.backup);
end

function T = changed_values(project, toolbox)
% Leaf values present in both configs that differ.
    keys = {}; pv = {}; tv = {};
    walk(project, toolbox, '');
    T = table(keys(:), pv(:), tv(:), 'VariableNames', {'key', 'project', 'toolbox'});

    function walk(p, t, prefix)
        if ~isstruct(p) || ~isstruct(t) || isempty(p) || isempty(t)
            return
        end
        fn = fieldnames(p(1));
        for k = 1:numel(fn)
            f = fn{k};
            if ~isfield(t, f)
                continue
            end
            path = f;
            if ~isempty(prefix); path = [prefix '.' f]; end
            a = p(1).(f);
            b = t(1).(f);
            if isstruct(a) && isstruct(b)
                walk(a, b, path);
            elseif ~isequaln(a, b) && ~(isempty(a) && isempty(b))
                keys{end+1} = path; %#ok<AGROW>
                pv{end+1}   = to_text(a); %#ok<AGROW>
                tv{end+1}   = to_text(b); %#ok<AGROW>
            end
        end
    end
end

function s = to_text(v)
    if ischar(v) || isstring(v)
        s = char(v);
    elseif isnumeric(v) || islogical(v)
        s = mat2str(v);
    else
        s = ['<' class(v) '>'];
    end
    if numel(s) > 40
        s = [s(1:37) '...'];
    end
end
