function [parameters, resolved] = resolve_config_paths(parameters, base_dir)
% RESOLVE_CONFIG_PATHS  Make relative file and folder settings absolute
%
% A relative value of one of the file and folder settings listed in
% CONFIG_PATH_KEYS (below) is taken relative to base_dir; PATH_LOG_SETUP
% passes the PRESTUS root. Absolute values, empty values, and values starting
% with ~ or $ are left unchanged. Path templates read relative to path.anat
% (path.*_pattern) and Python environment names are not listed. Text, string,
% and cell values (e.g. a list of library folders) are all handled; for a
% struct array (several transducers) each element is resolved.
%
% Resolving here, rather than changing MATLAB's working folder, keeps a run
% from moving the caller's session and makes every path in the parameter
% snapshot absolute.
%
% Use as:
%   [parameters, resolved] = resolve_config_paths(parameters, base_dir)
%
% Input:
%   parameters - PRESTUS config
%   base_dir   - absolute folder that relative values are taken from
%
% Output:
%   parameters - with the listed relative values made absolute
%   resolved   - cell of 'key: old -> new' lines, one per changed value
%
% See also: PATH_LOG_SETUP, LOAD_PARAMETERS

arguments
    parameters (1,1) struct
    base_dir {mustBeTextScalar}
end
    base_dir = char(base_dir);
    resolved = {};
    keys = config_path_keys();
    for k = 1:numel(keys)
        [parameters, lines] = resolve_key(parameters, strsplit(keys{k}, '.'), base_dir, keys{k});
        resolved = [resolved, lines]; %#ok<AGROW>
    end
end

function keys = config_path_keys()
% Settings that name a file or folder (overlay.yaml kind: dir/file, plus the
% Localite folders, the transducer library, and the tool install folders).
    keys = { ...
        'path.sim', 'path.anat', 'path.seg', 'path.placement', 'path.pct', ...
        'path.transducer_library', 'path.acoustic_input', ...
        'path.localite_raw', 'path.localite_post', 'path.localite_t1', ...
        'path.localite_export', 'path.localite', ...
        'transducer.library_path', ...
        'transducer.matrix.matrix_shape.extract_from_file.file_path', ...
        'analysis.target_roi_mni', 'placement.localite.file', ...
        'placement.plantus.script_path', 'placement.plantus.connectome_wb_path', ...
        'segmentation.settings_file', ...
        'startup.simnibs_bin_path', 'startup.matlab_bin_path', 'startup.fsl_bin_path', ...
        'startup.ants_bin_path', 'startup.ants_simg_path'};
end

function [s, lines] = resolve_key(s, parts, base, key)
    lines = {};
    if ~isstruct(s) || ~isfield(s, parts{1})
        return;
    end
    for i = 1:numel(s)
        if numel(parts) == 1
            old = s(i).(parts{1});
            [new, changed] = resolve_value(old, base);
            if changed
                s(i).(parts{1}) = new;
                lines{end+1} = sprintf('%s: %s -> %s', key, value_text(old), value_text(new)); %#ok<AGROW>
            end
        else
            [s(i).(parts{1}), sub] = resolve_key(s(i).(parts{1}), parts(2:end), base, key);
            lines = [lines, sub]; %#ok<AGROW>
        end
    end
end

function [v, changed] = resolve_value(v, base)
    changed = false;
    if ischar(v)
        [v, changed] = resolve_one(v, base);
    elseif isstring(v)
        for i = 1:numel(v)
            if ismissing(v(i)); continue; end
            [p, ch] = resolve_one(char(v(i)), base);
            if ch; v(i) = string(p); changed = true; end
        end
    elseif iscell(v)
        for i = 1:numel(v)
            if ischar(v{i}) || (isstring(v{i}) && isscalar(v{i}) && ~ismissing(v{i}))
                [p, ch] = resolve_one(char(v{i}), base);
                if ch
                    if isstring(v{i}); p = string(p); end
                    v{i} = p; changed = true;
                end
            end
        end
    end
end

function [p, changed] = resolve_one(p, base)
    changed = false;
    if isempty(strtrim(p)) || any(p(1) == '~$') || is_absolute(p)
        return;
    end
    p = fullfile(base, p);
    changed = true;
end

function tf = is_absolute(p)
    tf = startsWith(p, {'/', '\'}) || ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'));
end

function t = value_text(v)
    if ischar(v) || (isstring(v) && isscalar(v))
        t = char(v);
    else
        t = strjoin(cellstr(string(v)), ', ');
    end
end
