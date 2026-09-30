function parameters = apply_config_overrides(parameters, overrides)
% APPLY_CONFIG_OVERRIDES  Set config values by dotted key, e.g. from a study manifest row
%
% Each override is a (key, value) pair. The key is a dotted config path; a
% struct-array level can be indexed with parentheses:
%
%   'subject_id'                      top-level value
%   'placement.localite.target'       nested value
%   'transducer.target_isppa_wcm2'    transducer(1), allowed only while the
%                                     config defines a single transducer
%   'transducer(2).serial'            second transducer of a multi-transducer setup
%
% Values are assigned as given (already typed). Keys that are absent from
% the merged config are still set, with a warning (typo guard; some valid
% advanced keys are commented out in config_default.yaml).
%
% Use as:
%   parameters = apply_config_overrides(parameters, overrides)
%
% Input:
%   parameters - PRESTUS config struct (after merging default + study config)
%   overrides  - N×2 cell {key, value; ...}, or a struct array with fields
%                .key and .value
%
% Output:
%   parameters - updated config struct
%
% See also: LOAD_PARAMETERS, STUDY_READ_MANIFEST, STUDY_ROW_PARAMETERS

    if isempty(overrides)
        return
    end
    if isstruct(overrides)
        overrides = [{overrides.key}', {overrides.value}'];
    end
    if ~iscell(overrides) || size(overrides, 2) ~= 2
        error('apply_config_overrides:badInput', ...
            'overrides must be an N×2 cell array {key, value; ...}.');
    end

    for i = 1:size(overrides, 1)
        key = char(overrides{i, 1});
        toks = parse_key(key);
        if ~path_exists(parameters, toks)
            warning('apply_config_overrides:unknownKey', ...
                'Override "%s" is not a key of the loaded config; it is set anyway (check for typos).', key);
        end
        parameters = set_in(parameters, toks, overrides{i, 2}, key);
    end
end

% -------------------------------------------------------------------------
function toks = parse_key(key)
% Split 'a.b(2).c' into a struct array with fields name / idx.
    parts = strsplit(strtrim(key), '.');
    toks = struct('name', {}, 'idx', {});
    for p = 1:numel(parts)
        % Two patterns: MATLAB drops tokens nested in an optional group.
        m = regexp(parts{p}, '^([A-Za-z]\w*)\((\d+)\)$', 'tokens', 'once');
        idx = [];
        if ~isempty(m)
            idx = str2double(m{2});
        else
            m = regexp(parts{p}, '^([A-Za-z]\w*)$', 'tokens', 'once');
        end
        if isempty(m)
            error('apply_config_overrides:badKey', 'Invalid config key "%s".', key);
        end
        toks(end+1) = struct('name', m{1}, 'idx', idx); %#ok<AGROW>
    end
end

function tf = path_exists(s, toks)
    tf = true;
    for t = 1:numel(toks)
        if ~isstruct(s) || ~isfield(s, toks(t).name)
            tf = false;
            return
        end
        s = s(1).(toks(t).name);
        if ~isempty(toks(t).idx) && isstruct(s)
            if toks(t).idx > numel(s)
                tf = false;
                return
            end
            s = s(toks(t).idx);
        elseif isstruct(s) && ~isempty(s)
            s = s(1);
        end
    end
end

function s = set_in(s, toks, value, key)
    name = toks(1).name;
    idx  = toks(1).idx;
    if numel(toks) == 1
        if ~isempty(idx)
            error('apply_config_overrides:badKey', ...
                '"%s": only struct levels can be indexed.', key);
        end
        s.(name) = value;
        return
    end

    if ~isfield(s, name) || isempty(s.(name)) || isa(s.(name), 'yaml.Null')
        child = struct();
    else
        child = s.(name);
    end
    if ~isstruct(child)
        error('apply_config_overrides:notStruct', ...
            '"%s": "%s" is not a struct level in the loaded config.', key, name);
    end

    n = numel(child);
    if isempty(idx)
        if n > 1
            error('apply_config_overrides:ambiguous', ...
                ['"%s": "%s" has %d elements; index it, e.g. "%s(2).%s".'], ...
                key, name, n, name, strjoin({toks(2:end).name}, '.'));
        end
        k = 1;
    else
        k = idx;
        if k < 1 || k > n
            error('apply_config_overrides:indexRange', ...
                '"%s": %s(%d) does not exist; the config defines %d element(s).', key, name, k, n);
        end
    end

    elem = set_in(child(k), toks(2:end), value, key);
    new_fields = setdiff(fieldnames(elem), fieldnames(child));
    for f = 1:numel(new_fields)
        [child.(new_fields{f})] = deal([]);
    end
    child(k) = orderfields(elem, child(1));
    s.(name) = child;
end
