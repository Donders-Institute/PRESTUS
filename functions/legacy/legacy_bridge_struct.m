function parameters = legacy_bridge_struct(parameters, caller)
% LEGACY_BRIDGE_STRUCT  Convert old key names set on a parameters struct after LOAD_PARAMETERS
%
% LOAD_PARAMETERS converts old key names in YAML files (UPGRADE_CONFIG_KEYS).
% Scripts written for earlier versions also set old names directly on the
% struct, e.g.
%
%   parameters = load_parameters('study.yaml');
%   parameters.code_type = 'cpp_gpu';                % now simulation.code_type
%   parameters.transducer = yaml.loadFile(f).transducer;   % v0.6 transducer fields
%
% Such fields would otherwise be ignored. This function moves every field that
% has an entry in the deprecation list (CONFIG_DEPRECATIONS) into a separate
% layer, converts that layer with UPGRADE_CONFIG_KEYS (which warns,
% prestus:renamedConfigKey, and stops for keys that cannot be converted
% exactly), and merges it back. The script's value wins over the value
% LOAD_PARAMETERS set under the current name, as the script set it later.
% Transducer fields are converted per transducer. A struct without old names
% is returned unchanged.
%
% Use as:
%   parameters = legacy_bridge_struct(parameters)
%   parameters = legacy_bridge_struct(parameters, caller)
%
% Input:
%   parameters - parameters struct (LOAD_PARAMETERS output, possibly edited)
%   caller     - name used in messages (default 'prestus_pipeline_start')
%
% Output:
%   parameters - the struct with current key names
%
% See also: LEGACY_BRIDGE_ENTRY, LEGACY_BRIDGE_LAYER, UPGRADE_CONFIG_KEYS, CONFIG_DEPRECATIONS

    if nargin < 2 || isempty(caller)
        caller = 'prestus_pipeline_start';
    end
    if ~isstruct(parameters) || ~isscalar(parameters)
        return
    end
    reg = config_deprecations();
    reg = reg(ismember({reg.kind}, {'rename', 'convert', 'removed', 'error'}));
    if isempty(reg)
        return
    end
    source = sprintf('parameters set on the struct before %s', caller);
    parameters = legacy_bridge_layer(parameters, source);

    is_tr = startsWith({reg.key}, 'transducer.');
    layer = struct();
    for e = reg(~is_tr)
        parts = strsplit(e.key, '.');
        if has_path(parameters, parts)
            layer = set_path(layer, parts, get_path(parameters, parts));
            parameters = rm_path(parameters, parts);
        end
    end

    tr_layers = {};
    n_tr = 0;
    if isfield(parameters, 'transducer') && isstruct(parameters.transducer)
        n_tr = numel(parameters.transducer);
        trs = num2cell(parameters.transducer);
        tr_layers = repmat({struct()}, 1, n_tr);
        for k = 1:n_tr
            for e = reg(is_tr)
                parts = strsplit(e.key, '.');
                parts = parts(2:end);
                if has_path(trs{k}, parts)
                    tr_layers{k} = set_path(tr_layers{k}, parts, get_path(trs{k}, parts));
                    trs{k} = rm_path(trs{k}, parts);
                end
            end
        end
    end
    has_tr = any(cellfun(@(s) ~isempty(fieldnames(s)), tr_layers));
    if isempty(fieldnames(layer)) && ~has_tr
        return
    end
    if has_tr
        layer.transducer = tr_layers;
    end

    up = upgrade_config_keys(layer, source);

    if isfield(up, 'transducer')
        up_tr = up.transducer;
        if isstruct(up_tr); up_tr = num2cell(up_tr); end
        up = rmfield(up, 'transducer');
        for k = 1:n_tr
            if k > numel(up_tr) || ~isstruct(up_tr{k}) || isempty(fieldnames(up_tr{k}))
                continue
            end
            old_type = '';
            if isfield(trs{k}, 'type'); old_type = trs{k}.type; end
            trs{k} = mergestruct(trs{k}, up_tr{k});
            if ~isempty(old_type)
                trs{k}.type = old_type;    % keep a type the transducer already had
            end
        end
        parameters.transducer = cells_to_struct_array(trs);
    end
    if ~isempty(fieldnames(up))
        parameters = mergestruct(parameters, up);
    end
end

%% Struct paths ---------------------------------------------------------------

function tf = has_path(s, parts)
    tf = true;
    for i = 1:numel(parts)
        if ~isstruct(s) || ~isscalar(s) || ~isfield(s, parts{i})
            tf = false;
            return
        end
        s = s.(parts{i});
    end
end

function v = get_path(s, parts)
    v = s;
    for i = 1:numel(parts)
        v = v.(parts{i});
    end
end

function s = set_path(s, parts, v)
    if numel(parts) == 1
        s.(parts{1}) = v;
        return
    end
    sub = struct();
    if isfield(s, parts{1}) && isstruct(s.(parts{1}))
        sub = s.(parts{1});
    end
    s.(parts{1}) = set_path(sub, parts(2:end), v);
end

function s = rm_path(s, parts)
    if numel(parts) == 1
        s = rmfield(s, parts{1});
        return
    end
    s.(parts{1}) = rm_path(s.(parts{1}), parts(2:end));
end

function arr = cells_to_struct_array(c)
% Struct array from scalar structs with possibly different fields
    names = {};
    for k = 1:numel(c)
        names = union(names, fieldnames(c{k}), 'stable');
    end
    for k = 1:numel(c)
        missing = setdiff(names, fieldnames(c{k}));
        for m = 1:numel(missing)
            c{k}.(missing{m}) = [];
        end
        c{k} = orderfields(c{k}, names);
    end
    arr = [c{:}];
end
