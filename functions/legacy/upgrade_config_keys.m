function cfg = upgrade_config_keys(cfg, source)
% UPGRADE_CONFIG_KEYS  Rename old config keys (PRESTUS v0.5.0 onwards) in one config layer to their current names
%
% Applies the 'deprecated:' list of the Nexus overlay (CONFIG_DEPRECATIONS) to
% one config layer before it is merged: renamed keys are moved to their current
% name, exact value conversions are applied, keys without a replacement are
% dropped, and keys that cannot be converted exactly stop the run with a
% message naming the key to set instead. Converted keys are listed in one
% warning (prestus:renamedConfigKey). A layer that sets an old key and its
% current name to different values stops the run (prestus:configKeyConflict).
%
% transducer.* entries apply to every transducer of the layer. A transducer
% given with v0.5.0 fields (transducer.n_elements, ...) and without type or
% serial gets type 'annular'.
%
% Use as:
%   cfg = upgrade_config_keys(cfg)
%   cfg = upgrade_config_keys(cfg, source)
%   overrides = upgrade_config_keys(overrides, source)   % N x 2 {key, value} cell
%
% Input:
%   cfg    - one config layer (struct, as read from YAML), or an N x 2 cell of
%            dotted {key, value} overrides (LOAD_PARAMETERS)
%   source - name of the layer for messages (default 'config')
%
% Output:
%   cfg    - the layer with current key names
%
% See also: CONFIG_DEPRECATIONS, LOAD_PARAMETERS, PRESTUS_CONFIG_UPGRADE, CALIBRATION_TRANSDUCER

    if nargin < 2 || isempty(source)
        source = 'config';
    end
    source = char(string(source));
    reg = config_deprecations();
    reg = reg(ismember({reg.kind}, {'rename', 'convert', 'removed', 'error'}));
    if isempty(reg)
        return
    end

    if iscell(cfg)
        cfg = upgrade_overrides(cfg, reg, source);
        return
    end
    if ~isstruct(cfg) || ~isscalar(cfg) || isempty(fieldnames(cfg))
        return
    end
    cfg = legacy_bridge_layer(cfg, source);   % old keys the layer makes redundant (warns)

    log = struct('done', {{}}, 'errors', {{}});
    is_tr = startsWith({reg.key}, 'transducer.');
    [cfg, log] = apply_entries(cfg, reg(~is_tr), '', log);

    if isfield(cfg, 'transducer') && ~isempty(cfg.transducer) && any(is_tr)
        [trs, was_struct] = to_cell(cfg.transducer);
        for k = 1:numel(trs)
            if ~isstruct(trs{k})
                continue
            end
            n_before = numel(log.done);
            prefix = 'transducer';
            if numel(trs) > 1
                prefix = sprintf('transducer(%d)', k);
            end
            [trs{k}, log] = apply_entries(trs{k}, reg(is_tr), prefix, log);
            % v0.5.0 transducers had no type; the converted fields are annular
            added_annular = any(contains(log.done(n_before+1:end), '.annular.'));
            if added_annular && ~isfield(trs{k}, 'type') && ~isfield(trs{k}, 'serial')
                trs{k}.type = 'annular';
                log.done{end+1} = sprintf('%s.type set to ''annular'' (v0.5.0 transducers were annular)', prefix);
            end
        end
        cfg.transducer = from_cell(trs, was_struct);
    end

    report(log, source);
end

%% Entries --------------------------------------------------------------------

function [s, log] = apply_entries(s, reg, prefix, log)
% Renames first, so that converters see current names (e.g. data_path -> path.anat)
    order = {'rename', 'removed', 'convert', 'error'};
    for o = 1:numel(order)
        for e = reg(strcmp({reg.kind}, order{o}))
            parts = strsplit(e.key, '.');
            if ~isempty(prefix)
                parts = parts(2:end);      % 'transducer.' entries act on one transducer
            end
            if ~has_path(s, parts)
                continue
            end
            old = label(prefix, parts);
            val = get_path(s, parts);
            switch e.kind
                case 'rename'
                    new_parts = strip_tr(strsplit(e.use{1}, '.'), prefix);
                    [s, log] = move(s, parts, new_parts, val, e, prefix, log);
                case 'removed'
                    if e.has_keep && ~same_value(val, e.keep)
                        log.errors{end+1} = sprintf('%s = %s is no longer supported: %s', ...
                            old, value_str(val), e.note);
                    else
                        s = rm_path(s, parts);
                        log.done{end+1} = sprintf('%s removed (since %s): %s', old, e.since, e.note);
                    end
                case 'convert'
                    [s, log] = convert(e.convert, s, parts, val, e, prefix, log);
                case 'error'
                    log.errors{end+1} = sprintf('%s: %s', old, e.note);
            end
        end
    end
end

function [s, log] = move(s, parts, new_parts, val, e, prefix, log)
    s = rm_path(s, parts);
    [s, log] = put(s, new_parts, val, label(prefix, parts), e, prefix, log);
end

function [s, log] = put(s, new_parts, val, old, e, prefix, log)
% Write the converted value under its current name; an already-set different value is a conflict
    new = label(prefix, new_parts);
    if has_path(s, new_parts)
        cur = get_path(s, new_parts);
        if e.yield
            log.done{end+1} = sprintf('%s ignored: %s is set', old, new);
            return
        end
        if isstruct(cur) && isstruct(val)
            [merged, clash] = merge_checked(cur, val, new);
            if ~isempty(clash)
                log.errors{end+1} = sprintf(['%s and %s are both set, with different values for %s; ' ...
                    'keep only %s'], old, new, strjoin(clash, ', '), new);
                return
            end
            s = set_path(s, new_parts, merged);
        elseif ~same_value(cur, val)
            log.errors{end+1} = sprintf('%s and %s are both set, to different values; keep only %s', ...
                old, new, new);
            return
        end
    else
        s = set_path(s, new_parts, val);
    end
    msg = sprintf('%s -> %s', old, new);
    if ~isempty(e.since); msg = sprintf('%s (since %s)', msg, e.since); end
    log.done{end+1} = msg;
end

function [a, clash] = merge_checked(a, b, path)
    clash = {};
    for f = fieldnames(b)'
        name = f{1};
        if ~isfield(a, name)
            a.(name) = b.(name);
        elseif isstruct(a.(name)) && isstruct(b.(name))
            [a.(name), c] = merge_checked(a.(name), b.(name), [path '.' name]);
            clash = [clash, c]; %#ok<AGROW>
        elseif ~same_value(a.(name), b.(name))
            clash{end+1} = [path '.' name]; %#ok<AGROW>
        end
    end
end

%% Converters (named by 'convert:' in the overlay) ----------------------------

function [s, log] = convert(name, s, parts, val, e, prefix, log)
    old = label(prefix, parts);
    switch name
        case 'output_location'
            % v0.5.0: an absolute path, or a folder name inside data_path
            v = char(string(val));
            if is_absolute(v)
                target = v;
            elseif has_path(s, {'path', 'anat'}) && ~isempty(get_path(s, {'path', 'anat'}))
                target = fullfile(char(string(get_path(s, {'path', 'anat'}))), v);
            else
                log.errors{end+1} = sprintf(['%s = ''%s'' is relative to data_path, which this config ' ...
                    'does not set; set path.sim to the full output path'], old, v);
                return
            end
            [s, log] = move(s, parts, {'path', 'sim'}, target, e, prefix, log);

        case 'pct_mapping'
            [s, log] = convert_pct_mapping(s, e, log);

        case 'calibration_weights'
            % v0.5.0: 0 = Gaussian weighting, any other value = uniform; now 0 = uniform
            if isnumeric(val) && isscalar(val) && val ~= 0
                [s, log] = move(s, parts, {'calibration', 'opt_weights'}, 0, e, prefix, log);
            else
                log.errors{end+1} = sprintf('%s = %s: %s', old, value_str(val), e.note);
            end

        case 'seg_labels'
            % labels are fixed to the SimNIBS charm labels; custom values cannot be kept
            ref = charm_seg_labels();
            differ = {};
            if isstruct(val)
                for f = fieldnames(val)'
                    if ~isfield(ref, f{1}) || ~isequal(double(val.(f{1})(:)'), double(ref.(f{1})(:)'))
                        differ{end+1} = f{1}; %#ok<AGROW>
                    end
                end
            end
            if isempty(differ)
                s = rm_path(s, parts);
                log.done{end+1} = sprintf('%s removed (since %s): %s', old, e.since, e.note);
            else
                log.errors{end+1} = sprintf('%s.{%s} differ from the SimNIBS charm labels: %s', ...
                    old, strjoin(differ, ', '), e.note);
            end

        otherwise
            error('prestus:legacyConfigKey', 'No converter ''%s'' for %s (overlay.yaml).', name, old);
    end
end

function [s, log] = convert_pct_mapping(s, e, log)
% v0.5.0 medium_setup.m: pseudoCT_variant names a preset of (density, sound
% speed, attenuation) mappings; pct_skullmapping was read only for a variant
% that names no preset. 'yakuub' is now called 'mueller' (same formula).
    presets = containers.Map( ...
        {'kosciessa', 'k-plan', 'carpino', 'yakuub', 'marquet'}, ...
        {{'k-plan', 'marsac', 'yakuub'}, {'k-plan', 'k-plan', 'k-plan'}, ...
         {'k-wave', 'k-plan', 'yakuub'}, {'marsac', 'marsac', 'yakuub'}, ...
         {'marquet', 'marquet', 'k-plan'}});
    props = {'density', 'soundspeed', 'attenuation'};
    has_v = isfield(s, 'pseudoCT_variant');
    has_m = isfield(s, 'pct_skullmapping') && isstruct(s.pct_skullmapping);
    if has_v
        v = char(string(s.pseudoCT_variant));
        if isKey(presets, v)
            m = presets(v);
            src = sprintf('pseudoCT_variant ''%s''', v);
            if has_m
                log.done{end+1} = sprintf(['pct_skullmapping ignored: v0.5.0 read it only when ' ...
                    'pseudoCT_variant named no preset (%s does)'], v);
            end
        elseif has_m && all(isfield(s.pct_skullmapping, props))
            m = cellfun(@(p) char(string(s.pct_skullmapping.(p))), props, 'UniformOutput', false);
            src = 'pct_skullmapping';
        else
            log.errors{end+1} = sprintf(['pseudoCT_variant = ''%s'' names no preset and ' ...
                'pct_skullmapping.{density, soundspeed, attenuation} is not set: %s'], v, e.note);
            return
        end
    else
        log.errors{end+1} = ['pct_skullmapping without pseudoCT_variant: v0.5.0 then used the ' ...
            '''kosciessa'' preset and ignored pct_skullmapping; ' e.note];
        return
    end
    m = strrep(m, 'yakuub', 'mueller');
    if any(strcmp(m, 'marquet'))
        log.errors{end+1} = sprintf(['%s uses the ''marquet'' mapping, which PRESTUS no longer has: ' ...
            '%s'], src, e.note);
        return
    end
    if has_v; s = rmfield(s, 'pseudoCT_variant'); end
    if has_m; s = rmfield(s, 'pct_skullmapping'); end
    for p = 1:3
        [s, log] = put(s, {'pct', ['mapping_' props{p}]}, m{p}, src, e, '', log);
    end
end

%% Overrides ------------------------------------------------------------------

function ov = upgrade_overrides(ov, reg, source)
    done = {};
    errors = {};
    for i = 1:size(ov, 1)
        key = char(string(ov{i, 1}));
        idx = regexp(key, '^transducer(\(\d+\))', 'tokens', 'once');
        base = regexprep(key, '\(\d+\)', '');
        for e = reg
            if ~(strcmp(base, e.key) || startsWith(base, [e.key '.']))
                continue
            end
            if strcmp(e.kind, 'rename')
                new = [e.use{1} base(numel(e.key)+1:end)];
                if ~isempty(idx)
                    new = regexprep(new, '^transducer', ['transducer' idx{1}]);
                end
                ov{i, 1} = new;
                done{end+1} = sprintf('%s -> %s', key, new); %#ok<AGROW>
            else
                hint = strjoin(e.use, ', ');
                if isempty(hint); hint = e.note; end
                errors{end+1} = sprintf('%s: set %s instead', key, hint); %#ok<AGROW>
            end
            break
        end
    end
    report(struct('done', {done}, 'errors', {errors}), source);
end

%% Messages -------------------------------------------------------------------

function report(log, source)
    if ~isempty(log.errors)
        error('prestus:legacyConfigKey', ...
            ['Config %s has keys from an earlier PRESTUS version that cannot be converted:\n  %s\n' ...
             'See doc_configuration.md, "Configs from earlier PRESTUS versions".'], ...
            source, strjoin(log.errors, sprintf('\n  ')));
    end
    if ~isempty(log.done)
        warning('prestus:renamedConfigKey', ...
            ['Config %s has keys from an earlier PRESTUS version; they were converted for this run. ' ...
             'Update the config:\n  %s'], source, strjoin(log.done, sprintf('\n  ')));
    end
end

function s = label(prefix, parts)
    s = strjoin(parts, '.');
    if ~isempty(prefix)
        s = [prefix '.' s];
    end
end

function parts = strip_tr(parts, prefix)
    if ~isempty(prefix) && strcmp(parts{1}, 'transducer')
        parts = parts(2:end);
    end
end

function txt = value_str(v)
    if ischar(v) || isstring(v)
        txt = sprintf('''%s''', char(v));
    elseif isnumeric(v) || islogical(v)
        txt = mat2str(v);
    else
        txt = class(v);
    end
end

function tf = same_value(a, b)
    if (ischar(a) || isstring(a)) && (ischar(b) || isstring(b))
        tf = strcmp(char(a), char(b));
    elseif (isnumeric(a) || islogical(a)) && (isnumeric(b) || islogical(b))
        tf = isequal(double(a), double(b));
    else
        tf = isequal(a, b);
    end
end

function tf = is_absolute(p)
    tf = startsWith(p, {'/', '\', '~'}) || ~isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'));
end

%% Nested-struct helpers --------------------------------------------------------

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
    if isfield(s, parts{1}) && isstruct(s.(parts{1})) && isscalar(s.(parts{1}))
        sub = s.(parts{1});
    elseif isfield(s, parts{1}) && ~isempty(s.(parts{1}))
        error('prestus:configKeyConflict', ...
            'Cannot convert into %s: %s is set and is not a section.', strjoin(parts, '.'), parts{1});
    else
        sub = struct();
    end
    s.(parts{1}) = set_path(sub, parts(2:end), v);
end

function s = rm_path(s, parts)
    if numel(parts) == 1
        s = rmfield(s, parts{1});
    else
        s.(parts{1}) = rm_path(s.(parts{1}), parts(2:end));
    end
end

function [c, was_struct] = to_cell(t)
    was_struct = isstruct(t);
    if was_struct
        c = num2cell(t);
    else
        c = t;
    end
end

function t = from_cell(c, was_struct)
    if ~was_struct
        t = c;
        return
    end
    names = {};
    for k = 1:numel(c)
        names = union(names, fieldnames(c{k}), 'stable');
    end
    for k = 1:numel(c)
        missing = setdiff(names, fieldnames(c{k}));
        for j = 1:numel(missing)
            c{k}.(missing{j}) = [];
        end
        c{k} = orderfields(c{k}, names);
    end
    t = reshape([c{:}], size(c));
end
