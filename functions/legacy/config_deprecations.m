function reg = config_deprecations(prestus_root)
% CONFIG_DEPRECATIONS  Config keys renamed or removed since PRESTUS v0.5.0
%
% Reads the 'deprecated:' section of the PRESTUS Nexus overlay
% (nexus/prestus_nexus/overlay.yaml), the single list of old config keys that
% MATLAB (UPGRADE_CONFIG_KEYS) and the Nexus both use. Each entry is one line,
%
%   <old.dotted.key>: {use: <new key or [keys]>, kind: <kind>, since: <version>, ...}
%
% with kind
%   rename  - same meaning and units under a new name: the value is moved
%   convert - renamed with an exact value conversion (convert: <converter>,
%             a local function of UPGRADE_CONFIG_KEYS)
%   removed - no replacement: the key is dropped; with keep: <value>, any other
%             value stops the run, as it asked for behaviour PRESTUS no longer has
%   error   - no exact conversion: the run stops and note says what to set
% Entries without a kind are handled by the function named in 'source'.
% yield: true lets the old key give way to an already-set new key (aliases).
%
% Use as:
%   reg = config_deprecations()
%   reg = config_deprecations(prestus_root)
%
% Output:
%   reg - struct array: key, use (cellstr), kind, since, convert, keep,
%         has_keep, yield, note, source ([] when the overlay is not found)
%
% See also: UPGRADE_CONFIG_KEYS, CONFIG_KNOWN_EXTRA_KEYS, LOAD_PARAMETERS

    persistent cached cached_file cached_date

    if nargin < 1 || isempty(prestus_root)
        prestus_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    end
    reg = struct('key', {}, 'use', {}, 'kind', {}, 'since', {}, 'convert', {}, ...
        'keep', {}, 'has_keep', {}, 'yield', {}, 'note', {}, 'source', {});

    overlay = fullfile(prestus_root, 'nexus', 'prestus_nexus', 'overlay.yaml');
    if ~isfile(overlay)
        warning('prestus:noDeprecationList', ...
            ['%s not found: config keys from earlier PRESTUS versions are not converted ' ...
             'and are reported as unknown.'], overlay);
        return
    end
    info = dir(overlay);
    if ~isempty(cached) && strcmp(cached_file, overlay) && cached_date == info.datenum
        reg = cached;
        return
    end
    if isempty(which('yaml.load'))
        addpath(genpath(fullfile(prestus_root, 'external', 'yaml')));
    end

    lines = splitlines(fileread(overlay));
    in_section = false;
    for i = 1:numel(lines)
        ln = lines{i};
        top = regexp(ln, '^([A-Za-z_]\w*):', 'tokens', 'once');
        if ~isempty(top)
            in_section = strcmp(top{1}, 'deprecated');
            continue
        end
        if ~in_section
            continue
        end
        tok = regexp(ln, '^  ([A-Za-z_][\w.]*):\s*(\{.*\})\s*$', 'tokens', 'once');
        if isempty(tok)
            continue
        end
        e = yaml.load(tok{2});
        r = struct();
        r.key      = tok{1};
        r.use      = cellfun(@char, cellify(get_or(e, 'use', {})), 'UniformOutput', false);
        r.kind     = char(get_or(e, 'kind', ''));
        r.since    = char(get_or(e, 'since', ''));
        r.convert  = char(get_or(e, 'convert', ''));
        r.has_keep = isfield(e, 'keep');
        r.keep     = get_or(e, 'keep', []);
        if isstring(r.keep); r.keep = char(r.keep); end
        r.yield    = logical(get_or(e, 'yield', false));
        r.note     = char(get_or(e, 'note', ''));
        r.source   = char(get_or(e, 'source', ''));
        reg(end+1) = r; %#ok<AGROW>
    end

    cached = reg;
    cached_file = overlay;
    cached_date = info.datenum;
end

function v = get_or(s, f, default)
    if isfield(s, f)
        v = s.(f);
    else
        v = default;
    end
end

function c = cellify(v)
    if iscell(v)
        c = v;
    elseif isempty(v)
        c = {};
    else
        c = {v};
    end
end
