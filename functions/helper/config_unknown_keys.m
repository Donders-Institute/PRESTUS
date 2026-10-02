function [unknown, suggestions] = config_unknown_keys(user, reference, extra_keys)
% CONFIG_UNKNOWN_KEYS  Dotted keys of a user config that PRESTUS does not know
%
% Walks a user config struct against the merged default config. A key is known
% when the default config contains it, or when it is listed in extra_keys (keys
% the code reads that config_default.yaml omits, e.g. from the Nexus overlay).
% The walk descends only into groups that are structs in both the user config
% and the reference, so free-form values (maps such as neuronav.coil_map, or a
% struct given where the default holds a scalar) are not inspected. Struct
% arrays (e.g. transducer) are checked element by element against element 1 of
% the reference.
%
% Use as:
%   [unknown, suggestions] = config_unknown_keys(user, reference)
%   [unknown, suggestions] = config_unknown_keys(user, reference, extra_keys)
%
% Input:
%   user       - struct read from a study config YAML
%   reference  - struct of known keys (merged default config)
%   extra_keys - (optional) cellstr of additional known dotted keys; a key also
%                counts as known when it is a prefix of an extra key
%
% Output:
%   unknown     - cellstr of unknown dotted keys
%   suggestions - cellstr (same size); closest known sibling key, or ''
%
% See also: LOAD_PARAMETERS, APPLY_CONFIG_OVERRIDES

    if nargin < 3 || isempty(extra_keys)
        extra_keys = {};
    end
    extra_keys = cellstr(extra_keys);
    unknown = {};
    suggestions = {};
    if ~isstruct(user) || ~isstruct(reference)
        return
    end
    walk(user, reference, '');

    function walk(u, r, prefix)
        for ui = 1:numel(u)
            fn = fieldnames(u(ui));
            for k = 1:numel(fn)
                key = fn{k};
                path = join_key(prefix, key);
                if ~isempty(r) && isstruct(r) && isfield(r, key)
                    rv = r(1).(key);
                    uv = u(ui).(key);
                    if isstruct(uv) && isstruct(rv) && ~isempty(fieldnames(rv))
                        walk(uv, rv, path);
                    end
                elseif is_extra(path)
                    continue
                elseif any(strcmp(path, unknown))
                    continue
                else
                    unknown{end+1} = path; %#ok<AGROW>
                    siblings = {};
                    if isstruct(r) && ~isempty(r)
                        siblings = fieldnames(r);
                    end
                    suggestions{end+1} = closest(key, siblings); %#ok<AGROW>
                end
            end
        end
    end

    function tf = is_extra(path)
        tf = any(strcmp(path, extra_keys)) || ...
             any(startsWith(extra_keys, [path '.']));
    end
end

function s = join_key(prefix, key)
    if isempty(prefix)
        s = key;
    else
        s = [prefix '.' key];
    end
end

function best = closest(key, candidates)
% Closest candidate by edit distance, when it is plausibly a typo.
    best = '';
    best_d = inf;
    for i = 1:numel(candidates)
        d = edit_distance(lower(key), lower(candidates{i}));
        if d < best_d
            best_d = d;
            best = candidates{i};
        end
    end
    if best_d > max(2, floor(numel(key) / 3))
        best = '';
    end
end

function d = edit_distance(a, b)
    m = numel(a);
    n = numel(b);
    D = zeros(m + 1, n + 1);
    D(:, 1) = 0:m;
    D(1, :) = 0:n;
    for i = 1:m
        for j = 1:n
            D(i + 1, j + 1) = min([D(i, j + 1) + 1, D(i + 1, j) + 1, ...
                                   D(i, j) + (a(i) ~= b(j))]);
        end
    end
    d = D(m + 1, n + 1);
end
