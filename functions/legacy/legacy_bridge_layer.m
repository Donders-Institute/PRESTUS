function cfg = legacy_bridge_layer(cfg, source)
% LEGACY_BRIDGE_LAYER  Drop old config keys that the same layer makes redundant
%
% Some old keys cannot be converted exactly (UPGRADE_CONFIG_KEYS stops the run
% for them). One of them is harmless when the layer already holds what replaces
% it: the top-level expected_focal_distance_mm was, in v0.6, only a fallback for
% transducers without focal_distance_ep. When every transducer of the same layer
% sets focal_distance_ep or focal_distance_bowl, the key is dropped with a
% warning (prestus:legacyRedundantKey); otherwise it is left for
% UPGRADE_CONFIG_KEYS, which stops with instructions (v0.5.0 read this distance
% from the bowl apex in some code paths and from the exit plane in others).
%
% Use as:
%   cfg = legacy_bridge_layer(cfg, source)
%
% Input:
%   cfg    - one config layer (struct as read from YAML) or a parameters struct
%   source - name of the layer for messages
%
% Output:
%   cfg    - the layer without the redundant keys
%
% See also: UPGRADE_CONFIG_KEYS, LEGACY_BRIDGE_STRUCT

    if nargin < 2 || isempty(source)
        source = 'config';
    end
    if ~isstruct(cfg) || ~isscalar(cfg) || ~isfield(cfg, 'expected_focal_distance_mm')
        return
    end
    val = cfg.expected_focal_distance_mm;
    if isempty(val)
        cfg = rmfield(cfg, 'expected_focal_distance_mm');
        return
    end
    if ~isfield(cfg, 'transducer') || isempty(cfg.transducer)
        return
    end
    trs = cfg.transducer;
    if isstruct(trs); trs = num2cell(trs); end
    all_set = true;
    for k = 1:numel(trs)
        t = trs{k};
        all_set = all_set && isstruct(t) && ...
            (nonempty_field(t, 'focal_distance_ep') || nonempty_field(t, 'focal_distance_bowl'));
    end
    if ~all_set
        return
    end
    cfg = rmfield(cfg, 'expected_focal_distance_mm');
    warning('prestus:legacyRedundantKey', ...
        ['%s: top-level expected_focal_distance_mm (%g) is ignored: every transducer sets ' ...
         'focal_distance_ep or focal_distance_bowl, which took precedence in v0.6 as well. ' ...
         'Remove the key.'], char(string(source)), double(val(1)));
end

function tf = nonempty_field(s, f)
    tf = isfield(s, f) && ~isempty(s.(f));
end
