function parameters = medium_absorption_layer(parameters, layer, source)
% MEDIUM_ABSORPTION_LAYER  Keep one absorption field per tissue when merging a config layer
%
% Each tissue specifies its absorption either as absorption_fraction (of the
% attenuation) or as absorption_coeff [dB/(cm MHz^y)], not both. A config
% layer (YAML file, struct or key/value overrides) that sets one of them for
% a tissue clears the other one inherited from earlier layers, so a study
% config sets absorption_coeff without also clearing the default
% absorption_fraction. A layer that sets both for one tissue errors. Empty
% and NaN values count as unset.
%
% Use as:
%   parameters = mergestruct(parameters, layer);
%   parameters = medium_absorption_layer(parameters, layer, source)
%
% Input:
%   parameters - config struct after merging the layer
%   layer      - the merged layer (struct; its medium_properties are read)
%   source     - name of the layer for the error message
%
% Output:
%   parameters - config struct with the other absorption field emptied
%
% See also: LOAD_PARAMETERS, THERMAL_ABSORPTION_MODEL

if ~isstruct(layer) || ~isfield(layer, 'medium_properties') || ~isstruct(layer.medium_properties)
    return
end
for t = fieldnames(layer.medium_properties)'
    tissue = t{1};
    props = layer.medium_properties.(tissue);
    if ~isstruct(props)
        continue
    end
    sets_fraction = isfield(props, 'absorption_fraction') && is_set(props.absorption_fraction);
    sets_coeff    = isfield(props, 'absorption_coeff')    && is_set(props.absorption_coeff);
    if sets_fraction && sets_coeff
        error('medium_absorption_layer:bothSet', ...
            ['%s sets both medium_properties.%s.absorption_fraction and absorption_coeff; ' ...
             'set one of them (see doc_simulations-thermal.md#heat-deposition).'], source, tissue);
    elseif sets_coeff
        parameters.medium_properties.(tissue).absorption_fraction = [];
    elseif sets_fraction
        parameters.medium_properties.(tissue).absorption_coeff = [];
    end
end
end

function tf = is_set(v)
tf = ~isempty(v) && ~isa(v, 'yaml.Null') && ~(isnumeric(v) && all(isnan(v(:))));
end
