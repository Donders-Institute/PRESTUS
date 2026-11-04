function [model, forced] = thermal_absorption_model(parameters, tissue)
% THERMAL_ABSORPTION_MODEL  Resolve how acoustic heat deposition is specified for a tissue
%
% The absorption field that medium_properties.<tissue> sets selects the model:
%   'fraction'    - absorption = absorption_fraction x attenuation
%                   (alpha_coeff, alpha_power); also when neither is set
%                   (fraction 1)
%   'coefficient' - absorption = absorption_coeff [dB/(cm MHz^y)], with
%                   y = the tissue's alpha_power
% Setting both errors; empty and NaN count as unset. LOAD_PARAMETERS clears
% the inherited field when a config layer sets the other one
% (MEDIUM_ABSORPTION_LAYER); a script that sets one after loading must empty
% the other. A pseudo-CT skull with mapped attenuation (pct.enabled with
% pct.mapping_attenuation other than 'none') has no matching absorption map,
% so a skull coefficient is applied as the fraction absorption_coeff /
% alpha_coeff of the mapped attenuation (forced = true; MEDIUM_SETUP warns).
%
% Use as:
%   [model, forced] = thermal_absorption_model(parameters, tissue)
%
% Input:
%   parameters - PRESTUS config struct
%   tissue     - medium_properties label
%
% Output:
%   model  - 'fraction' | 'coefficient'
%   forced - true when a coefficient is applied as a fraction (pseudo-CT skull)
%
% See also: THERMAL_ABSORPTION_PROPERTIES, MEDIUM_ABSORPTION_LAYER, MEDIUM_SETUP

model = 'fraction';
forced = false;
if ~isfield(parameters, 'medium_properties') || ~isfield(parameters.medium_properties, tissue) || ...
        ~isstruct(parameters.medium_properties.(tissue))
    return
end
m = parameters.medium_properties.(tissue);
has_coeff    = isfield(m, 'absorption_coeff')    && is_set(m.absorption_coeff);
has_fraction = isfield(m, 'absorption_fraction') && is_set(m.absorption_fraction);
if has_coeff && has_fraction
    error('thermal_absorption_model:bothSet', ...
        ['medium_properties.%s sets both absorption_fraction and absorption_coeff; set one ' ...
         'and leave the other empty (see doc_simulations-thermal.md#heat-deposition).'], tissue);
end
if ~has_coeff
    return
end
model = 'coefficient';

% medium_setup maps pseudo-CT attenuation into 'skull' only
if strcmp(tissue, 'skull') && pct_skull_mapped(parameters)
    mapping = 'none';
    if isfield(parameters.pct, 'mapping_attenuation') && ~isempty(parameters.pct.mapping_attenuation)
        mapping = char(parameters.pct.mapping_attenuation);
    end
    if ~strcmp(mapping, 'none')
        model  = 'fraction';
        forced = true;
    end
end
end

function tf = is_set(v)
tf = ~isempty(v) && ~isa(v, 'yaml.Null') && ~(isnumeric(v) && all(isnan(v(:))));
end
