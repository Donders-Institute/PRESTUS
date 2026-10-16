function [model, forced] = thermal_absorption_model(parameters)
% THERMAL_ABSORPTION_MODEL  Resolve how acoustic heat deposition is specified
%
% Returns the absorption model the thermal simulation uses, from
% thermal.absorption_model:
%   'fraction'    - absorption = medium_properties.<tissue>.absorption_fraction
%                   x attenuation (alpha_coeff, alpha_power)
%   'coefficient' - absorption = medium_properties.<tissue>.absorption_coeff
%                   [dB/(cm MHz^y)], with y = the tissue's alpha_power
% Only one model applies to a run. A pseudo-CT skull with mapped attenuation
% (pct.enabled with pct.mapping_attenuation other than 'none') has no
% matching absorption map, so 'coefficient' falls back to 'fraction' there
% (forced = true; MEDIUM_SETUP warns).
%
% Use as:
%   [model, forced] = thermal_absorption_model(parameters)
%
% Input:
%   parameters - PRESTUS config struct
%
% Output:
%   model  - 'fraction' | 'coefficient'
%   forced - true when 'coefficient' was requested but 'fraction' is used
%
% See also: MEDIUM_SETUP, THERMAL_SIMULATION

model = 'fraction';
if isfield(parameters, 'thermal') && isfield(parameters.thermal, 'absorption_model') && ...
        ~isempty(parameters.thermal.absorption_model)
    model = lower(char(parameters.thermal.absorption_model));
end
if ~ismember(model, {'fraction', 'coefficient'})
    error('thermal_absorption_model:unknownModel', ...
        'thermal.absorption_model must be ''fraction'' or ''coefficient'' (got ''%s'').', model);
end

forced = false;
if strcmp(model, 'coefficient') && isfield(parameters, 'pct') && ...
        isfield(parameters.pct, 'enabled') && parameters.pct.enabled == 1
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
