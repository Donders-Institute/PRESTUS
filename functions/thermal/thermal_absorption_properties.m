function a = thermal_absorption_properties(parameters, tissue)
% THERMAL_ABSORPTION_PROPERTIES  Absorption model, coefficient and fraction of one tissue
%
% Resolves the tissue's absorption model (THERMAL_ABSORPTION_MODEL) and
% returns the specified quantity together with the one it implies at the
% tissue's configured attenuation (absorption shares the power law
% y = alpha_power):
%   absorption_fraction set - absorption_coeff    = absorption_fraction x alpha_coeff
%   absorption_coeff set    - absorption_fraction = absorption_coeff / alpha_coeff
% Neither set: fraction 1, as in MEDIUM_SETUP. A pseudo-CT skull with mapped
% attenuation applies its coefficient as the fraction absorption_coeff /
% alpha_coeff of the mapped attenuation (forced = true); its values here use
% the configured skull alpha_coeff (attenuation_mapped = true), and
% MEDIUM_SETUP logs the mapped range.
%
% Use as:
%   a = thermal_absorption_properties(parameters, tissue)
%
% Input:
%   parameters - PRESTUS config struct with medium_properties.<tissue>
%   tissue     - medium_properties label
%
% Output:
%   a - struct with fields:
%         model               - 'fraction' | 'coefficient'
%         forced              - true when a coefficient is applied as a fraction
%         absorption_coeff    - effective coefficient [dB/(cm MHz^y)]
%         absorption_fraction - effective fraction of attenuation [-]
%         attenuation_mapped  - true for a skull with pseudo-CT attenuation
%
% See also: THERMAL_ABSORPTION_MODEL, MEDIUM_SETUP

[a.model, a.forced] = thermal_absorption_model(parameters, tissue);
m = parameters.medium_properties.(tissue);
alpha = NaN;
if isfield(m, 'alpha_coeff') && ~isempty(m.alpha_coeff)
    alpha = m.alpha_coeff;
end

if strcmp(a.model, 'coefficient') || a.forced
    a.absorption_coeff = m.absorption_coeff;
    if alpha > 0
        a.absorption_fraction = a.absorption_coeff / alpha;
    else
        a.absorption_fraction = NaN;   % no attenuation: the fraction is undefined
    end
else
    a.absorption_fraction = 1;
    if isfield(m, 'absorption_fraction') && ~isempty(m.absorption_fraction) && ~any(isnan(m.absorption_fraction(:)))
        a.absorption_fraction = m.absorption_fraction;
    end
    a.absorption_coeff = a.absorption_fraction * alpha;
end

a.attenuation_mapped = false;
if strcmp(tissue, 'skull') && pct_skull_mapped(parameters) && isfield(parameters.pct, 'mapping_attenuation') && ...
        ~isempty(parameters.pct.mapping_attenuation) && ~strcmp(char(parameters.pct.mapping_attenuation), 'none')
    a.attenuation_mapped = true;
end
end
