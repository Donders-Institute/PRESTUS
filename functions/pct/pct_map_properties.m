function [props, mapping, info] = pct_map_properties(parameters, hu, skull_idx, hu_ref)
% PCT_MAP_PROPERTIES  Skull density, sound speed and attenuation from HU
%
% The one place where skull properties are mapped from a (pseudo-)CT, used by
% MEDIUM_SETUP, PCT_PREDICT_SKULL_PROPERTIES and PCT_MAPPING_CURVES:
%   1. HU limits (pct.limits.hu_*, PCT_APPLY_LIMITS)
%   2. density (pct.mapping_density, MEDIUM_PCT_DENSITY), then its limits
%   3. sound speed (pct.mapping_soundspeed, MEDIUM_PCT_SOUNDSPEED; may read
%      the density), then its limits
%   4. attenuation (pct.mapping_attenuation, MEDIUM_PCT_ATTENUATION), then the
%      limits of alpha_coeff
% Max-normalised mappings scale by hu_ref (PCT_HU_REFERENCE); without it, by
% the largest of the given skull HU. Voxels outside skull_idx keep the
% homogeneous skull values of medium_properties.skull.
%
% Use as:
%   [props, mapping, info] = pct_map_properties(parameters, hu, skull_idx)
%   [props, mapping, info] = pct_map_properties(parameters, hu, skull_idx, hu_ref)
%
% Input:
%   parameters - PRESTUS parameters (pct.mapping_*, pct.limits, medium_properties.skull/water)
%   hu         - Hounsfield values (any size)
%   skull_idx  - linear indices of the skull voxels into hu
%   hu_ref     - densest-bone HU for max-normalised mappings ([] = largest given skull HU)
%
% Output:
%   props   - struct: density [kg/m^3], sound_speed [m/s], alpha_coeff
%             [dB/(MHz^y cm)], alpha_power, each the size of hu; hu (after the limits)
%   mapping - struct: density, sound_speed, attenuation (algorithm names)
%   info    - struct: hu_reference, and per quantity (hu, density,
%             sound_speed, alpha_coeff) the number of skull voxels raised to
%             the minimum (n_low) and lowered to the maximum (n_high)
%
% See also: MEDIUM_SETUP, PCT_APPLY_LIMITS, PCT_HU_REFERENCE, PCT_PREDICT_SKULL_PROPERTIES

arguments
    parameters (1,1) struct
    hu         {mustBeNumeric}
    skull_idx  {mustBeNumericOrLogical}
    hu_ref     = []
end

    mapping = struct('density',     mapping_name(parameters, 'mapping_density'), ...
                     'sound_speed', mapping_name(parameters, 'mapping_soundspeed'), ...
                     'attenuation', mapping_name(parameters, 'mapping_attenuation'));
    p = parameters;
    if ~isfield(p, 'simulation') || ~isfield(p.simulation, 'debug'); p.simulation.debug = 0; end
    p.pct.hu_reference_value = hu_ref;      % read by MEDIUM_PCT_HU_MAX for this mapping only
    info = struct('hu_reference', hu_ref);

    [hu(skull_idx), info.hu.n_low, info.hu.n_high] = pct_apply_limits(p, 'hu', hu(skull_idx));

    skull = p.medium_properties.skull;
    density     = repmat(skull.density,     size(hu));
    sound_speed = repmat(skull.sound_speed, size(hu));
    alpha_coeff = repmat(skull.alpha_coeff, size(hu));
    alpha_power = repmat(skull.alpha_power, size(hu));

    density = medium_pct_density(p, density, hu, skull_idx, mapping.density);
    [density(skull_idx), info.density.n_low, info.density.n_high] = ...
        pct_apply_limits(p, 'density', density(skull_idx));
    sound_speed = medium_pct_soundspeed(p, sound_speed, density, hu, skull_idx, mapping.sound_speed);
    [sound_speed(skull_idx), info.sound_speed.n_low, info.sound_speed.n_high] = ...
        pct_apply_limits(p, 'sound_speed', sound_speed(skull_idx));
    [alpha_coeff, alpha_power] = medium_pct_attenuation(p, alpha_coeff, alpha_power, hu, skull_idx, ...
                                                        mapping.attenuation, density);
    [alpha_coeff(skull_idx), info.alpha_coeff.n_low, info.alpha_coeff.n_high] = ...
        pct_apply_limits(p, 'alpha_coeff', alpha_coeff(skull_idx));

    props = struct('density', density, 'sound_speed', sound_speed, ...
                   'alpha_coeff', alpha_coeff, 'alpha_power', alpha_power, 'hu', hu);
end

function name = mapping_name(parameters, key)
    name = 'none';
    if isfield(parameters, 'pct') && isfield(parameters.pct, key) && ~isempty(parameters.pct.(key))
        name = char(parameters.pct.(key));
    end
end
