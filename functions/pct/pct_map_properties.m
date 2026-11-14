function [props, mapping] = pct_map_properties(parameters, hu, skull_idx)
% PCT_MAP_PROPERTIES  Skull density, sound speed and attenuation from HU, as MEDIUM_SETUP maps them
%
% Applies pct.mapping_density / _soundspeed / _attenuation (MEDIUM_PCT_DENSITY,
% MEDIUM_PCT_SOUNDSPEED, MEDIUM_PCT_ATTENUATION) to the voxels skull_idx of
% hu, in the order MEDIUM_SETUP uses (sound speed and attenuation may read
% the mapped density). Voxels outside skull_idx keep the homogeneous skull
% values of medium_properties.skull. Shared by PCT_PREDICT_SKULL_PROPERTIES
% and PCT_MAPPING_CURVES.
%
% Use as:
%   [props, mapping] = pct_map_properties(parameters, hu, skull_idx)
%
% Input:
%   parameters - PRESTUS parameters (pct.mapping_*, medium_properties.skull/water)
%   hu         - Hounsfield values (any size)
%   skull_idx  - linear indices of the skull voxels into hu
%
% Output:
%   props   - struct: density [kg/m^3], sound_speed [m/s], alpha_coeff
%             [dB/(MHz^y cm)], alpha_power, each the size of hu
%   mapping - struct: density, sound_speed, attenuation (algorithm names)
%
% See also: MEDIUM_SETUP, PCT_PREDICT_SKULL_PROPERTIES, PCT_MAPPING_CURVES

arguments
    parameters (1,1) struct
    hu         {mustBeNumeric}
    skull_idx  {mustBeNumericOrLogical}
end

    mapping = struct('density',     mapping_name(parameters, 'mapping_density'), ...
                     'sound_speed', mapping_name(parameters, 'mapping_soundspeed'), ...
                     'attenuation', mapping_name(parameters, 'mapping_attenuation'));
    p = parameters;
    if ~isfield(p, 'simulation') || ~isfield(p.simulation, 'debug'); p.simulation.debug = 0; end
    skull = p.medium_properties.skull;
    density     = repmat(skull.density,     size(hu));
    sound_speed = repmat(skull.sound_speed, size(hu));
    alpha_coeff = repmat(skull.alpha_coeff, size(hu));
    alpha_power = repmat(skull.alpha_power, size(hu));
    density = medium_pct_density(p, density, hu, skull_idx, mapping.density);
    sound_speed = medium_pct_soundspeed(p, sound_speed, density, hu, skull_idx, mapping.sound_speed);
    [alpha_coeff, alpha_power] = medium_pct_attenuation(p, alpha_coeff, alpha_power, hu, skull_idx, ...
                                                        mapping.attenuation, density);
    props = struct('density', density, 'sound_speed', sound_speed, ...
                   'alpha_coeff', alpha_coeff, 'alpha_power', alpha_power);
end

function name = mapping_name(parameters, key)
    name = 'none';
    if isfield(parameters, 'pct') && isfield(parameters.pct, key) && ~isempty(parameters.pct.(key))
        name = char(parameters.pct.(key));
    end
end
