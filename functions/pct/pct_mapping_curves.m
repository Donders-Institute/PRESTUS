function curves = pct_mapping_curves(parameters, hu)
% PCT_MAPPING_CURVES  Skull properties as functions of HU, for the configured mappings
%
% Evaluates PCT_MAP_PROPERTIES on a range of Hounsfield values, every value
% treated as a skull voxel, so the Nexus can plot the HU -> density, sound
% speed and attenuation curves the run will apply (the same MATLAB code as
% MEDIUM_SETUP, no re-implementation).
%
% Use as:
%   curves = pct_mapping_curves(parameters)
%   curves = pct_mapping_curves(parameters, hu)
%
% Input:
%   parameters - PRESTUS parameters (pct.mapping_*, pct.density_calibration,
%                medium_properties, transducer.freq_hz)
%   hu         - HU values (default: -200 to 3000 HU in 10 HU steps)
%
% Output:
%   curves - struct: hu, density [kg/m^3], sound_speed [m/s], alpha_coeff
%            [dB/(MHz^y cm)], alpha_power (row vectors), mapping (algorithm
%            names, plus density_calibration for 'k-plan'), reference.skull
%            (the homogeneous skull values)
%
% See also: PCT_MAP_PROPERTIES, NEXUS_PCT

arguments
    parameters (1,1) struct
    hu (1,:) double = -200:10:3000
end

    [props, mapping] = pct_map_properties(parameters, hu, 1:numel(hu));
    if strcmp(mapping.density, 'k-plan')
        [~, ~, mapping.density_calibration] = medium_pct_density_calibration(parameters);
    end
    skull = parameters.medium_properties.skull;
    curves = struct('hu', hu, ...
        'density', double(props.density(:)'), 'sound_speed', double(props.sound_speed(:)'), ...
        'alpha_coeff', double(props.alpha_coeff(:)'), 'alpha_power', double(props.alpha_power(:)'), ...
        'mapping', mapping, ...
        'reference', struct('skull', struct('density', skull.density, 'sound_speed', skull.sound_speed, ...
                                            'alpha_coeff', skull.alpha_coeff)));
end
