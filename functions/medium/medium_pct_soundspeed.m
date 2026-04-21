function [sound_speed] = medium_pct_soundspeed(parameters, sound_speed, density, pseudoCT, skull_idx, algorithm)
% MEDIUM_PCT_SOUNDSPEED  Map pseudo-CT Hounsfield values to skull bone sound speed
%
% Applies one of four algorithms to fill sound_speed for skull voxels
% indexed by skull_idx. Algorithm 'k-plan' uses the linear
% density-to-speed relationship from the k-Plan pipeline; 'marsac' uses
% the linear mixture model from Marsac et al. 2017 (requires density to be
% already filled); 'aubry' uses the porosity mixture model from Aubry et
% al. 2003; 'none' uses the fixed scalar from
% parameters.medium_properties.skull. All algorithms clamp to at least
% water sound speed.
%
% Use as:
%   [sound_speed] = medium_pct_soundspeed(parameters, sound_speed, density, pseudoCT, skull_idx, algorithm)
%
% Input:
%   parameters  - PRESTUS config; must contain medium_properties.water.sound_speed [m/s],
%                 medium_properties.skull.sound_speed [m/s], and
%                 medium_properties.water.density [kg/m^3]
%   sound_speed - full-grid sound speed array to update [m/s]
%   density     - full-grid density array (already filled) [kg/m^3]
%   pseudoCT    - pseudo-CT Hounsfield values (full grid)
%   skull_idx   - linear indices of skull voxels into the grid
%   algorithm   - one of 'k-plan', 'marsac', 'aubry', 'none'
%
% Output:
%   sound_speed - updated sound speed with skull voxels filled [m/s]
%
% See also: MEDIUM_SETUP, MEDIUM_PCT_DENSITY, MEDIUM_PCT_ATTENUATION

arguments
    parameters  (1,1) struct
    sound_speed {mustBeNumeric}
    density     {mustBeNumeric}
    pseudoCT    {mustBeNumeric}
    skull_idx   {mustBeNumeric}
    algorithm   (1,:) char {mustBeMember(algorithm, {'k-plan','marsac','aubry','none'})}
end

switch algorithm
    case 'k-plan'

        sound_speed(skull_idx) = 1.33.*density(skull_idx) + 167;

        % regularize minimum to sound speed in water
        sound_speed(skull_idx) = max(parameters.medium_properties.water.sound_speed, sound_speed(skull_idx));

    case 'marsac'

        c_water       = parameters.medium_properties.water.sound_speed;     % sound speed [m/s]
        c_skull       = 3360;     % max. speed of sound in skull [m/s]
        rho_water     = parameters.medium_properties.water.density;      % density [kg/m^3]
        rho_bone      = 2100;     % max. skull density [kg/m3]

        sound_speed(skull_idx) = c_water + (c_skull - c_water) * ...
            (density(skull_idx) - rho_water) / (rho_bone - rho_water);
        
    case 'aubry'

        c_water = parameters.medium_properties.water.sound_speed;
        c_bone = parameters.medium_properties.skull.sound_speed;

        phi(skull_idx) = 1-(pseudoCT(skull_idx)/max(pseudoCT(skull_idx)));
        sound_speed(skull_idx) = c_water * phi(skull_idx) + ...
            c_bone * (1-phi(skull_idx));
        % regularize sound speed to a minimum of water
        sound_speed(skull_idx) = max(sound_speed(skull_idx),c_water);

    case 'none'

        sound_speed(skull_idx) = parameters.medium_properties.skull.sound_speed;

    otherwise
        error("Specified CT sound speed mapping is not supported.")
end
