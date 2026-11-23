function [sound_speed] = medium_pct_soundspeed(parameters, sound_speed, density, pseudoCT, skull_idx, algorithm)
% MEDIUM_PCT_SOUNDSPEED  Map pseudo-CT Hounsfield values to skull bone sound speed
%
% Applies one of seven algorithms to fill sound_speed for skull voxels
% indexed by skull_idx.
%   'k-plan'    linear ρ→c from k-Plan pipeline
%   'marsac'    linear ρ mixture model (Marsac et al. 2017)
%   'aubry'     porosity mixture model (Aubry et al. 2003)
%   'webb'      linear HU→c from Webb et al. 2018 (scanner-specific; requires
%               parameters.pct.webb_sos_slope and .webb_sos_intercept or uses
%               GE 120 kVp B-kernel defaults: slope=0.75, intercept=1.32)
%   'pichardo'  2D ρ–freq interpolation from MapPichardo.h5
%               (Pichardo et al.); requires external/BabelBrain submodule and
%               simulation frequency in parameters.transducer(1).freq_hz
%   'mcdannold' polynomial ρ→c from McDannold et al.
%   'none'      fixed scalar from parameters.medium_properties.skull
% Limits (by default at least water sound speed) are applied afterwards by
% PCT_MAP_PROPERTIES (pct.limits.sound_speed_*).
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
%   algorithm   - one of the seven algorithms above
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
    skull_idx   {mustBeNumericOrLogical}
    algorithm   (1,:) char {mustBeMember(algorithm, {'k-plan','marsac','aubry','webb','pichardo','mcdannold','none'})}
end

switch algorithm
    case 'k-plan'

        sound_speed(skull_idx) = 1.33.*density(skull_idx) + 167;


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

        phi(skull_idx) = medium_pct_porosity(parameters, pseudoCT(skull_idx));
        sound_speed(skull_idx) = c_water * phi(skull_idx) + ...
            c_bone * (1-phi(skull_idx));

    case 'webb'
        % Linear HU→SoS from Webb et al. 2018 (Tables I/II).
        % Formula: SoS = slope * HU + intercept * 1000  [m/s]
        % Default parameters (GE 120 kVp, B kernel, 0.5×0.6 mm):
        %   slope = 0.75, intercept = 1.32
        % Override via parameters.pct.webb_sos_slope and .webb_sos_intercept.
        % Reference: Webb et al. IEEE Trans Ultrason Ferroelectr Freq Control
        %   65(7):1111–1124, 2018. DOI: 10.1109/TUFFC.2018.2827899
        slope     = 0.75;
        intercept = 1.32;
        if isfield(parameters, 'pct')
            if isfield(parameters.pct, 'webb_sos_slope') && ~isempty(parameters.pct.webb_sos_slope),     slope     = parameters.pct.webb_sos_slope;     end
            if isfield(parameters.pct, 'webb_sos_intercept') && ~isempty(parameters.pct.webb_sos_intercept), intercept = parameters.pct.webb_sos_intercept; end
        end
        sound_speed(skull_idx) = slope * pseudoCT(skull_idx) + intercept * 1000;

    case 'pichardo'
        % Frequency-dependent ρ→SoS from Pichardo et al. 2D interpolation map.
        % Requires external/BabelBrain submodule and parameters.transducer(1).freq_hz.
        freq_hz = parameters.transducer(1).freq_hz;
        [F_sos, ~] = medium_pct_pichardo_load(parameters);
        sound_speed(skull_idx) = F_sos(density(skull_idx), repmat(freq_hz / 1e6, size(density(skull_idx))));

    case 'mcdannold'
        % Polynomial ρ→SoS from McDannold et al.
        % Formula: SoS = 1 / polyval(p, ρ)
        % Reference: McDannold et al.
        p = [-2.79e-18, 5.31e-16, 1.69e-10, -7.63e-7, 1.24e-3];
        sound_speed(skull_idx) = 1 ./ polyval(p, density(skull_idx));

    case 'none'

        sound_speed(skull_idx) = parameters.medium_properties.skull.sound_speed;

    otherwise
        error("Specified CT sound speed mapping is not supported.")
end
