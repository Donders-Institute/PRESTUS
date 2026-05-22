function [alpha_coeff, alpha_power] = medium_pct_attenuation(parameters, alpha_coeff, alpha_power, pseudoCT, skull_idx, algorithm)
% MEDIUM_PCT_ATTENUATION  Map pseudo-CT Hounsfield values to skull attenuation coefficients
%
% Applies one of four algorithms to fill alpha_coeff and alpha_power for
% skull voxels indexed by skull_idx. Algorithm 'k-plan' uses a fixed
% value; 'mueller' uses a square-root porosity model from Aubry et al.
% 2022 / Fry 1978; 'aubry' uses the porosity model from Aubry et al.
% 2003; 'none' uses the fixed scalar from parameters.medium_properties.skull.
%
% Use as:
%   [alpha_coeff, alpha_power] = medium_pct_attenuation(parameters, alpha_coeff, alpha_power, pseudoCT, skull_idx, algorithm)
%
% Input:
%   parameters  - PRESTUS config; must contain
%                 medium_properties.skull.alpha_coeff [dB/(MHz cm)] and alpha_power
%   alpha_coeff - full-grid attenuation coefficient array to update
%   alpha_power - full-grid attenuation power array to update
%   pseudoCT    - pseudo-CT Hounsfield values (full grid)
%   skull_idx   - linear indices of skull voxels into the grid
%   algorithm   - one of 'k-plan', 'mueller', 'aubry', 'none'
%
% Output:
%   alpha_coeff - updated alpha_coeff with skull voxels filled
%   alpha_power - updated alpha_power with skull voxels filled
%
% See also: MEDIUM_SETUP, MEDIUM_PCT_DENSITY, MEDIUM_PCT_SOUNDSPEED

arguments
    parameters  (1,1) struct
    alpha_coeff {mustBeNumeric}
    alpha_power {mustBeNumeric}
    pseudoCT    {mustBeNumeric}
    skull_idx   {mustBeNumericOrLogical}
    algorithm   (1,:) char {mustBeMember(algorithm, {'k-plan','mueller','aubry','none'})}
end

switch algorithm
    case 'k-plan'

        kPlan_alpha = 13.3; % https://dispatch.k-plan.io/static/docs/simulation-pipeline.html
        kPlan_alpha_power = 1;
        % Note that we allow different values to be specified in the config.
        % If replication of k-Wave is the goal, the above values should be specified.
        % Throw a warning in the case of deviations.
        if parameters.medium_properties.skull.alpha_coeff ~= kPlan_alpha || ...
            parameters.medium_properties.skull.alpha_power ~= kPlan_alpha_power
            warn('Specified attenuation varies from k-Plan setup.')
        end
        alpha_coeff(skull_idx) = parameters.medium_properties.skull.alpha_coeff;
        alpha_power(skull_idx) = parameters.medium_properties.skull.alpha_power;

    case 'mueller'

        alpha_min     = 4;        % cortical bone at 500 kHz [dB/cm] [Aubry et al., 2022] 
        alpha_max     = 8.7;      % bone at 500 kHz [dB/cm] [Fry 1978]

        % Finds maximum and minimum values
        HU_min = min(pseudoCT(skull_idx));
        HU_max = max(pseudoCT(skull_idx));

        % estimate attenuation based on (pseudo-)HU
        alpha_pseudoCT(skull_idx) = alpha_min + (alpha_max - alpha_min) * ...
            (1 - (pseudoCT(skull_idx) - HU_min) / (HU_max - HU_min)).^0.5;
        alpha_power(skull_idx) = parameters.medium_properties.skull.alpha_power;
        % convert alpha at 500 kHz into prefactor alpha0 (dB/MHz/cm) according to specified alpha_power
        % (definition of lower and upper attenuation bounds is derived from 500kHz)
        alpha_coeff(skull_idx) = alpha_pseudoCT(skull_idx)./(0.5^parameters.medium_properties.skull.alpha_power);

    case 'aubry'

        alpha_min     = 0.2;    % Aubry et al., 2003
        alpha_max     = 8;      % Aubry et al., 2003

        phi(skull_idx) = 1-(pseudoCT(skull_idx)/max(pseudoCT(skull_idx)));
        alpha_coeff(skull_idx) = alpha_min + (alpha_max - alpha_min) * (phi(skull_idx).^0.5);
        % regularize sound speed to a minimum of water
        alpha_power(skull_idx) = parameters.medium_properties.skull.alpha_power;

    case 'none'

        alpha_coeff(skull_idx) = parameters.medium_properties.skull.alpha_coeff;
        alpha_power(skull_idx) = parameters.medium_properties.skull.alpha_power;

    otherwise
        error("Specified pCT attenuation mapping is not supported.")
end