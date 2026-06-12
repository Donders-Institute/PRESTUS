function [alpha_coeff, alpha_power] = medium_pct_attenuation(parameters, alpha_coeff, alpha_power, pseudoCT, skull_idx, algorithm, density)
% MEDIUM_PCT_ATTENUATION  Map pseudo-CT Hounsfield values to skull attenuation coefficients
%
% Applies one of seven algorithms to fill alpha_coeff and alpha_power for
% skull voxels indexed by skull_idx.
%   'k-plan'    fixed value from config
%   'mueller'   square-root porosity model (Aubry et al. 2022 / Fry 1978)
%   'aubry'     porosity model (Aubry et al. 2003)
%   'webb'      exponential HU model (Webb et al. 2020, scanner-specific; uses
%               parameters.pct.webb_att_alpha0, .webb_att_beta, .webb_att_c or
%               GE 120 kVp B-kernel defaults: Alpha_0=26, Beta=1.3, c=-0.0016)
%   'pichardo'  2D ρ–freq interpolation from MapPichardo.h5
%               (requires external/BabelBrain submodule and density argument)
%   'mcdannold' polynomial ρ model (McDannold et al., linear in frequency)
%   'none'      fixed scalar from parameters.medium_properties.skull
% Webb, Pichardo, and McDannold results are converted to k-Wave alpha_coeff
% (dB/cm/MHz^alpha_power) by normalising out the frequency dependence at 1 MHz.
%
% Use as:
%   [alpha_coeff, alpha_power] = medium_pct_attenuation(parameters, alpha_coeff, alpha_power, pseudoCT, skull_idx, algorithm)
%   [alpha_coeff, alpha_power] = medium_pct_attenuation(..., density)  % required for pichardo/mcdannold
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
    algorithm   (1,:) char {mustBeMember(algorithm, {'k-plan','mueller','aubry','webb','pichardo','mcdannold','none'})}
    density     {mustBeNumeric} = []
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

    case 'webb'
        % Exponential HU-based attenuation from Webb et al. 2020 (Table IV).
        % Formula: Att[Np/m] = (Alpha_0 * f_MHz^Beta * exp(HU*c)) * 100
        % Converted to k-Wave format: alpha_coeff [dB/cm/MHz^Beta]
        %   alpha_coeff(x) = Alpha_0 * exp(HU(x) * c) * 100 * 0.08686
        %   alpha_power    = Beta
        % Default parameters (GE 120 kVp, B kernel, 0.5×0.6 mm):
        %   Alpha_0 = 26, Beta = 1.3, c = -0.0016
        % Override via parameters.pct.webb_att_alpha0, .webb_att_beta, .webb_att_c.
        % Reference: Webb et al. IEEE Trans Ultrason Ferroelectr Freq Control
        %   68(5):1532–1545, 2020. DOI: 10.1109/TUFFC.2020.3039743
        Alpha_0 = 26;
        Beta    = 1.3;
        c       = -0.0016;
        if isfield(parameters, 'pct')
            if isfield(parameters.pct, 'webb_att_alpha0'), Alpha_0 = parameters.pct.webb_att_alpha0; end
            if isfield(parameters.pct, 'webb_att_beta'),   Beta    = parameters.pct.webb_att_beta;   end
            if isfield(parameters.pct, 'webb_att_c'),      c       = parameters.pct.webb_att_c;      end
        end
        % 100 Np/m per (cm/m) × 0.08686 dB/Np = 8.686 dB/cm per unit of the formula
        alpha_coeff(skull_idx) = Alpha_0 .* exp(pseudoCT(skull_idx) * c) * 8.686;
        alpha_power(skull_idx) = Beta;

    case 'pichardo'
        % Frequency-dependent ρ→Att from Pichardo et al. 2D interpolation map.
        % Requires external/BabelBrain submodule, density argument, and parameters.transducer(1).freq_hz.
        % Map returns Att in Np/m; converted to k-Wave alpha_coeff normalised at 1 MHz
        % assuming linear frequency scaling (alpha_power = 1).
        if isempty(density)
            error(['medium_pct_attenuation: pichardo requires density argument. ' ...
                   'Call with 7 arguments: ..., algorithm, density']);
        end
        freq_hz  = parameters.transducer(1).freq_hz;
        freq_MHz = freq_hz / 1e6;
        [~, F_att] = medium_pct_pichardo_load(parameters);
        att_Npm = F_att(density(skull_idx), freq_MHz);
        att_dBcm = att_Npm * 0.08686;
        alpha_coeff(skull_idx) = att_dBcm ./ freq_MHz;  % normalise to 1 MHz, alpha_power = 1
        alpha_power(skull_idx) = 1;

    case 'mcdannold'
        % Polynomial ρ→Att from McDannold et al. with linear frequency scaling.
        % Att[Np/m] = polyval(p, ρ) * freq_Hz / 660e3
        % Converted to k-Wave alpha_coeff [dB/cm/MHz] (alpha_power = 1).
        if isempty(density)
            error(['medium_pct_attenuation: mcdannold requires density argument. ' ...
                   'Call with 7 arguments: ..., algorithm, density']);
        end
        p = [1.36e-10, -1.41e-6, 5.40e-3, -9.02, 5.71e3];
        % Att[Np/m] at f_MHz: polyval(p,ρ) * f_MHz / 0.66
        % alpha_coeff = Att_at_1MHz_Npm * 0.08686 [dB/cm/MHz]
        att_1MHz_Npm = polyval(p, density(skull_idx)) / 0.66;
        alpha_coeff(skull_idx) = att_1MHz_Npm * 0.08686;
        alpha_power(skull_idx) = 1;

    case 'none'

        alpha_coeff(skull_idx) = parameters.medium_properties.skull.alpha_coeff;
        alpha_power(skull_idx) = parameters.medium_properties.skull.alpha_power;

    otherwise
        error("Specified pCT attenuation mapping is not supported.")
end