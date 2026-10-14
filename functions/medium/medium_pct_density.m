function [density] = medium_pct_density(parameters, density, pseudoCT, skull_idx, algorithm)
% MEDIUM_PCT_DENSITY  Map pseudo-CT Hounsfield values to skull bone density
%
% Applies one of six algorithms to fill density for skull voxels indexed
% by skull_idx. Algorithm 'k-plan' uses a piece-wise linear HU-density
% table: k-Plan's default, or the scanner/protocol-specific table named by
% pct.density_calibration (MEDIUM_PCT_DENSITY_CALIBRATION); 'k-wave' uses hounsfield2density with
% a +1000 HU offset; 'marsac' uses the linear model from Marsac et al.
% 2017; 'aubry' uses the porosity mixture model from Aubry et al. 2003;
% 'air-tissue' uses the linear air-to-bone model (Pichardo and McDannold
% workflows); 'none' uses the
% fixed scalar from parameters.medium_properties.skull.
% All algorithms clamp density to at least water density.
%
% Use as:
%   [density] = medium_pct_density(parameters, density, pseudoCT, skull_idx, algorithm)
%
% Input:
%   parameters - PRESTUS config; must contain medium_properties.water.density [kg/m^3],
%                medium_properties.skull.density [kg/m^3], and simulation.debug
%   density    - full-grid density array to update [kg/m^3]
%   pseudoCT   - pseudo-CT Hounsfield values (full grid)
%   skull_idx  - linear indices of skull voxels into the grid
%   algorithm  - one of 'k-plan', 'k-wave', 'marsac', 'aubry', 'air-tissue', 'none'
%
% Output:
%   density - updated density with skull voxels filled [kg/m^3]
%
% See also: MEDIUM_SETUP, MEDIUM_PCT_SOUNDSPEED, MEDIUM_PCT_ATTENUATION

arguments
    parameters (1,1) struct
    density    {mustBeNumeric}
    pseudoCT   {mustBeNumeric}
    skull_idx  {mustBeNumericOrLogical}
    algorithm  (1,:) char {mustBeMember(algorithm, {'k-plan','k-wave','marsac','aubry','air-tissue','none'})}
end

% The k-plan / k-wave fits draw a figure; only when it is saved (debug mode)
show_fit = parameters.simulation.debug == 1;

% A UTE->HU line returns HU on the scale of the CT it was fitted to; density needs
% that CT's calibration (doc: Pairing a UTE->HU mapping with its CT calibration)
coupled = struct('treeby', 'ucl_lowdose_80kvp');       % skull_mapping -> its density_calibration
if isfield(parameters, 'pct') && isfield(parameters.pct, 'skull_mapping') && isfield(coupled, char(parameters.pct.skull_mapping))
    need = coupled.(char(parameters.pct.skull_mapping));
    have = '';
    if isfield(parameters.pct, 'density_calibration'); have = char(parameters.pct.density_calibration); end
    if ~strcmp(algorithm, 'k-plan') || ~strcmp(have, need)
        warn('prestus:pct:calibrationMismatch', ['pct.skull_mapping ''%s'' gives HU on the scale of its own CT; use ' ...
            'pct.mapping_density ''k-plan'' with pct.density_calibration ''%s'' (now: ''%s'', ''%s'').'], ...
            char(parameters.pct.skull_mapping), need, algorithm, have);
    end
end

% pct.density_calibration only acts on the k-plan table
if ~strcmp(algorithm, 'k-plan') && isfield(parameters, 'pct') && isfield(parameters.pct, 'density_calibration') ...
        && ~isempty(parameters.pct.density_calibration)
    warn(sprintf('pct.density_calibration (%s) is ignored: it applies to pct.mapping_density ''k-plan'', not ''%s''.', ...
        char(parameters.pct.density_calibration), algorithm));
end

switch algorithm
    case 'k-plan'

        % piece-wise linear HU -> mass density [kg/m^3]: k-Plan's default table
        % (-990, 60, 1000, 1950 HU -> 1.2, 1060, 1530, 2150 kg/m^3; see
        % https://dispatch.k-plan.io/static/docs/planning-images.html#ct-calibration)
        % or the table selected by pct.density_calibration
        [hounsfieldUnits, massDensity] = medium_pct_density_calibration(parameters);

        density(skull_idx) = fit_pairwiselinear(pseudoCT(skull_idx), hounsfieldUnits, massDensity, show_fit);

        % plot the mapping
        if show_fit
            output_plot = fullfile(parameters.io.dir_debug_medium,...
                sprintf('pCT_hounsfield-density_kplan.png'));
            exportgraphics(gcf, output_plot, 'Resolution', 150);
            close(gcf);
        end

        % regularize minimum to density in water
        density(skull_idx) = max(parameters.medium_properties.water.density, density(skull_idx));

    case 'k-wave'

        offset_HU   = 1000;
        rho_max     = 2100;     % max. density in skull [kg/m3]

        % Preprocess pCT values
        % Offset CT values to use housfield2density
        pseudoCT(skull_idx) = pseudoCT(skull_idx) + offset_HU;

        % set minimum to air tissue (pHU=-1000, pHU_scaled = 0)
        pseudoCT(skull_idx) = max(pseudoCT(skull_idx),0);

        % estimate density
        density(skull_idx) = hounsfield2density(pseudoCT(skull_idx), show_fit);
        
        % plot the mapping
        if show_fit
            output_plot = fullfile(parameters.io.dir_debug_medium,...
                sprintf('pCT_hounsfield-density_kwave.png'));
            exportgraphics(gcf, output_plot, 'Resolution', 150);
            close(gcf);
        end

        % regularize minimum density to water density
        density(skull_idx) = max(density(skull_idx),parameters.medium_properties.water.density);

        % regularize maximum density to rho_max
        density(skull_idx) = min(density(skull_idx),rho_max);

    case 'marsac'

        HU_min = 300;	  % minimum HU considered as skull
        HU_max = 2000;	  % maximum skull HU for regularization

        % cf. Marsac et al., 2017: do not exclude pHU < HU_min
        % as this is prone to create trabecular holes
        % skull_idx(pseudoCT(skull_idx) < HU_min) = [];

        % regularize minimum pHU to pHU_min
        pseudoCT(skull_idx) = max(pseudoCT(skull_idx),HU_min);

        % regularize maximum pHU to pHU_max
        pseudoCT(skull_idx) = min(pseudoCT(skull_idx),HU_max);

        rho_water     = parameters.medium_properties.water.density;      % density [kg/m^3]
        rho_bone      = 2100;     % max. skull density [kg/m3]

        % estimate density from CT HU based on Marsac et al., 2017 & Bancel et al., 2021
        % note: the original code hard-codes HU_min as 0, which may have been an error
        density(skull_idx) = rho_water + (rho_bone - rho_water) * ...
            (pseudoCT(skull_idx) - HU_min) / (HU_max - HU_min);

    case 'aubry'

        rho_water = parameters.medium_properties.water.density;
        rho_bone = parameters.medium_properties.skull.density;

        phi(skull_idx) = 1-(pseudoCT(skull_idx)/max(pseudoCT(skull_idx))); % [Aubry et al., 2003; Guo et al., 2019]
        density(skull_idx) = rho_water * phi(skull_idx) + ...
            rho_bone * (1-phi(skull_idx));

    case 'air-tissue'
        % Linear HU→density calibrated to air (HU≈-1000) and cortical bone (HU≈1000).
        % Coefficients: slope=1.01237293, intercept=1013.66593 kg/m³.
        density(skull_idx) = 1.01237293 * pseudoCT(skull_idx) + 1013.66593;
        density(skull_idx) = max(parameters.medium_properties.water.density, density(skull_idx));

    case 'none'

        density(skull_idx) = parameters.medium_properties.skull.density;

    otherwise
        error("Specified CT density mapping is not supported.")
end