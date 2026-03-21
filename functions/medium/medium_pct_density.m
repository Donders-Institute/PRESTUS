function [density] = medium_pct_density(parameters, density, pseudoCT, skull_idx, algorithm)

switch algorithm
    case 'k-plan'

        % define piece-wise linear mapping between HU and mass density in kg/m^3 hounsfieldUnits = [-990, 60, 1000, 1950]; 
        % see https://dispatch.k-plan.io/static/docs/planning-images.html#ct-calibration
        hounsfieldUnits = [-990, 60, 1000, 1950]; 
        massDensity = [1.2, 1060, 1530, 2150]; 

        density(skull_idx) = fit_pairwiselinear(pseudoCT(skull_idx), hounsfieldUnits, massDensity, 1);

        % plot the mapping
        if parameters.simulation.debug == 1
            output_plot = fullfile(parameters.io.debug_dir, ...
                sprintf('pCT_hounsfield-density_kplan.png'));
            exportgraphics(gcf, output_plot, 'Resolution', 150);
        end
        close(gcf);

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
        density(skull_idx) = hounsfield2density(pseudoCT(skull_idx), 1);
        
        % plot the mapping
        if parameters.simulation.debug == 1
            output_plot = fullfile(parameters.io.debug_dir, ...
                sprintf('pCT_hounsfield-density_kwave.png'));
            exportgraphics(gcf, output_plot, 'Resolution', 150);
        end
        close(gcf);

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

    case 'none'

        density(skull_idx) = parameters.medium_properties.skull.density;

    otherwise
        error("Specified CT density mapping is not supported.")
end