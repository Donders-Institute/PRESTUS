function kwave_medium = medium_setup(parameters, medium_masks, planimg, pseudoCT)

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                         Setup k-wave medium                       %
    %                                                                   %
    % This function sets up the medium for k-wave simulations to take   %
    % place in. The way this mask is set up is by starting with an      %
    % empty grid and adjust the acoustic and thermal parameters in the  %
    % shape of the provided mask. This can be repeated to incorporate   %
    % different kinds of tissue (see 'layered').                        %
    %                                                                   %
    % Some notes:                                                       %
    % - The alpha value is converted at the end of this script.         %
    % - pseudoCT is an optional input when using pseudoCTs              %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % Loads the medium settings from the config file
    medium = parameters.medium;
    
    % Create empty matrices for medium properties
    empty_grid = NaN(parameters.grid_dims);
    sound_speed = empty_grid; 
    density = empty_grid;
    alpha_coeff = empty_grid;
    alpha_power = empty_grid;
    thermal_conductivity = empty_grid;       % [W/(m.K)]
    specific_heat = empty_grid;            % [J/(kg.K)]
    perfusion = empty_grid;
    absorption_fraction = ones(parameters.grid_dims); % default: convert 100% of attenuation into absorption
    temp_0 = empty_grid;

    % Get layer and medium labels
    layer_labels = fieldnames(parameters.layers);
    medium_labels = fieldnames(parameters.medium);

    % Iterate through each layer & assign medium ID to create medium mask
    for label_i = 1:length(layer_labels)
        label_name = layer_labels{label_i};
        medium_i = find(strcmp(medium_labels, label_name));
        if parameters.usepseudoCT == 1 && strcmp(label_name, 'skull')
            skull_idx = find(ismember(medium_masks,medium_i));
            thermal_conductivity(skull_idx) = medium.(label_name).thermal_conductivity;                
            specific_heat(skull_idx) = medium.(label_name).specific_heat_capacity;
            if isfield(medium.(label_name), 'perfusion')
                perfusion(skull_idx) = medium.(label_name).perfusion;
            end
            if isfield(medium.(label_name), 'absorption_fraction')
                absorption_fraction(skull_idx) = medium.(label_name).absorption_fraction;
            end
            if isfield(parameters.thermal.temp_0, label_name)
                temp_0(skull_idx) = parameters.thermal.temp_0.(label_name);
            end

            % define default pCT-to-tissue conversion variant
            if isfield(parameters, "pseudoCT_variant")
                pCT_variant = parameters.pseudoCT_variant;
            else
                pCT_variant = "kosciessa";
            end

            switch pCT_variant
                case 'carpino'
                    pct_skullmapping.density = 'k-wave';
                    pct_skullmapping.soundspeed = 'k-plan';
                    pct_skullmapping.attenuation = 'yakuub';
                case 'yakuub'
                    pct_skullmapping.density = 'marsac';
                    pct_skullmapping.soundspeed = 'marsac';
                    pct_skullmapping.attenuation = 'yakuub';
                case 'k-plan'
                    pct_skullmapping.density = 'k-plan';
                    pct_skullmapping.soundspeed = 'k-plan';
                    pct_skullmapping.attenuation = 'k-plan';
                case 'marquet'
                    pct_skullmapping.density = 'marquet';
                    pct_skullmapping.soundspeed = 'marquet';
                    pct_skullmapping.attenuation = 'k-plan';
                case 'kosciessa'
                    pct_skullmapping.density = 'k-plan';
                    pct_skullmapping.soundspeed = 'marsac';
                    pct_skullmapping.attenuation = 'yakuub';
                otherwise
                    pct_skullmapping.density = parameters.pct_skullmapping.density;
                    pct_skullmapping.soundspeed = parameters.pct_skullmapping.soundspeed;
                    pct_skullmapping.attenuation = parameters.pct_skullmapping.attenuation;
            end

            switch pct_skullmapping.density
                case 'k-plan'

                    % define piece-wise linear mapping between HU and mass density in kg/m^3 hounsfieldUnits = [-990, 60, 1000, 1950]; 
                    % see https://dispatch.k-plan.io/static/docs/planning-images.html#ct-calibration
                    hounsfieldUnits = [-990, 60, 1000, 1950]; 
                    massDensity = [1.2, 1060, 1530, 2150]; 

                    density(skull_idx) = fit_pairwiselinear(pseudoCT(skull_idx), hounsfieldUnits, massDensity, 1);

                    % plot the mapping
                    if parameters.debug == 1
                        output_plot = fullfile(parameters.debug_dir, ...
                            sprintf('pCT_hounsfield-density_kplan.png'));
                        exportgraphics(gcf, output_plot, 'Resolution', 150);
                    end
                    close(gcf);

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
                    if parameters.debug == 1
                        output_plot = fullfile(parameters.debug_dir, ...
                            sprintf('pCT_hounsfield-density_kwave.png'));
                        exportgraphics(gcf, output_plot, 'Resolution', 150);
                    end
                    close(gcf);

                    % regularize minimum density to water density
                    density(skull_idx) = max(density(skull_idx),medium.water.density);
                    % regularize maximum density to rho_max
                    density(skull_idx) = min(density(skull_idx),rho_max);

                    % remove initial offset
                    pseudoCT(skull_idx) = pseudoCT(skull_idx)-offset_HU;

                case 'marsac'

                    HU_min 	      = 300;	  % minimum HU considered as skull
                    HU_max 	      = 2000;	  % maximum skull HU for regularization
                    rho_water     = 996;      % density [kg/m^3] [Baumgartner et al., 2024]
                    rho_bone      = 2100;     % max. skull density [kg/m3] [Baumgartner et al., 2024]

                    % if observed max HU is lower than threshold, set HU_max to actual max
                    HU_max = min(max(pseudoCT(skull_idx),[],'all'), HU_max);

                    % truncate CT HU (see Marsac et al., 2017)
                    % do not use pCT-based properties for presumed non-skull tissue
                    % BUT: this is prone to create trabecular holes
                    skull_idx(pseudoCT(skull_idx) < HU_min) = [];
                    % regularize maximum HU to HU_max
                    pseudoCT(skull_idx) = min(pseudoCT(skull_idx),HU_max);

                    % estimate density from CT HU based on Marsac et al., 2017 & Bancel et al., 2021
                    % note: the original code hard-codes HU_min as 0, which may have been an error
                    density(skull_idx) = rho_water + (rho_bone - rho_water) * ...
                        (pseudoCT(skull_idx) - HU_min) / (HU_max - HU_min);

                case 'marquet'

                    rho_water = 1000;
                    rho_bone = 2200;

                    phi(skull_idx) = 1-(pseudoCT(skull_idx)/1000);
                    density(skull_idx) = rho_water * phi(skull_idx) + ...
                        rho_bone * (1-phi(skull_idx));

                otherwise
                    error("Specified CT density mapping is not supported.")
            end

            switch pct_skullmapping.soundspeed
                case 'k-plan'

                    sound_speed(skull_idx) = 1.33.*density(skull_idx) + 167;

                case 'marsac'

                    c_water       = 1500;     % sound speed [m/s]
                    c_skull       = 3360;     % max. speed of sound in skull [m/s] [Baumgartner et al., 2024]
                    rho_water     = 996;      % density [kg/m^3] [Baumgartner et al., 2024]
                    rho_bone      = 2100;     % max. skull density [kg/m3] [Baumgartner et al., 2024]

                    sound_speed(skull_idx) = c_water + (c_skull - c_water) * ...
                        (density(skull_idx) - rho_water) / (rho_bone - rho_water);
                    
                case 'marquet'

                    c_water = 1500;
                    if strcmp(label_name, 'skull_cortical')
                        c_bone = 3100;
                    else
                        c_bone = 2200;
                    end

                    phi(skull_idx) = 1-(pseudoCT(skull_idx)/1000);
                    sound_speed(skull_idx) = c_water * phi(skull_idx) + ...
                        c_bone * (1-phi(skull_idx));
                    % regularize sound speed to a minimum of water
                    sound_speed(skull_idx) = max(sound_speed(skull_idx),c_water);

                otherwise
                    error("Specified CT sound speed mapping is not supported.")
            end

            switch pct_skullmapping.attenuation
                case 'k-plan'

                    kPlan_alpha = 13.3; % https://dispatch.k-plan.io/static/docs/simulation-pipeline.html
                    kPlan_alpha_power = 1;
                    % Note that we allow different values to be specified in the config.
                    % If replication of k-Wave is the goal, the above values should be specified.
                    % Throw a warning in the case of deviations.
                    if medium.(label_name).alpha_coeff ~= kPlan_alpha || ...
                        medium.(label_name).alpha_power ~= kPlan_alpha_power
                        warning('Specified attenuation varies from k-Plan setup.')
                    end
                    alpha_coeff(skull_idx) = medium.(label_name).alpha_coeff;
                    alpha_power(skull_idx) = medium.(label_name).alpha_power;

                case 'yakuub'

                    alpha_min     = 4;        % cortical bone at 500 kHz [dB/cm] [Aubry et al., 2022] 
                    alpha_max     = 8.7;      % bone at 500 kHz [dB/cm] [Fry 1978]

                    % Finds maximum and minimum values
                    HU_min = min(pseudoCT(skull_idx));
                    HU_max = max(pseudoCT(skull_idx));

                    % estimate attenuation based on (pseudo-)HU
                    alpha_pseudoCT(skull_idx) = alpha_min + (alpha_max - alpha_min) * ...
                        (1 - (pseudoCT(skull_idx) - HU_min) / (HU_max - HU_min)).^0.5;
                    alpha_power(skull_idx) = medium.(label_name).alpha_power;
                    % convert alpha at 500 kHz into prefactor alpha0 (dB/MHz/cm) according to specified alpha_power
                    % (definition of lower and upper attenuation bounds is derived from 500kHz)
                    alpha_coeff(skull_idx) = alpha_pseudoCT(skull_idx)./(0.5^medium.(label_name).alpha_power);

                otherwise
                    error("Specified pCT attenuation mapping is not supported.")
            end
            
            % [DEBUG] save pCT mapping overview
            if parameters.debug == 1
                h = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.25, 1]);
                subplot(4,1,1); 
                    hold on; histogram(pseudoCT(skull_idx)); xlabel("pseudo-HU")
                    title(['pseudoCT tissue property ranges: ', parameters.pseudoCT_variant]);
                subplot(4,1,2); 
                    hold on; histogram(density(skull_idx)); xlabel("Density [kg/m3]")
                    % add lines for the fixed parameters
                    xline(medium.skull_trabecular.density, 'r', 'LineWidth', 2);
                    xline(medium.skull_cortical.density, 'r', 'LineWidth', 2);
                subplot(4,1,3); 
                    hold on; histogram(sound_speed(skull_idx)); xlabel("Sound speed [m/s]")
                    xline(medium.skull_trabecular.sound_speed, 'r', 'LineWidth', 2);
                    xline(medium.skull_cortical.sound_speed, 'r', 'LineWidth', 2);
                subplot(4,1,4); 
                    hold on; histogram(alpha_coeff(skull_idx)); xlabel("Attenuation [dB/(cm.MHzy)]")
                    xline(medium.skull_trabecular.alpha_coeff, 'r', 'LineWidth', 2);
                    xline(medium.skull_cortical.alpha_coeff, 'r', 'LineWidth', 2);
                output_plot = fullfile(parameters.debug_dir, ...
                    sprintf('pCT_histograms_%s.png',parameters.pseudoCT_variant));
                exportgraphics(h, output_plot, 'Resolution', 150);
                close(h);
            end

            clear skull_idx

        elseif parameters.usepseudoCT == 1 && contains(label_name, 'skull') && skull_mapped == 1
            disp(['pCT: ', label_name, ' requested, but skull is already mapped. Skipping...']);
        else
            thermal_conductivity(medium_masks==medium_i) = medium.(label_name).thermal_conductivity;                 % [W/(m.K)]
            specific_heat(medium_masks==medium_i) = medium.(label_name).specific_heat_capacity;                      % [J/(kg.K)]
            sound_speed(medium_masks==medium_i) = medium.(label_name).sound_speed; 
            density(medium_masks==medium_i) = medium.(label_name).density;
            alpha_coeff(medium_masks==medium_i) = medium.(label_name).alpha_coeff; 
            alpha_power(medium_masks==medium_i) = medium.(label_name).alpha_power;
            if isfield(medium.(label_name), 'perfusion')
                perfusion(medium_masks==medium_i) = medium.(label_name).perfusion;
            end
            if isfield(medium.(label_name), 'absorption_fraction')
                absorption_fraction(medium_masks==medium_i) = medium.(label_name).absorption_fraction;
            end
            if isfield(parameters.thermal.temp_0, label_name)
                temp_0(medium_masks==medium_i) = parameters.thermal.temp_0.(label_name);                             % [degC]
            end
        end
    end

    % account for k-Wave's actual attenuation behaviour
    % limits discrepancies for high attenuation estimates (see https://doi.org/10.1121/1.4894790).
    % 'alpha_coeff' is rescaled to match the specified alpha_coeff and alpa_power_true for the center frequency

    if ~isfield(parameters, 'fit_alpha_power')
        parameters.fit_alpha_power = 1;
    end

    if parameters.fit_alpha_power == 1
        alpha_power_fixed = 2;
    
        if parameters.debug == 1
            plot_fit = true;
        else
            plot_fit = false;
        end
    
        if numel(unique([parameters.transducer.source_freq_hz])) > 1
            warning('Multiple source frequencies not yet supported. Using %i Hz from first transducer.', ...
                    parameters.transducer(1).source_freq_hz)
        end
        
        alpha_coeff_fixed = fitPowerLawParamsMulti(...
            alpha_coeff, ...
            alpha_power, ...
            sound_speed, ...
            parameters.transducer(1).source_freq_hz, ...
            alpha_power_fixed, ...
            plot_fit);
    
        % DEBUG mode: save plot of fitted attenuation values
        if parameters.debug == 1
            fig_path = fullfile(parameters.debug_dir, ...
            ['attenuation_fit', char(parameters.results_filename_affix), '.png']);
            saveas(gcf, fig_path);
            close(gcf);
        end
    else
        warning('Attenuation power-law is not remapped (as requested).')
        alpha_coeff_fixed = alpha_coeff;
        alpha_power_fixed = alpha_power;
    end

    % convert perfusion rate [mL/min/kg] into perfusion coefficient [1/s]
    perfusion_coeff = (perfusion ./ 60) .* density * 1e-6; % [1/s]

    % specify the medium as a kWave-compatible structure 
    % Note: absorption_fraction and temp_0 need to be removed later)
    kwave_medium = struct('sound_speed', sound_speed, ...
                          'density', density, ...
                          'alpha_coeff', alpha_coeff_fixed,...
                          'alpha_power', alpha_power_fixed, ...
                          'thermal_conductivity', thermal_conductivity,...
                          'specific_heat', specific_heat,...
                          'perfusion_coeff', perfusion_coeff,...
                          'absorption_fraction', absorption_fraction,...
                          'temp_0', temp_0);
    
    % save images for debugging
    if parameters.debug == 1
        try
            filename_density = fullfile(parameters.debug_dir, sprintf('matrix_density'));
            niftiwrite(density, filename_density, 'Compressed',true);
            pause(0.1);
            filename_sound_speed = fullfile(parameters.debug_dir, sprintf('matrix_sound_speed'));
            niftiwrite(sound_speed, filename_sound_speed, 'Compressed',true);
            pause(0.1);
            filename_alpha_coeff = fullfile(parameters.debug_dir, sprintf('matrix_alpha_coeff'));
            niftiwrite(alpha_coeff, filename_alpha_coeff, 'Compressed',true);
            pause(0.1);
            filename_alpha_power = fullfile(parameters.debug_dir, sprintf('matrix_alpha_power'));
            niftiwrite(alpha_power, filename_alpha_power, 'Compressed',true);
            pause(0.1);
            filename_alpha_coeff_fixed = fullfile(parameters.debug_dir, sprintf('matrix_alpha_coeff_fixed'));
            niftiwrite(alpha_coeff_fixed, filename_alpha_coeff_fixed, 'Compressed',true);
            pause(0.1);
            filename_perfusion = fullfile(parameters.debug_dir, sprintf('matrix_perfusion'));
            niftiwrite(perfusion_coeff, filename_perfusion, 'Compressed',true);
            pause(0.1);
            filename_absorption = fullfile(parameters.debug_dir, sprintf('matrix_absorption'));
            niftiwrite(absorption_fraction, filename_absorption, 'Compressed',true);
        catch
            warning("Error with saving debug images: medium mapping. May result from concurrent write attempts...")
        end
    end

    % save images of assigned medium properties
    if contains(parameters.simulation_medium, {'layered'}) && parameters.debug == 1 && (exist('planimg') & ~isempty(planimg))
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'sound_speed')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'density')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_coeff')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_power')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'thermal_conductivity')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'specific_heat')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'perfusion_coeff')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'absorption_fraction')
    end

    % save a pCT if used
    if parameters.usepseudoCT == 1 && parameters.debug == 1
        filename_pct = fullfile(parameters.debug_dir, sprintf('pct%s', parameters.results_filename_affix));
        niftiwrite(pseudoCT, filename_pct, 'Compressed',true);
    end

end