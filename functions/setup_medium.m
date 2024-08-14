function kwave_medium = setup_medium(parameters, medium_masks, pseudoCT)

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

    % Creates an empty grid according to the grid dimensions
    grid_of_ones = ones(parameters.grid_dims);
    
    % Fills this grid with a baseline medium
    if strcmp(parameters.simulation_medium, 'water') || strcmp(parameters.simulation_medium, 'water_and_skull') || strcmp(parameters.simulation_medium, 'layered') 
        baseline_medium_name = 'water';
        baseline_medium = medium.(baseline_medium_name);
    elseif strcmp(parameters.simulation_medium, 'brain') || strcmp(parameters.simulation_medium, 'brain_and_skull')
        baseline_medium_name = 'brain';
        baseline_medium = medium.(baseline_medium_name);
    end
    
    % 'sound_speed' and 'density' are 3D arrays determining the speed of sound and the medium density within the simulation grid
    sound_speed = baseline_medium.sound_speed*grid_of_ones; 
    density = baseline_medium.density*grid_of_ones;
    alpha_0_true = baseline_medium.alpha_0_true*grid_of_ones;
    alpha_power_true = baseline_medium.alpha_power_true*grid_of_ones;

    % 'thermal conductivity' and 'specific_heat' are constants that define
    % the transmission of heat and the amount of heat required to change
    % the temperature within a givven mass of the tissue.
    thermal_conductivity = baseline_medium.thermal_conductivity * grid_of_ones;                                    % [W/(m.K)]
    specific_heat = baseline_medium.specific_heat_capacity * grid_of_ones;                                         % [J/(kg.K)]
    if isfield(parameters.thermal.temp_0, baseline_medium_name)
        temp_0 = parameters.thermal.temp_0.(baseline_medium_name) * grid_of_ones;                                  % [degC]
    else % if a global temp0 is indicated
        temp_0 = parameters.thermal.temp_0 * grid_of_ones;
    end

    % Changes the values of the acoustic and thermal properties in the
    % baseline_medium in the shape of the labelled mask
    if strcmp(parameters.simulation_medium, 'layered')
        labels = fieldnames(parameters.layer_labels);
        % Loops through each labelled layer to create a new mask
        for label_i = 1:length(labels)
            label_name = labels{label_i};
            if strcmp(label_name, 'water')
               continue % Loops to next label since the baseline medium is set to water anyways
            end
            if parameters.usepseudoCT == 1 && (strcmp(label_name, 'skull_cortical') || strcmp(label_name, 'skull_trabecular'))
                skull_idx = find(medium_masks==label_i);
                thermal_conductivity(skull_idx) = medium.(label_name).thermal_conductivity;                
                specific_heat(skull_idx) = medium.(label_name).specific_heat_capacity;
                if isfield(parameters.thermal.temp_0, label_name)
                    temp_0(skull_idx) = parameters.thermal.temp_0.(label_name);
                end

                % define default pCT-to-tissue conversion variant
                if isfield(parameters, "pseudoCT_variant")
                    pCT_variant = parameters.pseudoCT_variant;
                else
                    pCT_variant = "carpino";
                end

                switch pCT_variant
                    case 'carpino'

                        alpha_min   = 4;
                        alpha_max   = 8.7;
                        offset_HU   = 1000;
                        c_max       = 3100;     % max. speed of sound in skull (F. A. Duck, 2013.) [m/s]
                        rho_max     = 2100;     % max. density in skull [kg/m3]

                        % Finds maximum and minimum values
                        HU_min = min(pseudoCT(:));
                        HU_max = max(pseudoCT(:));

                        % Preprocess CT values
                        % Offset CT values to use housfield2density
                        pseudoCT(skull_idx) = pseudoCT(skull_idx) + offset_HU;
                        % replace negative values with 1 [JQK: regularize via minimum a la Yaakub?]
                        % note: this does not account for the offset, but perhaps this is desired
                        pseudoCT(skull_idx) = max(pseudoCT(skull_idx),1);
                                        
                        % estimate density
                        density(skull_idx) = hounsfield2density(pseudoCT(skull_idx));
                        % regularize minimum density to water density
                        density(skull_idx) = max(density(skull_idx),medium.water.density);
                        % regularize maximum density to rho_max
                        density(skull_idx) = min(density(skull_idx),rho_max);
                                        
                        % estimate sound speed
                        sound_speed(skull_idx) = 1.33.*density(skull_idx) + 167;
                        % regularize minimum sound speed to water
                        sound_speed(skull_idx) = max(sound_speed(skull_idx),medium.water.sound_speed);
                        % regularize maximum sound speed to c_max
                        sound_speed(skull_idx) = min(sound_speed(skull_idx),c_max);
                                        
                        % remove initial offset
                        pseudoCT(skull_idx) = pseudoCT(skull_idx)-offset_HU;

                    case 'yaakub'
                        % see https://github.com/sitiny/BRIC_TUS_Simulation_Tools/blob/main/tussim_skull_3D.m

                        c_water       = 1500;     % sound speed [m/s]
                        c_skull       = 3100;     % max. speed of sound in skull (F. A. Duck, 2013.) [m/s]
                        rho_water     = 1000;     % density [kg/m^3]
                        rho_bone      = 1900;     % max. skull density [kg/m3]
                        alpha_min     = 4;        % Aubry, J.-F. et al., 2022 (cortical bone; α (dB/cm at 500 kHz))
                        alpha_max     = 8.7;      % Fry 1978 at 0.5MHz: 1 Np/cm (8.7 dB/cm)
                        HU_min 	      = 300;	  % minimum HU considered as skull
                        HU_max 	      = 2000;	  % maximum skull HU for regularization

                        % if observed max HU is lower than threshold, set HU_max to actual max
                        HU_max = min(max(pseudoCT(skull_idx),[],'all'), HU_max);

                        % truncate CT HU (see Marsac et al., 2017)
                        % do not use pCT-based properties for presumed non-skull tissue
                        skull_idx(pseudoCT(skull_idx) < HU_min) = [];
                        % regularize upper HU to HU_max
                        idx_HUtoohigh = pseudoCT(skull_idx) > HU_max;
                        pseudoCT(skull_idx(idx_HUtoohigh)) = HU_max;

                        % estimate density from CT HU based on Marsac et al., 2017 & Bancel et al., 2021
                        % note: the original code hard-codes HU_min as 0, which may have been an error
                        density(skull_idx) = rho_water + (rho_bone - rho_water) * ...
                            (pseudoCT(skull_idx) - HU_min) / (HU_max - HU_min);

                        % estimate sound speed
                        sound_speed(skull_idx) = c_water + (c_skull - c_water) * ...
                            (density(skull_idx) - rho_water) / (rho_bone - rho_water);
                end

                % estimate attenuation coefficients
                alpha_pseudoCT(skull_idx) = alpha_min + (alpha_max - alpha_min) * ...
                    (1 - (pseudoCT(skull_idx) - HU_min) / (HU_max - HU_min)).^0.5;
                alpha_power_true(skull_idx) = medium.(label_name).alpha_power_true;
                % convert alpha at 500 kHz into prefactor alpha0 (dB/Mhz/cm) according to specified alpha_power_true
                % (definition of lower and upper attenuation bounds is derived from 500kHz)
                alpha_0_true(skull_idx) = alpha_pseudoCT(skull_idx)./(0.5^medium.(label_name).alpha_power_true);
                clear skull_idx
            else
                % Sets the parameters in the shape of the mask
                thermal_conductivity(medium_masks==label_i) = medium.(label_name).thermal_conductivity;                 % [W/(m.K)]
                specific_heat(medium_masks==label_i) = medium.(label_name).specific_heat_capacity;                      % [J/(kg.K)]
                sound_speed(medium_masks==label_i) = medium.(label_name).sound_speed; 
                density(medium_masks==label_i) = medium.(label_name).density;
                alpha_0_true(medium_masks==label_i) =  medium.(label_name).alpha_0_true; 
                alpha_power_true(medium_masks==label_i) = medium.(label_name).alpha_power_true;
                if isfield(parameters.thermal.temp_0, label_name)
                    temp_0(medium_masks==label_i) = parameters.thermal.temp_0.(label_name);                             % [degC]
                end
            end
        end
    elseif ~isempty(medium_masks) % Use the medium_masks if one is specified and the simulation_medium is not layered
        % Sets the parameters in the shape of the medium_masks
        % identify tissue ids
        if isfield(parameters, 'layer_labels')
            labels = fieldnames(parameters.layer_labels);
            i_skull = find(contains(labels,  'skull')); i_skull = i_skull(1);
            % i_brain = find(contains(labels,  'brain'));
        else
            % warning("Trying to infer tissue labels; please check...")
            masklabels = unique(medium_masks);
            masklabels(masklabels==0) = [];
            % warning("Assuming brain layer is specified before skull...")
            i_skull = masklabels(1);
            % i_brain = masklabels(1);
        end
        if contains(parameters.simulation_medium, 'skull')
            thermal_conductivity(medium_masks==i_skull) = medium.skull.thermal_conductivity;
            specific_heat(medium_masks==i_skull) = medium.skull.specific_heat_capacity;
            sound_speed(medium_masks==i_skull) = medium.skull.sound_speed; 
            density(medium_masks==i_skull) = medium.skull.density;
            alpha_0_true(medium_masks==i_skull) =  medium.skull.alpha_0_true; 
            alpha_power_true(medium_masks==i_skull) = medium.skull.alpha_power_true;
            if isfield(parameters.thermal.temp_0, 'skull')
                temp_0(medium_masks==i_skull) = parameters.thermal.temp_0.skull;
            end
        end
    end

    % Account for actual absorption behaviour in k-Wave, which varies when high
    % absorption is used (see https://doi.org/10.1121/1.4894790).

    % 'alpha_coeff' is tweaked so any alpha_power could be used, 
    % as long as the alpha_0_true and alpa_power_true are correct for a
    % given frequency.

    alpha_power_fixed = 2;

    alpha_coeff = fitPowerLawParamsMulti(...
        alpha_0_true, ...
        alpha_power_true, ...
        sound_speed, ...
        parameters.transducer.source_freq_hz, ...
        alpha_power_fixed, ...
        false);
    
    % Outputs the medium as a structure
    kwave_medium = struct('sound_speed', sound_speed, ...
                          'density', density, ...
                          'alpha_coeff', alpha_coeff,...
                          'alpha_power', alpha_power_fixed, ...
                          'thermal_conductivity', thermal_conductivity,...
                          'specific_heat', specific_heat,...
                          'temp_0', temp_0);
    
    % save images for debugging
    if (contains(parameters.simulation_medium, 'skull') || contains(parameters.simulation_medium, 'layered'))
        if ~exist(fullfile(parameters.output_dir, 'debug')); mkdir(parameters.output_dir, 'debug'); end
        % filename_density = fullfile(parameters.output_dir, 'debug', sprintf('matrix_density'));
        % niftiwrite(density, filename_density, 'Compressed',true);
        % filename_sound_speed = fullfile(parameters.output_dir, 'debug', sprintf('matrix_sound_speed'));
        % niftiwrite(sound_speed, filename_sound_speed, 'Compressed',true);
        % filename_alpha_0_true = fullfile(parameters.output_dir, 'debug', sprintf('matrix_alpha_0_true'));
        % niftiwrite(alpha_0_true, filename_alpha_0_true, 'Compressed',true);
        % filename_alpha_power_true = fullfile(parameters.output_dir, 'debug', sprintf('matrix_alpha_power'));
        % niftiwrite(alpha_power_true, filename_alpha_power_true, 'Compressed',true);
        % filename_alpha_coeff = fullfile(parameters.output_dir, 'debug', sprintf('matrix_alpha_coeff'));
        % niftiwrite(alpha_coeff, filename_alpha_coeff, 'Compressed',true);
    end
end