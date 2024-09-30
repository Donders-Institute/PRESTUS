function kwave_medium = setup_medium(parameters, medium_mask)

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
    % - The alpha value can be tweaked at the end of this script.       %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % Loads the medium settings from the config file
    medium = parameters.medium;

    % Creates an empty grid according to the grid dimensions
    grid_of_ones = ones(parameters.grid_dims);
    
    % Fills this grid with a baseline medium
    if strcmp(parameters.simulation_medium, 'water') || strcmp(parameters.simulation_medium, 'water_and_skull') || strcmp(parameters.simulation_medium, 'layered') 
        baseline_medium = medium.('water');
    elseif strcmp(parameters.simulation_medium, 'brain') || strcmp(parameters.simulation_medium, 'brain_and_skull')
        baseline_medium = medium.('brain');
    end
    
    % 'sound_speed' and 'density' are 3D arrays determining the speed of sound and the medium density within the simulation grid
    sound_speed = baseline_medium.sound_speed*grid_of_ones; 
    density = baseline_medium.density*grid_of_ones;    %waterDensity(temp_0) * ones(Nx,Ny,Nz);                     % [kg/m^3]

    % Assume a homogeneous attenuation across the skull of 0.6 dB/MHz/cm in water and 7.4 dB/MHz/cm in the skull
    alpha_0_true =  baseline_medium.alpha_0_true*grid_of_ones; %0.6 * ones(Nx,Ny,Nz);                              % [dB/MHz/cm] (Kremkau et al., 1981)
    alpha_power_true = baseline_medium.alpha_power_true*grid_of_ones;

    % 'thermal conductivity' and 'specific_heat' are constants that define
    % the transmission of heat and the amount of heat required to change
    % the temperature within a givven mass of the tissue.
    thermal_conductivity = baseline_medium.thermal_conductivity * grid_of_ones;                                    % [W/(m.K)]
    specific_heat = baseline_medium.specific_heat_capacity * grid_of_ones;                                         % [J/(kg.K)]
    
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
            % Sets the parameters in the shape of the mask
            thermal_conductivity(medium_mask==label_i) = medium.(label_name).thermal_conductivity;                 % [W/(m.K)]
            specific_heat(medium_mask==label_i) = medium.(label_name).specific_heat_capacity;                      % [J/(kg.K)]

            sound_speed(medium_mask==label_i) = medium.(label_name).sound_speed; 
            density(medium_mask==label_i) = medium.(label_name).density;
            alpha_0_true(medium_mask==label_i) =  medium.(label_name).alpha_power_true; 
            alpha_power_true(medium_mask==label_i) = medium.(label_name).alpha_power_true;
        end
    elseif ~isempty(medium_mask) % Use the medium_mask if one is specified and the simulation_medium is not layered
        % Sets the parameters in the shape of the medium_mask
        thermal_conductivity(medium_mask) = medium.skull.thermal_conductivity;                                     % [W/(m.K)]
        specific_heat(medium_mask) = medium.skull.specific_heat_capacity;                                          % [J/(kg.K)]

        sound_speed(medium_mask) = medium.skull.sound_speed; 
        density(medium_mask) = medium.skull.density;
        alpha_0_true(medium_mask) =  medium.skull.alpha_power_true; 
        alpha_power_true(medium_mask) = medium.skull.alpha_power_true;
    end
    
    % Account for actual absorption behaviour in k-Wave, which varies when high
    % absorption is used (see https://doi.org/10.1121/1.4894790).

    % 'alpha_coeff' is tweaked so any alpha_power could be used, 
    % as long as the alpha_0_true and alpa_power_true are correct for a
    % given frequency.
    % Here, I use 2.
    alpha_power_fixed = 2;
    
    alpha_coeff = fitPowerLawParamsMulti(alpha_0_true, alpha_power_true, sound_speed, parameters.transducer.source_freq_hz, alpha_power_fixed );

    % Outputs the medium as a structure
    kwave_medium = struct('sound_speed', sound_speed, ...
                          'density', density, ...
                          'alpha_coeff',alpha_coeff,...
                          'alpha_power', alpha_power_fixed , ...
                          'thermal_conductivity', thermal_conductivity,...
                          'specific_heat', specific_heat);

end