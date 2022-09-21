function kwave_medium = setup_medium(parameters, medium_mask)

    medium = parameters.medium;

    grid_of_ones = ones(parameters.grid_dims);
    
    if strcmp(parameters.simulation_medium, 'water') || strcmp(parameters.simulation_medium, 'water_and_skull') || strcmp(parameters.simulation_medium, 'layered') 
        baseline_medium = medium.('water');
    elseif strcmp(parameters.simulation_medium, 'brain') || strcmp(parameters.simulation_medium, 'brain_and_skull')
        baseline_medium = medium.('brain');
    end
    
    % sound_speed and density are 3D arrays determining the speed of sound and the medium density within the simulation grid
    sound_speed = baseline_medium.sound_speed*grid_of_ones; 
    density = baseline_medium.density*grid_of_ones;    %waterDensity(temp_0) * ones(Nx,Ny,Nz);             % [kg/m^3]

    % assume a homogeneous attenuation across the skull with 0.6 dB/MHz/cm in water and 7.4 dB/MHz/cm in the skull
    alpha_0_true =  baseline_medium.alpha_0_true*grid_of_ones; %0.6 * ones(Nx,Ny,Nz);                            % [dB/MHz/cm] (Kremkau et al., 1981)
    alpha_power_true = baseline_medium.alpha_power_true*grid_of_ones;

    thermal_conductivity = baseline_medium.thermal_conductivity * grid_of_ones;          % [W/(m.K)]
    specific_heat = baseline_medium.specific_heat_capacity * grid_of_ones;               % [J/(kg.K)]
    
    if strcmp(parameters.simulation_medium, 'layered')
        labels = fieldnames(parameters.layer_labels);
        for label_i = 1:length(labels)
            label_name = labels{label_i};
            if strcmp(label_name, 'water')
               continue
            end
            thermal_conductivity(medium_mask==label_i) = medium.(label_name).thermal_conductivity;             % [W/(m.K)]
            specific_heat(medium_mask==label_i) = medium.(label_name).specific_heat_capacity;                  % [J/(kg.K)]

            sound_speed(medium_mask==label_i) = medium.(label_name).sound_speed; 
            density(medium_mask==label_i) = medium.(label_name).density;
            alpha_0_true(medium_mask==label_i) =  medium.(label_name).alpha_power_true; 
            alpha_power_true(medium_mask==label_i) = medium.(label_name).alpha_power_true;
        end
    elseif ~isempty(medium_mask)
        % thermal parameters
        
        thermal_conductivity(medium_mask) = medium.skull.thermal_conductivity;             % [W/(m.K)]
        specific_heat(medium_mask) = medium.skull.specific_heat_capacity;                  % [J/(kg.K)]

        sound_speed(medium_mask) = medium.skull.sound_speed; 
        density(medium_mask) = medium.skull.density;
        alpha_0_true(medium_mask) =  medium.skull.alpha_power_true; 
        alpha_power_true(medium_mask) = medium.skull.alpha_power_true;
    end
    
    % account for actual absorption behaviour in k-Wave, which varies when high
    % absorption is used (see https://doi.org/10.1121/1.4894790)
    % alpha_coeff is tweaked so any alpha_power could be used, 
    % as long as the alpha_0_true and alpa_power_true are correct for a
    % given frequency.
    % here, I use 2
    
    alpha_power_fixed = 2;
    
    alpha_coeff = fitPowerLawParamsMulti(alpha_0_true, alpha_power_true, sound_speed, parameters.transducer.source_freq_hz, alpha_power_fixed );

    kwave_medium = struct('sound_speed', sound_speed, ...
                          'density', density, ...
                          'alpha_coeff',alpha_coeff,...
                          'alpha_power', alpha_power_fixed , ...
                          'thermal_conductivity', thermal_conductivity,...
                          'specific_heat', specific_heat);
end