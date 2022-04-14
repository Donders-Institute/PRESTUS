function kwave_medium = setup_medium(parameters, skull_mask)

    medium = parameters.medium;

    grid_of_ones = ones(parameters.grid_dims);
    % sound_speed and density are 3D arrays determining the speed of sound and the medium density within the simulation grid
    sound_speed = medium.water_sound_speed*grid_of_ones; 
    density = medium.water_density*grid_of_ones;    %waterDensity(temp_0) * ones(Nx,Ny,Nz);             % [kg/m^3]

    % assume a homogeneous attenuation across the skull with 0.6 dB/MHz/cm in water and 7.4 dB/MHz/cm in the skull
    alpha_coeff_true =  medium.water_alpha_coef_true*grid_of_ones; %0.6 * ones(Nx,Ny,Nz);                            % [dB/MHz/cm] (Kremkau et al., 1981)
    % alpha_coeff_true(air_mask) = 12;
    alpha_power_true = medium.water_alpha_power_true*grid_of_ones;
    if ~isempty(skull_mask)
        sound_speed(skull_mask) = medium.skull_sound_speed; 
        density(skull_mask) = medium.skull_density;
        alpha_coeff_true(skull_mask) =  medium.skull_alpha_power_true; %7.4; % [dB/MHz/cm] (FUN21 Conference tutorial, Treeby, 2021)
        alpha_power_true(skull_mask) = medium.skull_alpha_power_true;
    end
    % account for actual absorption behaviour in k-Wave, which varies when high
    % absorption is used (see https://doi.org/10.1121/1.4894790)
    alpha_coeff = fitPowerLawParamsMulti(alpha_coeff_true, alpha_power_true, sound_speed, parameters.transducer.source_freq_hz, medium.alpha_power);

    kwave_medium = struct('sound_speed', sound_speed, ...
                          'density', density, ...
                          'alpha_coeff',alpha_coeff,...
                          'alpha_power', parameters.medium.alpha_power);
end