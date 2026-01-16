function adjusted_profile_focus = scale_real_intensity_profile(parameters, desired_intensity, profile_focus)
    % Scale the real intensity profile to match the desired maximum intensity.
    %
    % Arguments:
    % - parameters: Simulation parameters, including transducer properties.
    %       .medium.water.density: Density of water [kg/m^3].
    %       .medium.water.sound_speed: Speed of sound in water [m/s].
    %       .transducer: transducer information
    % - desired_intensity: Desired maximum intensity for the profile [W/cm^2].
    % - profile_focus: Measured intensity profile at the focus [W/cm^2].
    %
    % Returns:
    % - adjusted_profile_focus: Adjusted intensity profile to match the desired intensity [W/cm^2].

    % Calculate maximum intensity and log adjustment details
    max_intens = max(profile_focus);
    fprintf('Current maximum intensity: %.2f \nDesired maximum intensity: %.2f \nAdjust profile to match desired maximum intensity.\n', max_intens, desired_intensity);

    % Calculate the corresponding pressure amplitude for the desired intensity

    DENSITY_WATER = parameters.medium.water.density;
    SOUND_SPEED_WATER = parameters.medium.water.sound_speed;
    p_pa = sqrt(2 * desired_intensity * 1e4 * DENSITY_WATER * SOUND_SPEED_WATER);

    % Update source amplitude in simulation parameters
    parameters.transducer.source_amp = repmat(p_pa, 1, parameters.transducer.n_elements);

    % Linearly scale the intensity profile
    adjustment_factor_intensity = max_intens / desired_intensity;
    adjusted_profile_focus = profile_focus' ./ adjustment_factor_intensity;

end