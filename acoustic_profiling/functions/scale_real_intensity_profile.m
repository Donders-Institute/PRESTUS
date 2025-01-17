function adjusted_profile_focus = scale_real_intensity_profile(sim_param, desired_intensity, DENSITY_WATER, SOUND_SPEED_WATER, profile_focus)
    % Scale the real intensity profile to match the desired maximum intensity.
    %
    % Arguments:
    % - sim_param: Simulation parameters, including transducer properties.
    % - desired_intensity: Desired maximum intensity for the profile [W/cm^2].
    % - DENSITY_WATER: Density of water [kg/m^3].
    % - SOUND_SPEED_WATER: Speed of sound in water [m/s].
    % - profile_focus: Measured intensity profile at the focus [W/cm^2].
    %
    % Returns:
    % - adjusted_profile_focus: Adjusted intensity profile to match the desired intensity [W/cm^2].

    % Calculate maximum intensity and log adjustment details
    max_intens = max(profile_focus);
    fprintf('Current maximum intensity: %.2f \nDesired maximum intensity: %.2f \nAdjust profile to match desired maximum intensity.\n', max_intens, desired_intensity);

    % Calculate the corresponding pressure amplitude for the desired intensity
    p_pa = sqrt(2 * desired_intensity * 1e4 * DENSITY_WATER * SOUND_SPEED_WATER);

    % Update source amplitude in simulation parameters
    sim_param.transducer.source_amp = repmat(p_pa, 1, sim_param.transducer.n_elements);

    % Scale the intensity profile
    adjustment_factor_intensity = max_intens / desired_intensity;
    adjusted_profile_focus = profile_focus' ./ adjustment_factor_intensity;

end