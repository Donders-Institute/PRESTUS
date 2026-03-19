function [profile_target, parameters] = scale_real_intensity_profile(parameters, profile_empirical, desired_intensity)
    % Scale the real intensity profile to match the desired maximum intensity.
    % Match the length of the target profile to the simulation axis length.
    %
    % Arguments:
    % - parameters: Simulation parameters, including transducer properties.
    %       .medium.water.density: Density of water [kg/m^3].
    %       .medium.water.sound_speed: Speed of sound in water [m/s].
    %       .transducer: transducer information
    % - profile_empirical.axial_intensity: Measured or theoretical intensity profile along beam axis [W/cm^2]
    % - profile_empirical.axial_distance_bowl: Distance from transducer reference point (mm)
    % - desired_intensity: Desired maximum intensity for the profile [W/cm^2].
    %
    % Returns:
    % - profile_target.axial_intensity: Adjusted intensity profile to match the desired intensity [W/cm^2].
    % - profile_target.axial_distance_bowl: Distance from transducer reference point (mm)
    % - parameters: updated with scaled parameters.transducer.source_amp

    % Calculate maximum intensity and log adjustment details
    max_intens = max(profile_empirical.axial_intensity);
    fprintf('Current maximum intensity: %.2f \nDesired maximum intensity: %.2f \nAdjust profile to match desired maximum intensity.\n', max_intens, desired_intensity);

    % Calculate the corresponding pressure amplitude for the desired intensity

    DENSITY_WATER = parameters.medium.water.density;
    SOUND_SPEED_WATER = parameters.medium.water.sound_speed;
    p_pa = sqrt(2 * desired_intensity * 1e4 * DENSITY_WATER * SOUND_SPEED_WATER);

    % Update source amplitude in simulation parameters
    parameters.transducer.source_amp = repmat(p_pa, 1, parameters.transducer.n_elements);

    % Linearly scale the intensity profile
    adjustment_factor_intensity = max_intens / desired_intensity;
    adjusted_axial_intensity = profile_empirical.axial_intensity' ./ adjustment_factor_intensity;

    profile_target.axial_intensity = adjusted_axial_intensity;
    profile_target.axial_distance_bowl = profile_empirical.axial_distance_bowl';

    % Calculate size of the simulation domain
    sim_axis_mm = (parameters.grid.default_dims(end)-parameters.grid.pml_size-1)*parameters.grid.resolution_mm;

    % Dynamically truncate or pad to match simulation domain
    if profile_empirical.axial_distance_bowl(end) > sim_axis_mm
        % Case 1: Remove values beyond simulation length
        valid_idx = profile_target.axial_distance_bowl <= sim_axis_mm;
        profile_target.axial_distance_bowl = profile_target.axial_distance_bowl(valid_idx);
        profile_target.axial_intensity = profile_target.axial_intensity(valid_idx);
    else
        % Case 2: Pad with NaN to simulation length
        pad_length = round((sim_axis_mm - profile_target.axial_distance_bowl(end)) / parameters.grid.resolution_mm);
        if pad_length > 0
            pad_dist = profile_target.axial_distance_bowl(end) + ...
                       (1:pad_length) * parameters.grid.resolution_mm;
            profile_target.axial_distance_bowl = [profile_target.axial_distance_bowl, pad_dist];
            profile_target.axial_intensity = [profile_target.axial_intensity, NaN(1, pad_length)];
        end
    end    

end