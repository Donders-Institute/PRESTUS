function [profile_target, parameters] = scale_real_intensity_profile(parameters, profile_empirical, desired_intensity)
% SCALE_REAL_INTENSITY_PROFILE  Scale an empirical intensity profile to a desired peak intensity
%
% Linearly scales the empirical profile so its maximum equals desired_intensity,
% updates parameters.transducer.annular.elem_amp from the corresponding pressure,
% and truncates or pads the profile to match the simulation axis length.
%
% Use as:
%   [profile_target, parameters] = scale_real_intensity_profile(parameters, profile_empirical, desired_intensity)
%
% Input:
%   parameters        - PRESTUS config; uses medium_properties.water, grid.default_dims,
%                       grid.resolution_mm [mm], transducer.annular.elem_n
%   profile_empirical - struct with axial_intensity [W/cm²] and axial_distance_bowl [mm]
%   desired_intensity - desired maximum intensity [W/cm²]
%
% Output:
%   profile_target - struct with axial_intensity [W/cm²] and axial_distance_bowl [mm]
%   parameters     - updated with scaled transducer.annular.elem_amp
%
% See also: EXTRACT_REAL_INTENSITY_PROFILE, CALIBRATION_TRANSDUCER

arguments
    parameters        (1,1) struct
    profile_empirical (1,1) struct
    desired_intensity (1,1) {mustBeNumeric}
end

    % Calculate maximum intensity and log adjustment details
    max_intens = max(profile_empirical.axial_intensity);
    fprintf('Current maximum intensity: %.2f \nDesired maximum intensity: %.2f \nAdjust profile to match desired maximum intensity.\n', max_intens, desired_intensity);

    % Calculate the corresponding pressure amplitude for the desired intensity

    DENSITY_WATER = parameters.medium_properties.water.density;
    SOUND_SPEED_WATER = parameters.medium_properties.water.sound_speed;
    p_pa = sqrt(2 * desired_intensity * 1e4 * DENSITY_WATER * SOUND_SPEED_WATER);

    % Update source amplitude in simulation parameters
    parameters.transducer.annular.elem_amp = repmat(p_pa, 1, parameters.transducer.annular.elem_n);

    % Linearly scale the intensity profile
    adjustment_factor_intensity = max_intens / desired_intensity;
    adjusted_axial_intensity = profile_empirical.axial_intensity' ./ adjustment_factor_intensity;

    profile_target.axial_intensity = adjusted_axial_intensity;
    profile_target.axial_distance_bowl = profile_empirical.axial_distance_bowl';

    % Calculate size of the simulation domain.
    % The PML is added outside the grid (PMLInside=false), so the full axial extent
    % is active. The default transducer sits at position 2, leaving default_dims(end)-2
    % voxels of propagation distance.
    sim_axis_mm = (parameters.grid.default_dims(end) - 2) * parameters.grid.resolution_mm;

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