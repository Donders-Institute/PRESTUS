function error = phase_optimization_annulus(phase, parameters, velocity, axial_position, desired_focal_dist_mm)

% PHASE_OPTIMIZATION_ANNULUS  Compute focal-distance error for an annular transducer
%
% Evaluates the O'Neil axial pressure model, locates the focal distance via
% the FWHM centre of the pressure profile, and returns the absolute error
% relative to the desired focal distance. Used as an objective function in
% single-point focal-depth optimisation.
%
% Use as:
%   error = phase_optimization_annulus(phase, parameters, velocity, axial_position, desired_focal_dist_mm)
%
% Input:
%   phase                 - phase shifts per transducer element [rad]
%   parameters            - PRESTUS config with transducer.annular geometry and
%                           medium_properties.water
%   velocity              - particle velocities per element [m/s]
%   axial_position        - axial positions along the beam axis [mm]
%   desired_focal_dist_mm - desired focal distance [mm]
%
% Output:
%   error - absolute error between actual and desired focal distances [mm]
%
% See also: PHASE_OPTIMIZATION_ANNULUS_FULL_CURVE, PERFORM_GLOBAL_SEARCH

arguments
    phase                 (1,:) {mustBeNumeric}
    parameters            (1,1) struct
    velocity              (1,1) {mustBeNumeric}
    axial_position        (1,:) {mustBeNumeric}
    desired_focal_dist_mm (1,1) {mustBeNumeric}
end

    %% Compute acoustic pressure profile using O'Neil's model
    % Focused annulus model based on transducer geometry and medium properties
    p_axial_oneil = focusedAnnulusONeil(parameters.transducer.annular.curv_radius_mm / 1e3, ...
        [parameters.transducer.annular.elem_id_mm; parameters.transducer.annular.elem_od_mm] / 1e3, ...
        repmat(velocity, 1, parameters.transducer.annular.elem_n), ...
        [0 phase], parameters.transducer.freq_hz, ...
        parameters.medium_properties.water.sound_speed, ...
        parameters.medium_properties.water.density, ...
        (axial_position - 0.5) * 1e-3);

    %% Identify focal distance based on FWHM
    % Find maximum amplitude and its position
    [max_amp, max_pos] = max(p_axial_oneil);

    % Compute half-max value
    halfMax = (min(p_axial_oneil) + max_amp) / 2;

    % Find index where pressure first drops below half-max on the left side of the peak
    index1 = find(p_axial_oneil <= halfMax & axial_position < axial_position(max_pos), 1, 'last') + 1;

    % Find index where pressure last rises above half-max on the right side of the peak
    index2 = find(p_axial_oneil <= halfMax & axial_position > axial_position(max_pos), 1, 'first') - 1;

    % Calculate actual focal distance in mm based on FWHM center
    actual_focal_dist_mm = (axial_position(index2) - axial_position(index1)) / 2 + axial_position(index1);

    %% Compute error between actual and desired focal distances
    error = abs(desired_focal_dist_mm - actual_focal_dist_mm); % Absolute error

end