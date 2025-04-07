function error = phase_optimization_annulus(phase, parameters, velocity, axial_position, desired_focal_dist_mm)

% PHASE_OPTIMIZATION_ANNULUS Computes the focal distance error for a focused annular transducer.
%
% This function calculates the acoustic pressure profile along the axial direction 
% for a focused annular transducer using O'Neil's model. It identifies the focal 
% distance based on the full-width half-maximum (FWHM) of the pressure profile and 
% computes the error relative to a desired focal distance.
%
% Input:
%   phase                 - Array specifying the phase shifts for each transducer element.
%   parameters            - Struct containing transducer and medium properties.
%   velocity              - Array specifying particle velocities for each transducer element.
%   axial_position        - Array specifying axial positions (0 = transducer surface).
%   desired_focal_dist_mm - Scalar specifying the desired focal distance in mm.
%
% Output:
%   error                 - Scalar value representing the absolute error between 
%                           the actual and desired focal distances.

    %% Compute acoustic pressure profile using O'Neil's model
    % Focused annulus model based on transducer geometry and medium properties
    p_axial_oneil = focusedAnnulusONeil(parameters.transducer.curv_radius_mm / 1e3, ...
        [parameters.transducer.Elements_ID_mm; parameters.transducer.Elements_OD_mm] / 1e3, ...
        repmat(velocity, 1, parameters.transducer.n_elements), ...
        [0 phase], parameters.transducer.source_freq_hz, ...
        parameters.medium.water.sound_speed, ...
        parameters.medium.water.density, ...
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