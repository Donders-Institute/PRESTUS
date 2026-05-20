function phases = compute_phases(SOUND_SPEED_WATER, tran, desired_focal_distance_ep, tran_ini_data)
% COMPUTE_PHASES  Compute element steering phases for a target focal depth
%
% Calculates the per-element phase delays required to focus an annular
% transducer at a specified depth by computing geometric path lengths from
% element positions to the target point.
%
% Use as:
%   phases = compute_phases(SOUND_SPEED_WATER, tran, desired_focal_distance_ep, tran_ini_data)
%
% Input:
%   SOUND_SPEED_WATER     - sound speed in water [m/s]
%   tran                  - struct with transducer parameters and element data
%   desired_focal_distance_ep - focal depth relative to the exit plane [mm]
%   tran_ini_data         - ini data struct with element positions under tran_ini_data.elements
%
% Output:
%   phases - [1xN] computed phase values per element [°]
%
% See also: PHASE_OPTIMIZATION_ANNULUS, CALIBRATION_TRANSDUCER

arguments
    SOUND_SPEED_WATER     (1,1) {mustBeNumeric}
    tran                  (1,1) struct
    desired_focal_distance_ep (1,1) {mustBeNumeric}
    tran_ini_data         (1,1) struct
end

    % Calculate wavelength based on sound speed and source frequency
    wavelen = SOUND_SPEED_WATER / tran.transducer.freq_hz; % [m]
    
    % Initialize phases array for all transducer elements
    phases = zeros(1, tran.n_elem);
    
    % Convert focal distance relative to the exit plane to relative to the mid-bowl
    dist_to_mid_bowl = tran.transducer.annular.curv_radius_mm - ...
        tran.transducer.annular.dist_geom_ep_mm;
    focus_wrt_mid_bowl = desired_focal_distance_ep + dist_to_mid_bowl;

    % Compute the target point in relation to the natural focus [mm]
    % Positive z = beyond natural focus; negative z = between bowl and natural focus.
    % Bowl apex sits at z = -curv_radius in this frame, so a point at
    % focus_wrt_mid_bowl mm from the apex is at z = focus_wrt_mid_bowl - curv_radius.
    aim_wrt_natural_focus = focus_wrt_mid_bowl - tran.transducer.annular.curv_radius_mm;
    point_mm = [0, 0, aim_wrt_natural_focus]; % Target point coordinates [x, y, z] in mm

    % Target point in meters
    x = point_mm(1) / 1000.0;
    y = point_mm(2) / 1000.0;
    z = point_mm(3) / 1000.0;
    
    % Compute phase for each transducer element
    for i = 1:tran.n_elem
        % Retrieve element position as a string from tran_ini_data
        field_name = sprintf('x%d', i); 
        str_elem = tran_ini_data.elements.(field_name);
        
        % Parse the string and convert element position to numeric array [m]
        elem_cell = strsplit(str_elem, '|');
        elem = str2double(elem_cell) / 1000; % Convert from mm to meters

        % Calculate distance from the element to the target point
        dist = sqrt((elem(1) - x)^2 + (elem(2) - y)^2 + (elem(3) - z)^2);
        
        % Compute the fractional wavelength component of the distance
        rem = mod(dist / wavelen, 1);
        
        % Calculate the corresponding phase in degrees
        phases(i) = rem * 360.0;
    end

end