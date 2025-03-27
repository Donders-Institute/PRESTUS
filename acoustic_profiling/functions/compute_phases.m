function phases = compute_phases(SOUND_SPEED_WATER, tran, focus_wrt_exit_plane, tran_ini_data)
    % Computes the phases necessary to aim at the specified focal point.
    %
    % Arguments:
    % - SOUND_SPEED_WATER: Speed of sound in water [m/s].
    % - tran: Structure containing transducer parameters and element data.
    % - focus_wrt_exit_plane: Focal depth relative to the transducer exit plane [mm].
    % - tran_ini_data: Ini data structure containing transducer element positions.
    %
    % Returns:
    % - phases: Computed phase values [degrees] for each transducer element.

    % Calculate wavelength based on sound speed and source frequency
    wavelen = SOUND_SPEED_WATER / tran.prestus.transducer.source_freq_hz; % [m]
    
    % Initialize phases array for all transducer elements
    phases = zeros(1, tran.n_elem);
    
    % Convert focal distance relative to the exit plane to relative to the mid-bowl
    dist_to_mid_bowl = tran.prestus.transducer.curv_radius_mm - tran.prestus.transducer.dist_to_plane_mm;
    focus_wrt_mid_bowl = focus_wrt_exit_plane + dist_to_mid_bowl;

    % Compute the target point in relation to the natural focus [mm]
    aim_wrt_natural_focus = tran.prestus.transducer.curv_radius_mm - focus_wrt_mid_bowl;
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