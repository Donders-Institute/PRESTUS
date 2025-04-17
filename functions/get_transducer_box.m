function [transducer_box, ex_plane_pos_trig, geom_focus_pos, dist_to_ep_mm] = get_transducer_box(trans_pos, focus_pos, grid_step, parameters, plot)

% GET_TRANSDUCER_BOX Computes the transducer box dimensions and positions.
%
% This function calculates the bounding box of a transducer based on its position 
% (`trans_pos`), focus position (`focus_pos`), and transducer parameters. It also 
% computes the geometric focus position, exit plane position, and distance to the 
% exit plane. Optionally, it can visualize the transducer box on a plot.
%
% Input:
%   trans_pos   - [1x2] array specifying the transducer position in grid coordinates.
%   focus_pos   - [1x2] array specifying the focus position in grid coordinates.
%   grid_step   - Scalar specifying the grid step size (in mm).
%   parameters  - Struct containing transducer properties (e.g., curvature radius, element diameters).
%   plot        - Boolean flag to enable/disable visualization of the transducer box (default: 1).
%
% Output:
%   transducer_box     - [4x2] matrix specifying the coordinates of the bounding box corners.
%   ex_plane_pos_trig  - [1x2] array specifying the exit plane position in grid coordinates.
%   geom_focus_pos     - [1x2] array specifying the geometric focus position in grid coordinates.
%   dist_to_ep_mm      - Scalar specifying the distance to the exit plane (in mm).

    arguments
        trans_pos (1, 2) % Transducer position in grid coordinates
        focus_pos (1, 2) % Focus position in grid coordinates
        grid_step (1, 1) % Grid step size in mm
        parameters struct % Struct containing transducer properties
        plot = 1 % Enable/disable visualization (default: enabled)
    end

    if parameters.unique_tran_design == 1

        %% Compute focal slope and angle
        % Calculate unit vector pointing from focus to transducer and its angle
        focal_slope = (trans_pos - focus_pos) / norm(trans_pos - focus_pos);
        focal_angle = atan2(focal_slope(2), focal_slope(1));
    
        %% Compute geometric focus position
        % Calculate geometric focus position based on curvature radius and focal angle
        geom_focus_pos = trans_pos - (parameters.transducer.curv_radius_mm) / grid_step * [cos(focal_angle), sin(focal_angle)];
    
        %% Compute distance to exit plane
        % Maximum outer diameter of transducer elements
        max_od = max(parameters.transducer.Elements_OD_mm);
    
        % Distance from geometric focus to exit plane in mm
        dist_to_ep_mm = 0.5 * sqrt(4 * parameters.transducer.curv_radius_mm^2 - max_od^2);
    
        % Convert distance to exit plane from mm to grid units
        dist_to_ep_grid = dist_to_ep_mm / grid_step;
    
        %% Compute exit plane position
        % Calculate exit plane position based on focal angle and distance to exit plane
        ex_plane_pos_trig = geom_focus_pos + dist_to_ep_grid * [cos(focal_angle), sin(focal_angle)];
    
        %% Compute orthogonal angle for bounding box calculation
        % Orthogonal angle perpendicular to focal slope
        ort_angle = atan(-focal_slope(1) / focal_slope(2));
    
        %% Compute bounding box dimensions
        % Radius of bounding box based on maximum outer diameter of elements
        r = max(parameters.transducer.Elements_OD_mm) / 2 / grid_step;
    
        % Depth of transducer in grid units
        trans_full_depth = 16 / grid_step;
    
        % Back end position of transducer based on depth and focal slope
        trans_back = ex_plane_pos_trig + trans_full_depth * focal_slope;
    
        % Coordinates of bounding box corners
        transducer_box = [[trans_back(2) - r * sin(ort_angle), trans_back(1) - r * cos(ort_angle)], ...
                          [trans_back(2) + r * sin(ort_angle), trans_back(1) + r * cos(ort_angle)], ...
                          [ex_plane_pos_trig(2) - r * sin(ort_angle), ex_plane_pos_trig(1) - r * cos(ort_angle)], ...
                          [ex_plane_pos_trig(2) + r * sin(ort_angle), ex_plane_pos_trig(1) + r * cos(ort_angle)]];

    elseif parameters.unique_tran_design == 2
        % Convert physical dimensions to grid points
        elem_width_grid = round(parameters.transducer.elem_width / parameters.grid_step_m);  
        elem_height_grid = round(parameters.transducer.elem_height / parameters.grid_step_m);
        elem_spacing_width_grid = round(parameters.transducer.elem_spacing_width / parameters.grid_step_m);
        elem_spacing_height_grid = round(parameters.transducer.elem_spacing_height / parameters.grid_step_m);
        
        % Define geometry for the mask
        [X, Y] = meshgrid(1:parameters.transducer.n_elem_col, 1:parameters.transducer.n_elem_row);
        centerX = (parameters.transducer.n_elem_col + 1) / 2;
        centerY = (parameters.transducer.n_elem_row + 1) / 2;

        % Create circular mask
        mask = (X - centerX).^2 + (Y - centerY).^2 <= parameters.transducer.modular_radius^2;
        mask(:,1) = [];
        totalMask = [fliplr(mask), mask];
        
        % Find the boundaries of the active elements
        [rows, cols] = find(totalMask);
        min_row = min(rows);
        max_row = max(rows);
        min_col = min(cols);
        max_col = max(cols);
        
        % Calculate the grid extents of the array
        half_width = (max_col - min_col + 1) * (elem_width_grid + elem_spacing_width_grid) / 2;
        half_height = (max_row - min_row + 1) * (elem_height_grid + elem_spacing_height_grid) / 2;
        
        % Calculate the focal slope similar to design 1
        focal_slope = (trans_pos - focus_pos) / norm(trans_pos - focus_pos);
        
        % Calculate the orthogonal angle
        ort_angle = atan(-focal_slope(1) / focal_slope(2));
        
        % Depth of transducer in grid units
        trans_full_depth = 16 / grid_step;
        
        % Front face center position
        front_center = trans_pos;
        
        % Back face center position
        back_center = front_center + trans_full_depth * focal_slope;
        
        % Calculate corners of the bounding box
        % Front face corners
        front_corner1 = [front_center(1) - half_width * cos(ort_angle), front_center(2) - half_width * sin(ort_angle)];
        front_corner2 = [front_center(1) + half_width * cos(ort_angle), front_center(2) + half_width * sin(ort_angle)];
        front_corner3 = [front_center(1) - half_height * sin(ort_angle), front_center(2) + half_height * cos(ort_angle)];
        front_corner4 = [front_center(1) + half_height * sin(ort_angle), front_center(2) - half_height * cos(ort_angle)];
        
        % Back face corners
        back_corner1 = [back_center(1) - half_width * cos(ort_angle), back_center(2) - half_width * sin(ort_angle)];
        back_corner2 = [back_center(1) + half_width * cos(ort_angle), back_center(2) + half_width * sin(ort_angle)];
        back_corner3 = [back_center(1) - half_height * sin(ort_angle), back_center(2) + half_height * cos(ort_angle)];
        back_corner4 = [back_center(1) + half_height * sin(ort_angle), back_center(2) - half_height * cos(ort_angle)];
        
        back_y_values = [back_corner1(2), back_corner2(2), back_corner3(2), back_corner4(2)];
        back_x_values = [back_corner1(1), back_corner2(1), back_corner3(1), back_corner4(1)];
        front_y_values = [front_corner1(2), front_corner2(2), front_corner3(2), front_corner4(2)];
        front_x_values = [front_corner1(1), front_corner2(1), front_corner3(1), front_corner4(1)];
        
        transducer_box = [
            [min(back_y_values), min(back_x_values)],
            [max(back_y_values), max(back_x_values)],
            [min(front_y_values), min(front_x_values)],
            [max(front_y_values), max(front_x_values)]
        ];
    else
        fprintf("Transducer design number %i hasn't been implemented. \n", parameters.unique_tran_design)
    end


    %% Visualization (optional)
    if plot 
        overlay_weight = 0; % Weight for overlay color blending
        overlay_color = [0, 0.2, 0.7]; % Overlay color (blue)
        lineWidth = 1; % Line width for visualization

        boxColor = [235, 185, 47] / 255 * (1 - overlay_weight) + overlay_color * overlay_weight; % Box color blending
        LineSmoothing = 'on'; % Enable line smoothing

        if parameters.unique_tran_design == 1
            % Draw bounding box lines for visualization
            line([trans_back(2) - r * sin(ort_angle), trans_back(2) + r * sin(ort_angle)], ...
                 [trans_back(1) - r * cos(ort_angle), trans_back(1) + r * cos(ort_angle)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
    
            line([ex_plane_pos_trig(2), trans_back(2)] - r * sin(ort_angle), ...
                 [ex_plane_pos_trig(1), trans_back(1)] - r * cos(ort_angle), ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
    
            line([ex_plane_pos_trig(2), trans_back(2)] + r * sin(ort_angle), ...
                 [ex_plane_pos_trig(1), trans_back(1)] + r * cos(ort_angle), ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
        elseif parameters.unique_tran_design == 2
            % Draw back face
            line([back_corner1(2), back_corner2(2)], [back_corner1(1), back_corner2(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            line([back_corner2(2), back_corner4(2)], [back_corner2(1), back_corner4(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            line([back_corner4(2), back_corner3(2)], [back_corner4(1), back_corner3(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            line([back_corner3(2), back_corner1(2)], [back_corner3(1), back_corner1(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            
            % Draw front face
            line([front_corner1(2), front_corner2(2)], [front_corner1(1), front_corner2(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            line([front_corner2(2), front_corner4(2)], [front_corner2(1), front_corner4(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            line([front_corner4(2), front_corner3(2)], [front_corner4(1), front_corner3(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            line([front_corner3(2), front_corner1(2)], [front_corner3(1), front_corner1(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            
            % Connect front to back
            line([front_corner1(2), back_corner1(2)], [front_corner1(1), back_corner1(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            line([front_corner2(2), back_corner2(2)], [front_corner2(1), back_corner2(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            line([front_corner3(2), back_corner3(2)], [front_corner3(1), back_corner3(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
            line([front_corner4(2), back_corner4(2)], [front_corner4(1), back_corner4(1)], ...
                 'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);
        else
            fprintf("Transducer design number %i hasn't been implemented. \n", parameters.unique_tran_design)
        end
    end

end
