function [transducer_box, ort_angle, front_center, max_y_projection] = setup_dolphin_transducer_box(totalMask, transducer_pars, trans_pos, focus_pos, grid_step_mm)
% This function computes the transducer box dimensions and positions.
% This function calculates the bounding box of a transducer based on its position 
% (`trans_pos`), focus position (`focus_pos`), and transducer parameters. 
% Input:
%   totalMask         - binary mask representation of Dolphin transducer.
%   transducer_pars   - Struct containing geometric parameters of the transducer:
%                       * n_elements: Number of elements in the transducer.
%                       * Elements_OD_mm: Outer diameter of each element (in mm).
%                       * Elements_ID_mm: Inner diameter of each element (in mm).
%                       * curv_radius_mm: Curvature radius of the elements (in mm).
%   trans_pos         - [1x2] array specifying the transducer position in grid coordinates.
%   focus_pos         - [1x2] array specifying the focus position in grid coordinates.
%   grid_step_mm      - Scalar specifying the step size of the computational grid (in mm).
% Output:
%   transducer_box     - [4x2] matrix specifying the coordinates of the bounding box corners.
%   ort_angle          -
%   front_center       -
%   max_y_projection   -

    % Convert physical dimensions to grid points
    grid_step_m = grid_step_mm / 1000;
    elem_width_grid = round(transducer_pars.elem_width / grid_step_m);  
    elem_height_grid = round(transducer_pars.elem_height / grid_step_m);
    elem_spacing_width_grid = round(transducer_pars.elem_spacing_width / grid_step_m);
    elem_spacing_height_grid = round(transducer_pars.elem_spacing_height / grid_step_m);  
    
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
    trans_full_depth = transducer_pars.depth_mm / grid_step_mm;
    
    % Front face center position
    front_center = trans_pos;
    
    % Back face center position
    back_center = front_center + trans_full_depth * focal_slope;
    
    % Calculate the focal angle
    focal_angle = atan2(focal_slope(2), focal_slope(1));
    
    % The projection of the width onto the viewing plane depends on the focal angle
    % As the transducer rotates, the width component in the viewing plane increases
    width_projection = half_width * abs(sin(focal_angle));
    
    % The projection of the height is always fully visible on the y-axis
    height_projection = half_height * abs(cos(focal_angle));
    
    % For non-zero angles, we need the maximum of both projections for each direction
    max_y_projection = max(width_projection, height_projection);
    
    % Define the bounding box with these projections
    transducer_box = [[back_center(2) - max_y_projection * sin(ort_angle), back_center(1) - max_y_projection * cos(ort_angle)], ...
                      [back_center(2) + max_y_projection * sin(ort_angle), back_center(1) + max_y_projection * cos(ort_angle)], ...
                      [front_center(2) - max_y_projection * sin(ort_angle), front_center(1) - max_y_projection * cos(ort_angle)], ...
                      [front_center(2) + max_y_projection * sin(ort_angle), front_center(1) + max_y_projection * cos(ort_angle)]];
end