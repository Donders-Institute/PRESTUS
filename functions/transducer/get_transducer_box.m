function [transducer_box, ex_plane_pos_trig, geom_focus_pos, dist_to_ep_mm] = ...
    get_transducer_box(trans_pos, focus_pos, natural_foc, grid_step, parameters, is_plot)

% GET_TRANSDUCER_BOX  Compute transducer bounding box, exit plane, and geometric focus
%
% Calculates the 2-D bounding box of a transducer from its grid position,
% focus position, and transducer parameters. Also returns the geometric
% focus position, exit plane position, and distance to the exit plane.
%
% Use as:
%   [transducer_box, ex_plane_pos_trig, geom_focus_pos, dist_to_ep_mm] = ...
%       get_transducer_box(trans_pos, focus_pos, natural_foc, grid_step, parameters)
%   [transducer_box, ex_plane_pos_trig, geom_focus_pos, dist_to_ep_mm] = ...
%       get_transducer_box(trans_pos, focus_pos, natural_foc, grid_step, parameters, is_plot)
%
% Input:
%   trans_pos   - [1x2] transducer position in grid coordinates
%   focus_pos   - [1x2] focus position in grid coordinates
%   natural_foc - [1x2] natural focus position in grid coordinates
%   grid_step   - grid step size [mm]
%   parameters  - (1,1) simulation parameters struct with transducer properties
%   is_plot     - enable/disable bounding box visualisation (default: 1)
%
% Output:
%   transducer_box    - [4x2] bounding box corner coordinates
%   ex_plane_pos_trig - [1x2] exit plane position in grid coordinates
%   geom_focus_pos    - [1x2] geometric focus position in grid coordinates
%   dist_to_ep_mm     - distance to exit plane [mm]
%
% See also: TRANSDUCER_SETUP, GET_ARC

    arguments
        trans_pos   (1,2) double
        focus_pos   (1,2) double
        natural_foc (:,:) double
        grid_step   (1,1) double
        parameters  (1,1) struct
        is_plot     (1,1) double = 1
    end

    %% Compute focal slope0
    % Calculate unit vector pointing from focus to transducer
    if isfield(parameters.transducer(1), 'align_to_focus') && parameters.transducer(1).align_to_focus == 0 && ~isempty(natural_foc)
        focal_slope = (trans_pos - natural_foc) / norm(trans_pos - natural_foc);
    else
        focal_slope = (trans_pos - focus_pos) / norm(trans_pos - focus_pos);
    end
    
    focal_angle = atan2(focal_slope(2),focal_slope(1));

    %% Compute geometric focus position
    % Calculate geometric focus position based on curvature radius and focal angle
    tr = parameters.transducer(1);
    geom_focus_pos = trans_pos - (tr.(tr.type).curv_radius_mm) / grid_step * focal_slope;

    %% Compute distance to exit plane
    % Maximum outer diameter of transducer elements
    switch tr.type
        case 'annular'
            max_od = max(tr.annular.elem_od_mm);
        case 'matrix'
            max_od = tr.matrix.outer_diameter_mm;
        otherwise
            error('Array type %s is unknown or not implemented.', tr.type)
    end

    % Distance from geometric focus to exit plane in mm
    dist_to_ep_mm = 0.5 * sqrt(4 * tr.(tr.type).curv_radius_mm^2 - max_od^2);

    % Convert distance to exit plane from mm to grid units
    dist_to_ep_grid = dist_to_ep_mm / grid_step;

    %% Compute exit plane position
    % Calculate exit plane position based on focal angle and distance to exit plane
    ex_plane_pos_trig = geom_focus_pos + dist_to_ep_grid * focal_slope;

    %% Compute orthogonal angle for bounding box calculation
    % Orthogonal angle perpendicular to focal slope
    ort_angle = atan(-focal_slope(1) / focal_slope(2));

    arc_halfangle = atan(max_od/2/dist_to_ep_grid/grid_step);

    %% Compute bounding box dimensions
    % Radius of bounding box based on maximum outer diameter of elements
    r = max_od / 2 / grid_step;

    % Depth of transducer in grid units
    trans_full_depth = tr.(tr.type).depth_mm / grid_step;

    % Back end position of transducer based on depth and focal slope
    trans_back = ex_plane_pos_trig + trans_full_depth * focal_slope;

    % Coordinates of bounding box corners
    transducer_box = [[trans_back(2) - r * sin(ort_angle), trans_back(1) - r * cos(ort_angle)], ...
                      [trans_back(2) + r * sin(ort_angle), trans_back(1) + r * cos(ort_angle)], ...
                      [ex_plane_pos_trig(2) - r * sin(ort_angle), ex_plane_pos_trig(1) - r * cos(ort_angle)], ...
                      [ex_plane_pos_trig(2) + r * sin(ort_angle), ex_plane_pos_trig(1) + r * cos(ort_angle)]];

    %% Visualization (optional)
    if is_plot 
        hold on;
        overlay_weight = 0; % Weight for overlay color blending
        overlay_color = [0, 0.2, 0.7]; % Overlay color (blue)
        lineWidth = 1; % Line width for visualization

        boxColor = [235, 185, 47] / 255 * (1 - overlay_weight) + overlay_color * overlay_weight; % Box color blending

        % Draw bounding box lines for visualization
        line([trans_back(2) - r * sin(ort_angle), trans_back(2) + r * sin(ort_angle)], ...
             [trans_back(1) - r * cos(ort_angle), trans_back(1) + r * cos(ort_angle)], ...
             'LineWidth', lineWidth, 'Color', boxColor);

        line([ex_plane_pos_trig(2), trans_back(2)] - r * sin(ort_angle), ...
             [ex_plane_pos_trig(1), trans_back(1)] - r * cos(ort_angle), ...
             'LineWidth', lineWidth, 'Color', boxColor);

        line([ex_plane_pos_trig(2), trans_back(2)] + r * sin(ort_angle), ...
             [ex_plane_pos_trig(1), trans_back(1)] + r * cos(ort_angle), ...
             'LineWidth', lineWidth, 'Color', boxColor);

        % exit plane
        line([ex_plane_pos_trig(2), ex_plane_pos_trig(2)]+r*sin(ort_angle), ...
            [trans_back(1)-r*cos(ort_angle), trans_back(1) + r*cos(ort_angle)], ...
            'LineWidth', lineWidth,'Color', boxColor,'LineStyle', ':');

        % transducer curvature
        [arc_x, arc_y] = get_arc(geom_focus_pos, ...
            tr.(tr.type).curv_radius_mm/grid_step, ...
            focal_angle-arc_halfangle, ...
            focal_angle+arc_halfangle );
        plot(arc_y, arc_x, 'Color', boxColor, 'LineWidth', lineWidth)
    end

end
