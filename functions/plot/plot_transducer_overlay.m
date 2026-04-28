function plot_transducer_overlay(parameters, trans_pos, focus_pos, natural_foc, max_data_pos, options, overlay_weight, overlay_color)
% PLOT_TRANSDUCER_OVERLAY  Draw transducer geometry with labelled rectangles on the current axes
%
% Draws the transducer bowl, exit plane, and focus position on the current
% axes using the computed geometry, with optional weighted colour blending.
%
% Use as:
%   plot_transducer_overlay(parameters, trans_pos, focus_pos, natural_foc, ...
%       max_data_pos, options, overlay_weight, overlay_color)
%
% Input:
%   parameters     - (1,1) simulation parameters struct with transducer geometry
%   trans_pos      - [2x1] transducer centre position [grid units]
%   focus_pos      - [2x1] acoustic focus position [grid units]
%   natural_foc    - [2x1] natural focus position [grid units]
%   max_data_pos   - [2x1] maximum data position for visualisation context
%   options        - struct with .grid_step, .show_rectangles, .rect_size
%   overlay_weight - blend weight in [0, 1] for overlay_color
%   overlay_color  - [1x3] RGB colour for blending
%
% See also: PLOT_OVERLAY_2D, GET_TRANSDUCER_BOX

if isempty(trans_pos)
    return;
end

tr = parameters.transducer(1);

%% 1. GEOMETRY CALCULATIONS (curved transducer details)
%% Compute focal slope
% Calculate unit vector pointing from focus to transducer
if ~isfield(tr, 'align_to_focus') || tr.align_to_focus
    focal_slope = (trans_pos(:) - focus_pos(:)) / norm(trans_pos(:) - focus_pos(:));

else
    focal_slope = (trans_pos(:) - natural_foc(:)) / norm(trans_pos(:) - natural_foc(:));
end

focal_angle = atan2(focal_slope(2), focal_slope(1));

grid_step = options.grid_step;
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

r = max_od / 2 / grid_step;

% Exit plane positions (two methods for consistency)
dist_to_exit_plane_mm = tr.(tr.type).curv_radius_mm - tr.(tr.type).dist_geom_ep_mm;
dist_to_exit_plane_vox = round(dist_to_exit_plane_mm / parameters.grid.resolution_mm);
ex_plane_x_simple = trans_pos(2) - 1 + dist_to_exit_plane_vox;  % Rectangle method

% Full geometry (for curved visualization)
geom_focus_pos = trans_pos(:)' - (tr.(tr.type).curv_radius_mm / grid_step) * [cos(focal_angle), sin(focal_angle)];
dist_to_ep = 0.5 * sqrt(4*tr.(tr.type).curv_radius_mm^2 - max_od^2) / grid_step;
ex_plane_pos_trig = geom_focus_pos + dist_to_ep * [cos(focal_angle), sin(focal_angle)];
ort_angle = atan2(-focal_slope(1), focal_slope(2));

%% 2. VISUAL PROPERTIES
hold on;
lineWidth = 1.5;
boxColor = [235, 185, 47]/255 * (1-overlay_weight) + overlay_color * overlay_weight;
LineSmoothing = 'on';
rect_size = options.rect_size;

%% 3. RECTANGLE OVERLAYS (if enabled)
if options.show_rectangles
    % TRANSDUCER CENTER
    rectangle('Position', [trans_pos(2)-1-rect_size/2, trans_pos(1)-rect_size/2, ...
                          rect_size*2+1, rect_size*2+1], ...
              'EdgeColor', boxColor, 'LineWidth', 2, 'LineStyle', '-');
    text(trans_pos(2)-1, trans_pos(1)+10, ...
         num2str(round((trans_pos(2)-1)*parameters.grid.resolution_mm)), ...
         'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
    
    % GRID ONSET MARKERS
    try
        rectangle('Position', [trans_pos(2)-tr.trans_pos(3)-rect_size/2, ...
                              trans_pos(1)-rect_size/2, rect_size*2+1, rect_size*2+1], ...
                  'EdgeColor', 'w', 'LineWidth', 1, 'LineStyle', ':');
        rectangle('Position', [parameters.grid.dims(3)-rect_size/2, trans_pos(1)-rect_size/2, ...
                              rect_size*2+1, rect_size*2+1], ...
                  'EdgeColor', 'w', 'LineWidth', 1, 'LineStyle', ':');
    catch ME
        disp(ME);
    end
    
    % EXIT PLANE (simple method matching original)
    rectangle('Position', [ex_plane_x_simple-rect_size/2, trans_pos(1)-rect_size/2, ...
                          rect_size*2+1, rect_size*2+1], ...
              'EdgeColor', boxColor, 'LineWidth', 1, 'LineStyle', ':');
    text(ex_plane_x_simple, trans_pos(1)+10, ...
         sprintf('%.0f mm', round(ex_plane_x_simple*parameters.grid.resolution_mm)), ...
         'Color', 'w', 'FontSize', 10, 'FontWeight', 'bold');
end

%% 4. FOCUS AND DATA POSITION
if options.show_rectangles
    rectangle('Position', [focus_pos(2)-rect_size/2, focus_pos(1)-rect_size/2, ...
                          rect_size*2+1, rect_size*2+1], ...
              'EdgeColor', 'r', 'LineWidth', 1, 'LineStyle', '-');
    rectangle('Position', [max_data_pos(2)-rect_size/2, max_data_pos(1)-rect_size/2, ...
                          rect_size*2+1, rect_size*2+1], ...
              'EdgeColor', 'b', 'LineWidth', 1, 'LineStyle', '-');
    text(max_data_pos(2), max_data_pos(1)+10, ...
         sprintf('%.0f mm', round(max_data_pos(2)*parameters.grid.resolution_mm)), ...
         'Color', 'w', 'FontSize', 10);
end

%% 5. DETAILED TRANSDUCER GEOMETRY (always shown)
trans_full_depth = tr.(tr.type).depth_mm / grid_step;
trans_back = ex_plane_pos_trig + trans_full_depth * focal_slope(:)';

% BACK PLANE, SIDE WALLS, EXIT CHORD (as before)
line([trans_back(2)-r*sin(ort_angle), trans_back(2)+r*sin(ort_angle)], ...
     [trans_back(1)-r*cos(ort_angle), trans_back(1)+r*cos(ort_angle)], ...
     'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);

line([ex_plane_pos_trig(2), trans_back(2)] - r*sin(ort_angle), ...
     [ex_plane_pos_trig(1), trans_back(1)] - r*cos(ort_angle), ...
     'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);

line([ex_plane_pos_trig(2), trans_back(2)] + r*sin(ort_angle), ...
     [ex_plane_pos_trig(1), trans_back(1)] + r*cos(ort_angle), ...
     'LineWidth', lineWidth, 'Color', boxColor, 'LineSmoothing', LineSmoothing);

% Exit plane chord
left_edge_x  = ex_plane_pos_trig(2) - r * sin(ort_angle);
left_edge_y  = ex_plane_pos_trig(1) - r * cos(ort_angle);
right_edge_x = ex_plane_pos_trig(2) + r * sin(ort_angle);
right_edge_y = ex_plane_pos_trig(1) + r * cos(ort_angle);
line([left_edge_x, right_edge_x], [left_edge_y, right_edge_y], ...
     'LineWidth', lineWidth, 'Color', boxColor, 'LineStyle', ':', 'LineSmoothing', LineSmoothing);

% Curved front surface
arc_halfangle = atan((max_od/2) / (dist_to_ep * grid_step));
[arc_x, arc_y] = get_arc(geom_focus_pos, tr.(tr.type).curv_radius_mm/grid_step, ...
                        focal_angle-arc_halfangle, focal_angle+arc_halfangle);
plot(arc_y, arc_x, 'Color', boxColor, 'LineWidth', lineWidth, 'LineSmoothing', LineSmoothing);

hold off;
end