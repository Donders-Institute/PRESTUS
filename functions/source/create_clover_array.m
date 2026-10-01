function [elem_pos_m, aim_pos_m] = create_clover_array(parameters, tr_matrix, elem_pos_m, trans_pos_m, focus_pos_m)
% CREATE_CLOVER_ARRAY  Replicate a single matrix sub-array into a multi-leaf clover configuration
%
% Each of the N leaves (default 3, 360/N degrees apart) is placed on a parent
% spherical bowl of radius ROC_parent, rotated azimuthally around Z, then
% tilted in elevation toward the parent bowl center so that all leaves
% converge on a common natural focus. By default the elevation is the
% smallest tilt at which adjacent leaf centres are one aperture diameter
% apart (no overlap on the parent sphere); a single leaf stays on axis.
% Optional debug figures show the 3-D leaf arrangement.
%
% Use as:
%   [elem_pos_m, aim_pos_m] = create_clover_array(parameters, tr_matrix, elem_pos_m, trans_pos_m, focus_pos_m)
%
% Input:
%   parameters  - PRESTUS config; parameters.simulation.debug controls figure generation
%   tr_matrix   - matrix + clover geometry fields: tr_matrix.clover.n_leaves,
%                 tr_matrix.clover.ROC_parent [mm], tr_matrix.outer_diameter_mm [mm],
%                 tr_matrix.curv_radius_mm [mm], tr_matrix.elem_shape;
%                 optional tr_matrix.clover.az_sep_deg (azimuthal separation
%                 between adjacent leaves [deg], default 360/n_leaves) and
%                 tr_matrix.clover.elevation_deg (leaf tilt from the parent
%                 axis [deg], default derived from the aperture as above)
%   elem_pos_m  - [3xN] element positions of one sub-array [m]
%   trans_pos_m - [3x1] transducer origin in the simulation grid [m]
%   focus_pos_m - [3x1] acoustic focus position in the simulation grid [m]
%
% Output:
%   elem_pos_m  - [3xN*n_leaves] full clover element positions for all leaves [m]
%   aim_pos_m   - [3xN*n_leaves] natural focus of each element's leaf [m]; element
%                 normals point here (equals the parent centre when
%                 curv_radius_mm == ROC_parent)
%
% See also: CREATE_MATRIX_KARRAY, SOURCE_CREATE

arguments
    parameters  (1,1) struct
    tr_matrix   (1,1) struct
    elem_pos_m  (3,:) {mustBeNumeric}
    trans_pos_m (3,1) {mustBeNumeric}
    focus_pos_m (3,1) {mustBeNumeric}
end

    % --------------------------------------------------------------------
    % Clover geometry parameters
    % --------------------------------------------------------------------
    n_leaves = tr_matrix.clover.n_leaves;
    ROC_parent_mm = tr_matrix.clover.ROC_parent;

    if isfield(tr_matrix.clover, 'az_sep_deg') && ~isempty(tr_matrix.clover.az_sep_deg)
        theta_az = deg2rad(tr_matrix.clover.az_sep_deg);
    else
        theta_az = 2*pi / n_leaves;
    end

    if isfield(tr_matrix.clover, 'elevation_deg') && ~isempty(tr_matrix.clover.elevation_deg)
        elevation_angle = deg2rad(tr_matrix.clover.elevation_deg);
    elseif n_leaves == 1
        elevation_angle = 0;
    else
        % Smallest elevation at which adjacent sub-apertures do not overlap
        % on the parent sphere: chord between adjacent leaf centres
        % 2*ROC_parent*sin(el)*sin(theta_az/2) equals the aperture diameter.
        aperture_diam = tr_matrix.outer_diameter_mm;
        radius_circle = aperture_diam / (2 * ROC_parent_mm * sin(theta_az / 2));
        if radius_circle > 1
            error(['create_clover_array: %d leaves of %.1f mm aperture %.1f deg apart do not fit ' ...
                'on a parent bowl of ROC %.1f mm. Increase ROC_parent or az_sep_deg, or set elevation_deg.'], ...
                n_leaves, aperture_diam, rad2deg(theta_az), ROC_parent_mm);
        end
        elevation_angle = asin(radius_circle);
    end

    % --------------------------------------------------------------------
    % Extract base geometry (convert to mm for geometric operations)
    % --------------------------------------------------------------------
    base_pos_mm = elem_pos_m' * 1e3;   % [N x 3]

    % Parent bowl center (global reference)
    parent_center_m = trans_pos_m;
    parent_center_m(3) = trans_pos_m(3) + ROC_parent_mm * 1e-3;
    parent_center_mm = parent_center_m * 1e3;

    % Apex and natural focus of the base sub-array (beam along +Z)
    ROC_leaf = tr_matrix.curv_radius_mm;
    base_apex_mm = trans_pos_m' * 1e3;
    base_nat_focus_mm = base_apex_mm + [0, 0, ROC_leaf];

    % --------------------------------------------------------------------
    % Allocate storage
    % --------------------------------------------------------------------
    n_base = size(base_pos_mm, 1);
    elem_all = zeros(n_base * n_leaves, 3);
    aim_all  = zeros(n_base * n_leaves, 3);
    h_leaves = gobjects(1, n_leaves);  % handles for each leaf


    % [DEBUG] visualize leaf orientation
    if parameters.simulation.debug == 1
        h = figure;
        hold on;
        axis equal;
        colors = lines(n_leaves);
        legend_entries = strings(1, n_leaves);
    end

    % --------------------------------------------------------------------
    % Generate each clover leaf
    % --------------------------------------------------------------------
    for i = 1:n_leaves

        % --- Rotation matrices --- 
        angle_z = (i-1) * theta_az;

        % Rotate around Z to spread evenly
        Rz = [cos(angle_z), -sin(angle_z), 0;
              sin(angle_z),  cos(angle_z), 0;
              0,             0,            1];

        % Rotate around Y to tilt downward (toward focus)
        Ry = [cos(elevation_angle), 0, sin(elevation_angle);
              0,                    1, 0;
             -sin(elevation_angle), 0, cos(elevation_angle)];

        R = Rz * Ry;

        % --- Rotate the leaf rigidly about the parent bowl centre ---
        rot = @(p) (R * (p - parent_center_mm')')' + parent_center_mm';
        elems_rot = rot(base_pos_mm);
        apex      = rot(base_apex_mm);
        nat_focus = rot(base_nat_focus_mm);

        % Offset of the leaf's natural focus from the parent centre (sanity
        % check; zero when curv_radius_mm == ROC_parent)
        dist_focus = norm(parent_center_mm' - nat_focus);

        % [DEBUG] visualize leaf orientation
        if parameters.simulation.debug == 1
            legend_entries(i) = sprintf('Leaf %d (focus offset: %.2f mm)', i, dist_focus);

            scatter3(elems_rot(:,1), elems_rot(:,2), elems_rot(:,3), ...
                15, colors(i,:), 'filled');

            h_leaves(i) = plot3([apex(1), parent_center_mm(1)], ...
                [apex(2), parent_center_mm(2)], ...
                [apex(3), parent_center_mm(3)], 'k--');
        end

        % --- Store ---
        idx = (i-1)*n_base + (1:n_base);
        elem_all(idx, :) = elems_rot;
        aim_all(idx, :)  = repmat(nat_focus, n_base, 1);
    end
    
    % --------------------------------------------------------------------
    % Finalize outputs
    % --------------------------------------------------------------------
    elem_pos_m = elem_all' / 1e3;  % back to meters
    aim_pos_m  = aim_all'  / 1e3;

    % --------------------------------------------------------------------
    % Debug plot
    % --------------------------------------------------------------------
    % [DEBUG] visualize leaf orientation
    if parameters.simulation.debug == 1
        h_parent = scatter3(parent_center_mm(1), parent_center_mm(2), ...
            parent_center_mm(3), 100, 'r', 'filled');
        h_focus  = scatter3(focus_pos_m(1)*1e3, focus_pos_m(2)*1e3, focus_pos_m(3)*1e3, 100, 'b', 'filled');
        h_center = scatter3(trans_pos_m(1)*1e3, trans_pos_m(2)*1e3, trans_pos_m(3)*1e3, 100, 'g', 'filled');

        % Combine handles for legend
        legend_handles = [h_leaves, h_parent, h_focus, h_center];

        % Labels
        par_bowl_label = "Middle of parent bowl, ROC " + sprintf('%.2f', ...
            ROC_parent_mm) + " mm";
        legend_labels = [legend_entries, ...
            par_bowl_label, "Focus", "Transducer center"];

        legend(legend_handles, legend_labels);

        xlabel('X [mm]');
        ylabel('Y [mm]');
        zlabel('Z [mm]');
        title(sprintf('Clover Array (%d leaves) ROC sub-arrray %.1f mm', n_leaves, ROC_leaf));

        grid on;
        view([20 25 30]);

        output_file = fullfile(parameters.io.dir_debug_source, ...
            sprintf('sub-%03d_%s_clover%s.png',...
            parameters.subject_id, parameters.simulation.medium, ...
            parameters.io.output_affix));

        saveas(h, output_file);
        close(h);
    end

end