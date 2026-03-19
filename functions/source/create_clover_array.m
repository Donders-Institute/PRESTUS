function elem_pos_m = create_clover_array(parameters, matrix_tp, elem_pos_m, trans_pos_m)    
% CREATE_CLOVER_ARRAY Generate a multi-leaf clover transducer configuration.
%
% This function replicates a single matrix transducer layout into a clover
% configuration (multi-aperture arrangement), where multiple identical
% sub-arrays are distributed on a spherical surface (parent bowl).
%
% Each sub-array is:
%   1) Positioned on a sphere with radius ROC_parent
%   2) Rotated in azimuth (Z-axis)
%   3) Tilted toward the natural focus of the multi-leaf configuration
%      (Y-axis rotation)
%
% Additionally, a sphere fit is performed per leaf to determine its apex
% and validate its distance to the intended focus.
%
% INPUTS:
%   parameters   - global parameter struct (used for debug output)
%   matrix_tp    - struct containing matrix + clover configuration
%   elem_pos_m   - [3 x N] element positions (single sub-array) [m]
%   trans_pos_m  - [3 x 1] transducer origin in simulation grid [m]
%
% OUTPUTS:
%   elem_pos_m   - [3 x N_total] full clover element positions [m]

    % --------------------------------------------------------------------
    % Clover geometry parameters
    % --------------------------------------------------------------------
    n_leaves = matrix_tp.clover.n_leaves;
    ROC_parent_mm = matrix_tp.clover.ROC_parent;

    theta_az = 2*pi / 3;  % 120° spacing (fixed geometry)

    % Estimate elevation angle to ensure sub-apertures do not overlap on
    % the parent sphere
    aperture_diam = matrix_tp.Elements_OD_mm(end);
    radius_circle = aperture_diam / (2 * ROC_parent_mm * sin(theta_az / 2));
    elevation_angle = asin(radius_circle);

    % --------------------------------------------------------------------
    % Extract base geometry (convert to mm for geometric operations)
    % --------------------------------------------------------------------
    base_pos_mm = elem_pos_m' * 1e3;   % [N x 3]

    % Parent bowl center (global reference)
    parent_center_m = trans_pos_m;
    parent_center_m(3) = trans_pos_m(3) + ROC_parent_mm * 1e-3;
    parent_center_mm = parent_center_m * 1e3;

    % --------------------------------------------------------------------
    % Define reference (first leaf center)
    % Place first transducer at desired elevation on ROC sphere 
    % (along X-axis azimuthally)
    % --------------------------------------------------------------------
    center0 = parent_center_mm + [ ...
        ROC_parent * cos(elevation_angle);
        0;
        ROC_parent * sin(elevation_angle)
    ];

    % Express base positions relative to first center
    positions_local = base_pos_mm - center0';

    % --------------------------------------------------------------------
    % Allocate storage
    % --------------------------------------------------------------------
    elem_all = [];
    h_leaves = gobjects(1, matrix_tp.clover.n_leaves);  % handles for each leaf

    h = figure;
    hold on; 
    axis equal;
    colors = lines(n_leaves);
    legend_entries = strings(1, n_leaves);

    % --------------------------------------------------------------------
    % Generate each clover leaf
    % --------------------------------------------------------------------
    for i = 0:n_leaves-1

        % --- Rotation matrices --- 
        angle_z = i * theta_az;

        % Rotate around Z to spread evenly
        Rz = [cos(angle_z), -sin(angle_z), 0;
              sin(angle_z),  cos(angle_z), 0;
              0,             0,            1];

        % Rotate around Y to tilt downward (toward focus)
        Ry = [cos(elevation_angle), 0, sin(elevation_angle);
              0,                    1, 0;
             -sin(elevation_angle), 0, cos(elevation_angle)];

        R = Rz * Ry;

        % --- Rotate center ---
        center_rot = (R * (center0 - parent_center_mm)) + parent_center_mm;

        % --- Rotate elements ---
        elems_rot = (R * positions_local')' + center_rot';

        % --- Fit sphere to find apex (validation) ---
        ROC_leaf = matrix_tp.curved.curv_radius_mm;

        residuals = @(c) vecnorm(elems_rot - c, 2, 2) - ROC_leaf;

        % Least squares optimization
        center_fit = lsqnonlin(residuals, mean(elems_rot)); 

        % Compute apex (bowl center point)
        dir_vec = mean(elems_rot) - center_fit;
        dir_vec = dir_vec / norm(dir_vec);
        apex = center_fit + ROC_leaf * dir_vec;

        % Distance to parent center (sanity check)
        dist_focus = norm(parent_center_mm' - apex);

        legend_entries(i+1) = sprintf('Leaf %d (dist: %.2f mm)', i+1, dist_focus);

        scatter3(elems_rot(:,1), elems_rot(:,2), elems_rot(:,3), ...
            15, colors(i+1,:), 'filled');

        h_leaves(i) = plot3([apex(1), parent_center_mm(1)], ...
            [apex(2), parent_center_mm(2)], ...
            [apex(3), parent_center_mm(3)], 'k--');

        % --- Store ---
        elem_all = [elem_all; elems_rot];
    end
    
    % --------------------------------------------------------------------
    % Finalize outputs
    % --------------------------------------------------------------------
    elem_pos_m = elem_all' / 1e3;  % back to meters

    % --------------------------------------------------------------------
    % Debug plot
    % --------------------------------------------------------------------
    h_parent = scatter3(parent_center_mm(1), parent_center_mm(2), ...
        parent_center_mm(3), 100, 'r', 'filled');
    h_focus  = scatter3(focus_pos_m(1)*1e3, focus_pos_m(2)*1e3, focus_pos_m(3)*1e3, 100, 'b', 'filled');
    h_center = scatter3(trans_pos_m(1)*1e3, trans_pos_m(2)*1e3, trans_pos_m(3)*1e3, 100, 'g', 'filled');

    % Combine handles for legend
    legend_handles = [h_leaves, h_parent, h_focus, h_center];

    % Labels
    par_bowl_label = "Middle of parent bowl, ROC " + sprintf('%.2f', ...
        ROC_parent) + " mm";
    legend_labels = [legend_entries, ...
       par_bowl_label, "Focus", "Transducer center"];

    legend(legend_handles, legend_labels);

    xlabel('X [mm]');
    ylabel('Y [mm]');
    zlabel('Z [mm]');
    title(sprintf('Clover Array (%d leaves) ROC sub-arrray %.1f mm', n_leaves, ROC_leaf));

    grid on;
    view([20 25 30]);

    output_file = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_%s_clover%s.png', ...
        parameters.subject_id, parameters.simulation_medium, ...
        parameters.results_filename_affix));

    saveas(h, output_file);
    close(h);

end