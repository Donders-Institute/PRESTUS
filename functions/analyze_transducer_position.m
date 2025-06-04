function [dist_to_target, dist_to_target_mm, prop_intersect, mean_dts, var_dts] = analyze_transducer_position(trans_pos, pos_shift_mm, target, pixel_size, parameters, coord_mesh, full_skull_mask, skin_boundary_coords)

% ANALYZE_TRANSDUCER_POSITION Computes geometric and statistical measures for a transducer position.
%
% This function analyzes the position of a transducer relative to a target and 
% computes various measures such as the distance to the target, proportion of 
% intersection with the skull mask, and statistical measures (mean and variance) 
% of distances to the skin boundary.
%
% Input:
%   trans_pos           - [1x3] array specifying the transducer position in grid coordinates.
%   pos_shift_mm        - Scalar specifying the positional shift applied to the transducer (in mm).
%   target              - [1x3] array specifying the target position in grid coordinates.
%   pixel_size          - Scalar specifying the voxel size (in mm).
%   parameters          - Struct containing transducer properties (e.g., curvature radius, diameters).
%   coord_mesh          - Struct containing x, y, z coordinates of the computational grid.
%   full_skull_mask     - Binary mask indicating voxels belonging to the skull.
%   skin_boundary_coords- [Mx3] array of coordinates defining the skin boundary.
%
% Output:
%   dist_to_target      - Scalar distance from transducer to target (in grid units).
%   dist_to_target_mm   - Distance from transducer to target (in mm).
%   prop_intersect      - Proportion of intersection between the orthogonal plane and the skull mask.
%   mean_dts            - Mean distance from non-intersecting voxels to the skin boundary.
%   var_dts             - Variance of distances from non-intersecting voxels to the skin boundary.

    %% Step 1: Compute maximum outer diameter and exit plane distance
    % Calculate maximum outer diameter of transducer elements
    max_od = max(parameters.transducer.Elements_OD_mm);

    % Compute distance from geometric focus to exit plane in mm
    dist_to_ep_mm = 0.5 * sqrt(4 * parameters.transducer.curv_radius_mm^2 - max_od^2);

    %% Step 2: Compute normal vector and adjust transducer position
    % Calculate normal vector pointing from transducer to target
    norm_v = (trans_pos - target) / norm(target - trans_pos);

    % Adjust transducer position based on positional shift (converted to grid units)
    trans_pos = trans_pos + norm_v * pos_shift_mm / pixel_size;

    %% Step 3: Compute geometric focus and exit plane positions
    % Compute geometric focus position based on curvature radius
    geom_focus_pos = trans_pos - norm_v * (parameters.transducer.curv_radius_mm) / pixel_size;

    % Compute exit plane position based on distance to exit plane
    ex_plane_pos = geom_focus_pos + norm_v * dist_to_ep_mm / pixel_size;

    %% Step 4: Define orthogonal plane for analysis
    % Compute distance `d` for the plane equation using dot product of normal vector and exit plane position
    d = sum(norm_v .* ex_plane_pos);

    % Identify voxels near this orthogonal plane within a threshold distance (0.5 grid units)
    orth_plane_disk = find(abs(coord_mesh.x * norm_v(1) + ...
                               coord_mesh.y * norm_v(2) + ...
                               coord_mesh.z * norm_v(3) - d) < 0.5);

    % Further restrict to voxels within a disk centered at `ex_plane_pos` with radius `(max_od/2)`
    orth_plane_disk = orth_plane_disk(sqrt((coord_mesh.x(orth_plane_disk) - ex_plane_pos(1)).^2 + ...
                                           (coord_mesh.y(orth_plane_disk) - ex_plane_pos(2)).^2 + ...
                                           (coord_mesh.z(orth_plane_disk) - ex_plane_pos(3)).^2) < max_od / 2);

    %% Step 5: Compute proportion of intersection with skull mask
    % Calculate proportion of voxels in `orth_plane_disk` that intersect with `full_skull_mask`
    prop_intersect = sum(full_skull_mask(orth_plane_disk)) / length(orth_plane_disk);

    %% Step 6: Identify non-intersecting voxels
    % Extract indices of non-intersecting voxels in `orth_plane_disk`
    non_intersecting_vox_idx = intersect(orth_plane_disk, find(~full_skull_mask));

    % Extract coordinates of non-intersecting voxels
    non_intersecting_vox_coords = [coord_mesh.x(non_intersecting_vox_idx), ...
                                   coord_mesh.y(non_intersecting_vox_idx), ...
                                   coord_mesh.z(non_intersecting_vox_idx)];

    %% Step 7: Filter relevant skin boundary coordinates
    % Extract skin boundary coordinates close to `ex_plane_pos` within a threshold distance (40 mm)
    skin_boundary_coords = skin_boundary_coords(pdist2(skin_boundary_coords, ex_plane_pos) < 40,:);

    %% Step 8: Compute distances to skin boundary
    % Calculate pairwise distances between non-intersecting voxels and skin boundary
    dist_to_skin = pdist2(non_intersecting_vox_coords, skin_boundary_coords);

    % Find minimum distance for each voxel (closest point on the skin boundary)
    dist_to_skin = min(dist_to_skin, [], 2);

    %% Step 9: Compute statistical measures for distances
    % Calculate mean and variance for distances to skin boundary
    mean_dts = mean(dist_to_skin(:));
    var_dts = var(dist_to_skin(:));

end
