function [dist_to_target, dist_to_target_mm, prop_intersect, mean_dts, var_dts] = transducer_analyze_position(trans_pos, pos_shift_mm, target, pixel_size, parameters, coord_mesh, full_skull_mask, skin_boundary_coords)

% TRANSDUCER_ANALYZE_POSITION  Compute geometric and statistical measures for a transducer position
%
% Analyses a transducer position relative to a target, computing distance
% to target, proportion of exit-plane intersection with the skull mask,
% and mean/variance of distances from non-intersecting voxels to the skin
% boundary.
%
% Use as:
%   [dist_to_target, dist_to_target_mm, prop_intersect, mean_dts, var_dts] = ...
%       transducer_analyze_position(trans_pos, pos_shift_mm, target, pixel_size, ...
%           parameters, coord_mesh, full_skull_mask, skin_boundary_coords)
%
% Input:
%   trans_pos            - [1x3] transducer position in grid coordinates
%   pos_shift_mm         - positional shift applied to the transducer [mm]
%   target               - [1x3] target position in grid coordinates
%   pixel_size           - voxel size [mm]
%   parameters           - (1,1) simulation parameters struct
%   coord_mesh           - struct with x, y, z grid coordinate fields
%   full_skull_mask      - [Nx x Ny x Nz] binary skull mask
%   skin_boundary_coords - [Mx3] skin boundary coordinates
%
% Output:
%   dist_to_target    - transducer-to-target distance [grid units]
%   dist_to_target_mm - transducer-to-target distance [mm]
%   prop_intersect    - fraction of exit-plane area intersecting skull mask
%   mean_dts          - mean distance from non-intersecting voxels to skin boundary
%   var_dts           - variance of distances to skin boundary
%
% See also: TRANSDUCER_ANALYZE_POSITION_FAST, TP_EVALUATE_CANDIDATE_POSITIONS

    %% Step 1: Compute maximum outer diameter and exit plane distance
    % Calculate maximum outer diameter of transducer elements
    max_od = max(parameters.transducer(1).(parameters.transducer(1).type).elem_od_mm);

    % Compute distance from geometric focus to exit plane in mm
    dist_to_ep_mm = 0.5 * sqrt(4 * parameters.transducer(1).(parameters.transducer(1).type).curv_radius_mm^2 - max_od^2);

    %% Step 2: Compute normal vector and adjust transducer position
    % Calculate normal vector pointing from transducer to target
    norm_v = (trans_pos - target) / norm(target - trans_pos);

    % Adjust transducer position based on positional shift (converted to grid units)
    trans_pos = trans_pos + norm_v * pos_shift_mm / pixel_size;

    %% Step 3: Compute geometric focus and exit plane positions
    % Compute geometric focus position based on curvature radius
    geom_focus_pos = trans_pos - norm_v * (parameters.transducer(1).(parameters.transducer(1).type).curv_radius_mm) / pixel_size;

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
