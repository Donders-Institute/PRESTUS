function [prop_intersect, mean_dist_skin, var_dist_skin, mean_dist_skull, var_dist_skull] = analyze_transducer_position_fast(i, norm_v, ex_plane_pos, coord_mesh_gpu, full_skull_mask_idx, skin_boundary_coords, skull_boundary_coords, max_od)

% ANALYZE_TRANSDUCER_POSITION_FAST Analyzes the transducer position relative to the skull and skin.
%
% This function computes the proportion of intersection between a transducer's 
% orthogonal plane and the skull mask. It also calculates statistical measures 
% (mean and variance) of distances from non-intersecting voxels to the skin and 
% skull boundaries.
%
% Input:
%   i                     - Index of the current transducer element.
%   norm_v                - [Nx3] matrix of normal vectors for each transducer element.
%   ex_plane_pos          - [Nx3] matrix of positions for each transducer element's exit plane.
%   coord_mesh_gpu        - Struct containing x, y, z coordinates of the computational grid on GPU.
%   full_skull_mask_idx   - Indices of voxels belonging to the skull mask.
%   skin_boundary_coords  - [Mx3] array of coordinates defining the skin boundary.
%   skull_boundary_coords - [Px3] array of coordinates defining the skull boundary.
%   max_od                - Maximum outer diameter of the transducer elements (in grid units).
%
% Output:
%   prop_intersect        - Proportion of intersection between the orthogonal plane and the skull mask.
%   mean_dist_skin        - Mean distance from non-intersecting voxels to the skin boundary.
%   var_dist_skin         - Variance of distances from non-intersecting voxels to the skin boundary.
%   mean_dist_skull       - Mean distance from intersecting voxels to the skull boundary.
%   var_dist_skull        - Variance of distances from intersecting voxels to the skull boundary.

    %% Step 1: Define orthogonal plane for the current transducer element
    % Compute distance `d` for the plane equation using dot product of normal vector and position
    d = sum(norm_v(i,:) .* ex_plane_pos(i,:));

    % Identify voxels near this orthogonal plane within a threshold distance (0.5 grid units)
    orth_plane_disk = find(abs(coord_mesh_gpu.x * norm_v(i,1) + ...
                               coord_mesh_gpu.y * norm_v(i,2) + ...
                               coord_mesh_gpu.z * norm_v(i,3) - d) < 0.5);

    % Further restrict to voxels within a disk centered at `ex_plane_pos` with radius `(max_od/2 + 4)`
    orth_plane_disk = orth_plane_disk(sqrt((coord_mesh_gpu.x(orth_plane_disk) - ex_plane_pos(i,1)).^2 + ...
                                           (coord_mesh_gpu.y(orth_plane_disk) - ex_plane_pos(i,2)).^2 + ...
                                           (coord_mesh_gpu.z(orth_plane_disk) - ex_plane_pos(i,3)).^2) < (max_od / 2 + 4));

    %% Step 2: Compute proportion of intersection with skull mask
    % Calculate proportion of voxels in `orth_plane_disk` that intersect with `full_skull_mask_idx`
    prop_intersect = length(intersect(full_skull_mask_idx, orth_plane_disk)) / length(orth_plane_disk);

    %% Step 3: Handle cases with high intersection or empty planes
    % If there are no intersecting or non-intersecting voxels or if `prop_intersect > 0.3`, return early
    non_intersecting_vox_idx = setdiff(orth_plane_disk, full_skull_mask_idx);
    if isempty(orth_plane_disk) || isempty(non_intersecting_vox_idx) || prop_intersect > 0.3
        if isempty(orth_plane_disk) || isempty(non_intersecting_vox_idx)
            prop_intersect = 1; % Set full intersection if no valid voxels exist
        end
        mean_dist_skin = gpuArray(0);
        var_dist_skin = gpuArray(0);
        mean_dist_skull = gpuArray(0);
        var_dist_skull = gpuArray(0);
        return;
    end

    %% Step 4: Filter relevant boundary coordinates
    % Extract coordinates close to `ex_plane_pos` for both skin and skull boundaries
    non_intersecting_vox_coords = coord_mesh_gpu.xyz(non_intersecting_vox_idx,:);
    skin_boundary_coords = skin_boundary_coords(pdist2(skin_boundary_coords, ex_plane_pos(i,:)) < 50,:);
    skull_boundary_coords = skull_boundary_coords(pdist2(skull_boundary_coords, ex_plane_pos(i,:)) < 100,:);

    %% Step 5: Compute distances to skin boundary
    % Calculate pairwise distances between non-intersecting voxels and skin boundary
    dist_to_skin = pdist2(non_intersecting_vox_coords, skin_boundary_coords);
    
    % Find minimum distance for each voxel (closest point on the skin boundary)
    dist_to_skin = min(dist_to_skin, [], 2);

    %% Step 6: Compute distances to skull boundary
    % Calculate pairwise distances between all voxels in `orth_plane_disk` and skull boundary
    dist_to_skull = pdist2(coord_mesh_gpu.xyz(orth_plane_disk,:), skull_boundary_coords);
    
    % Find minimum distance for each voxel (closest point on the skull boundary)
    dist_to_skull = min(dist_to_skull, [], 2);

    %% Step 7: Compute statistical measures for distances
    % Calculate mean and variance for distances to skin boundary
    mean_dist_skin = mean(dist_to_skin(:));
    var_dist_skin = var(dist_to_skin(:));

    % Calculate mean and variance for distances to skull boundary
    mean_dist_skull = mean(dist_to_skull(:));
    var_dist_skull = var(dist_to_skull(:));

end
