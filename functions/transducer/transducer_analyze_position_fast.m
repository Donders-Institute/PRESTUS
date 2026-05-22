function [prop_intersect, mean_dist_skin, var_dist_skin, mean_dist_skull, var_dist_skull] = transducer_analyze_position_fast(i, norm_v, ex_plane_pos, coord_mesh_gpu, full_skull_mask_idx, skin_boundary_coords, skull_boundary_coords, max_od)

% TRANSDUCER_ANALYZE_POSITION_FAST  GPU-accelerated transducer position analysis
%
% Computes the fraction of exit-plane area intersecting the skull mask and
% mean/variance distances to skin and skull boundaries. Designed to be
% called via arrayfun for parallel evaluation of all candidates.
%
% Use as:
%   [prop_intersect, mean_dist_skin, var_dist_skin, mean_dist_skull, var_dist_skull] = ...
%       transducer_analyze_position_fast(i, norm_v, ex_plane_pos, coord_mesh_gpu, ...
%           full_skull_mask_idx, skin_boundary_coords, skull_boundary_coords, max_od)
%
% Input:
%   i                     - candidate index into norm_v and ex_plane_pos rows
%   norm_v                - [Nx3] unit normal vectors for each candidate (gpuArray)
%   ex_plane_pos          - [Nx3] exit plane centre positions (gpuArray)
%   coord_mesh_gpu        - struct with x, y, z grid coordinates (gpuArray)
%   full_skull_mask_idx   - linear indices of skull mask voxels (gpuArray)
%   skin_boundary_coords  - [Mx3] skin boundary coordinates (gpuArray)
%   skull_boundary_coords - [Px3] skull boundary coordinates (gpuArray)
%   max_od                - aperture diameter in voxels
%
% Output:
%   prop_intersect   - fraction of exit-plane area intersecting skull mask
%   mean_dist_skin   - mean distance to skin boundary [voxels]
%   var_dist_skin    - variance of distances to skin boundary
%   mean_dist_skull  - mean distance to skull boundary [voxels]
%   var_dist_skull   - variance of distances to skull boundary
%
% See also: TRANSDUCER_ANALYZE_POSITION, TP_EVALUATE_CANDIDATE_POSITIONS

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
