function mesh = tp_candidate_mesh(img, target, parameters, pixel_size)
% TP_CANDIDATE_MESH  Build coordinate mesh and candidate transducer geometry
%
% Constructs the full 3-D voxel coordinate grid, skin and skull boundary
% point clouds, and candidate outer-surface positions within tp_dist_close
% of the target. Computes transducer axis, geometric focus, and exit plane
% for all candidate positions.
%
% Use as:
%   mesh = tp_candidate_mesh(img, target, parameters, pixel_size)
%
% Input:
%   img        - [Nx x Ny x Nz] segmented head volume
%   target     - [1x3] target position in voxel space
%   parameters - (1,1) simulation parameters struct
%   pixel_size - voxel size [mm]
%
% Output:
%   mesh - struct with fields: coord_mesh, skin_coords, skull_coords,
%          outer_boundary, close_enough_idx, trans_pos, norm_v,
%          geom_focus, ex_plane, max_od_grid, all_masks_idx
%
% See also: TP_EVALUATE_CANDIDATE_POSITIONS, TP_SELECT_HEURISTIC_POSITION

    %--- 3D voxel grid (CPU + GPU) ---------------------------------------
    sz = size(img);
    [t1_x, t1_y, t1_z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));

    coord_mesh.x   = gpuArray(t1_x);
    coord_mesh.y   = gpuArray(t1_y);
    coord_mesh.z   = gpuArray(t1_z);
    coord_mesh.xyz = gpuArray([reshape(t1_x,[],1), ...
                               reshape(t1_y,[],1), ...
                               reshape(t1_z,[],1)]);  % [Nx*Ny*Nz x 3]

    %--- Tissue masks and boundaries -------------------------------------
    all_masks     = img > 0;  % Brain+skull+skin
    outer_boundary = imdilate(all_masks, strel('sphere',1)) - all_masks;  % Outer surface

    skin_boundary = all_masks - imerode(all_masks, strel('sphere',1));
    skin_coords   = gpuArray(coord_mesh.xyz(find(skin_boundary), :));

    img_cp = img;
    img_cp(img_cp == 7 | img_cp == 8) = 4;  % Merge trabecular/CSF → compact bone
    skull      = img_cp == 4;
    skull_fill = imfill(img_cp > 0 & img_cp <= 4, 'holes');
    skull_fill = imerode(skull_fill, strel('sphere',1));
    skull_boundary = skull & ~skull_fill;
    skull_coords   = gpuArray(coord_mesh.xyz(find(skull_boundary), :));

    %--- Candidate positions near target ---------------------------------
    outer_idx        = find(outer_boundary);
    outer_coords     = coord_mesh.xyz(outer_idx, :);
    coord_rel_to_targ = outer_coords - target;
    distances_to_target = sqrt(sum(coord_rel_to_targ.^2, 2));

    % Default dist_close to focal distance + 20 mm if not set
    if ~isfield(parameters.placement.heuristic, 'dist_close') || isempty(parameters.placement.heuristic.dist_close)
        parameters.placement.heuristic.dist_close = parameters.transducer(1).focal_distance_bowl + 20;
    end
    close_enough_idx = outer_idx(distances_to_target < (parameters.placement.heuristic.dist_close / pixel_size));
    trans_pos_coords = coord_mesh.xyz(close_enough_idx, :);  % [N_cand x 3]

    %--- Transducer axis + geometry --------------------------------------
    norm_v = gpuArray((trans_pos_coords - target) ./ ...
             repmat(sqrt(sum((trans_pos_coords - target).^2, 2)), [1, 3]));

    max_od_mm      = max(parameters.transducer(1).(parameters.transducer(1).type).elem_od_mm);
    dist_gf_to_ep_mm = 0.5 * sqrt(4*parameters.transducer(1).(parameters.transducer(1).type).curv_radius_mm^2 - max_od_mm^2);
    dist_tp_to_ep_mm = parameters.transducer(1).(parameters.transducer(1).type).curv_radius_mm - dist_gf_to_ep_mm;

    pos_shift_mm   = 5 + dist_tp_to_ep_mm;
    shifted_trans_pos_coords = trans_pos_coords + norm_v * (pos_shift_mm / pixel_size);

    geom_focus_pos_all = shifted_trans_pos_coords - norm_v * (parameters.transducer(1).(parameters.transducer(1).type).curv_radius_mm / pixel_size);
    ex_plane_pos_all   = geom_focus_pos_all + norm_v * (dist_gf_to_ep_mm / pixel_size);

    all_masks_indx     = find(img > 0);
    max_od_grid        = max_od_mm / pixel_size;

    %--- Pack outputs ----------------------------------------------------
    mesh.coord_mesh       = coord_mesh;
    mesh.skin_coords      = skin_coords;
    mesh.skull_coords     = skull_coords;
    mesh.outer_boundary   = outer_boundary;
    mesh.close_enough_idx = close_enough_idx;
    mesh.trans_pos        = shifted_trans_pos_coords;
    mesh.norm_v           = norm_v;
    mesh.geom_focus       = geom_focus_pos_all;
    mesh.ex_plane         = ex_plane_pos_all;
    mesh.max_od_grid      = max_od_grid;
    mesh.all_masks_idx    = all_masks_indx;
end