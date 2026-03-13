function mesh = tp_candidate_mesh(img, target, parameters, pixel_size)
% TP_CANDIDATE_MESH Build coordinate mesh and candidate transducer geometry
%
% Constructs:
%   - Full 3D voxel coordinate grid (CPU + GPU)
%   - Skin and skull boundary point clouds (GPU)
%   - Candidate outer-surface positions within tp_dist_close of target
%   - Geometry for transducer axis, geometric focus, and exit plane (all candidates)
%
% INPUT
%   img        - Segmented head volume [Nx Ny Nz]
%   target     - 1x3 target position [x y z] in voxel space
%   parameters - Struct with fields:
%                  .tp_dist_close (mm), .transducer.curv_radius_mm,
%                  .transducer.Elements_OD_mm (vector of element ODs in mm)
%   pixel_size - Scalar voxel size in mm
%
% OUTPUT (struct mesh)
%   mesh.coord_mesh      - struct with fields x,y,z,xyz (all gpuArray)
%   mesh.skin_coords     - [N_skin x 3] gpuArray skin boundary points
%   mesh.skull_coords    - [N_skull x 3] gpuArray skull boundary points
%   mesh.outer_boundary  - logical volume (CPU) outer surface mask
%   mesh.close_enough_idx- linear indices of candidate positions (CPU)
%   mesh.trans_pos       - [N_cand x 3] gpuArray shifted transducer positions
%   mesh.norm_v          - [N_cand x 3] gpuArray unit direction vectors
%   mesh.geom_focus      - [N_cand x 3] gpuArray geometric focus positions
%   mesh.ex_plane        - [N_cand x 3] gpuArray exit plane centers
%   mesh.max_od_grid     - scalar, aperture diameter in voxels
%   mesh.all_masks_idx   - gpuArray linear indices of tissue voxels

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

    close_enough_idx = outer_idx(distances_to_target < (parameters.tp_dist_close / pixel_size));
    trans_pos_coords = coord_mesh.xyz(close_enough_idx, :);  % [N_cand x 3]

    %--- Transducer axis + geometry --------------------------------------
    norm_v = gpuArray((trans_pos_coords - target) ./ ...
             repmat(sqrt(sum((trans_pos_coords - target).^2, 2)), [1, 3]));

    max_od_mm      = max(parameters.transducer.Elements_OD_mm);
    dist_gf_to_ep_mm = 0.5 * sqrt(4*parameters.transducer.curv_radius_mm^2 - max_od_mm^2);
    dist_tp_to_ep_mm = parameters.transducer.curv_radius_mm - dist_gf_to_ep_mm;

    pos_shift_mm   = 5 + dist_tp_to_ep_mm;
    shifted_trans_pos_coords = trans_pos_coords + norm_v * (pos_shift_mm / pixel_size);

    geom_focus_pos_all = shifted_trans_pos_coords - norm_v * (parameters.transducer.curv_radius_mm / pixel_size);
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