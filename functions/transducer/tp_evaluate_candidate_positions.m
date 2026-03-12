function tpos_pars = tp_evaluate_candidate_positions(img, target, parameters, pixel_size)
%% TP_EVALUATE_CANDIDATE_POSITIONS Evaluate candidate transducer positions
%  SEARCH CRITERIA:
%  1. Skull exterior surface points (img==0 → imdilate gradient)
%  2. Within parameters.tp_dist_close mm Euclidean distance from target
%  3. Exit aperture intersects skin surface (>0% overlap required)
%  4. Scores: skin_distance uniformity, skull_distance uniformity, intersection fraction
%
%  Exhaustively evaluates all candidate transducer positions within distance threshold using GPU.
%  Computes intersection fraction, skin/skull distances, and variance metrics for position optimization.
%
% INPUT
%   img        - Segmented head image (nifti final_tissues.nii.gz)  
%   target     - 1x3 target coordinates [x,y,z] in voxel space
%   parameters - Struct with .tp_dist_close, .transducer properties
%   pixel_size - Scalar voxel size in mm
%
% OUTPUT
%   tpos_pars  - Table with columns: idx, trans_x/y/z, targ_x/y/z, dist_to_target, 
%                prop_intersect, mean_dist_skin, var_dist_skin, mean_dist_skull, var_dist_skull

disp("[TP] Evaluating candidate positions (fraction, skin/skull distances, and variance) ...")

sz = size(img);  
[t1_x, t1_y, t1_z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));  % Full volume voxel coordinates

% GPU coordinate mesh (vectorized for parallel evaluation)
coord_mesh.x = gpuArray(t1_x); 
coord_mesh.y = gpuArray(t1_y); 
coord_mesh.z = gpuArray(t1_z);
coord_mesh.xyz = gpuArray([reshape(t1_x,[],1) reshape(t1_y,[],1) reshape(t1_z,[],1)]);  % Nx*Ny*Nz x 3

%% Extract tissue boundaries for distance calculations
all_masks = img > 0;  % Brain+skull+skin
outer_boundary = imdilate(all_masks, strel('sphere',1)) - all_masks;  % Exterior surface (air interface)    
skin_boundary = all_masks - imerode(all_masks, strel('sphere',1));  % Skin surface gradient
skin_coords = gpuArray(coord_mesh.xyz(find(skin_boundary),:));  % Skin surface points

% Skull boundary (compact bone exterior)
img_cp = img; img_cp(img_cp==7|img_cp==8) = 4;  % Merge trabecular/CSF → compact bone
skull = img_cp == 4;
skull_fill = imfill(img_cp>0 & img_cp<=4, 'holes');  % Fill inner skull cavity
skull_fill = imerode(skull_fill, strel('sphere',1));
skull_boundary = skull & ~skull_fill;  % Outer compact bone surface
skull_coords = gpuArray(coord_mesh.xyz(find(skull_boundary),:));

%% Filter candidate positions: outer surface within tp_dist_close of target

outer_idx = find(outer_boundary);
outer_coords = coord_mesh.xyz(outer_idx, :);
coord_rel_to_targ = outer_coords - target;
distances_to_target = sqrt(sum(coord_rel_to_targ.^2, 2));
close_enough_idx = outer_idx(distances_to_target < (parameters.tp_dist_close/pixel_size));

%% GPU vectorized transducer geometry for all candidates
trans_pos_coords = coord_mesh.xyz(close_enough_idx,:);  % N_candidates x 3
norm_v = gpuArray((trans_pos_coords-target) ./ repmat(sqrt(sum((trans_pos_coords-target).^2,2)), [1, 3]));  % Unit vectors

% Physical → voxel geometry conversion
max_od_mm = max(parameters.transducer.Elements_OD_mm);
dist_gf_to_ep_mm = 0.5*sqrt(4*parameters.transducer.curv_radius_mm^2 - max_od_mm^2);  % Sagitta/2
dist_tp_to_ep_mm = parameters.transducer.curv_radius_mm - dist_gf_to_ep_mm;  % Transducer-to-exit distance

% Shift positions forward along ray by 5mm + transducer-to-exit distance
pos_shift_mm = 5 + dist_tp_to_ep_mm;  
shifted_trans_pos_coords = trans_pos_coords + norm_v*(pos_shift_mm/pixel_size);

% Compute geometric focus and exit plane for all candidates
geom_focus_pos_all = shifted_trans_pos_coords - norm_v*(parameters.transducer.curv_radius_mm/pixel_size);
ex_plane_pos_all = gpuArray(geom_focus_pos_all + norm_v*(dist_gf_to_ep_mm/pixel_size));
all_masks_indx_gpu = gpuArray(find(img>0)); 
max_od_grid = max_od_mm / pixel_size;

%% Parallel evaluation of all candidate positions
[prop_intersect, mean_dts, var_dts, mean_dist_skull, var_dist_skull] = ...
    arrayfun(@(x) transducer_analyze_position_fast(x, norm_v, ex_plane_pos_all, coord_mesh, ...
    all_masks_indx_gpu, skin_coords, skull_coords, max_od_grid), 1:length(close_enough_idx));

%% Assemble results table (CPU)
tpos_pars = array2table(gather([close_enough_idx, shifted_trans_pos_coords, repmat(target,[length(close_enough_idx),1]), ...
    pdist2(shifted_trans_pos_coords, target), [prop_intersect; mean_dts; var_dts; mean_dist_skull; var_dist_skull]']), ...
    'VariableNames', ["idx","trans_x","trans_y","trans_z","targ_x","targ_y","targ_z","dist_to_target",...
    "prop_intersect","mean_dist_skin","var_dist_skin","mean_dist_skull","var_dist_skull"]);

end
