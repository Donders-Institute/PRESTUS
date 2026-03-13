function [tpos_pars, mesh] = tp_evaluate_candidate_positions(img, target, parameters, pixel_size)
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
%   mesh       - Structure with mesh information

disp("[TP] Evaluating candidate position via intersection of skull with expanding sphere ...")

%--- Build candidate mesh and candidate geometry ---------------------------------
mesh = tp_candidate_mesh(img, target, parameters, pixel_size);

%--- Extract required variables for position analysis ----------------------
coord_mesh         = mesh.coord_mesh;
skin_coords        = mesh.skin_coords;
skull_coords       = mesh.skull_coords;
close_enough_idx   = mesh.close_enough_idx;
norm_v             = mesh.norm_v;
ex_plane           = mesh.ex_plane;
all_masks_idx      = mesh.all_masks_idx;
max_od_grid        = mesh.max_od_grid;

%% Parallel evaluation of all candidate positions
N_cand = length(close_enough_idx);
[prop_intersect, mean_dts, var_dts, mean_dist_skull, var_dist_skull] = ...
    arrayfun(@(x) transducer_analyze_position_fast(x, norm_v, ex_plane, coord_mesh, ...
    all_masks_idx, skin_coords, skull_coords, max_od_grid), 1:N_cand);

shifted_trans_pos_coords = gather(mesh.trans_pos);  % [N_cand x 3] CPU for table

%% Assemble results table (CPU)
metrics = [prop_intersect; mean_dts; var_dts; mean_dist_skull; var_dist_skull]';
tpos_pars = array2table([close_enough_idx, shifted_trans_pos_coords, ...
    repmat(target, N_cand, 1), ...
    pdist2(shifted_trans_pos_coords, target), ...
    metrics], ...
    'VariableNames', ["idx","trans_x","trans_y","trans_z","targ_x","targ_y","targ_z",...
    "dist_to_target","prop_intersect","mean_dist_skin","var_dist_skin","mean_dist_skull","var_dist_skull"]);

end