function [tpos_pars, mesh] = tp_evaluate_candidate_positions(img, target, parameters, pixel_size)
% TP_EVALUATE_CANDIDATE_POSITIONS  Evaluate candidate transducer positions on skull surface
%
% Exhaustively evaluates all candidate positions within a distance
% threshold of the target using GPU acceleration. Computes intersection
% fraction, skin/skull distance metrics, and variance scores for position
% optimisation.
%
% Use as:
%   [tpos_pars, mesh] = tp_evaluate_candidate_positions(img, target, parameters, pixel_size)
%
% Input:
%   img        - [Nx x Ny x Nz] segmented head volume
%   target     - [1x3] target coordinates in voxel space
%   parameters - (1,1) simulation parameters struct
%   pixel_size - voxel size [mm]
%
% Output:
%   tpos_pars - table with columns: idx, trans_x/y/z, targ_x/y/z,
%               dist_to_target, prop_intersect, mean_dist_skin,
%               var_dist_skin, mean_dist_skull, var_dist_skull
%   mesh      - struct with mesh information (see TP_CANDIDATE_MESH)
%
% See also: TP_CANDIDATE_MESH, TP_SELECT_HEURISTIC_POSITION

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
trans_pos          = mesh.trans_pos;

%% Parallel evaluation of all candidate positions

N_cand = length(close_enough_idx);
[prop_intersect, mean_dts, var_dts, mean_dist_skull, var_dist_skull] = ...
    arrayfun(@(x) transducer_analyze_position_fast(x, norm_v, ex_plane, coord_mesh, ...
    all_masks_idx, skin_coords, skull_coords, max_od_grid), 1:N_cand);

%% Assemble results table (CPU)

% Ensure all components are on CPU before table assembly
close_enough_idx = gather(close_enough_idx);
trans_pos = gather(trans_pos);
metrics = gather([prop_intersect; mean_dts; var_dts; mean_dist_skull; var_dist_skull]');
dist_to_target = gather(pdist2(trans_pos, target));

tpos_pars = array2table([close_enough_idx, trans_pos, ...
    repmat(target, N_cand, 1), ...
    dist_to_target, ...
    metrics], ...
    'VariableNames', ["idx","trans_x","trans_y","trans_z","targ_x","targ_y","targ_z",...
    "dist_to_target","prop_intersect","mean_dist_skin","var_dist_skin","mean_dist_skull","var_dist_skull"]);

end