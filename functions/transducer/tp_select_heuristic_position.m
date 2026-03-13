function [trans_pos, target_pos, best_trans_pos] = ...
    tp_select_heuristic_position(locs, pixel_size, subject_id, target_name, parameters)
%% TP_SELECT_HEURISTIC_POSITION Select optimal transducer position and export Localite coordinates
%
% Multi-step HEURISTIC selection:
% 1. Filter by intersection fraction < criterion_intersection (default 0.05)
% 2. Auto-expand criterion if no candidates found
% 3. Select minimum distance to target from valid candidates
%
% INPUT
%   locs        - Table with transducer candidate metrics (from tp_evaluate_candidate_positions)
%   pixel_size  - Scalar voxel size (mm/voxel)
%   subject_id  - Scalar subject ID
%   target_name - String target identifier (e.g. 'right_PUL')
%   parameters  - [OPTIONAL] Struct with .tp_criterion_intersection, .output_dir, .localite_path
%
% OUTPUT
%   trans_pos        - [1x3] integer voxel coordinates [x y z]
%   target_pos       - [1x3] integer voxel coordinates [x y z]  
%   best_trans_pos   - Full table row of HEURISTIC position with added Localite fields

%% Select optimal position (multi-criteria heuristic)
%tppf = tpos_pars(tpos_pars.prop_intersect<0.05,:);
%tppf = tppf(tppf.mean_dist_skull <= quantile(tppf.mean_dist_skull, 0.5) & ...
%           tppf.var_dist_skull <= quantile(tppf.var_dist_skull, 0.1),:);
%tppf = tppf(tppf.var_dist_skin==min(tppf.var_dist_skin),:);

%% Check for configurable intersection criterion
if nargin < 5 || ~isfield(parameters, 'tp_criterion_intersection') || isempty(parameters.tp_criterion_intersection)
    criterion_intersection = 0.05;  % Default
    fprintf('[TP_HEURISTIC] Using default intersection criterion: %.2f\n', criterion_intersection);
else
    criterion_intersection = parameters.tp_criterion_intersection;
    fprintf('[TP_HEURISTIC] Using specified intersection criterion: %.2f\n', criterion_intersection);
end

%% Multi-step HEURISTIC selection with adaptive criterion
tppf = locs(locs.prop_intersect < criterion_intersection, :);

% Auto-expand criterion if no candidates found
if isempty(tppf)
    criterion_intersection = criterion_intersection + 0.01;
    warning('Expanding intersection criterion to %.2f (no candidates at %.2f)', ...
            criterion_intersection, criterion_intersection - 0.01);
    tppf = locs(locs.prop_intersect < criterion_intersection, :);
end

% Select minimum distance to target
tppf = tppf(tppf.dist_to_target == min(tppf.dist_to_target), :);
i = find(locs.idx == tppf.idx(1));  % Global index in original locs table

%% Extract HEURISTIC position
best_trans_pos = locs(i, :);
trans_pos = floor(table2array(best_trans_pos(1, ["trans_x", "trans_y", "trans_z"])));
target_pos = floor(table2array(best_trans_pos(1, ["targ_x", "targ_y", "targ_z"])));

%% Convert to Localite physical coordinates (mm)
trans_pos_localite = round(trans_pos * pixel_size);
target_pos_localite = round(target_pos * pixel_size);

%% Add Localite coordinates to output table
best_trans_pos.trans_x = round(best_trans_pos.trans_x);
best_trans_pos.trans_y = round(best_trans_pos.trans_y); 
best_trans_pos.trans_z = round(best_trans_pos.trans_z);

best_trans_pos.loc_trans_x = trans_pos_localite(1);
best_trans_pos.loc_trans_y = trans_pos_localite(2);
best_trans_pos.loc_trans_z = trans_pos_localite(3);

best_trans_pos.loc_targ_x = target_pos_localite(1);
best_trans_pos.loc_targ_y = target_pos_localite(2);
best_trans_pos.loc_targ_z = target_pos_localite(3);

%% Export Localite-compatible table
table_path = fullfile(parameters.output_dir, sprintf('sub-%03.0f_%s.txt', subject_id, target_name));
writetable(best_trans_pos, table_path, 'Delimiter', 'tab');
fprintf('[TP_HEURISTIC] Saved heuristic position: %s\n', table_path);

%% [Optional] Save copy in localite directory
if isfield(parameters, 'localite_path') && ~isempty(parameters.localite_path)
    localite_path = fullfile(parameters.localite_path,...
        sprintf('sub-%03.0f_%s.txt', subject_id, target_name));
    writetable(best_trans_pos, localite_path, 'Delimiter', 'tab');
    fprintf('[TP_HEURISTIC] Additional Localite copy saved: %s\n', localite_path);
end

end