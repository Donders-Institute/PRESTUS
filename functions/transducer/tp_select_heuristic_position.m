function [best_trans_pos, trans_pos, target_pos, table_path] = ...
    tp_select_heuristic_position(locs, pixel_size, subject_id, target_name, rootpath)
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
%   rootpath    - Root output directory
%
% OUTPUT
%   best_trans_pos   - Full table row of HEURISTIC position with added Localite fields
%   trans_pos        - [1x3] integer voxel coordinates [x y z]
%   target_pos       - [1x3] integer voxel coordinates [x y z]  
%   table_path       - Full path to saved Localite .txt file

%% Multi-step HEURISTIC selection with adaptive intersection criterion
criterion_intersection = 0.05;
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
table_path = fullfile(rootpath, 'data', 'localite', sprintf('sub-%03.0f', subject_id), ...
                     sprintf('sub-%03.0f_%s.txt', subject_id, target_name));

writetable(best_trans_pos, table_path, 'Delimiter', 'tab');
fprintf('[TP_HEURISTIC] Saved heuristic position: %s\n', table_path);

end