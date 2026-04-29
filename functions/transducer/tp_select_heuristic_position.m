function [trans_pos_grid, target_pos_grid, best_trans_pos] = ...
    tp_select_heuristic_position(locs, subject_id, target_name, parameters, img_header)
% TP_SELECT_HEURISTIC_POSITION  Select optimal transducer position and export Localite coordinates
%
% Applies a multi-step heuristic: filters candidates by intersection
% fraction, auto-expands the criterion if no candidates survive, then
% selects the position with minimum distance to the target. Exports
% Localite-compatible coordinates.
%
% Use as:
%   [trans_pos_grid, target_pos_grid, best_trans_pos] = ...
%       tp_select_heuristic_position(locs, subject_id, target_name, parameters, img_header)
%
% Input:
%   locs        - table of candidate metrics (from TP_EVALUATE_CANDIDATE_POSITIONS)
%   subject_id  - subject identifier used in output files
%   target_name - target label string
%   parameters  - (1,1) simulation parameters struct
%   img_header  - NIfTI header struct (from niftiinfo)
%
% Output:
%   trans_pos_grid  - [1x3] transducer position in voxel coordinates
%   target_pos_grid - [1x3] target position in voxel coordinates
%   best_trans_pos  - table row of the selected position
%
% See also: TP_EVALUATE_CANDIDATE_POSITIONS, TP_REMOVE_EAR_LOCATIONS

%% Select heuristic position (multi-criteria heuristic)

% Default parameter values if not provided
defaults = struct( ...
    'criterion_intersection',    0.05, ...
    'criterion_skin_mean',    NaN,  ...
    'criterion_skull_mean',   NaN,  ...
    'criterion_skin_var',     NaN,   ...
    'criterion_skull_var',    NaN,  ...
    'expand_step',            0.01 ...
);

% Fill in missing parameter fields
if ~isfield(parameters, 'placement') || ~isfield(parameters.placement, 'heuristic')
    parameters.placement.heuristic = struct();
end
fields = fieldnames(defaults);
for k = 1:numel(fields)
    if ~isfield(parameters.placement.heuristic, fields{k}) || isempty(parameters.placement.heuristic.(fields{k}))
        parameters.placement.heuristic.(fields{k}) = defaults.(fields{k});
        fprintf('[TP_HEURISTIC] Using default %s = %.3f\n', ...
                fields{k}, defaults.(fields{k}));
    end
end

%% STEP 1 — Intersection (adaptive expansion)

if ~isnan(parameters.placement.heuristic.criterion_intersection)
    criterion = parameters.placement.heuristic.criterion_intersection;
    tppf = locs(locs.prop_intersect < criterion, :);

    while isempty(tppf) && criterion < 0.3
        criterion = criterion + parameters.placement.heuristic.expand_step;
        tppf = locs(locs.prop_intersect < criterion, :);
    end
else
    tppf = locs;
end

%% STEP 2 — Skin proximity
if ~isnan(parameters.placement.heuristic.criterion_skin_mean) && ~isempty(tppf)
    % Use quantile
    skin_lim = quantile(tppf.mean_dist_skin, parameters.placement.heuristic.criterion_skin_mean);
    tmp = tppf(tppf.mean_dist_skin <= skin_lim, :);
    if ~isempty(tmp), tppf = tmp;
    else, fprintf('[TP_HEURISTIC] Skipping skin_mean criterion.\n'); end
end

%% STEP 3 — Skull variance

if ~isnan(parameters.placement.heuristic.criterion_skull_var) && ~isempty(tppf)
    % Use quantile
    skull_var_lim = quantile(tppf.var_dist_skull, parameters.placement.heuristic.criterion_skull_var);
    tmp = tppf(tppf.var_dist_skull <= skull_var_lim, :);
    if ~isempty(tmp), tppf = tmp;
    else, fprintf('[TP_HEURISTIC] Skipping skull_var criterion.\n'); end
end

%% STEP 4 — Skull mean distance (always applied)
if ~isnan(parameters.placement.heuristic.criterion_skull_mean) && ~isempty(tppf)
    % Use quantile
    skull_mean_lim = quantile(tppf.mean_dist_skull, parameters.placement.heuristic.criterion_skull_mean);
    tmp = tppf(tppf.mean_dist_skull <= skull_mean_lim, :);
    if ~isempty(tmp), tppf = tmp;
    else, fprintf('[TP_HEURISTIC] Skipping skull_mean criterion.\n'); end
end

%% STEP 5 — Skin variance
if ~isnan(parameters.placement.heuristic.criterion_skin_var) && ~isempty(tppf)
    % Use quantile
    skin_var_lim = quantile(tppf.var_dist_skin, parameters.placement.heuristic.criterion_skin_var);
    tmp = tppf(tppf.var_dist_skin <= skin_var_lim, :);
    if ~isempty(tmp)
        tppf = tmp;
    else
        fprintf('[TP_HEURISTIC] Skipping skin_var criterion.\n');
    end
end

%% STEP 6 — Minimum distance to target

if ~isempty(tppf)
    tppf = tppf(tppf.dist_to_target == min(tppf.dist_to_target), :);
    i = find(locs.idx == tppf.idx(1));
    fprintf('[TP_HEURISTIC] Selected %d final candidate(s).\n', size(tppf,1));
else
    warn('[TP_HEURISTIC] No candidates remain.');
    i = [];
end

%% Extract HEURISTIC position

best_trans_pos = locs(i, :);
trans_pos_grid = floor(table2array(best_trans_pos(1, ["trans_x", "trans_y", "trans_z"])));
target_pos_grid = floor(table2array(best_trans_pos(1, ["targ_x", "targ_y", "targ_z"])));

%% Convert to Localite physical coordinates (mm)
% IMPORTANT: Assumes that localite planning image has
% canonical_affine_transform applied! Always visually check coordinates ...

trans_pos_localite = round(trans_pos_grid .* img_header.PixelDimensions(1:3), 2);
target_pos_localite = round(target_pos_grid .* img_header.PixelDimensions(1:3), 2);

%% Convert to RAS+ physical coordinates (mm)

trans_pos_ras = round(transform_coordinates(parameters, trans_pos_grid, 'grid', 'ras_plus', img_header),2);
target_pos_ras = round(transform_coordinates(parameters, target_pos_grid, 'grid', 'ras_plus', img_header),2);

%% Add Localite coordinates to output table

best_trans_pos.trans_x = round(best_trans_pos.trans_x);
best_trans_pos.trans_y = round(best_trans_pos.trans_y); 
best_trans_pos.trans_z = round(best_trans_pos.trans_z);

best_trans_pos.LOCALITE_Transducer_mm = Inf;

best_trans_pos.loc_trans_x_mm = trans_pos_localite(1);
best_trans_pos.loc_trans_y_mm = trans_pos_localite(2);
best_trans_pos.loc_trans_z_mm = trans_pos_localite(3);

best_trans_pos.LOCALITE_Target_mm = Inf;

best_trans_pos.loc_targ_x_mm = target_pos_localite(1);
best_trans_pos.loc_targ_y_mm = target_pos_localite(2);
best_trans_pos.loc_targ_z_mm = target_pos_localite(3);

best_trans_pos.RAS_Transducer_mm = Inf;

best_trans_pos.ras_trans_x_mm = trans_pos_ras(1);
best_trans_pos.ras_trans_y_mm = trans_pos_ras(2);
best_trans_pos.ras_trans_z_mm = trans_pos_ras(3);

best_trans_pos.RAS_Target_mm = Inf;

best_trans_pos.ras_targ_x_mm = target_pos_ras(1);
best_trans_pos.ras_targ_y_mm = target_pos_ras(2);
best_trans_pos.ras_targ_z_mm = target_pos_ras(3);

%% Export Localite-compatible table

table_formatted = rows2vars(best_trans_pos);

table_path = fullfile(parameters.io.dir_output, sprintf('sub-%03.0f_%s.txt', subject_id, target_name));
writetable(table_formatted, table_path, 'Delimiter', 'tab', 'WriteVariableNames', false);
fprintf('[TP_HEURISTIC] Saved heuristic position: %s\n', table_path);

%% [Optional] Save copy in localite directory
if isfield(parameters.path, 'localite') && ~isempty(parameters.path.localite)
    localite_path = fullfile(parameters.path.localite,...
        sprintf('sub-%03.0f_%s.txt', subject_id, target_name));
    writetable(table_formatted, localite_path, 'Delimiter', 'tab', 'WriteVariableNames', false);
    fprintf('[TP_HEURISTIC] Additional Localite copy saved: %s\n', localite_path);
end

end