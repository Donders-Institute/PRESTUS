function [geom_corr_m, apex_to_target_mm, source_tag] = babelbrain_placement_correction( ...
        tr, focus_pos, trans_pos, res_mm, preset_focal_length_m)
% BABELBRAIN_PLACEMENT_CORRECTION  Mechanical Tx offset that puts BabelBrain's apex at trans_pos
%
% BabelBrain places the transducer apex at its preset FocalLength from the
% target voxel. PRESTUS places the apex (trans_pos = bowl rear centre) at
% |focus_pos - trans_pos| from the target. The difference is applied as
% TxMechanicalAdjustmentZ so both solvers see the same source position:
%
%   geom_corr = preset_FL - |focus_pos - trans_pos| * res_mm
%
% That distance equals focal_distance_bowl only when the target sits at the
% nominal focal depth; planned / pFOV placements generally differ (e.g. 71 mm
% vs 92 mm), and using focal_distance_bowl there shifts the whole BabelBrain
% field axially relative to the head. focal_distance_bowl (or
% focal_distance_ep + bowl offset) is used only as a fallback when trans_pos is
% unavailable.
%
% Use as:
%   [geom_corr_m, apex_to_target_mm, source_tag] = ...
%       babelbrain_placement_correction(tr, focus_pos, trans_pos, res_mm, preset_FL_m)
%
% Input:
%   tr                    - parameters.transducer(i) struct
%   focus_pos, trans_pos  - [1x3] voxel positions in the same grid ([] if unknown)
%   res_mm                - grid step [mm] of that grid
%   preset_focal_length_m - BabelBrain preset FocalLength [m] ([] if unknown)
%
% Output:
%   geom_corr_m        - placement correction [m] (0 when it cannot be derived)
%   apex_to_target_mm  - apex-to-target distance used [mm] ([] if none)
%   source_tag         - 'trans_pos', 'focal_distance_bowl', or 'none'
%
% See also: BABELBRAIN_WRAPPER, EXPORT_BABELBRAIN_INPUTS

geom_corr_m       = 0.0;
apex_to_target_mm = [];
source_tag        = 'none';

if isempty(preset_focal_length_m)
    return
end

if ~isempty(trans_pos) && numel(trans_pos) >= 3 && ~isempty(focus_pos) && numel(focus_pos) >= 3 && ...
        all(isfinite(trans_pos(1:3))) && all(isfinite(focus_pos(1:3)))
    apex_to_target_mm = norm(double(focus_pos(1:3)) - double(trans_pos(1:3))) * res_mm;
    source_tag = 'trans_pos';
else
    if isfield(tr, 'focal_distance_offset') && ~isempty(tr.focal_distance_offset)
        bowl_off_mm = tr.focal_distance_offset;
    elseif isfield(tr, 'annular') && all(isfield(tr.annular, {'curv_radius_mm','dist_geom_ep_mm'}))
        bowl_off_mm = tr.annular.curv_radius_mm - tr.annular.dist_geom_ep_mm;
    else
        bowl_off_mm = 0;
    end
    if isfield(tr, 'focal_distance_bowl') && ~isempty(tr.focal_distance_bowl)
        apex_to_target_mm = tr.focal_distance_bowl;
    elseif isfield(tr, 'focal_distance_ep') && ~isempty(tr.focal_distance_ep)
        apex_to_target_mm = tr.focal_distance_ep + bowl_off_mm;
    end
    if ~isempty(apex_to_target_mm)
        source_tag = 'focal_distance_bowl';
    end
end

if ~isempty(apex_to_target_mm)
    geom_corr_m = preset_focal_length_m - apex_to_target_mm / 1000;
end
end
