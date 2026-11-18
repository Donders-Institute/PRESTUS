function [plane_keep, plane_stride] = thermal_plane_frames(total_timepoints, plane_numel, budget_gb)
% THERMAL_PLANE_FRAMES  Timepoints whose focal-plane frames the thermal simulation keeps.
%
%   [plane_keep, plane_stride] = thermal_plane_frames(total_timepoints, plane_numel, budget_gb)
%
% The thermal simulation stores three double frames (T, CEM43, CEM43 ISO) of
% plane_numel voxels per timepoint. When all timepoints would exceed budget_gb
% (thermal.max_focal_plane_gb), every plane_stride-th timepoint is kept, plus
% the first (initial state) and the last (final state). budget_gb = 0 keeps
% every timepoint.
%
% Inputs:
%   total_timepoints - number of timepoints (initial state + one per segment)
%   plane_numel      - voxels per frame
%   budget_gb        - memory budget for the three frame series [GB]; 0 = no limit
%
% Outputs:
%   plane_keep   - 1 x total_timepoints logical, true for kept timepoints
%   plane_stride - step between kept timepoints (1 = all kept)

arguments
    total_timepoints (1,1) double {mustBeInteger, mustBePositive}
    plane_numel      (1,1) double {mustBeInteger, mustBePositive}
    budget_gb        (1,1) double {mustBeNonnegative}
end

if budget_gb == 0
    plane_stride = 1;
else
    plane_stride = max(1, ceil(total_timepoints * 3 * 8 * plane_numel / (budget_gb * 1e9)));
end
n_seg = total_timepoints - 1;
plane_keep = false(1, total_timepoints);
plane_keep([1, 1 + (plane_stride:plane_stride:n_seg), total_timepoints]) = true;

end
