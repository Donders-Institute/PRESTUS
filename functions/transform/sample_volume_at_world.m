function vals = sample_volume_at_world(vol, nii_header, world_pts)
% SAMPLE_VOLUME_AT_WORLD  Nearest-neighbour values of a NIfTI volume at world (mm) positions
%
% Maps world coordinates (RAS mm for subject images, MNI mm for MNI-space
% images) through the inverse NIfTI affine to 1-based voxel indices and
% returns the nearest voxel's value; points outside the volume give NaN.
% Used to compare volumes stored on different grids (e.g. an MNI ROI mask
% at 2 mm against 1 mm PRESTUS MNI outputs) without resampling whole volumes.
%
% Use as:
%   vals = sample_volume_at_world(vol, nii_header, world_pts)
%
% Input:
%   vol        - 3-D array as returned by niftiread
%   nii_header - niftiinfo header of vol (Transform.T, 0-based indices)
%   world_pts  - [N×3] positions in the volume's world space [mm]
%
% Output:
%   vals       - [N×1] double; NaN outside the volume
%
% See also: GRID_TO_RAS, RAS_TO_GRID

arguments
    vol
    nii_header (1,1) struct
    world_pts (:,3) {mustBeNumeric}
end

    n = size(world_pts, 1);
    vals = nan(n, 1);
    if n == 0; return; end
    % Row-vector convention: [ras 1] = [vox0 1] * T  ->  [vox0 1] = [ras 1] / T
    vox0 = [double(world_pts), ones(n, 1)] / nii_header.Transform.T;
    idx = round(vox0(:, 1:3)) + 1;
    sz = size(vol);
    sz(end+1:3) = 1;
    inside = all(idx >= 1, 2) & idx(:, 1) <= sz(1) & idx(:, 2) <= sz(2) & idx(:, 3) <= sz(3);
    lin = sub2ind(sz(1:3), idx(inside, 1), idx(inside, 2), idx(inside, 3));
    vals(inside) = double(vol(lin));
end
