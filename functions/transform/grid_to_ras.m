function ras_pos = grid_to_ras(grid_pos, nii_header)
% GRID_TO_RAS  Convert 1-based voxel indices to RAS world coordinates
%
% Inverse of RAS_TO_GRID: maps 1-based MATLAB array indices (niftiread
% order; fractional indices allowed) to RAS mm. niftiinfo's Transform.T is
% the NIfTI affine, which takes 0-based indices, so 1 is subtracted first.
%
% Use as:
%   ras_pos = grid_to_ras(grid_pos, nii_header)
%
% Input:
%   grid_pos   - [Nx3] (or [3x1] for one point) 1-based voxel indices
%   nii_header - niftiinfo header containing Transform.T (4x4 affine)
%
% Output:
%   ras_pos    - [Nx3] positions in RAS world coordinates [mm]
%
% See also: RAS_TO_GRID, TRANSFORM_COORDINATES, NIFTIINFO

arguments
    grid_pos   (:,:) {mustBeNumeric}
    nii_header (1,1) struct
end

    if size(grid_pos, 2) ~= 3 && size(grid_pos, 1) == 3
        grid_pos = grid_pos';  % Convert 3x1 → 1x3
    end

    % Row-vector convention: [ras 1] = [vox0 1] * Transform.T
    pts = [grid_pos - 1, ones(size(grid_pos, 1), 1)] * nii_header.Transform.T;
    ras_pos = pts(:, 1:3);
end
