function [grid_pos] = ras_to_grid(ras_pos, nii_header)
% RAS_TO_GRID  Convert RAS world coordinates to integer voxel grid coordinates
%
% Applies the inverse of the NIfTI affine transformation to map a physical
% position in RAS (Right-Anterior-Superior) space to the nearest voxel index.
% Accepts either a [3x1] or [1x3] position vector and always returns [1x3].
%
% Use as:
%   grid_pos = ras_to_grid(ras_pos, nii_header)
%
% Input:
%   ras_pos    - [3x1] or [1x3] position in RAS world coordinates [mm]
%   nii_header - niftiinfo header containing Transform.T (4x4 affine)
%
% Output:
%   grid_pos   - [1x3] nearest voxel indices [voxels]
%
% See also: TRANSFORM_COORDINATES, CANONICAL_AFFINE_TRANSFORM, NIFTIINFO

arguments
    ras_pos    (:,:) {mustBeNumeric}
    nii_header (1,1) struct
end

    if size(ras_pos, 2) == 3 && size(ras_pos, 1) == 1
        ras_pos = ras_pos';  % Convert 1x3 → 3x1
    end

    % Apply the inverse transformation matrix from the NIfTI header
    grid_pos = round(nii_header.Transform.T' \ [ras_pos; 1]);

    % Extract x, y, z components (ignore homogeneous coordinate)
    grid_pos = grid_pos(1:3)';
end