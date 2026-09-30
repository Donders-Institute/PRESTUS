function [grid_pos] = ras_to_grid(ras_pos, nii_header)
% RAS_TO_GRID  Convert RAS world coordinates to 1-based voxel indices
%
% Maps a physical position in RAS (Right-Anterior-Superior) space to the
% nearest voxel, as a 1-based MATLAB index into the image array returned by
% niftiread. niftiinfo's Transform.T is the NIfTI affine itself, which maps
% 0-based voxel indices to RAS, so the inverse affine gives a 0-based index
% and 1 is added. Accepts either a [3x1] or [1x3] position vector and always
% returns [1x3]. The inverse is GRID_TO_RAS.
%
% Use as:
%   grid_pos = ras_to_grid(ras_pos, nii_header)
%
% Input:
%   ras_pos    - [3x1] or [1x3] position in RAS world coordinates [mm]
%   nii_header - niftiinfo header containing Transform.T (4x4 affine)
%
% Output:
%   grid_pos   - [1x3] nearest voxel, 1-based array index [voxels]
%
% See also: GRID_TO_RAS, TRANSFORM_COORDINATES, CANONICAL_AFFINE_TRANSFORM, NIFTIINFO

arguments
    ras_pos    (:,:) {mustBeNumeric}
    nii_header (1,1) struct
end

    if size(ras_pos, 2) == 3 && size(ras_pos, 1) == 1
        ras_pos = ras_pos';  % Convert 1x3 → 3x1
    end

    % Inverse NIfTI affine gives the 0-based voxel; MATLAB arrays are 1-based
    grid_pos = round(nii_header.Transform.T' \ [ras_pos; 1]) + 1;

    % Extract x, y, z components (ignore homogeneous coordinate)
    grid_pos = grid_pos(1:3)';
end
