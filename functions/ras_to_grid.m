function [grid_pos] = ras_to_grid(ras_pos, nii_header)

% RAS_TO_GRID Converts RAS coordinates to grid (voxel) coordinates.
%
% This function converts a position in RAS (Right-Anterior-Superior) coordinates 
% to grid (voxel) coordinates using the transformation matrix from a NIfTI header.
%
% Input:
%   ras_pos    - [1x3] array specifying the position in RAS coordinates.
%   nii_header - Struct containing the NIfTI header, including the transformation matrix.
%
% Output:
%   grid_pos   - [1x3] array specifying the position in grid (voxel) coordinates.

    % Apply the inverse transformation matrix from the NIfTI header
    grid_pos = round(nii_header.Transform.T' \ [ras_pos; 1]);

    % Extract x, y, z components (ignore homogeneous coordinate)
    grid_pos = grid_pos(1:3);
end