function grid_pos = fov_to_grid_pos(fov_pos, fov_offset_ac)
% FOV_TO_GRID_POS  Convert FOV-grid voxel indices to full sim-grid indices.
%
% Applies the inverse of the acoustic FOV crop: a position in the cropped
% FOV grid (1-based) maps to the full simulation grid via
%
%   grid_pos = fov_pos + fov_offset_ac - 1
%
% The `- 1` correction arises because fov_offset_ac is the 1-based start
% index of the FOV box in the full grid, so the first FOV voxel (index 1)
% corresponds to full-grid index fov_offset_ac, not fov_offset_ac + 1.
%
% Use as:
%   grid_pos = fov_to_grid_pos(fov_pos, fov_offset_ac)
%
% Input:
%   fov_pos      - [1×N] position in FOV-grid voxels (1-based)
%   fov_offset_ac - [1×N] 1-based start voxel of the FOV box in the full grid
%                  (returned by acoustic_grid_fov)
%
% Output:
%   grid_pos - [1×N] position in full sim-grid voxels (1-based)
%
% See also: ACOUSTIC_GRID_FOV, FOV_BACKPROJECT_VOLUME

arguments
    fov_pos       (1,:) {mustBeNumeric}
    fov_offset_ac (1,:) {mustBeNumeric}
end

grid_pos = fov_pos + fov_offset_ac - 1;
end
