function vol_full = fov_backproject_volume(vol_fov, fov_offset, full_dims)
% FOV_BACKPROJECT_VOLUME  Zero-pad a FOV sub-grid volume back to the full acoustic grid.
%
% Inserts vol_fov into a zero-initialised array of size full_dims at the
% 1-based position fov_offset.  Values outside the inserted region are zero.
%
% Use as:
%   vol_full = fov_backproject_volume(vol_fov, fov_offset, full_dims)
%
% Input:
%   vol_fov    - numeric array of any class with the FOV sub-grid dimensions
%   fov_offset - [1×3] 1-based start voxel index in the full grid
%   full_dims  - [1×3] full acoustic grid dimensions
%
% Output:
%   vol_full   - array of class matching vol_fov, size full_dims
%
% See also: ACOUSTIC_GRID_FOV, ASSEMBLE_LIMITED_FOV_FIELDS

arguments
    vol_fov    {mustBeNumeric}
    fov_offset (1,3) {mustBeNumeric}
    full_dims  (1,3) {mustBeNumeric}
end

vol_full = zeros(full_dims, 'like', vol_fov);
fov_sz   = size(vol_fov);
fov_end  = fov_offset + fov_sz - 1;
vol_full(fov_offset(1):fov_end(1), ...
         fov_offset(2):fov_end(2), ...
         fov_offset(3):fov_end(3)) = vol_fov;
end
