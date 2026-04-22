function [val, Ix, Iy, Iz] = masked_max_3d(array_3d, mask)
% MASKED_MAX_3D  Find the maximum value and its 3D subscript within a masked region
%
% Sets voxels outside the mask to NaN, then returns the max and its
% [x, y, z] subscript indices.
%
% Use as:
%   [val, Ix, Iy, Iz] = masked_max_3d(array_3d, mask)
%
% Input:
%   array_3d - [Nx x Ny x Nz] 3D value array
%   mask     - [Nx x Ny x Nz] binary mask (1 = included, 0 = excluded)
%
% Output:
%   val - maximum value within the masked region
%   Ix  - x-index of the maximum
%   Iy  - y-index of the maximum
%   Iz  - z-index of the maximum
%
% See also: ACOUSTIC_ANALYSIS, THERMAL_ANALYSIS

arguments
    array_3d {mustBeNumeric}
    mask     {mustBeNumericOrLogical}
end

    % Exclude values outside the mask by setting them to NaN
    array_3d(~mask) = nan;

    % Find the maximum value and its linear index within the masked region
    [val, I] = max(array_3d(:));

    % Convert linear index to subscripts (x, y, z coordinates)
    [Ix, Iy, Iz] = ind2sub(size(array_3d), I);
end
