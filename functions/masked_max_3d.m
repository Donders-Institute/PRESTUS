function [val, Ix, Iy, Iz] = masked_max_3d(array_3d, mask)

% MASKED_MAX_3D Computes the maximum intensity within a masked 3D region.
%
% This function takes a 3D intensity matrix (`array_3d`) and a binary mask (`mask`) 
% to calculate the maximum intensity value within the masked region. It also returns 
% the corresponding indices along the x, y, and z axes.
%
% Input:
%   array_3d - [Nx x Ny x Nz] matrix representing the 3D intensity values.
%   mask     - [Nx x Ny x Nz] binary matrix defining the region of interest (1 = included, 0 = excluded).
%
% Output:
%   val      - Maximum intensity value within the masked region.
%   Ix       - Index along the x-axis corresponding to the maximum intensity value.
%   Iy       - Index along the y-axis corresponding to the maximum intensity value.
%   Iz       - Index along the z-axis corresponding to the maximum intensity value.

    % Exclude values outside the mask by setting them to NaN
    array_3d(~mask) = nan;

    % Find the maximum value and its linear index within the masked region
    [val, I] = max(array_3d(:));

    % Convert linear index to subscripts (x, y, z coordinates)
    [Ix, Iy, Iz] = ind2sub(size(array_3d), I);
end
