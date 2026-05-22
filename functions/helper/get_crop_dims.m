function [min_dims, max_dims, grid_size] = get_crop_dims(image, margin)
% GET_CROP_DIMS  Compute bounding-box crop dimensions for a 3D binary/label image
%
%   Finds the smallest axis-aligned bounding box that contains all non-zero
%   voxels in image, then expands the box symmetrically by margin voxels.
%
% Use as:
%   [min_dims, max_dims, grid_size] = get_crop_dims(image, margin)
%
% Input:
%   image  - [Nx x Ny x Nz] 3D label or binary volume
%   margin - [1x1] voxels to add around the bounding box [voxels]
%
% Output:
%   min_dims  - [1x3] lower-bound voxel indices [x y z] (may be < 1 if margin is large)
%   max_dims  - [1x3] upper-bound voxel indices [x y z]
%   grid_size - [1x3] size of the cropped region [voxels]
%
% See also: PREPROC_CROP_GRID, FIND_MIN_FACTOR

arguments
    image  (:,:,:) {mustBeNumeric}
    margin (1,1) {mustBeNonnegative, mustBeInteger}
end


    % Find indices of non-zero elements in the image
    I = find(image);

    % Convert linear indices to subscripts (x, y, z coordinates)
    [x, y, z] = ind2sub(size(image), I);

    % Compute minimum dimensions and apply margin
    min_dims = [min(x), min(y), min(z)] - margin;

    % Compute maximum dimensions and apply margin
    max_dims = [max(x), max(y), max(z)] + margin;

    % Compute grid size based on min and max dimensions
    grid_size = max_dims - min_dims + 1;
end
