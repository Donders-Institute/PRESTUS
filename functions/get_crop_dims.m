function [min_dims, max_dims, grid_size] = get_crop_dims(image, margin)

% GET_CROP_DIMS Computes the cropping dimensions for an image with a margin.
%
% This function calculates the minimum and maximum dimensions required to crop 
% a 3D image (`image`) based on its non-zero elements. A specified `margin` is 
% added around the bounding box of the non-zero elements to ensure sufficient padding.
%
% Input:
%   image  - [Nx x Ny x Nz] matrix representing the 3D image.
%   margin - Scalar specifying the margin to add around the bounding box.
%
% Output:
%   min_dims - [1x3] array specifying the minimum dimensions of the cropped region.
%   max_dims - [1x3] array specifying the maximum dimensions of the cropped region.
%   grid_size - [1x3] array specifying the size of the cropped region (in voxels).


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
