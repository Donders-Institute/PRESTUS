function median_montage(t1_image)

% MEDIAN_MONTAGE Creates a montage of central slices from a 3D T1 image.
%
% This function extracts the central slices along the x, y, and z axes from 
% a 3D T1-weighted image (`t1_image`) and creates a montage for visualization. 
% The slices are normalized for display using `mat2gray`.
%
% Input:
%   t1_image - [Nx x Ny x Nz] matrix representing the 3D T1-weighted image.
%
% Output:
%   None. The function displays a montage of the central slices.

    % Compute the center position of the image
    t1_center = round((size(t1_image) + 1) / 2);

    % Extract and normalize the central slices along each axis
    im1 = mat2gray(squeeze(t1_image(t1_center(1), :, :))); % Slice along x-axis
    im2 = mat2gray(squeeze(t1_image(:, t1_center(2), :))); % Slice along y-axis
    im3 = mat2gray(squeeze(t1_image(:, :, t1_center(3)))); % Slice along z-axis

    % Create and display a montage of the slices
    montage({im1, im2, im3}, 'Size', [1 nan]);
end
