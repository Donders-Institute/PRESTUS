function plot_median_montage(t1_image)

% PLOT_MEDIAN_MONTAGE  Display a montage of central orthogonal slices from a 3D image
%
% Extracts the central slices along x, y, and z from a 3D image, normalises
% each slice with mat2gray, and displays them as a montage.
%
% Use as:
%   plot_median_montage(t1_image)
%
% Input:
%   t1_image - [Nx x Ny x Nz] 3D image volume
%
% See also: PLOT_CORONAL_SLICES, PLOT_OVERLAY

arguments
    t1_image (:,:,:) {mustBeNumeric}
end

    % Compute the center position of the image
    t1_center = round((size(t1_image) + 1) / 2);

    % Extract and normalize the central slices along each axis
    im1 = mat2gray(squeeze(t1_image(t1_center(1), :, :))); % Slice along x-axis
    im2 = mat2gray(squeeze(t1_image(:, t1_center(2), :))); % Slice along y-axis
    im3 = mat2gray(squeeze(t1_image(:, :, t1_center(3)))); % Slice along z-axis

    % Create and display a montage of the slices
    montage({im1, im2, im3}, 'Size', [1 nan]);
end
