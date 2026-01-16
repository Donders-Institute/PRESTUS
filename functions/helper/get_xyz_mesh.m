function coord_mesh_xyz = get_xyz_mesh(img)

% GET_XYZ_MESH Generates a mesh of 3D coordinates for a given image.
%
% This function creates an Nx3 array of point coordinates for a 3D image (`img`).
% Each row in the output corresponds to the (x, y, z) coordinates of a voxel in 
% the image, where N is the total number of voxels.
%
% Input:
%   img - [Nx x Ny x Nz] matrix representing the 3D image.
%
% Output:
%   coord_mesh_xyz - [N x 3] array of voxel coordinates, where each row contains 
%                    the (x, y, z) coordinates of a voxel in the image.

    %% Generate coordinate grids for each dimension
    % Create 3D grids for x, y, and z dimensions using `ndgrid`
    [x, y, z] = ndgrid(1:size(img, 1), 1:size(img, 2), 1:size(img, 3));

    %% Reshape coordinate grids into Nx1 vectors
    % Flatten each grid into a column vector and combine them into an Nx3 array
    coord_mesh_xyz = [reshape(x, [], 1), reshape(y, [], 1), reshape(z, [], 1)];
end
