function coord_mesh_xyz = get_xyz_mesh(img)
% GET_XYZ_MESH  Generate an [Nx3] coordinate mesh for all voxels in a 3D image
%
% Creates an Nx3 array where each row contains the (x, y, z) integer grid
% indices of one voxel, with N = numel(img).
%
% Use as:
%   coord_mesh_xyz = get_xyz_mesh(img)
%
% Input:
%   img - [Nx x Ny x Nz] 3D volume
%
% Output:
%   coord_mesh_xyz - [Nx3] voxel coordinate array [voxels]
%
% See also: NDGRID, RESHAPE

arguments
    img (:,:,:) {mustBeNumeric}
end

    %% Generate coordinate grids for each dimension
    % Create 3D grids for x, y, and z dimensions using `ndgrid`
    [x, y, z] = ndgrid(1:size(img, 1), 1:size(img, 2), 1:size(img, 3));

    %% Reshape coordinate grids into Nx1 vectors
    % Flatten each grid into a column vector and combine them into an Nx3 array
    coord_mesh_xyz = [reshape(x, [], 1), reshape(y, [], 1), reshape(z, [], 1)];
end
