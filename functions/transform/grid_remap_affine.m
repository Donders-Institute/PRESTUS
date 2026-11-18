function [M, info] = grid_remap_affine(A_src, A_dst, offset_src, offset_dst)
% GRID_REMAP_AFFINE  Voxel mapping from one simulation grid (or FOV box) to another
%
% Composes the forward affine that takes voxel indices of an array stored on
% the source grid to voxel indices of the destination grid, through world
% space: [vox_dst 1] = [vox_src 1] * M. Each grid is given by its
% grid-to-world affine (GRID_WORLD_AFFINE, full simulation grid) and,
% for an array that holds only an acoustic FOV box, the box's offset in that
% grid (1-based full-grid index of the box's first voxel). Use M with
% AFFINE_RESAMPLE_3D or REMAP_ACOUSTIC_FIELD.
%
% Use as:
%   [M, info] = grid_remap_affine(A_src, A_dst)
%   [M, info] = grid_remap_affine(A_src, A_dst, offset_src, offset_dst)
%
% Input:
%   A_src, A_dst         - 4x4 grid-to-world affines (row-vector convention)
%   offset_src, offset_dst - [1x3] FOV box offsets ([1 1 1] or [] = the full grid)
%
% Output:
%   M    - 4x4 forward affine, source array voxel -> destination array voxel
%   info - .scale: singular values of the linear part (destination voxels per
%          source voxel; all 1 for two grids of one resolution)
%          .rigid: true when the linear part is a rotation times one scale
%          (no shear, isotropic; tolerance 1e-3)
%          .shift_vox: destination index of source voxel [1 1 1] minus [1 1 1]
%
% See also: GRID_WORLD_AFFINE, REMAP_ACOUSTIC_FIELD, AFFINE_RESAMPLE_3D

    if nargin < 3 || isempty(offset_src); offset_src = [1 1 1]; end
    if nargin < 4 || isempty(offset_dst); offset_dst = [1 1 1]; end
    % box voxel -> full-grid voxel: [full 1] = [box 1] * T(offset - 1)
    T_src = box_to_full(offset_src);
    T_dst = box_to_full(offset_dst);
    M = (T_src * A_src) / (T_dst * A_dst);

    s = svd(M(1:3, 1:3));
    info.scale = s(:)';
    info.rigid = (max(s) - min(s)) <= 1e-3 * max(s);
    p1 = [1 1 1 1] * M;
    info.shift_vox = p1(1:3) - [1 1 1];
end

function T = box_to_full(offset)
    T = eye(4);
    T(4, 1:3) = double(offset(:)') - 1;
end
