function [coords, mask] = tp_scalp_boundary(img)
% TP_SCALP_BOUNDARY  Outer (scalp) boundary voxels of a segmented head volume
%
% Returns the one-voxel shell immediately OUTSIDE the head (the set of
% background voxels adjacent to any tissue), which is where candidate
% transducer positions sit. Defined as a one-voxel dilation of (img>0)
% minus (img>0), using a 6-connected (face-adjacent) structuring element —
% identical to the scalp definition used by the heuristic placement
% (see TP_CANDIDATE_MESH). CPU-only (no gpuArray) so it is reusable from the
% non-GPU 'mni' placement path.
%
% Use as:
%   coords         = tp_scalp_boundary(img)
%   [coords, mask] = tp_scalp_boundary(img)
%
% Input:
%   img    - [Nx x Ny x Nz] segmented head volume (tissue labels; 0 = background)
%
% Output:
%   coords - [M x 3] voxel coordinates (1-based) of the outer-boundary shell,
%            in the same column-major order as ndgrid(1:Nx,1:Ny,1:Nz)
%   mask   - [Nx x Ny x Nz] double array, 1 on the outer-boundary shell, else 0
%
% See also: TP_CANDIDATE_MESH, TRANSDUCER_SCALP_GEOMETRY

arguments
    img
end

    all_masks     = img > 0;                       % brain + skull + skin
    [kx1,ky1,kz1] = ndgrid(-1:1,-1:1,-1:1);
    se1           = (kx1.^2 + ky1.^2 + kz1.^2) <= 1;
    mask          = (convn(logical(all_masks), se1, 'same') > 0) - logical(all_masks);

    idx           = find(mask);
    [ix, iy, iz]  = ind2sub(size(img), idx);
    coords        = [ix, iy, iz];
end
