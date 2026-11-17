function SKULL_BALLON = skull_rubber_wrap_reference(BW, BRAIN, SKIN, CSF, r)
% SKULL_RUBBER_WRAP_REFERENCE  Full-grid rubber wrap, the algorithm SKULL_RUBBER_WRAP replaced (test reference)
%
% Closes the largest skull component with a sphere of radius r on the whole
% grid, fills holes, keeps the largest component of the result, and adds
% the wrapped voxels that touch both brain and skin, CSF voxels at the skin
% next to them, and single-voxel z gaps. Used by the tests and the
% comparison job to check that SKULL_RUBBER_WRAP gives the same mask.
%
% Use as:
%   SKULL_BALLON = skull_rubber_wrap_reference(BW, BRAIN, SKIN, CSF, r)
%
% Input:
%   BW    - [Nx x Ny x Nz] logical skull mask
%   BRAIN - [Nx x Ny x Nz] logical brain mask
%   SKIN  - [Nx x Ny x Nz] logical skin mask
%   CSF   - [Nx x Ny x Nz] logical CSF mask (SimNIBS labels)
%   r     - wrap radius [grid voxels]
%
% Output:
%   SKULL_BALLON - [Nx x Ny x Nz] logical filled skull mask
%
% See also: SKULL_RUBBER_WRAP

    sz = size(BW);
    CC = bwconncomp(BW, 26);
    [~, iBig] = max(cellfun(@numel, CC.PixelIdxList));
    BW1 = false(sz);
    BW1(CC.PixelIdxList{iBig}) = true;

    dilated = bwdist(BW1) <= r;
    BWwrap  = ~(bwdist(~dilated) <= r);
    BWwrap  = imfill(BWwrap, 'holes');
    CCw = bwconncomp(BWwrap, 26);
    [~, iw] = max(cellfun(@numel, CCw.PixelIdxList));
    W = false(sz);
    W(CCw.PixelIdxList{iw}) = true;

    se1 = sphere_kernel(1);
    B_out = W & (convn(BRAIN, se1, 'same') > 0) & (convn(SKIN, se1, 'same') > 0) & ~BRAIN;
    near = convn(B_out, sphere_kernel(2), 'same') > 0;
    B_out2 = (B_out | (CSF & (convn(SKIN, se1, 'same') > 0) & near)) & ~BRAIN;
    B0 = BW1 | B_out2;

    BWm1 = circshift(B0, [0 0  1]); BWm1(:, :, 1)   = false;
    BWp1 = circshift(B0, [0 0 -1]); BWp1(:, :, end) = false;
    SKULL_BALLON = B0 | (~B0 & BWm1 & BWp1);
end

function se = sphere_kernel(r)
    [kx, ky, kz] = ndgrid(-r:r);
    se = double(kx.^2 + ky.^2 + kz.^2 <= r^2);
end
