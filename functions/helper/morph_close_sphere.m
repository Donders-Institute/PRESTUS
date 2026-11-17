function W = morph_close_sphere(S, r, use_ipt)
% MORPH_CLOSE_SPHERE  Morphological closing of a 3-D mask by a sphere of radius r
%
% Dilates S by a sphere of radius r and erodes the result by the same
% sphere, using Euclidean distance transforms: with the Image Processing
% Toolbox bwdist, otherwise an exact base-MATLAB distance transform capped
% at r. Both give the same mask. Voxels outside the array are not part of
% the mask or its background (as bwdist treats them).
%
% Use as:
%   W = morph_close_sphere(S, r)
%   W = morph_close_sphere(S, r, use_ipt)
%
% Input:
%   S       - [Nx x Ny x Nz] logical mask
%   r       - sphere radius [voxels], > 0
%   use_ipt - true to use bwdist (default: when the toolbox is on the path
%             and licensed)
%
% Output:
%   W - [Nx x Ny x Nz] logical closed mask
%
% See also: SKULL_RUBBER_WRAP

    if nargin < 3 || isempty(use_ipt)
        use_ipt = license('test', 'image_toolbox') && ~isempty(which('bwdist'));
    end
    S = logical(S);
    if use_ipt
        D = bwdist(S) <= r;
        W = ~(bwdist(~D) <= r);
    else
        D = sqdist_capped(S, r) <= r^2;
        W = ~(sqdist_capped(~D, r) <= r^2);
    end
end

function D2 = sqdist_capped(S, r)
% Squared Euclidean distance to the true voxels of S, exact where it is <= r^2
%
% One min-plus pass per axis over offsets |k| <= floor(r): every voxel
% within distance r differs by at most floor(r) along each axis, so the
% three passes find it; larger values may be overestimates.
    R = floor(r);
    D2 = inf(size(S), 'single');
    D2(S) = 0;
    for d = 1:3
        order = [d, setdiff(1:3, d)];
        P = permute(D2, order);
        G = P;
        n = size(P, 1);
        for k = 1:min(R, n - 1)
            w = single(k^2);
            G(1+k:n, :, :) = min(G(1+k:n, :, :), P(1:n-k, :, :) + w);
            G(1:n-k, :, :) = min(G(1:n-k, :, :), P(1+k:n, :, :) + w);
        end
        D2 = ipermute(G, order);
    end
end
