function Minv = affine_inverse(M)
% AFFINE_INVERSE  Inverse of a 4x4 affine, in the convention it was given
%
% Inverts the linear block and the translation separately, so the fixed row
% or column of the result is exactly [0 0 0 1] (no floating-point residue
% for affine3d or later products to trip over). Both conventions are read:
%   row-vector     [x y z 1] * M    last column [0;0;0;1] (NIfTI Transform.T,
%                                   maketform tdata.T, PRESTUS grid transforms)
%   column-vector  M * [x;y;z;1]    last row [0 0 0 1] (SimNIBS/FSL text
%                                   matrices, Transform.T')
% A matrix that is both (no translation) gives the same inverse either way.
% Every 4x4 affine inversion in PRESTUS goes through here (RAS_TO_VOXEL for
% RAS mm to voxel indices).
%
% Use as:
%   Minv = affine_inverse(M)
%
% Input:
%   M    - (4,4) affine; a maketform struct (tdata.T) is also accepted
%
% Output:
%   Minv - (4,4) double inverse, same convention as M
%
% Errors:
%   PRESTUS:transform:badAffine       neither last row nor last column is [0 0 0 1]
%   PRESTUS:transform:singularAffine  the linear block cannot be inverted reliably
%
% See also: RAS_TO_VOXEL, AFFINE_APPLY_PTS, AFFINE_RESAMPLE_3D

arguments
    M
end
    if isstruct(M) && isfield(M, 'tdata') && isfield(M.tdata, 'T')
        M = M.tdata.T;
    end
    if ~(isnumeric(M) && isequal(size(M), [4 4]) && all(isfinite(M(:))))
        error('PRESTUS:transform:badAffine', 'Expected a finite 4x4 affine matrix.');
    end
    M = double(M);
    tol = 1e-9 * max(1, max(abs(M(:))));
    row_vector = all(abs(M(:, 4) - [0; 0; 0; 1]) <= tol);
    col_vector = all(abs(M(4, :) - [0 0 0 1]) <= tol);
    if ~row_vector && ~col_vector
        error('PRESTUS:transform:badAffine', ...
            'Not an affine: neither the last column nor the last row is [0 0 0 1].');
    end

    if row_vector
        A = M(1:3, 1:3); t = M(4, 1:3);           % [x 1]*M = [x*A + t, 1]
    else
        A = M(1:3, 1:3)'; t = M(1:3, 4)';         % same map, written row-wise
    end
    if rcond(A) < 1e3 * eps
        error('PRESTUS:transform:singularAffine', ...
            'The affine''s linear part is singular or near singular (rcond %.1e).', rcond(A));
    end
    Ainv = A \ eye(3);
    Minv = eye(4);
    Minv(1:3, 1:3) = Ainv;
    Minv(4, 1:3) = -t * Ainv;
    if ~row_vector
        Minv = Minv';
    end
end
