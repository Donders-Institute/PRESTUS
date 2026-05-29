function bounds = affine_bounds(M, corners)
% AFFINE_BOUNDS  Compute output bounding box after applying an affine transform
%
% Drop-in replacement for:
%   findbounds(maketform('affine', M), corners)
%
% Transforms all 8 corners of the input bounding box and returns the
% axis-aligned bounding box in the output space.
%
% Use as:
%   bounds = affine_bounds(M, corners)
%
% Input:
%   M       - 4×4 affine (row-vector convention) or maketform struct
%   corners - [2×3] = [min_corner; max_corner] of the input bounding box
%
% Output:
%   bounds - [2×3] = [min_out; max_out] bounding box in output space
%
% See also: AFFINE_RESAMPLE_3D, AFFINE_APPLY_PTS

    if isstruct(M) && isfield(M, 'tdata') && isfield(M.tdata, 'T')
        M = M.tdata.T;
    end

    lo = corners(1,:);
    hi = corners(2,:);
    eight = [lo(1) lo(2) lo(3);
             hi(1) lo(2) lo(3);
             lo(1) hi(2) lo(3);
             hi(1) hi(2) lo(3);
             lo(1) lo(2) hi(3);
             hi(1) lo(2) hi(3);
             lo(1) hi(2) hi(3);
             hi(1) hi(2) hi(3)];

    out = affine_apply_pts(eight, M);
    bounds = [min(out, [], 1); max(out, [], 1)];
end
