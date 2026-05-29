function out = affine_apply_pts(pts, M)
% AFFINE_APPLY_PTS  Apply a 4×4 affine to a set of 3-D points (row-vector convention)
%
% Drop-in replacement for:
%   tformfwd(pts, maketform('affine', M))
%
% M is the 4×4 forward transform in MATLAB's row-vector convention:
%   [out_pt 1] = [in_pt 1] * M
%
% Old maketform tform structs are accepted for backward-compatibility.
%
% Use as:
%   out = affine_apply_pts(pts, M)
%
% Input:
%   pts - [N×3] input points
%   M   - 4×4 affine (row-vector convention) or maketform struct
%
% Output:
%   out - [N×3] transformed points
%
% See also: AFFINE_RESAMPLE_3D, AFFINE_BOUNDS

    if isstruct(M) && isfield(M, 'tdata') && isfield(M.tdata, 'T')
        M = M.tdata.T;
    end

    pts_h = [pts, ones(size(pts, 1), 1)];   % N×4
    out_h = pts_h * M;                        % N×4
    out   = out_h(:, 1:3);
end
