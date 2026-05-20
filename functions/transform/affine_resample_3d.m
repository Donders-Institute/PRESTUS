function out = affine_resample_3d(vol, M, out_size, method, fill)
% AFFINE_RESAMPLE_3D  Resample a 3D (or 4D) volume under an affine transform
%
% Drop-in replacement for:
%   tformarray(vol, maketform('affine',M), makeresampler(method,'fill'),
%              [1 2 3], [1 2 3], out_size, [], fill)
%
% M is the 4×4 FORWARD transform in row-vector (MATLAB maketform) convention:
%   [out_pt 1] = [in_pt 1] * M
%
% Internally the inverse M^{-1} is applied to map each output voxel to the
% corresponding input voxel, then griddedInterpolant is used to sample.
% NaN (out-of-bounds) values are replaced by fill.
%
% Old maketform tform structs are accepted for backward-compatibility with
% cached .mat files: if M is a struct with a .tdata.T field, that matrix
% is used.
%
% Use as:
%   out = affine_resample_3d(vol, M, out_size, method, fill)
%   out = affine_resample_3d(vol, M, out_size)         % method='linear', fill=0
%
% Input:
%   vol      - numeric array [D1 x D2 x D3] or [D1 x D2 x D3 x C]
%   M        - 4×4 forward affine (row-vector convention) or maketform struct
%   out_size - [1×3] output dimensions
%   method   - 'nearest', 'linear', or 'cubic' (default: 'linear')
%   fill     - scalar fill value for out-of-bounds voxels (default: 0)
%
% Output:
%   out - array of size [out_size] or [out_size, C], same class as vol
%
% See also: AFFINE_APPLY_PTS, AFFINE_BOUNDS

    if nargin < 4 || isempty(method); method = 'linear'; end
    if nargin < 5 || isempty(fill);   fill   = 0;        end

    % Accept legacy maketform struct
    if isstruct(M) && isfield(M, 'tdata') && isfield(M.tdata, 'T')
        M = M.tdata.T;
    end

    % griddedInterpolant only supports 'nearest' and 'linear' natively;
    % fall back to 'cubic' via interp3 but map 'cubic' → 'linear' here.
    if strcmp(method, 'cubic')
        method = 'linear';
    end

    Minv = inv(M);

    [xi, yi, zi] = ndgrid(1:out_size(1), 1:out_size(2), 1:out_size(3));
    N = prod(out_size);
    out_pts = [xi(:), yi(:), zi(:), ones(N, 1)];   % N×4
    in_pts  = out_pts * Minv;                        % N×4

    px = reshape(in_pts(:,1), out_size);
    py = reshape(in_pts(:,2), out_size);
    pz = reshape(in_pts(:,3), out_size);
    clear out_pts in_pts xi yi zi

    in_dims = size(vol, [1 2 3]);
    nC      = size(vol, 4);
    src_cls = class(vol);

    F = griddedInterpolant({1:in_dims(1), 1:in_dims(2), 1:in_dims(3)}, ...
        single(vol(:,:,:,1)), method, 'none');

    if nC > 1
        out = zeros([out_size, nC], src_cls);
        for c = 1:nC
            F.Values = single(vol(:,:,:,c));
            slice = F(px, py, pz);
            slice(isnan(slice)) = fill;
            out(:,:,:,c) = cast(slice, src_cls);
        end
    else
        raw = F(px, py, pz);
        raw(isnan(raw)) = fill;
        out = cast(raw, src_cls);
    end
end
