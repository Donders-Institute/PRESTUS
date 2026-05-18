function [img_out, hdr_out] = ensure_ras_plus(img_in, hdr_in)
% ENSURE_RAS_PLUS  Reorient a NIfTI volume to RAS+ if needed
%
% Inspects the NIfTI affine (hdr_in.Transform.T) and, when the image is not
% already in RAS+ orientation, applies the minimal combination of axis
% permutations and flips to bring the voxel array into RAS+ and updates the
% affine accordingly. If the image is already RAS+ the inputs are returned
% unchanged.
%
% "RAS+" means:
%   dim 1 increases toward Right   (+X world)
%   dim 2 increases toward Anterior (+Y world)
%   dim 3 increases toward Superior (+Z world)
%
% The function does NOT interpolate — it only permutes axes and flips them
% (equivalent to fslreorient2std for data whose affine has no off-diagonal
% shear terms beyond axis permutation).
%
% Use as:
%   [img_out, hdr_out] = ensure_ras_plus(img_in, hdr_in)
%
% Input:
%   img_in  - [Nx x Ny x Nz] or [Nx x Ny x Nz x ...] numeric/logical array
%   hdr_in  - niftiinfo header struct with Transform.T (4x4 affine, MATLAB
%             convention: T maps voxel row-vectors to RAS mm via v*T)
%
% Output:
%   img_out - reoriented image (same data, possibly different axis order/direction)
%   hdr_out - updated header: Transform.T and ImageSize updated to match img_out;
%             PixelDimensions updated if axes were permuted
%
% See also: NIFTIINFO, NIFTIREAD, RAS_TO_GRID

arguments
    img_in  {mustBeNumericOrLogical}
    hdr_in  (1,1) struct
end

    % MATLAB niftiinfo stores the affine as T such that:
    %   [x_ras, y_ras, z_ras, 1] = [i, j, k, 1] * T
    % i.e. T is the TRANSPOSE of the conventional column-vector affine.
    % The rotation sub-matrix is T(1:3, 1:3)'.
    T = hdr_in.Transform.T;         % [4x4], MATLAB convention
    R = T(1:3, 1:3)';               % column-vector affine rotation [3x3]

    % Determine which world axis each voxel axis most closely aligns with.
    % abs(R): row = world axis (X/Y/Z), col = voxel axis (i/j/k)
    [~, vox_for_world] = max(abs(R), [], 2);  % [3x1]: for each world axis, which voxel axis
    [~, world_for_vox] = max(abs(R), [], 1);  % [1x3]: for each voxel axis, which world axis

    % Desired permutation: voxel axis that most aligns with R(right), A(ant), S(sup)
    perm = vox_for_world';   % [1x3]: new dim order so that dim1→R, dim2→A, dim3→S

    if isequal(perm, [1 2 3]) && all(diag(R) > 0)
        % Already RAS+, nothing to do
        img_out = img_in;
        hdr_out = hdr_in;
        return
    end

    % --- Step 1: permute axes so R/A/S land on dims 1/2/3 ---
    ndims_extra = ndims(img_in) - 3;
    if ndims_extra > 0
        full_perm = [perm, 4:(3 + ndims_extra)];
    else
        full_perm = perm;
    end

    if ~isequal(perm, [1 2 3])
        img_out = permute(img_in, full_perm);
        % Permute columns of R accordingly (columns = voxel axes)
        R = R(:, perm);
    else
        img_out = img_in;
    end

    % --- Step 2: flip any axis where the world direction is negative ---
    sz = size(img_out);
    flip_dims = find(diag(R) < 0);
    for d = flip_dims(:)'
        img_out = flip(img_out, d);
    end

    % --- Step 3: rebuild affine ---
    % After permutation + flip the rotation is the identity (pure isotropic RAS+).
    % Pixel dimensions in the new axis order.
    new_pixdims = hdr_in.PixelDimensions(perm);

    % Build new column-vector affine: diagonal scale + origin at voxel (1,1,1)
    % We cannot recover the original world origin faithfully here (it depended
    % on the full rotation), so we set origin to zero — same convention as
    % canonical_affine_transform.
    new_R = diag(new_pixdims);   % 3x3 diagonal, always positive after flips

    % Translate origin: old voxel (1,1,1) → RAS mm.
    % After flips the corner that maps to RAS origin may have moved; compute it.
    orig_corner = ones(1,3);
    for d = flip_dims(:)'
        orig_corner(d) = sz(d);   % flipped axis: voxel 1 was formerly voxel N
    end
    % Map original corner through original R (before permutation reorder)
    R_orig = T(1:3, 1:3)';
    origin_ras = R_orig * orig_corner(perm)' + T(4, 1:3)';

    new_T        = eye(4);
    new_T(1:3,1:3) = new_R';    % MATLAB convention: row-vector * T
    new_T(4, 1:3)  = origin_ras';

    hdr_out = hdr_in;
    hdr_out.Transform        = affineTransform3d(new_T);
    hdr_out.ImageSize        = sz(1:3);
    hdr_out.PixelDimensions  = new_pixdims;

    if ~isequal(perm, [1 2 3]) || ~isempty(flip_dims)
        fprintf('[ensure_ras_plus] Reoriented: perm=[%d %d %d], flipped dims=[%s]\n', ...
            perm(1), perm(2), perm(3), num2str(flip_dims(:)'));
    end
end
