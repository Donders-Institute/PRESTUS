function [img_out, hdr_out, vox_map] = ensure_ras_plus(img_in, hdr_in)
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
% (as fslreorient2std). Every voxel keeps its world position: an oblique
% affine keeps its rotation, in the new axis order.
%
% Use as:
%   [img_out, hdr_out, vox_map] = ensure_ras_plus(img_in, hdr_in)
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
%   vox_map - [4x4] affine (row-vector convention, see AFFINE_APPLY_PTS)
%             mapping 1-based voxel indices of img_in to the corresponding
%             1-based voxel indices of img_out; eye(4) if nothing changed.
%             Use it to carry voxel coordinates defined on the input grid
%             (e.g. configured transducer positions) into the RAS+ grid.
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
        vox_map = eye(4);
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
    % Pixel dimensions in the new axis order.
    new_pixdims = hdr_in.PixelDimensions(perm);

    % Build new column-vector affine: the old axis vectors, permuted and with
    % flipped axes negated, so an oblique image keeps its rotation (diagonal,
    % positive for an image without one); origin from the old affine
    new_R = R;
    new_R(:, flip_dims) = -new_R(:, flip_dims);

    % Origin: RAS of the new array's first voxel. The NIfTI affine takes
    % 0-based indices, so find that voxel's 0-based index in the old array.
    corner_new = zeros(1,3);
    for d = flip_dims(:)'
        corner_new(d) = sz(d) - 1;   % flipped axis: first voxel was the last one
    end
    corner_old = zeros(1,3);
    corner_old(perm) = corner_new;   % new axis n is old axis perm(n)
    R_orig = T(1:3, 1:3)';
    origin_ras = R_orig * corner_old' + T(4, 1:3)';

    new_T        = eye(4);
    new_T(1:3,1:3) = new_R';    % MATLAB convention: row-vector * T
    new_T(4, 1:3)  = origin_ras';

    % Voxel-index map (1-based): new axis d reads old axis perm(d), reversed
    % when flipped, i.e. new_d = old_perm(d) or sz(d) + 1 - old_perm(d).
    vox_map = zeros(4);
    vox_map(4,4) = 1;
    for d = 1:3
        if ismember(d, flip_dims)
            vox_map(perm(d), d) = -1;
            vox_map(4, d)       = size(img_out, d) + 1;
        else
            vox_map(perm(d), d) = 1;
        end
    end

    hdr_out = hdr_in;
    hdr_out.Transform        = affine3d(new_T);
    hdr_out.ImageSize        = sz(1:3);
    hdr_out.PixelDimensions  = new_pixdims;

    if ~isequal(perm, [1 2 3]) || ~isempty(flip_dims)
        fprintf('[ensure_ras_plus] Reoriented: perm=[%d %d %d], flipped dims=[%s]\n', ...
            perm(1), perm(2), perm(3), num2str(flip_dims(:)'));
    end
end
