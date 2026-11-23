function vox = ras_to_voxel(ras_pos, affine, opts)
% RAS_TO_VOXEL  RAS world coordinates to 1-based voxel indices, fractional or rounded
%
% The one place where a NIfTI affine is inverted to go from RAS mm to array
% indices (RAS_TO_GRID calls it). niftiinfo's Transform.T maps 0-based voxel
% indices to RAS in MATLAB's row-vector convention, [ras 1] = [vox0 1] * T,
% so the inverse of T (AFFINE_INVERSE) gives the 0-based index and 1 is added
% for MATLAB arrays.
% The inverse is GRID_TO_RAS.
%
% Use as:
%   vox = ras_to_voxel(ras_pos, nii_header)
%   vox = ras_to_voxel(ras_pos, T)
%   vox = ras_to_voxel(ras_pos, T, Round=true)
%
% Input:
%   ras_pos - [Nx3] positions in RAS world coordinates [mm] (one per row)
%   affine  - niftiinfo header with Transform.T, or the [4x4] Transform.T itself
%   Round   - true: nearest voxel; false (default): fractional index
%
% Output:
%   vox     - [Nx3] 1-based voxel indices
%
% See also: RAS_TO_GRID, GRID_TO_RAS, AFFINE_INVERSE, TRANSFORM_COORDINATES

arguments
    ras_pos (:,3) {mustBeNumeric}
    affine
    opts.Round (1,1) logical = false
end

    if isstruct(affine)
        T = affine.Transform.T;
    else
        T = affine;
    end
    if ~isequal(size(T), [4 4]) || any(abs(T(:, 4)' - [0 0 0 1]) > 1e-9)
        error('PRESTUS:transform:badAffine', ...
            'Expected a 4x4 affine in row-vector convention (last column [0 0 0 1]).');
    end

    % [vox0 1] = [ras 1] * inv(T); AFFINE_INVERSE also rejects a singular T
    vox0 = [double(ras_pos), ones(size(ras_pos, 1), 1)] * affine_inverse(T);
    vox0 = vox0(:, 1:3);
    if opts.Round
        vox0 = round(vox0);   % on the 0-based index, as RAS_TO_GRID always did
    end
    vox = vox0 + 1;
end
