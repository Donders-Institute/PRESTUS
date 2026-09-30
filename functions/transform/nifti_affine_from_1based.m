function T_nifti = nifti_affine_from_1based(T_1based)
% NIFTI_AFFINE_FROM_1BASED  NIfTI header affine from a 1-based-index affine
%
% PRESTUS composes grid transforms on 1-based MATLAB array indices, while a
% NIfTI header (niftiinfo/niftiwrite Transform) maps 0-based indices to RAS.
% Prepends the 0-based -> 1-based index shift so the result can be written
% as hdr.Transform = affine3d(T_nifti).
%
% Use as:
%   T_nifti = nifti_affine_from_1based(T_1based)
%
% Input:
%   T_1based - [4x4] affine, row-vector convention: [ras 1] = [vox1 1] * T
%
% Output:
%   T_nifti  - [4x4] affine, row-vector convention: [ras 1] = [vox0 1] * T
%
% See also: GRID_TO_RAS_AFFINE, GRID_TO_RAS, RAS_TO_GRID

arguments
    T_1based (4,4) double
end

    M_0based_to_1based = eye(4);
    M_0based_to_1based(4, 1:3) = 1;
    T_nifti = M_0based_to_1based * T_1based;
    T_nifti(1:3, 4) = 0;
    T_nifti(4, 4) = 1;
end
