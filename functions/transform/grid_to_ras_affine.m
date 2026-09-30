function T = grid_to_ras_affine(planimg, fov_offset_1based)
% GRID_TO_RAS_AFFINE  Build the acoustic-grid-voxel -> RAS-mm affine
%
% Composes the head-alignment transform (planimg.inv_transf: full simulation
% grid voxel -> T1 voxel) with the T1 NIfTI's own voxel -> RAS affine
% (planimg.t1_header), and optionally an acoustic-FOV-crop voxel offset —
% needed whenever the array being exported is a FOV-cropped sub-grid rather
% than the full simulation grid (e.g. acoustic_fov_crop.m's fov_offset_ac).
%
% Row-vector convention throughout (matching affine_resample_3d /
% affine_apply_pts): [ras_mm 1] = [grid_vox 1] * T, with grid_vox 1-based.
% A NIfTI header for an array on this grid needs 0-based indices: use
% NIFTI_AFFINE_FROM_1BASED(T).
%
% Use as:
%   T = grid_to_ras_affine(planimg)                    % full simulation grid
%   T = grid_to_ras_affine(planimg, fov_offset_1based)  % FOV-cropped sub-grid
%
% Input:
%   planimg            - struct with .inv_transf [4x4] (full-grid-vox -> T1-vox,
%                         1-based) and .t1_header (affine3d, T1-vox -> RAS mm)
%   fov_offset_1based   - (optional) [1x3] 1-based start index of the cropped
%                         sub-grid within the full simulation grid
%                         (acoustic_grid_fov.m's fov_offset_ac). Omit or pass
%                         [1 1 1] when the array covers the full grid.
%
% Output:
%   T - [4x4] affine, row-vector convention: [ras_mm 1] = [grid_vox 1] * T
%
% See also: WRITE_BABELBRAIN_INPUT, ACOUSTIC_GRID_FOV, PREPROC_HEAD

arguments
    planimg            (1,1) struct
    fov_offset_1based   (1,3) double = [1 1 1]
end

M_crop_to_full = eye(4);
M_crop_to_full(4, 1:3) = fov_offset_1based(:)' - 1;

% inv_transf yields 1-based T1 voxels; the NIfTI affine takes 0-based ones
M_1based_to_0based = eye(4);
M_1based_to_0based(4, 1:3) = -1;

T = M_crop_to_full * planimg.inv_transf * M_1based_to_0based * planimg.t1_header.Transform.T;

% Chaining matrix products (incl. an inv() inside planimg.inv_transf) leaves
% floating-point noise (~1e-14) in what must be an exact [0;0;0;1] last
% column for affine3d to accept the matrix — clean it up.
T(1:3, 4) = 0;
T(4, 4) = 1;

end
