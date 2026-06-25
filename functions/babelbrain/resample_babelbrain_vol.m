function vol_out = resample_babelbrain_vol(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm, focus_pos, focal_length_mm)
% RESAMPLE_BABELBRAIN_VOL  Resample a BabelBrain volume onto the PRESTUS acoustic grid.
%
% Builds the BabelBrain-voxel → PRESTUS-grid affine in row-vector convention
% (as required by affine_resample_3d: [out 1] = [in 1] * M) and resamples
% using linear interpolation with zero fill for out-of-bounds voxels.
%
% Callers (read_babelbrain_output / read_babelbrain_thermal) permute the raw
% h5read array from [N_z, N_y, N_x] to PRESTUS axis order [N_x, N_y, N_z]
% BEFORE calling this function, so the dimension mapping is:
%   MATLAB dim 1 = BabelBrain x = PRESTUS lateral X (0 at focal spot)
%   MATLAB dim 2 = BabelBrain y = PRESTUS lateral Y (0 at focal spot)
%   MATLAB dim 3 = BabelBrain z = PRESTUS axial/depth (FocalLength at focal spot)
%
% The direct BB-voxel → PRESTUS-voxel mapping:
%   prestus_dim1 = x_mm(i) / res_mm + focus_pos(1)
%   prestus_dim2 = y_mm(j) / res_mm + focus_pos(2)
%   prestus_dim3 = (z_mm(k) - focal_length_mm) / res_mm + focus_pos(3)
%
% Input:
%   vol_bb          - [Nx x Ny x Nz] volume, already permuted to PRESTUS axis order
%   x_mm            - [1 x Nx] lateral-X coords (PRESTUS dim 1), [mm], 0 at focal spot
%   y_mm            - [1 x Ny] lateral-Y coords (PRESTUS dim 2), [mm], 0 at focal spot
%   z_mm            - [1 x Nz] depth coords    (PRESTUS dim 3), [mm], FocalLength at focal spot
%                     Decreasing (flip-corrected BabelBrain z_vec).
%   planimg         - planimg struct (kept for signature compatibility; unused)
%   prestus_dims    - [1 x 3] output grid dimensions
%   res_mm          - PRESTUS acoustic grid resolution [mm]
%   focus_pos       - [1 x 3] focal-spot position in PRESTUS grid (1-based voxels)
%   focal_length_mm - transducer preset FocalLength [mm]; z_mm at focal-spot voxel = this.
%
% Output:
%   vol_out - [prestus_dims] double, resampled volume
%
% See also: READ_BABELBRAIN_OUTPUT, READ_BABELBRAIN_THERMAL, AFFINE_RESAMPLE_3D

    dx_bb = mean(diff(x_mm));   % lateral X step, positive
    dy_bb = mean(diff(y_mm));   % lateral Y step, positive
    dz_bb = mean(diff(z_mm));   % depth step, negative (flip-corrected z_vec decreases)

    % Row-vector affine [in_vox 1] * M = [out_vox 1].
    % Anchor: focal-spot voxel (x_mm=0, y_mm=0, z_mm=focal_length_mm) → focus_pos.

    s1 = dx_bb / res_mm;   % scale MATLAB dim 1 (lateral X, positive)
    s2 = dy_bb / res_mm;   % scale MATLAB dim 2 (lateral Y, positive)
    s3 = dz_bb / res_mm;   % scale MATLAB dim 3 (depth, negative)

    t1 =  x_mm(1)                    / res_mm + focus_pos(1) - s1;
    t2 =  y_mm(1)                    / res_mm + focus_pos(2) - s2;
    t3 = (z_mm(1) - focal_length_mm) / res_mm + focus_pos(3) - s3;

    M_bb_to_acgrid = [s1,  0,   0,  0;
                       0,  s2,  0,  0;
                       0,   0, s3,  0;
                      t1,  t2, t3,  1];

    vol_out = double(affine_resample_3d(single(vol_bb), M_bb_to_acgrid, prestus_dims, 'linear', 0));
end
