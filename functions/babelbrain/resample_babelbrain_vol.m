function vol_out = resample_babelbrain_vol(vol_bb, x_mm, y_mm, z_mm, planimg, prestus_dims, res_mm, focus_pos, focal_length_mm, method, beam_sign)
% RESAMPLE_BABELBRAIN_VOL  Resample a BabelBrain volume onto the PRESTUS acoustic grid
%
% Builds the BabelBrain-voxel → PRESTUS-grid affine in row-vector convention
% (as required by affine_resample_3d: [out 1] = [in 1] * M) and resamples
% with zero fill for out-of-bounds voxels. Interpolation defaults to 'linear'
% (for continuous fields like pressure).
%
% For integer LABEL volumes (e.g. BabelBrain's MaterialMap) pass:
%   'nearest' - point-sampled nearest-neighbour (fast; fine for UPsampling).
%   'mode'    - majority vote: each BabelBrain-native voxel is forward-scattered
%               into the PRESTUS cell it lands in, and each output cell takes the
%               most common label among the native voxels that fell into it.
%               Use this when DOWNsampling (BabelBrain-native finer than the
%               PRESTUS grid), where point-sampled 'nearest' can drop thin
%               structures (e.g. sub-voxel cortical bone) that happen to miss
%               the single sample point even though they occupy real volume in
%               the coarse cell.
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
% Use as:
%   vol_out = resample_babelbrain_vol(vol_bb, x_mm, y_mm, z_mm, planimg, ...
%                 prestus_dims, res_mm, focus_pos, focal_length_mm)
%   vol_out = resample_babelbrain_vol(..., method, beam_sign)
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
%   focal_length_mm - z anchor [mm]: the z_mm value mapped onto focus_pos(3)
%                     (z of the /TargetLocation marker, or preset FocalLength as
%                     fallback)
%   method          - (optional) 'linear' (default) | 'nearest' | 'mode'
%   beam_sign       - (optional, default +1) direction of the beam axis in the
%                     PRESTUS grid: +1 when physical depth (into tissue past the
%                     focus) increases with dim 3 (focus_pos(3) > trans_pos(3)),
%                     −1 when it decreases (focus_pos(3) < trans_pos(3)). Only the
%                     depth mapping (dim 3) is affected; the focal anchor is fixed
%                     regardless. Must match the input-side depth orientation
%                     chosen in write_babelbrain_input (same focus_pos vs
%                     trans_pos comparison), or the returned field is mirrored
%                     along the beam axis.
%
% Output:
%   vol_out - [prestus_dims] double, resampled volume
%
% See also: READ_BABELBRAIN_OUTPUT, READ_BABELBRAIN_THERMAL, AFFINE_RESAMPLE_3D

    if nargin < 10 || isempty(method)
        method = 'linear';
    end
    if nargin < 11 || isempty(beam_sign)
        beam_sign = 1;
    end

    dx_bb = mean(diff(x_mm));   % lateral X step, positive
    dy_bb = mean(diff(y_mm));   % lateral Y step, positive
    dz_bb = mean(diff(z_mm));   % depth step, negative (flip-corrected z_vec decreases)

    % Row-vector affine [in_vox 1] * M = [out_vox 1].
    % Anchor: focal-spot voxel (x_mm=0, y_mm=0, z_mm=focal_length_mm) → focus_pos.
    % beam_sign flips only the depth (dim 3) mapping about the focal anchor, so a
    % −z beam is placed the right way round; lateral axes are unaffected.

    s1 = dx_bb / res_mm;               % scale MATLAB dim 1 (lateral X, positive)
    s2 = dy_bb / res_mm;               % scale MATLAB dim 2 (lateral Y, positive)
    s3 = beam_sign * dz_bb / res_mm;   % scale MATLAB dim 3 (depth), sign per beam direction

    t1 =  x_mm(1)                                 / res_mm + focus_pos(1) - s1;
    t2 =  y_mm(1)                                 / res_mm + focus_pos(2) - s2;
    t3 = beam_sign * (z_mm(1) - focal_length_mm)  / res_mm + focus_pos(3) - s3;

    M_bb_to_acgrid = [s1,  0,   0,  0;
                       0,  s2,  0,  0;
                       0,   0, s3,  0;
                      t1,  t2, t3,  1];

    if strcmpi(method, 'mode')
        vol_out = mode_downsample(vol_bb, s1, s2, s3, t1, t2, t3, prestus_dims);
    else
        vol_out = double(affine_resample_3d(single(vol_bb), M_bb_to_acgrid, prestus_dims, method, 0));
    end
end

% -------------------------------------------------------------------------
function vol_out = mode_downsample(vol_bb, s1, s2, s3, t1, t2, t3, out_dims)
% Majority-vote label downsample by forward-scattering each input voxel into
% the output cell it maps to, then taking the most common label per output
% cell. The BB→PRESTUS map is the diagonal scale+translation encoded by
% (s1,s2,s3,t1,t2,t3): out_dim_d = in_idx_d * s_d + t_d (1-based).
%
% Out-of-range input voxels (mapping outside [1, out_dims]) are dropped, so
% output cells that no input voxel reaches stay 0 — matching the zero-fill of
% the interpolating path.
    [n1, n2, n3] = size(vol_bb);
    [I1, I2, I3] = ndgrid(1:n1, 1:n2, 1:n3);

    o1 = round(I1(:) * s1 + t1);
    o2 = round(I2(:) * s2 + t2);
    o3 = round(I3(:) * s3 + t3);
    lab = double(vol_bb(:));

    in_bounds = o1 >= 1 & o1 <= out_dims(1) & ...
                o2 >= 1 & o2 <= out_dims(2) & ...
                o3 >= 1 & o3 <= out_dims(3);
    o1 = o1(in_bounds); o2 = o2(in_bounds); o3 = o3(in_bounds);
    lab = lab(in_bounds);

    vol_out = zeros(out_dims);
    if isempty(lab)
        return;   % no overlap between BabelBrain domain and PRESTUS grid
    end

    % Unique output cells and a per-vote group index into them.
    out_lin = sub2ind(out_dims, o1, o2, o3);
    [cell_lin, ~, grp] = unique(out_lin);

    % Per (cell, label) vote counts via accumarray, then argmax label per cell.
    labels = unique(lab);
    [~, lab_idx] = ismember(lab, labels);
    counts = accumarray([grp, lab_idx], 1, [numel(cell_lin), numel(labels)]);
    [~, win_idx] = max(counts, [], 2);   % ties: lowest label index (stable)
    win_label = labels(win_idx);

    vol_out(cell_lin) = win_label;
end
