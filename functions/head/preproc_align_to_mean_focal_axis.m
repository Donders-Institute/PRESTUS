function [rotated_img, trans_pos_new, focus_pos_new, transformation_matrix, rotation_matrix, angle_x_rad, angle_y_rad, montage_img] = ...
    preproc_align_to_mean_focal_axis(nii_image, nii_header, trans_pos_grid, focus_pos_grid, scale_factor, parameters)

% PREPROC_ALIGN_TO_MEAN_FOCAL_AXIS  Rotate a 3D image so the mean focal axis aligns with z
%
% Computes a rotation that aligns the mean focal unit vector of all
% configured transducers with the grid z-axis. This ensures z corresponds
% to the beam convergence direction for any transducer geometry:
%
%   N=1  : z = focal axis of transducer(1)  (honours align_to_focus flag)
%   N=2+ : z = normalize( mean of the focal unit vectors, each taken as an
%          axis: a vector pointing against transducer(1)'s is flipped first )
%
% Taking the beams as axes (lines, not directions) handles both geometries
% with one rule: beams converging from the same side average to their
% bisector, and opposed (bilateral) beams average to the bisector of
% transducer(1)'s axis and the reversed axis of the opposite transducer,
% so each beam is tilted from z by about half the angle between the axes.
% Plain direction averaging would cancel opposed beams. z points along
% transducer(1)'s beam.
%
% Guard: if norm(mean axis vector) < parameters.grid.focal_axis_min_norm
% (default 0.5; for two transducers the norm is always >= cos(45 deg) =
% 0.71, so this only triggers for three or more widely spread beams), a
% warning is issued and the function falls back to the focal axis of
% transducer(1) alone.
%
% Use as:
%   [rotated_img, trans_pos_new, focus_pos_new, transformation_matrix, ...
%    rotation_matrix, angle_x_rad, angle_y_rad, montage_img] = ...
%       preproc_align_to_mean_focal_axis(nii_image, nii_header, ...
%           trans_pos_grid, focus_pos_grid, scale_factor, parameters)
%
% Input:
%   nii_image      - [Nx x Ny x Nz] 3D image volume
%   nii_header     - NIfTI header struct (from niftiinfo)
%   trans_pos_grid - [1x3] transducer(1) position in T1-voxel coordinates
%   focus_pos_grid - [1x3] focus(1) position in T1-voxel coordinates
%   scale_factor   - image scaling factor (T1_res / grid_res; 1 = no scaling)
%   parameters     - (1,1) simulation parameters struct
%
% Output:
%   rotated_img           - rotated and scaled image volume
%   trans_pos_new         - [1x3] transducer(1) position after transformation
%   focus_pos_new         - [1x3] focus(1) position after transformation
%   transformation_matrix - [4x4] combined affine (rotate + scale + recenter)
%   rotation_matrix       - [4x4] rotation-only matrix
%   angle_x_rad           - rotation angle around x-axis [rad]
%   angle_y_rad           - rotation angle around y-axis [rad]
%   montage_img           - side-by-side comparison montage
%
% See also: PREPROC_ALIGN_TO_FOCAL_AXIS, PREPROC_ALIGN_TO_PRINCIPAL_AXIS, PREPROC_HEAD

    arguments
        nii_image (:,:,:) {mustBeNumericOrLogical}
        nii_header     (1,1) struct
        trans_pos_grid (1,3) double
        focus_pos_grid (1,3) double
        scale_factor   (1,1) double
        parameters     (1,1) struct
    end

    nii_image      = gather(nii_image);
    trans_pos_grid = gather(trans_pos_grid);
    focus_pos_grid = gather(focus_pos_grid);

    %% Step 1: Compute the mean focal axis across all transducers (beams as axes)

    min_norm = 0.5; % guard threshold (see MEAN_BEAM_AXIS)
    if isfield(parameters, 'grid') && isfield(parameters.grid, 'focal_axis_min_norm')
        min_norm = parameters.grid.focal_axis_min_norm;
    end

    n_trans = numel(parameters.transducer);

    % Focal unit vector for transducer(1) (honour align_to_focus flag)
    if ~isfield(parameters.transducer(1), 'align_to_focus') || parameters.transducer(1).align_to_focus
        raw1 = focus_pos_grid - trans_pos_grid;
    else
        curv_radius = parameters.transducer(1).(parameters.transducer(1).type).curv_radius_mm;
        natural_focus = trans_pos_grid + [0, 0, curv_radius / parameters.grid.resolution_mm];
        raw1 = natural_focus - trans_pos_grid;
    end
    u1 = raw1 / norm(raw1);

    % Beam unit vectors of transducers 2..N (planned positions, T1 voxels)
    focal_units = u1;
    for ti = 2:n_trans
        u_i = parameters.transducer(ti).focus_pos(:)' - parameters.transducer(ti).trans_pos(:)';
        focal_units(end+1,:) = u_i / norm(u_i); %#ok<AGROW>
    end
    align_vec = mean_beam_axis(focal_units, min_norm);

    focal_axis = [align_vec, 1]';

    %% Steps 2-8: rotation / scale / resample

    % Rotate around y-axis to bring x-component to zero
    angle_y_rad = atan(focal_axis(1) / focal_axis(3));
    if angle_y_rad < 0
        angle_y_rad = angle_y_rad + pi;
    end
    vox_transf_Y = makehgtform('yrotate', -angle_y_rad);

    focal_axis_new = vox_transf_Y * focal_axis;

    % Rotate around x-axis to bring y-component to zero
    angle_x_rad = atan(focal_axis_new(2) / focal_axis_new(3));
    vox_transf_X = makehgtform('xrotate', angle_x_rad);

    rotation_matrix           = vox_transf_X * vox_transf_Y;
    scale_matrix              = makehgtform('scale', scale_factor);
    rotation_and_scale_matrix = rotation_matrix * scale_matrix;

    rotation_center = (size(nii_image) + 1) / 2;
    forward         = makehgtform('translate', -rotation_center);
    TF_centered     = forward' * rotation_and_scale_matrix';

    img_bounds = affine_bounds(TF_centered, [0 0 0; size(nii_image)]);
    tp_trans   = affine_apply_pts(trans_pos_grid, TF_centered);
    fp_trans   = affine_apply_pts(focus_pos_grid, TF_centered);

    % Include all transducer and focus positions so none ends up outside the grid
    extra_pts = zeros(0, 3);
    if n_trans >= 2
        for ti = 2:n_trans
            tp_i = parameters.transducer(ti).trans_pos(:)';
            fp_i = parameters.transducer(ti).focus_pos(:)';
            extra_pts(end+1,:) = affine_apply_pts(tp_i, TF_centered); %#ok<AGROW>
            extra_pts(end+1,:) = affine_apply_pts(fp_i, TF_centered); %#ok<AGROW>
        end
    end
    all_points = [img_bounds; tp_trans; fp_trans; extra_pts];

    max_abs = max(abs(all_points), [], 1);
    newdims  = ceil(2 * max_abs) + 2;

    backward = makehgtform('translate', (newdims + 1) / 2);
    transformation_matrix = forward' * rotation_and_scale_matrix' * backward';

    if numel(unique(nii_image)) < 20
        parameters.interpolation = 'nearest';
    else
        parameters.interpolation = 'linear';
    end
    rotated_img = affine_resample_3d(nii_image, transformation_matrix, newdims, parameters.interpolation, 0);

    out_mat       = round(affine_apply_pts([trans_pos_grid; focus_pos_grid], transformation_matrix));
    trans_pos_new = out_mat(1,:);
    focus_pos_new = out_mat(2,:);

    orig_with_transducer_img    = plot_t1_with_transducer(nii_image, nii_header.PixelDimensions(1), ...
        trans_pos_grid, focus_pos_grid, parameters);
    rotated_with_transducer_img = plot_t1_with_transducer(rotated_img, ...
        nii_header.PixelDimensions(1) / scale_factor, trans_pos_new, focus_pos_new, parameters);
    montage_img = imtile({orig_with_transducer_img, rotated_with_transducer_img}, ...
        'GridSize', [1 nan], 'BackgroundColor', [0.5 0.5 0.5]);
end
