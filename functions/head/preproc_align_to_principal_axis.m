function [rotated_img, trans_pos_new, focus_pos_new, transformation_matrix, rotation_matrix, angle_x_rad, angle_y_rad, montage_img] = ...
    preproc_align_to_principal_axis_axis(nii_image, nii_header, trans_pos_grid, focus_pos_grid, scale_factor, parameters)

% PREPROC_ALIGN_TO_PRINCIPAL_AXIS  Rotate a 3D image so the principal axis of transducer positions aligns with z
%
% For multi-transducer (typically bilateral) setups this function computes
% a rotation that aligns the axis from transducer(1) to transducer(2)
% with the z-axis, giving a grid that is symmetric with respect to both
% bowls. The rotation is derived from the unit vector
%
%   alignment_axis = normalize( trans_pos_2 - trans_pos_1 )
%
% where trans_pos_1 is the (T1-space voxel) position passed in as
% trans_pos_grid and trans_pos_2 is read from
% parameters.transducer(2).trans_pos. If fewer than two transducers are
% configured a warning is issued and the function falls back to aligning
% the focal axis of transducer(1) with z (identical to
% PREPROC_ALIGN_TO_FOCAL_AXIS).
%
% The returned transformation matrix, trans_pos_new and focus_pos_new
% always refer to transducer(1) so that the downstream pipeline
% (preproc_head → head_smooth_and_crop → preproc_crop_grid) is unaffected;
% additional transducer positions are mapped separately via
% affine_apply_pts in preproc_head.
%
% Use as:
%   [rotated_img, trans_pos_new, focus_pos_new, transformation_matrix, ...
%    rotation_matrix, angle_x_rad, angle_y_rad, montage_img] = ...
%       preproc_align_to_principal_axis_axis(nii_image, nii_header, ...
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
% See also: PREPROC_ALIGN_TO_FOCAL_AXIS, PREPROC_HEAD

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

    %% Step 1: Compute alignment axis via PCA of all transducer positions
    %
    % N=1: degenerate (zero variance) — silently fall back to focal axis of trans(1).
    % N=2: PCA reduces to normalize(trans_pos_2 - trans_pos_1), identical to principal_axis.
    % N≥3: first principal component of the transducer position cloud.
    %
    % If all transducer positions are coincident (degenerate regardless of N)
    % a warning is issued and the focal axis of trans(1) is used instead.

    n_trans = numel(parameters.transducer);

    if n_trans == 1
        % Single transducer — no inter-transducer axis; use focal axis silently.
        align_vec = focus_pos_grid - trans_pos_grid;
    else
        % Collect all transducer positions (rows = transducers, cols = x/y/z).
        all_trans_pos = zeros(n_trans, 3);
        all_trans_pos(1,:) = trans_pos_grid;
        for ti = 2:n_trans
            all_trans_pos(ti,:) = parameters.transducer(ti).trans_pos(:)';
        end

        % PCA: first principal component = direction of maximum spread.
        centred = all_trans_pos - mean(all_trans_pos, 1);
        [~, ~, V] = svd(centred, 'econ');   % V columns are principal directions
        pc1 = V(:,1)';                       % [1x3] first principal component

        if norm(pc1) < 1e-6
            warn('All transducer positions are coincident — falling back to focal-axis alignment.');
            align_vec = focus_pos_grid - trans_pos_grid;
        else
            % Orient pc1 consistently: ensure it points from trans(1) toward the
            % centroid of the remaining transducers (so z-min = trans(1) side).
            other_centroid = mean(all_trans_pos(2:end,:), 1);
            if dot(pc1, other_centroid - trans_pos_grid) < 0
                pc1 = -pc1;
            end
            align_vec = pc1;
        end
    end

    % Express as homogeneous column vector [x y z 1]'
    focal_axis = [align_vec / norm(align_vec), 1]';

    %% Steps 2–8: identical rotation/scale/recenter logic as preproc_align_to_focal_axis

    % Rotate around y-axis to bring x-component to zero
    angle_y_rad = atan(focal_axis(1) / focal_axis(3));
    if angle_y_rad < 0
        angle_y_rad = angle_y_rad + pi;
    end
    vox_transf_Y = makehgtform('yrotate', -angle_y_rad);

    % Compute rotated axis after y-rotation
    focal_axis_new = vox_transf_Y * focal_axis;

    % Rotate around x-axis to bring y-component to zero
    angle_x_rad = atan(focal_axis_new(2) / focal_axis_new(3));
    vox_transf_X = makehgtform('xrotate', angle_x_rad);

    % Combined rotation and scaling
    rotation_matrix         = vox_transf_X * vox_transf_Y;
    scale_matrix            = makehgtform('scale', scale_factor);
    rotation_and_scale_matrix = rotation_matrix * scale_matrix;

    % Compute new image dimensions
    rotation_center = (size(nii_image) + 1) / 2;
    forward         = makehgtform('translate', -rotation_center);
    T_centered      = forward' * rotation_and_scale_matrix';
    TF_centered     = T_centered;

    img_bounds = affine_bounds(TF_centered, [0 0 0; size(nii_image)]);
    tp_trans   = affine_apply_pts(trans_pos_grid, TF_centered);
    fp_trans   = affine_apply_pts(focus_pos_grid, TF_centered);

    % Include all transducer positions so none ends up outside the padded grid.
    extra_pts = zeros(0, 3);
    if n_trans >= 2
        for ti = 2:n_trans
            tp_i = parameters.transducer(ti).trans_pos(:)';
            extra_pts(end+1,:) = affine_apply_pts(tp_i, TF_centered); %#ok<AGROW>
        end
    end
    all_points = [img_bounds; tp_trans; fp_trans; extra_pts];

    max_abs = max(abs(all_points), [], 1);
    newdims = ceil(2 * max_abs) + 2;

    backward = makehgtform('translate', (newdims + 1) / 2);
    transformation_matrix = forward' * rotation_and_scale_matrix' * backward';

    % Resample image
    if numel(unique(nii_image)) < 20
        parameters.interpolation = 'nearest';
    else
        parameters.interpolation = 'linear';
    end
    rotated_img = affine_resample_3d(nii_image, transformation_matrix, newdims, parameters.interpolation, 0);

    % Update transducer(1) and focus(1) positions
    out_mat = round(affine_apply_pts([trans_pos_grid; focus_pos_grid], transformation_matrix));
    trans_pos_new = out_mat(1,:);
    focus_pos_new = out_mat(2,:);

    % Diagnostic montage
    orig_with_transducer_img = plot_t1_with_transducer(nii_image, nii_header.PixelDimensions(1), ...
        trans_pos_grid, focus_pos_grid, parameters);
    rotated_with_transducer_img = plot_t1_with_transducer(rotated_img, ...
        nii_header.PixelDimensions(1) / scale_factor, trans_pos_new, focus_pos_new, parameters);
    montage_img = imtile({orig_with_transducer_img, rotated_with_transducer_img}, ...
        'GridSize', [1 nan], 'BackgroundColor', [0.5 0.5 0.5]);
end
