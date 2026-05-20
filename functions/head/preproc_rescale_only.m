function [rescaled_img, trans_pos_new, focus_pos_new, transformation_matrix, rotation_matrix, montage_img] = ...
    preproc_rescale_only(nii_image, nii_header, trans_pos_grid, focus_pos_grid, scale_factor, parameters)

% PREPROC_RESCALE_ONLY  Rescale a 3D image to grid resolution without rotation
%
% Isotropically rescales a 3D image to the target grid resolution while
% keeping the volume in its original (RAS+) orientation. Used as the
% grid.mode='ras_plus' alternative to PREPROC_ALIGN_TO_FOCAL_AXIS, which
% additionally rotates the volume to align the focal axis with z.
%
% The transformation matrix has the same structure as the one returned by
% PREPROC_ALIGN_TO_FOCAL_AXIS so that downstream code (preproc_head,
% final_transformation_matrix composition) needs no special-casing.
%
% Use as:
%   [rescaled_img, trans_pos_new, focus_pos_new, transformation_matrix, ...
%    rotation_matrix, montage_img] = ...
%       preproc_rescale_only(nii_image, nii_header, trans_pos_grid, ...
%           focus_pos_grid, scale_factor, parameters)
%
% Input:
%   nii_image      - [Nx x Ny x Nz] 3D image volume
%   nii_header     - NIfTI header struct (from niftiinfo)
%   trans_pos_grid - [1x3] transducer position in voxel coordinates
%   focus_pos_grid - [1x3] focus position in voxel coordinates
%   scale_factor   - image scaling factor (T1_res / grid_res; 1 = no scaling)
%   parameters     - (1,1) simulation parameters struct
%
% Output:
%   rescaled_img       - rescaled image volume (same orientation as input)
%   trans_pos_new      - [1x3] transducer position after rescaling
%   focus_pos_new      - [1x3] focus position after rescaling
%   transformation_matrix - [4x4] scale-and-recenter affine (no rotation)
%   rotation_matrix    - [4x4] identity (returned for API compatibility)
%   montage_img        - side-by-side comparison montage
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

    %% Build scale-only affine (no rotation)

    scale_matrix = makehgtform('scale', scale_factor);

    % Translate image centre to origin, scale, translate to new centre
    rotation_center = (size(nii_image) + 1) / 2;
    forward  = makehgtform('translate', -rotation_center);

    T_centered = forward' * scale_matrix';
    TF_centered = T_centered;

    img_bounds = affine_bounds(TF_centered, [0 0 0; size(nii_image)]);
    tp_trans   = affine_apply_pts(trans_pos_grid, TF_centered);
    fp_trans   = affine_apply_pts(focus_pos_grid, TF_centered);

    all_points = [img_bounds; tp_trans; fp_trans];
    max_abs    = max(abs(all_points), [], 1);
    newdims    = ceil(2 * max_abs) + 2;

    backward = makehgtform('translate', (newdims + 1) / 2);
    transformation_matrix = forward' * scale_matrix' * backward';

    % Identity rotation — no orientation change in ras_plus mode
    rotation_matrix = eye(4);

    %% Resample image

    if numel(unique(nii_image)) < 20
        parameters.interpolation = 'nearest';
    else
        parameters.interpolation = 'linear';
    end

    rescaled_img = affine_resample_3d(nii_image, transformation_matrix, newdims, parameters.interpolation, 0);

    %% Update transducer and focus positions

    out_mat = round(affine_apply_pts([trans_pos_grid; focus_pos_grid], transformation_matrix));
    trans_pos_new = out_mat(1,:);
    focus_pos_new = out_mat(2,:);

    %% Diagnostic montage

    orig_with_transducer_img = plot_t1_with_transducer(nii_image, nii_header.PixelDimensions(1), ...
        trans_pos_grid, focus_pos_grid, parameters);

    rescaled_with_transducer_img = plot_t1_with_transducer(rescaled_img, ...
        nii_header.PixelDimensions(1) / scale_factor, trans_pos_new, focus_pos_new, parameters);

    montage_img = imtile({orig_with_transducer_img, rescaled_with_transducer_img}, ...
        'GridSize', [1 nan], 'BackgroundColor', [0.5 0.5 0.5]);
end
