function [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, trans_pos_final, focus_pos_final, ...
        final_transformation_matrix, crop_translation_matrix] = head_model_place(parameters, hm, trans_pos_t1, focus_pos_t1, preproc_fn, grid_mode)
% HEAD_MODEL_PLACE  Simulation grid of one transducer placement from the shared head model
%
% Rotates the pose-independent head model of HEAD_MODEL_LOAD_OR_BUILD to the
% beam (grid.mode 'focal_axis' / 'transducer_axis', with the same function as
% the unshared path, at scale 1: labels nearest-neighbour, pseudo-CT linear)
% or keeps it unrotated ('ras_plus': no resampling at all), and crops it
% around head and transducers (PREPROC_CROP_GRID). This replaces the rotation,
% smoothing and hole filling that PREPROC_HEAD otherwise repeats for every
% placement; in the rotated modes the smoothing then precedes the rotation
% instead of following it.
%
% Use as:
%   [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, trans_pos_final, focus_pos_final, ...
%    final_transformation_matrix, crop_translation_matrix] = ...
%       head_model_place(parameters, hm, trans_pos_t1, focus_pos_t1, preproc_fn, grid_mode)
%
% Input:
%   parameters   - PRESTUS config; transducer positions in the RAS+ T1 grid
%   hm           - head model (HEAD_MODEL_LOAD_OR_BUILD)
%   trans_pos_t1, focus_pos_t1 - [1x3] transducer(1) positions [RAS+ T1 voxels]
%   preproc_fn   - rotation function of grid_mode (PREPROC_HEAD)
%   grid_mode    - grid.mode
%
% Output: as PREPROC_HEAD / HEAD_SMOOTH_AND_CROP; final_transformation_matrix
%   maps RAS+ T1 voxels to the simulation grid (row-vector convention)
%
% See also: HEAD_MODEL_LOAD_OR_BUILD, PREPROC_HEAD, PREPROC_CROP_GRID

arguments
    parameters   (1,1) struct
    hm           (1,1) struct
    trans_pos_t1 (1,3) double
    focus_pos_t1 (1,3) double
    preproc_fn   (1,1) function_handle
    grid_mode    (1,:) char
end

    n_trans = numel(parameters.transducer);
    tp_h = affine_apply_pts(trans_pos_t1, hm.M);
    fp_h = affine_apply_pts(focus_pos_t1, hm.M);
    % Transducers 2..N in head-model voxels (the rotation functions read them)
    p_h = parameters;
    for ti = 2:n_trans
        p_h.transducer(ti).trans_pos = affine_apply_pts(parameters.transducer(ti).trans_pos(:)', hm.M);
        p_h.transducer(ti).focus_pos = affine_apply_pts(parameters.transducer(ti).focus_pos(:)', hm.M);
    end

    if strcmp(grid_mode, 'ras_plus')
        medium_masks = hm.medium_masks; segmentation = hm.segmentation;
        bone_mask = hm.bone_mask; pseudoCT = hm.pseudoCT;
        M_place = eye(4);
        trans_pos = round(tp_h); focus_pos = round(fp_h);
    else
        hdr = struct('PixelDimensions', repmat(hm.resolution_mm, 1, 3), 'ImageSize', size(hm.medium_masks));
        [medium_masks, trans_pos, focus_pos, M_place] = preproc_fn(hm.medium_masks, hdr, tp_h, fp_h, 1, p_h);
        segmentation = preproc_fn(hm.segmentation, hdr, tp_h, fp_h, 1, p_h);
        bone_mask    = preproc_fn(hm.bone_mask, hdr, tp_h, fp_h, 1, p_h);
        pseudoCT = [];
        if ~isempty(hm.pseudoCT)
            pseudoCT = preproc_fn(hm.pseudoCT, hdr, tp_h, fp_h, 1, p_h);
        end
    end
    scale_rotate_recenter_matrix = hm.M * M_place;      % RAS+ T1 voxel -> rotated grid

    all_trans = zeros(n_trans, 3); all_focus = zeros(n_trans, 3);
    all_trans(1, :) = trans_pos; all_focus(1, :) = focus_pos;
    for ti = 2:n_trans
        all_trans(ti, :) = round(affine_apply_pts(parameters.transducer(ti).trans_pos(:)', scale_rotate_recenter_matrix));
        all_focus(ti, :) = round(affine_apply_pts(parameters.transducer(ti).focus_pos(:)', scale_rotate_recenter_matrix));
    end

    [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, ~, trans_pos_final, focus_pos_final, ...
        crop_translation_matrix] = preproc_crop_grid(parameters, medium_masks, segmentation, bone_mask, ...
        trans_pos, focus_pos, pseudoCT, all_trans, all_focus);
    final_transformation_matrix = scale_rotate_recenter_matrix * crop_translation_matrix';
end
