function [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, trans_pos_final, ...
    focus_pos_final, translation_matrix] = head_smooth_and_crop(parameters, segmentation, ...
    bone_mask_img, trans_pos_grid, focus_pos_grid, pseudoCT_img, all_trans_pos_grid, all_focus_pos_grid)
% HEAD_SMOOTH_AND_CROP  Convert segmentation to medium masks and crop to simulation grid
%
% Converts layered tissue segmentations into medium masks indexed by
% parameters.medium_properties, smooths layer transitions, fills skull
% gaps, and crops the volume to the simulation region (PREPROC_CROP_GRID).
%
% Use as:
%   [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, trans_pos_final, ...
%    focus_pos_final, translation_matrix] = head_smooth_and_crop( ...
%       parameters, segmentation, bone_mask_img, trans_pos_grid, focus_pos_grid)
%   [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, trans_pos_final, ...
%    focus_pos_final, translation_matrix] = head_smooth_and_crop( ...
%       parameters, segmentation, bone_mask_img, trans_pos_grid, focus_pos_grid, pseudoCT_img)
%   [...] = head_smooth_and_crop(..., pseudoCT_img, all_trans_pos_grid, all_focus_pos_grid)
%
% Input:
%   parameters     - (1,1) simulation configuration struct
%   segmentation   - [Nx x Ny x Nz] SimNIBS tissue label volume
%   bone_mask_img  - [Nx x Ny x Nz] binary skull mask
%   trans_pos_grid - [1x3] transducer position in voxel coordinates
%   focus_pos_grid - [1x3] focus position in voxel coordinates
%   pseudoCT_img   - [Nx x Ny x Nz] Hounsfield-unit skull image (optional, [] if unused)
%   all_trans_pos_grid - [N x 3] positions of all transducers, kept inside the
%                        crop (optional, default: zeros(0,3))
%   all_focus_pos_grid - [N x 3] focus positions of all transducers (optional)
%
% Output:
%   medium_masks      - [Nx x Ny x Nz] medium label array (cropped)
%   segmentation_crop - [Nx x Ny x Nz] tissue label array (cropped)
%   bone_mask_crop    - [Nx x Ny x Nz] binary skull mask (cropped)
%   pseudoCT_crop     - [Nx x Ny x Nz] pseudoCT image (cropped; [] when not provided)
%   trans_pos_final   - [1x3] transducer position in cropped grid
%   focus_pos_final   - [1x3] focus position in cropped grid
%   translation_matrix - [4x4] homogeneous translation matrix
%
% See also: PREPROC_CROP_GRID, PREPROC_MEDIUM_MASK

    arguments
        parameters          (1,1) struct
        segmentation        (:,:,:) {mustBeNumericOrLogical}
        bone_mask_img       (:,:,:) {mustBeNumericOrLogical}
        trans_pos_grid      (1,:)   double
        focus_pos_grid      (1,:)   double
        pseudoCT_img                {mustBeNumericOrLogical} = []
        all_trans_pos_grid  (:,3)   double = zeros(0,3)
        all_focus_pos_grid  (:,3)   double = zeros(0,3)
    end

    % This function turns the original `layered` segmentations into medium masks such
    % that MEDIUM_SETUP can fill in the tissue-dependent parameters.
    % Tissue masks will contain IDs according to the order of tissues in parameters.medium_properties.

    grid.resolution_mm = parameters.grid.resolution_mm;
   
    % Segmentations will be postprocessed. 
    % Incl. smoothing layer transitions & filling potential skull segmentation gaps. 

    % create "medium_masks" that contains indices according to the label order in parameters.medium_properties
    % each mask will be smoothed in the process
    log_timer('start','preproc_medium_mask', parameters.io.dir_output);
    [medium_masks] = preproc_medium_mask(segmentation, parameters);
    log_timer('stop','preproc_medium_mask');
    seg_qc('add_mapping', parameters, segmentation, medium_masks, focus_pos_grid, ...
        seg_qc('grid_note', parameters, 'placement'));

    % Fill gaps in skull mask
    requested_layers = fieldnames(parameters.layers);
    if any(contains(requested_layers, 'skull'))        
        [medium_masks, ~] = skull_fill_holes(parameters, ...
            medium_masks, focus_pos_grid, segmentation);
    end

    % [DEBUG] Plot segmentation and smoothed medium mask
    if debug_on(parameters, 'preproc', 'figures')
        h = figure;
        imshowpair(label2rgb(squeeze(segmentation(:,trans_pos_grid(2),:))), ...
            label2rgb(squeeze(medium_masks(:,trans_pos_grid(2),:))), 'montage')
        title('Original segmentation (left) and smoothed medium mask (right)')
        debug_save_figure(parameters, 'preproc', h, ...
            sprintf('sub-%03d_%s_segmented_img_smoothing_changes%s', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix), ...
            'Segmentation (left) and smoothed, hole-filled medium mask (right), transducer y-slice');
    end

    % Crop the simulation grid outside layered medium + transducer + PML for efficiency
    [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, parameters, trans_pos_final, focus_pos_final, translation_matrix] = ...
         preproc_crop_grid(parameters, medium_masks, segmentation, bone_mask_img, trans_pos_grid, focus_pos_grid, pseudoCT_img, all_trans_pos_grid, all_focus_pos_grid);

    % [DEBUG] plot the smoothed and unsmoothed skull segmentation with transducer and focus locations
    if debug_on(parameters, 'preproc', 'figures')
        h = figure;
        imshowpair(plot_t1_with_transducer(segmentation, grid.resolution_mm, trans_pos_grid, focus_pos_grid, parameters), ...
            plot_t1_with_transducer(medium_masks, grid.resolution_mm, trans_pos_final, focus_pos_final, parameters),...
            'montage')
        title('Original segmentation (left) and Cropped, padded, & smoothed medium mask (right)')
        debug_save_figure(parameters, 'preproc', h, ...
            sprintf('sub-%03d_%s_seg_smoothing_and_cropping%s', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix), ...
            'Segmentation (left) and cropped, padded, smoothed medium mask (right) with transducer and focus');
    end
    
end

