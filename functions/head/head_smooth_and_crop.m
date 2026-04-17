function [medium_masks, segmentation_crop, bone_crop, trans_pos_final, focus_pos_final, ...
    translation_matrix] = head_smooth_and_crop(parameters, segmentation, bone_img, ...
    trans_pos_grid, focus_pos_grid)
    arguments
        parameters struct
        segmentation (:,:,:) % SimNIBS tissue segmentation
        bone_img (:,:,:) % binary bone mask or pseudoCT
        trans_pos_grid
        focus_pos_grid
    end

    % This function turns the original `layered` segmentations into medium masks such
    % that the setup_medium.m function can fill in the tissue-dependent parameters.
    % Tissue masks will contain IDs according to the order of tissues in parameters.medium_properties.

    grid.resolution_mm = parameters.grid.resolution_mm;
   
    % Segmentations will be postprocessed. 
    % Incl. smoothing layer transitions & filling potential skull segmentation gaps. 

    % create "medium_masks" that contains indices according to the label order in parameters.layers
    % each mask will be smoothed in the process
    log_timer('start','preproc_medium_mask', parameters.io.output_dir);
    [medium_masks] = preproc_medium_mask(segmentation, parameters);
    log_timer('stop','preproc_medium_mask');

    % Fill gaps in skull mask
    requested_layers = fieldnames(parameters.layers);
    if any(contains(requested_layers, 'skull'))        
        [medium_masks, ~] = skull_fill_holes(parameters, ...
            medium_masks, focus_pos_grid, segmentation);
    end

    % [DEBUG] Plot segmentation and smoothed medium mask
    if parameters.simulation.debug == 1
        h = figure;
        imshowpair(label2rgb(squeeze(segmentation(:,trans_pos_grid(2),:))), ...
            label2rgb(squeeze(medium_masks(:,trans_pos_grid(2),:))), 'montage')
        title('Original segmentation (left) and smoothed medium mask (right)')
        output_plot_filename = fullfile(parameters.io.debug_dir_preproc, ...
            sprintf('sub-%03d_%s_segmented_img_smoothing_changes%s.png', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);
    end

    % Crop the simulation grid outside layered medium + transducer + PML for efficiency
    [medium_masks, segmentation_crop, bone_crop, parameters, trans_pos_final, focus_pos_final, translation_matrix] = ...
         preproc_crop_grid(parameters, medium_masks, segmentation, bone_img, trans_pos_grid, focus_pos_grid);

    % [DEBUG] plot the smoothed and unsmoothed skull segmentation with transducer and focus locations
    if parameters.simulation.debug == 1
        h = figure;
        imshowpair(plot_t1_with_transducer(segmentation, grid.resolution_mm, trans_pos_grid, focus_pos_grid, parameters), ...
            plot_t1_with_transducer(medium_masks, grid.resolution_mm, trans_pos_final, focus_pos_final, parameters),...
            'montage')
        title('Original segmentation (left) and Cropped, padded, & smoothed medium mask (right)')
        output_plot_filename = fullfile(parameters.io.debug_dir_preproc, ...
            sprintf('sub-%03d_%s_seg_smoothing_and_cropping_%s.png', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);
    end
    
end

