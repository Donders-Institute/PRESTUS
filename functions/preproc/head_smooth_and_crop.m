function [medium_masks, segmented_image_cropped, trans_pos_final, focus_pos_final, ...
    translation_matrix] = head_smooth_and_crop(parameters, segmented_img, bone_img, ...
    trans_pos_grid, focus_pos_grid)
    arguments
        parameters struct
        segmented_img (:,:,:) % tissue segmentation from SimNibs [alt: pseudoCT tissue mask]
        bone_img (:,:,:) % infilled bone mask [alt: pseudoCT]
        trans_pos_grid
        focus_pos_grid
    end

    % This function turns the original `layered` segmentations into medium masks such
    % that the setup_medium.m function can fill in the tissue-dependent parameters.

    % Note that tissue masks will assume the labels specified below.

    % [pCT] If bone_img is a pseudoCT, segmented_image_cropped will contain the
    % continuous values that will be used by setup_medium.m.

    grid_step_mm = parameters.grid_step_mm;

    labels = fieldnames(parameters.layers);
   
    % If a standard layered solution is used, the segmented images have to
    % be postprocessed. This includes smoothing layer transitions and
    % filling in potential gaps in the skull segmentation. 
    
    % These operations are expected to already have occured during pseudoCT creation.
    % But: other tissues and mask to fill may still benefit from this...
    
    if parameters.usepseudoCT==0

        % create "medium_masks" that contains indices according to the label order in parameters.layers
        % each mask will be smoothed in the process
        log_timer('start','preproc_medium_mask', parameters.output_dir);
        [medium_masks] = preproc_medium_mask(segmented_img, parameters);
        log_timer('stop','preproc_medium_mask');

        % Fill gaps in skull (and skull-skin) mask
        if any(contains(labels, 'skull'))        
            [medium_masks, skull_i] = skull_fill_holes(parameters, ...
                medium_masks, labels, focus_pos_grid, segmented_img);
        end

        % [DEBUG] Plot segmentation and smoothed medium mask
        if parameters.debug == 1
            h = figure;
            imshowpair(label2rgb(squeeze(segmented_img(:,trans_pos_grid(2),:))), ...
                label2rgb(squeeze(medium_masks(:,trans_pos_grid(2),:))), 'montage')
            title('Original segmentation (left) and smoothed medium mask (right)')
            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_segmented_img_smoothing_changes%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
        end
        
    elseif parameters.usepseudoCT==1
        % pseudoCT already has smoothing etc. applied
        % we still need to map the segmentation labels to medium masks
        % Creates an empty grid the size of the segmented image
        medium_masks = zeros(size(segmented_img));
        % add a smoothing threshold to bone and other non-water tissue.
        for label_i = 1:length(labels)
            if strcmp(labels{label_i}, 'water')
               continue
            end
            sim_nibs_layers = parameters.layers.(labels{label_i});
            layer_mask = ismember(segmented_img, sim_nibs_layers);
            medium_masks(layer_mask) = label_i;
        end
        % overwrite the segmentation with pseudoCT values
        segmented_img = bone_img;
        % use cortical skull for edge definition
        skull_i = find(strcmp(labels, 'skull_cortical'));
    end

    % Crop the simulation grid outside layered medium + transducer + PML for efficiency
    [medium_masks, segmented_image_cropped, parameters, trans_pos_final, focus_pos_final, translation_matrix] = ...
         preproc_crop_grid(parameters, medium_masks, segmented_img, trans_pos_grid, focus_pos_grid);
    
    % [DEBUG] plot the smoothed and unsmoothed skull segmentation with transducer and focus locations
    if parameters.debug == 1
        h = figure;
        imshowpair(plot_t1_with_transducer(segmented_img, grid_step_mm, trans_pos_grid, focus_pos_grid, parameters), ...
            plot_t1_with_transducer(medium_masks, grid_step_mm, trans_pos_final, focus_pos_final, parameters),...
            'montage')
        title('Original segmentation (left) and Cropped, padded, & smoothed medium mask (right)')
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_seg_smoothing_and_cropping_%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);
    end
    
end

