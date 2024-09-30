function [medium_masks, skull_edge, segmented_image_cropped, trans_pos_final, focus_pos_final, min_dims, max_dims, ...
    new_grid_dims, translation_matrix] = smooth_and_crop(segmented_img, bone_img,  voxel_size_mm, ...
    trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters)
    arguments
        segmented_img (:,:,:) % tissue segmentation from SimNibs [alt: pseudoCT tissue mask]
        bone_img (:,:,:) % infilled bone mask [alt: pseudoCT]
        voxel_size_mm double % parameters.grid_step_mm
        trans_pos_upsampled_grid 
        focus_pos_upsampled_grid
        parameters struct
    end
    
    % This function turns the original segmentations into medium masks such
    % that the setup_medium.m function can fill in the tissue-dependent parameters.
    % Note that tissue masks will assume the labels specified below.
    % If bone_img is a pseudoCT, segmented_image_cropped will contain the
    % continuous values that will be used by setup_medium.m.

    labels = fieldnames(parameters.layer_labels);

    % create CSF mask for later
    if isfield(parameters, 'seg_labels') && any(strcmp(fieldnames(parameters.seg_labels), 'csf'))
        csf_mask = segmented_img==getidx(parameters.seg_labels,'csf');
    else
        csf_mask = segmented_img==3;
        warning("Choosing default csf segmentation label = 3");
    end

   %% skull postprocessing
   
    % If a standard layered solution is used, the segmented images have to
    % be postprocessed. This includes smoothing layer transitions and
    % filling in potential gaps in the skull segmentation. These are redundant
    % when using pseudoCTs.
    
    if strcmp(parameters.simulation_medium, 'layered') && parameters.usepseudoCT==0
        % Note: here, we create a new image of "medium_masks" that is
        % structured according to the indices in parameters.layer_labels.

        % smoothing window size
        windowSize = 4;
        % Creates an empty grid the size of the segmented image
        medium_masks = zeros(size(segmented_img));
        % add a smoothing threshold to bone and other non-water tissue.
        for label_i = 1:length(labels)
            if strcmp(labels{label_i}, 'water')
               continue
            end
            sim_nibs_layers = parameters.layer_labels.(labels{label_i});
            layer_mask = ismember(segmented_img, sim_nibs_layers);
            if contains(labels{label_i}, 'skull')
                smooth_threshold = parameters.skull_smooth_threshold;    
                if any(contains(labels, 'skull_cortical')) % two bone types are smoothed together later
                    continue
                end
            else
                smooth_threshold = parameters.other_smooth_threshold;
            end
            layer_mask_smoothed = smooth_img(layer_mask, windowSize, smooth_threshold);
            medium_masks(layer_mask_smoothed~=0) = label_i;
        end
        % fill gaps in the skull by using the boundary of a bone image
        if any(contains(labels, 'skull_cortical'))  
            skull_i = find(strcmp(labels, 'skull_cortical')); % gives to skull_i the index of skull_cortical in labels array
            % combine all skull masks for smoothing
            layer_mask = ismember(segmented_img, getidx(parameters.layer_labels, 'skull'));
            smooth_threshold = parameters.skull_smooth_threshold;
            layer_mask_smoothed = smooth_img(layer_mask, windowSize, smooth_threshold);
            medium_masks(layer_mask_smoothed~=0) = skull_i; % add it to image as cortical bone
            trabecular_i = find(strcmp(labels,  'skull_trabecular'));
            trabecular_mask = ismember(segmented_img, getidx(parameters.layer_labels, 'skull_trabecular'));
            trabecular_mask_smoothed = smooth_img(trabecular_mask, windowSize, smooth_threshold);
            medium_masks(trabecular_mask_smoothed~=0) = trabecular_i;
        else
            skull_i = find(strcmp(labels,  'skull')); % gives to skull_i the index of skull in labels array, because find finds the index of each non zero element --> where labels == 'skull'
        end

        smoothed_segmented_img_with_gaps = medium_masks;

        h = figure;
        imshowpair(label2rgb(squeeze(segmented_img(:,trans_pos_upsampled_grid(2),:))), label2rgb(squeeze(medium_masks(:,trans_pos_upsampled_grid(2),:))), 'montage')
        title('Original (left) and smoothed (right) segmented images')
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_segmented_img_smoothing_changes%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);

        % increase skull perimeter
        if any(contains(labels, 'skull'))        
            skull = medium_masks==skull_i;
            % skull is 1 when smoothed segm img == skull_i (=4 by default in headreco) and it is 0 otherwise
            smoothed_bone_img = smooth_img(bone_img, windowSize, parameters.skull_smooth_threshold);
            bone_perimeter = smoothed_bone_img - imerode(smoothed_bone_img, strel('cube',3));
            new_skull = skull | bone_perimeter;
            medium_masks(new_skull) = skull_i;

            % plot skull expansion
            h = figure;
            montage({squeeze(skull(:,focus_pos_upsampled_grid(2),:))*255, ...
                squeeze(bone_perimeter(:,focus_pos_upsampled_grid(2),:))*255, ...
                squeeze(new_skull(:,focus_pos_upsampled_grid(2),:))*255},gray, 'Size',[1 3])
            title('Skull (left), bone perimeter (center), and expanded skull (right)')
            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_skull_expansion%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
        end
        % remove gaps between skull & skin
        if any(contains(labels, 'skull')) && any(strcmp(labels, 'skin'))
            skin_i = find(strcmp(labels, 'skin'));
            skin = medium_masks==skin_i;
            skull = medium_masks==skull_i;
            skin_skull = skin+skull;
            skin_skull_filled = imfill(skin_skull);

            [labeledImage, ~] = bwlabeln(~skin_skull);
            blobMeasurements = regionprops(labeledImage, 'area');
            allAreas = [blobMeasurements.Area];

            [~, sortIndexes] = sort(allAreas, 'descend');
            biggestBlob = ismember(labeledImage, sortIndexes(2));

            medium_masks((skin_skull_filled-skin_skull - biggestBlob)>0) = skull_i;
            if any(contains(labels, 'skull_cortical'))  
                medium_masks(trabecular_mask_smoothed~=0) = trabecular_i; % makes sure that the trabecular mask is not affected
            end
            % make sure that eyes are not labeled as bone
            if isfield(parameters.seg_labels, 'eye')
                eye_i = parameters.seg_labels.eye;
                eye = segmented_img==eye_i;
                medium_masks(eye~=0) = 0; % default to water
            end
        end

        % Shows the segmented figures before and after closing the gaps side-by-side
        h = figure;
        imshowpair(label2rgb(squeeze(smoothed_segmented_img_with_gaps(:,trans_pos_upsampled_grid(2),:))), label2rgb(squeeze(medium_masks(:,trans_pos_upsampled_grid(2),:))), 'montage')
        title('Smoothed (left) and closed off (right) segmented images')
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_segmented_img_closing_gaps_changes%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);
        
    elseif contains(parameters.simulation_medium, 'skull') && parameters.usepseudoCT==0
        
        % Smoothes the skull, meaning that it removes the differentiation
        % between gray and white matter for the simulations.
        windowSize = 4;
        skull_mask_unsmoothed = ismember(segmented_img, getidx(parameters.layer_labels, 'skull'));
        skull_mask_smoothed = smooth_img(skull_mask_unsmoothed, windowSize, parameters.skull_smooth_threshold);
        % add a boundary of the bone mask to fill in potential gaps
        bone_img_smoothed = smooth_img(bone_img, windowSize, parameters.skull_smooth_threshold);
        bone_perimeter = bone_img_smoothed - imerode(bone_img_smoothed, strel('cube',3));
        % create binary medium musk (skull only)
        medium_masks = skull_mask_smoothed | bone_perimeter;
        % find specified skull index & adjuts medium mask accordingly
        skull_i = find(contains(labels,  'skull')); skull_i = skull_i(1); 
        medium_masks = medium_masks.*skull_i;
        % add brain segmentation if requested
        % if contains(parameters.simulation_medium, 'brain') &&...
        %         contains(parameters.simulation_medium, 'skull')
        %     brain_mask_unsmoothed = ismember(segmented_img, getidx(parameters.layer_labels, 'brain'));
        %     brain_mask_smoothed = smooth_img(brain_mask_unsmoothed, windowSize, parameters.skull_smooth_threshold);
        %     brain_i = find(contains(labels,  'brain'));
        %     % add to medium mask
        %     medium_masks = medium_masks + brain_i.*brain_mask_smoothed;
        %     % for slight overlap, assign to skull
        %     medium_masks(~ismember(medium_masks,[0,skull_i,brain_i])) = skull_i;
        % end

        % Shows the unsmoothed and smoothed image side-by-side
        h = figure;
        imshowpair(mat2gray(squeeze(skull_mask_unsmoothed(:,trans_pos_upsampled_grid(2),:))),...
            mat2gray(squeeze(medium_masks(:,trans_pos_upsampled_grid(2),:))), 'montage');
        title('Original (left) and smoothed (right) skull masks')
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_skull_smoothing_changes%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);
        
    elseif strcmp(parameters.simulation_medium, 'layered') && parameters.usepseudoCT==1
        
        % pseudoCT already has smoothing etc. applied
        % we only need to apply reassign labels
        
        % Creates an empty grid the size of the segmented image
        medium_masks = zeros(size(segmented_img));
        % add a smoothing threshold to bone and other non-water tissue.
        for label_i = 1:length(labels)
            if strcmp(labels{label_i}, 'water')
               continue
            end
            sim_nibs_layers = parameters.layer_labels.(labels{label_i});
            layer_mask = ismember(segmented_img, sim_nibs_layers);
            medium_masks(layer_mask) = label_i;
        end
        
        % overwrite the segmentation with pseudoCT values
        segmented_img = bone_img;

        % use cortical skull for edge definition
        skull_i = find(strcmp(labels, 'skull_cortical'));
    end

    %% Cropping the bone tissue
    % Some bone tissue is not close to the neural tissue, and thus not essential for simulations.
    % The csf will be expanded to be used as a guide for what bone tissue should be left in.
    
    % To do this, the csf volume will first be expanded by a lot.
    SE = strel('cube', parameters.csf_mask_expansion_factor/voxel_size_mm);
    csf_mask_expanded = imdilate(csf_mask, SE);
    
    % Plots the result of the CSF expansion
    h = figure;
    imshowpair(squeeze(medium_masks(:,trans_pos_upsampled_grid(2),:)), squeeze(csf_mask_expanded(:,trans_pos_upsampled_grid(2),:)), 'falsecolor')
    title('CSF mask (pink) and the segmented image (white)')
    output_plot_filename = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_%s_skull_csf_mask%s.png', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    saveas(h, output_plot_filename, 'png')
    close(h);

    % The margin around the edge of the figure where the CSF will be cropped
    crop_margin = parameters.pml_size+1;
    
    % get_crop_dims ensures that the CSF mask is not expanded outside the
    % bounds of the image dimensions
    medium_masks(~csf_mask_expanded) = 0;
    % calculate size if the transducer bowl
    transducer_bowl = transducer_setup(parameters.transducer, trans_pos_upsampled_grid, focus_pos_upsampled_grid, size(segmented_img), voxel_size_mm);
    [min_dims, max_dims, new_grid_dims] = get_crop_dims(medium_masks+transducer_bowl, crop_margin);
    
    % In case the minimum dimensions are too small, 
    if any(min_dims < 1)
        pad_amount = abs(min(min_dims, [1 1 1]));
        segmented_img = padarray(segmented_img, pad_amount,0,'pre');
        medium_masks = padarray(medium_masks, pad_amount,0,'pre');
        min_dims = max(min_dims, [1 1 1]);
        max_dims = max_dims+pad_amount;
        new_grid_dims = max_dims - min_dims + 1;
    end

    new_grid_dims(1) = find_min_factor(new_grid_dims(1),new_grid_dims(1)+parameters.prime_factor_max_grid_expansion);
    new_grid_dims(2) = find_min_factor(new_grid_dims(2),new_grid_dims(2)+parameters.prime_factor_max_grid_expansion);
    new_grid_dims(3) = find_min_factor(new_grid_dims(3),new_grid_dims(3)+parameters.prime_factor_max_grid_expansion);
    max_dims = min_dims+new_grid_dims - 1;
    
    % Ensures that the maximum size of the new grid matches the skull mask
    if any( max_dims > size(medium_masks) )
        pad_amount = max(max_dims, size(segmented_img))-size(segmented_img);
        segmented_img = padarray(segmented_img, pad_amount,0,'post');
        medium_masks = padarray(medium_masks, pad_amount,0,'post');
    end

    % crop the processed medium mask
    medium_masks = medium_masks(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
    % crop the original segmentation/pseudoCT
    segmented_image_cropped = segmented_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
    % define the edge of the skull
    skull_edge = edge3(medium_masks==skull_i, 'approxcanny',0.1);

    % save the skull mask in the parameters structure as the new grid dimensions
    parameters.grid_dims = size(medium_masks);

    % moves the transducer and focus positions to the correct location in the cropped grid
    trans_pos_final = trans_pos_upsampled_grid-min_dims;
    focus_pos_final = focus_pos_upsampled_grid-min_dims;
    
    % create a translation matrix for later reference
    translation_matrix = makehgtform('translate', -min_dims);
    
    disp('New grid dimensions:')
    disp(new_grid_dims)
    disp('New transducer position:')
    disp(trans_pos_final)
    
    % plot the smoothed and unsmoothed skull segmentation with transducer and focus locations
    h = figure;
    imshowpair(plot_t1_with_transducer(medium_masks, voxel_size_mm, trans_pos_final, focus_pos_final, parameters),...
        plot_t1_with_transducer(segmented_img, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters),...
        'montage')
    title('Cropped, padded, & smoothed (left) and original (right) segmented mask')
    output_plot_filename = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_%s_seg_smoothing_and_cropping_changes%s.png', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    saveas(h, output_plot_filename, 'png')
    close(h);

    % save medium mask
    segmented_file = fullfile(parameters.debug_dir, sprintf('sub-%03d_medium_masks_final', parameters.subject_id));
    niftiwrite(uint8(medium_masks), segmented_file  ,'Compressed', 1);
    
    % plot the smoothed and unsmoothed skull segmentation with transducer and focus locations
    h = figure;
    imshowpair(plot_t1_with_transducer(segmented_image_cropped, voxel_size_mm, trans_pos_final, focus_pos_final, parameters),...
        plot_t1_with_transducer(bone_img, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters),...
        'montage')
    title('Cropped, padded, & smoothed (left) and original (right) bone image')
    output_plot_filename = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_%s_skull_smoothing_and_cropping_changes%s.png', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    saveas(h, output_plot_filename, 'png')
    close(h);
    
    % save skull mask/pseudoCT
    skull_mask_file = fullfile(parameters.debug_dir, sprintf('sub-%03d_skull_final', parameters.subject_id));
    niftiwrite(uint8(segmented_image_cropped), skull_mask_file ,'Compressed', 1);
end

