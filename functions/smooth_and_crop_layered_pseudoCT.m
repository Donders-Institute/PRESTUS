function [tissues_mask_img_cropped, skull_edge, segmented_image_cropped, trans_pos_final, focus_pos_final, min_dims, max_dims, ...
    new_grid_dims, translation_matrix] = smooth_and_crop_layered_pseudoCT(segmented_img,  pseudoCT_img, voxel_size_mm, ...
    trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters)
    arguments
        segmented_img (:,:,:)
        pseudoCT_img (:,:,:)
        voxel_size_mm double
        trans_pos_upsampled_grid 
        focus_pos_upsampled_grid
        parameters struct
    end
    %segmented_img is segmented_img_rr in the preprocess_brain function,
    %bone_img is bone_img_rr, pseudoCT_img is pCT_img_rr, segmented_img is mask_img_rr and voxel_size_mm is parameters.grid_step_mm
    %% split into layers of interest
   
    labels = fieldnames(parameters.layer_labels);
    
    if any(strcmp(labels, 'csf'))
        csf_mask = segmented_img==find(strcmp(labels, 'csf'));
    else
        csf_mask = segmented_img==6;
    end
        
    % Creates a transducer for overlay in a figure
    transducer_bowl = transducer_setup(parameters.transducer, trans_pos_upsampled_grid, focus_pos_upsampled_grid, ...
                                                            size(segmented_img), voxel_size_mm);

    %% Cropping the bone tissue
    % Some bone tissue is not close to the neural tissue, and thus not essential for simulations.
    % The csf will be expanded to be used as a guide for what bone tissue should be left in.
    
    % To do this, the csf volume will first be expanded by a lot.
    SE = strel('cube', parameters.csf_mask_expansion_factor/voxel_size_mm);
    csf_mask_expanded = imdilate(csf_mask, SE);
    
    % Plots the result of the CSF expansion
    imshowpair(squeeze(segmented_img(:,trans_pos_upsampled_grid(2),:)), squeeze(csf_mask_expanded(:,trans_pos_upsampled_grid(2),:)), 'falsecolor')
    title('CSF mask (pink) and the segmented image (white)')
    
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_csf_mask_pCT%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    export_fig(output_plot, '-native')

    % The margin around the edge of the figure where the CSF will be cropped
    crop_margin = parameters.pml_size+1;
    
    % get_crop_dims ensures that the CSF mask is not expanded outside the
    % bounds of the image dimensions
    segmented_img(~csf_mask_expanded) = 0;
    %imshow(squeeze(skull_mask_smoothed(trans_pos_upsampled_grid(1),:,:)))
    %montage({label2rgb(squeeze(segmented_img(trans_pos_upsampled_grid(1),:,:)),'jet','k'), label2rgb(squeeze(segmented_img(:,trans_pos_upsampled_grid(2),:)),'jet','k'),label2rgb(squeeze(segmented_img(:,:,3)),'jet','k')},'Size',[1 3])
    [min_dims, max_dims, new_grid_dims] = get_crop_dims(segmented_img+transducer_bowl, crop_margin);
    
    % In case the minimum dimensions are too small, 
    if any(min_dims < 1)
        pad_amount = abs(min(min_dims, [1 1 1]));
        segmented_img = padarray(segmented_img, pad_amount,0,'pre');
        pseudoCT_img = padarray(pseudoCT_img, pad_amount,0,'pre');
        min_dims = max(min_dims, [1 1 1]);
        max_dims = max_dims+pad_amount;
        new_grid_dims = max_dims - min_dims + 1;
    end

    new_grid_dims(1) = find_min_factor(new_grid_dims(1),new_grid_dims(1)+parameters.prime_factor_max_grid_expansion);
    new_grid_dims(2) = find_min_factor(new_grid_dims(2),new_grid_dims(2)+parameters.prime_factor_max_grid_expansion);
    new_grid_dims(3) = find_min_factor(new_grid_dims(3),new_grid_dims(3)+parameters.prime_factor_max_grid_expansion);
    max_dims = min_dims+new_grid_dims - 1;
    
    % Ensures that the maximum size of the new grid matches the skull mask
    if any( max_dims > size(segmented_img) )
        pad_amount = max( max_dims, size(segmented_img))-size(segmented_img);
        segmented_img = padarray(segmented_img, pad_amount,0,'post');
        pseudoCT_img = padarray(pseudoCT_img, pad_amount,0,'post');
    end

    % Crops the processed segmented image
    tissues_mask_img_cropped = segmented_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3)); % Crop image
    
    % Creates a mask around the edge of the figure to define the edge of the skull
    skull_i = ismember(tissues_mask_img_cropped, find(strcmp(labels, 'skull_cortical') | strcmp(labels, 'skull_trabecular')));
    skull_edge = edge3(tissues_mask_img_cropped==skull_i, 'approxcanny',0.1);
    
    % Crops the original pseudoCT image in the same wave as the processed tissues mask one
    segmented_image_cropped = pseudoCT_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));

    % Saves the skull mask in the parameters structure as the new grid dimensions
    parameters.grid_dims = size(tissues_mask_img_cropped);

    % Moves the transducer and focus positions to the correct location in
    % the cropped grid
    trans_pos_final = trans_pos_upsampled_grid-min_dims;
    focus_pos_final = focus_pos_upsampled_grid-min_dims;
    
    % Creates a translation matrix for later reference
    translation_matrix = makehgtform('translate', -min_dims);
    
    disp('New grid dimensions:')
    disp(new_grid_dims)
    disp('New transducer position:')
    disp(trans_pos_final)
    
     % Shows and saves figures of the smoothed and unsmoothed skull mask with the transducer and focus locations
    imshowpair(plot_t1_with_transducer(tissues_mask_img_cropped, voxel_size_mm, trans_pos_final, focus_pos_final, parameters),...
        plot_t1_with_transducer(segmented_img, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters),...
    'montage')
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_final_pCT%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    title('Cropped (left) and original (right) tissue mask')
    export_fig(output_plot, '-native')

    skull_mask_file = fullfile(parameters.output_dir, sprintf('sub-%03d_skull_mask_final_pCT', parameters.subject_id));
    skull_mask = zeros(size(segmented_image_cropped));
    for layer_i = find(contains(labels,  'skull'))'
        i = segmented_image_cropped==layer_i; % which voxels contain skull
        skull_mask(i) = segmented_image_cropped(i);
    end
    niftiwrite(uint8(skull_mask), skull_mask_file ,'Compressed', 1);
   
    segmented_file = fullfile(parameters.output_dir, sprintf('sub-%03d_medium_masks_final_pCT', parameters.subject_id));
    niftiwrite(uint8(tissues_mask_img_cropped), segmented_file  ,'Compressed', 1);

    % niftiwrite creates a file with the name segment_file (so
    % sub-009_medium_masks_final in the specified directory) using the
    % volume V uint8(smoothed_segmented_img)
    % uint8( ) transforms the values of smoothed_segmented_img in type
    % uint8 (8 bit = 1 byte unsigned integers)

    imshowpair(plot_t1_with_transducer(segmented_image_cropped, voxel_size_mm, trans_pos_final, focus_pos_final, parameters),...
        plot_t1_with_transducer(pseudoCT_img, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters),...
    'montage')
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_cropped_pCT%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    title('Cropped (left) and original (right) pseudoCT')
    export_fig(output_plot, '-native')
    
end

