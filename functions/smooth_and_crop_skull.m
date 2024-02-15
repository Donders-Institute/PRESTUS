function [skull_mask, skull_edge, segmented_image_cropped, trans_pos_final, focus_pos_final, min_dims, max_dims, new_grid_dims, translation_matrix] = smooth_and_crop_skull(segmented_img, bone_img, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters)
    
    %% Smoothes the skull
    % Renames images from segmented_img to make script more readible
    % Unsmoothed skull & skin masks

    
    skull_mask_unsmoothed = segmented_img==4;
    csf_mask = segmented_img==3;
    % csf stands for cerebrospinal fluid, we use it to crop the image along
    % with the skull mask (because skull mask includes other bones)
    
    %skin_mask_unsmoothed = segmented_image==5;

    % Smoothes the skull, meaning that it removes the differentiation
    % between gray and white matter for the simulations.
    windowSize = 4;
    skull_mask_smoothed = smooth_img(skull_mask_unsmoothed, windowSize, parameters.skull_smooth_threshold);
    
    % add a boundary of the bone mask to fill in potential gaps
    bone_img_smoothed = smooth_img(bone_img, windowSize, parameters.skull_smooth_threshold);
    
    bone_perimeter = bone_img_smoothed - imerode(bone_img_smoothed, strel('cube',3));
    skull_mask_smoothed = skull_mask_smoothed | bone_perimeter;


    % Shows the unsmoothed and smoothed image side-by-side
    imshowpair(mat2gray(squeeze(skull_mask_unsmoothed(:,trans_pos_upsampled_grid(2),:))),...
        mat2gray(squeeze(skull_mask_smoothed(:,trans_pos_upsampled_grid(2),:))), 'montage');
    
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_before_after_smoothing%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    title('Original (left) and smoothed (right) skull mask')
    export_fig(output_plot, '-native')

    % This creates a transducer overlay for use in a figure
    transducer_bowl = transducer_setup(parameters.transducer, trans_pos_upsampled_grid, focus_pos_upsampled_grid, ...
                                                            size(segmented_img), voxel_size_mm);

    %% Cropping the bone tissue
    % Some bone tissue is not close to the neural tissue, and thus not essential for simulations.
    % The csf will be expanded to be used as a guide for what bone tissue should be left in.
    
    % To do this, the CSF volume will first be expanded by a lot.
    SE = strel('cube', parameters.csf_mask_expansion_factor/voxel_size_mm);
    csf_mask_expanded = imdilate(csf_mask, SE);

%    imshow(plot_t1_with_transducer(segmented_image, segmented_hdr, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters))
    
    % Plots the result of the CSF expansion
    imshowpair(squeeze(skull_mask_smoothed(:,trans_pos_upsampled_grid(2),:)), squeeze(csf_mask_expanded(:,trans_pos_upsampled_grid(2),:)), 'falsecolor')
    title('CSF mask (pink) and the skull (white)')
    
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_csf_mask%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    export_fig(output_plot, '-native')

    % The margin around the edge of the figure where the CSF will be cropped
    crop_margin = parameters.pml_size+1;
    
    % get_crop_dims ensures that the CSF mask is not expanded outside the
    % bounds of the image dimensions
    skull_mask_smoothed = skull_mask_smoothed & csf_mask_expanded;
    %imshow(squeeze(skull_mask_smoothed(trans_pos_upsampled_grid(1),:,:)))
    %montage({label2rgb(squeeze(segmented_img(trans_pos_upsampled_grid(1),:,:)),'jet','k'), label2rgb(squeeze(segmented_img(:,trans_pos_upsampled_grid(2),:)),'jet','k'),label2rgb(squeeze(segmented_img(:,:,3)),'jet','k')},'Size',[1 3])
    [min_dims, max_dims, new_grid_dims] = get_crop_dims(skull_mask_smoothed+transducer_bowl, crop_margin);
    
    % In case the minimum dimensions are too small, 
    if any(min_dims < 1)
        pad_amount = abs(min(min_dims, [1 1 1]));
        
        segmented_img = padarray(segmented_img, pad_amount,0,'pre');
        skull_mask_smoothed = padarray(skull_mask_smoothed, pad_amount,0,'pre');
        min_dims = max(min_dims, [1 1 1]);
        max_dims = max_dims+pad_amount;
        new_grid_dims = max_dims - min_dims + 1;
        
    end

    new_grid_dims(1) = find_min_factor(new_grid_dims(1),new_grid_dims(1)+parameters.prime_factor_max_grid_expansion);
    new_grid_dims(2) = find_min_factor(new_grid_dims(2),new_grid_dims(2)+parameters.prime_factor_max_grid_expansion);
    new_grid_dims(3) = find_min_factor(new_grid_dims(3),new_grid_dims(3)+parameters.prime_factor_max_grid_expansion);
    max_dims = min_dims+new_grid_dims - 1;
    
    % Ensures that the maximum size of the new grid matches the skull mask
    if any( max_dims > size(skull_mask_smoothed) )
        pad_amount = max( max_dims, size(segmented_img))-size(segmented_img);
        segmented_img = padarray(segmented_img, pad_amount,0,'post');
        skull_mask_smoothed = padarray(skull_mask_smoothed, pad_amount,0,'post');        
    end

    % Crops the figure
    skull_mask = skull_mask_smoothed(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
    % Creates a mask around the edge of the figure to define the edge of the skull
    skull_edge = edge3(skull_mask, 'approxcanny',0.1);

    % Crops the segmented figure
    segmented_image_cropped = segmented_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));

    % Saves the skull mask in the parameters structure as the new grid dimensions
    parameters.grid_dims = size(skull_mask);

    % Moves the transducer and focus positions to the correct location in
    % the cropped grid
    trans_pos_final = trans_pos_upsampled_grid-min_dims';
    focus_pos_final = focus_pos_upsampled_grid-min_dims';
    
    % Creates a translation matrix for later reference
    translation_matrix = makehgtform('translate', -min_dims);
    
    disp('New grid dimensions:')
    disp(new_grid_dims)
    disp('New transducer position:')
    disp(trans_pos_final)
    
    % Shows and saves figures of the smoothed and unsmoothed skull mask with the transducer and focus locations
    imshowpair(plot_t1_with_transducer(skull_mask, voxel_size_mm, trans_pos_final, focus_pos_final, parameters),...
        plot_t1_with_transducer(skull_mask_unsmoothed, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters),...
    'montage')
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_final%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    title('Cropped, padded, & smoothed (left) and original (right) skull mask')
    export_fig(output_plot, '-native')

end