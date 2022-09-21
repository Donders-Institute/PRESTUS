function [skull_mask, skull_edge, segmented_image_cropped, trans_pos_final, focus_pos_final, min_dims, max_dims, new_grid_dims, translation_matrix] = smooth_and_crop_skull(segmented_img, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters)
    
    % unsmoothed skull & skin masks
    skull_mask_unsmoothed = segmented_img==4;
    % csf stands for cerebrospinal fluid, we use it to crop the image along
    % with the skull mask (because skull mask includes other bones)
    csf_mask = segmented_img==3;
    %skin_mask_unsmoothed = segmented_image==5;

    windowSize = 4;
    kernel = ones(windowSize,windowSize,windowSize)/ windowSize ^ 3;
    blurryImage = convn(double(skull_mask_unsmoothed), kernel, 'same');
    skull_mask_smoothed = blurryImage > parameters.skull_smooth_threshold; % Rethreshold
    
    imshowpair(mat2gray(squeeze(skull_mask_unsmoothed(:,trans_pos_upsampled_grid(2),:))),...
        mat2gray(squeeze(skull_mask_smoothed(:,trans_pos_upsampled_grid(2),:))), 'montage');
    
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_before_after_smoothing%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    title('Original (left) and smoothed (right) skull mask')
    export_fig(output_plot, '-native')

    transducer_bowl = transducer_setup(parameters.transducer, trans_pos_upsampled_grid, focus_pos_upsampled_grid, ...
                                                            size(segmented_img), voxel_size_mm);

    % restrict skull based on CSF 
    % expand CSF by a large amount
    
    SE = strel('cube', parameters.csf_mask_expansion_factor/voxel_size_mm);
    csf_mask_expanded = imdilate(csf_mask, SE);

%     imshow(plot_t1_with_transducer(segmented_image, segmented_hdr, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters))
    imshowpair(squeeze(skull_mask_smoothed(:,trans_pos_upsampled_grid(2),:)), squeeze(csf_mask_expanded(:,trans_pos_upsampled_grid(2),:)), 'falsecolor')
    title('CSF mask (pink) and the skull (white)')
    
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_csf_mask%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    export_fig(output_plot, '-native')

    crop_margin = parameters.pml_size+1; 
    
    % combine expanded CSF mask and skull mask, so that only the parts of
    % the skull close to the brain remain
    
    skull_mask_smoothed = skull_mask_smoothed & csf_mask_expanded;
    %imshow(squeeze(skull_mask_smoothed(trans_pos_upsampled_grid(1),:,:)))
    %montage({label2rgb(squeeze(segmented_img(trans_pos_upsampled_grid(1),:,:)),'jet','k'), label2rgb(squeeze(segmented_img(:,trans_pos_upsampled_grid(2),:)),'jet','k'),label2rgb(squeeze(segmented_img(:,:,3)),'jet','k')},'Size',[1 3])
    [min_dims, max_dims, new_grid_dims] = get_crop_dims(skull_mask_smoothed+transducer_bowl, crop_margin);
    
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
    
    
    if any( max_dims > size(skull_mask_smoothed) )
        pad_amount = max( max_dims, size(segmented_img))-size(segmented_img);
        segmented_img = padarray(segmented_img, pad_amount,0,'post');
        skull_mask_smoothed = padarray(skull_mask_smoothed, pad_amount,0,'post');        
    end

    skull_mask = skull_mask_smoothed(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3)); % Crop image
    skull_edge = edge3(skull_mask, 'approxcanny',0.1);

    segmented_image_cropped = segmented_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3)); % Crop image

    parameters.grid_dims = size(skull_mask);

    trans_pos_final = trans_pos_upsampled_grid-min_dims';
    focus_pos_final = focus_pos_upsampled_grid-min_dims';
    
    translation_matrix = makehgtform('translate', -min_dims);
    
    disp('New grid dimensions:')
    disp(new_grid_dims)
    disp('New transducer position:')
    disp(trans_pos_final)
    
    imshowpair(plot_t1_with_transducer(skull_mask, voxel_size_mm, trans_pos_final, focus_pos_final, parameters),...
        plot_t1_with_transducer(skull_mask_unsmoothed, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters),...
    'montage')
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_skull_final%s.png', parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    title('Cropped, padded, & smoothed (left) and original (right) skull mask')
    export_fig(output_plot, '-native')


end