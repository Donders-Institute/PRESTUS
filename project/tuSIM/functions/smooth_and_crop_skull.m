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
    skull_mask_smoothed = blurryImage > 0.5; % Rethreshold

    imshowpair(mat2gray(squeeze(skull_mask_unsmoothed(:,200,:))),...
        mat2gray(squeeze(skull_mask_smoothed(:,200,:))), 'montage')


    transducer_bowl = transducer_setup(parameters.transducer, trans_pos_upsampled_grid, focus_pos_upsampled_grid, ...
                                                            size(segmented_img), voxel_size_mm);

    % restrict skull based on CSF 
    % expand CSF by a large amount
    
    SE = strel('cube', parameters.csf_mask_expansion_factor/voxel_size_mm);
    csf_mask_expanded = imdilate(csf_mask, SE);

%     imshow(plot_t1_with_transducer(segmented_image, segmented_hdr, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters))
    imshowpair(squeeze(skull_mask_unsmoothed(:,trans_pos_upsampled_grid(2),:)), squeeze(csf_mask_expanded(:,trans_pos_upsampled_grid(2),:)), 'falsecolor')

    crop_margin = 6; % make sure margin is at least 1

    [min_dims, max_dims, new_grid_dims] = get_crop_dims(csf_mask_expanded+transducer_bowl, crop_margin);
    if any( max_dims > size(skull_mask_smoothed) )
        pad_amount = max( max_dims, size(segmented_img))-size(segmented_img);
        segmented_img = padarray(segmented_img, pad_amount,0,'post');
        skull_mask_smoothed = padarray(skull_mask_smoothed, pad_amount,0,'post');        
    end
    if any(min_dims < 0)
        pad_amount = abs(min(min_dims, [0 0 0]));
        
        segmented_img = padarray(segmented_img, pad_amount,0,'pre');
        skull_mask_smoothed = padarray(skull_mask_smoothed, pad_amount,0,'pre');
        min_dims = max(min_dims, [0 0 0]);
        max_dims = max_dims+pad_amount;
        
    end

    skull_mask = skull_mask_smoothed(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3)); % Crop image
    skull_edge = edge3(skull_mask, 'approxcanny',0.1);

    segmented_image_cropped = segmented_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3)); % Crop image

    parameters.grid_dims = size(skull_mask);

    trans_pos_final = trans_pos_upsampled_grid-min_dims';
    focus_pos_final = focus_pos_upsampled_grid-min_dims';
    
    translation_matrix = makehgtform('translate', -min_dims);
    
    imshowpair(plot_t1_with_transducer(skull_mask, voxel_size_mm, trans_pos_final, focus_pos_final, parameters),...
        plot_t1_with_transducer(skull_mask_smoothed, voxel_size_mm, trans_pos_upsampled_grid, focus_pos_upsampled_grid, parameters),...
    'montage')

end