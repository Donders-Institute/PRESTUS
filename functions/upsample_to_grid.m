function [upsampled_image, transformation_matrix, trans_pos_new, focus_pos_new] = ...
    upsample_to_grid(nii_image, previous_transformation_matrix, resample_factor, trans_pos_grid, focus_pos_grid)

    transformation_matrix = diag([repmat(resample_factor, [1 3]), 1])*previous_transformation_matrix;
    % Compute new dimensions
    new_dims = ceil(diff(findbounds(maketform('affine', transformation_matrix), [0 0 0; size(nii_image)])));
    previous_dims = ceil(diff(findbounds(maketform('affine', previous_transformation_matrix), [0 0 0; size(nii_image)])));

    % Move to the new centre
    translation_matrix = makehgtform('translate',(new_dims - previous_dims)/2);   
    
    TF = maketform('affine', transformation_matrix*translation_matrix');    
    
    upsampled_image = tformarray(nii_image, TF, makeresampler('nearest', 'fill'), ...
        [1 2 3], [1 2 3], new_dims, [], 0) ;

    % and the new positions for the transducer and the focus can be computed
    out_mat = round(tformfwd([trans_pos_grid focus_pos_grid]', TF));

    trans_pos_new = out_mat(1,:)';
    focus_pos_new = out_mat(2,:)';
end