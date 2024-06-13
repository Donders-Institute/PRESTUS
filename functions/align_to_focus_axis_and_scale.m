function [rotated_img, trans_pos_new, focus_pos_new, transformation_matrix, rotation_matrix, angle_x_rad, angle_y_rad, montage_img] = ...
    align_to_focus_axis_and_scale(nii_image, nii_header, trans_pos_grid, focus_pos_grid, scale_factor, parameters)
    arguments
        nii_image (:,:,:)
        nii_header struct
        trans_pos_grid (1, 3)
        focus_pos_grid (1, 3)
        scale_factor double
        parameters struct
    end

    % What we need to do is to align the focal axis with the coordinate system.
    % This is done by rotating the image so that the first two coordinates of
    % the transducer and the focus are the same, the focal axis becomes the
    % z-axis. 

    % First we create the transformation (affine) matrix for the voxel array.

    focal_axis = [focus_pos_grid - trans_pos_grid, 1]';

    % First, the transformation of matrix (3d voxel array) indices.
    % The transformation is done by two sequential rotations to align the focal
    % axis with the z axis. 

    % Rotate in Y first
    % Compute the axis angle
    angle_y_rad  = atan(focal_axis(1) / focal_axis(3)); 
    if angle_y_rad < 0
        angle_y_rad = angle_y_rad + pi;
    end
    % Create an affine matrix for Y
    vox_transf_Y  = makehgtform('yrotate', -angle_y_rad); 

    % Compute rotated axis
    focal_axis_new = vox_transf_Y*focal_axis;

    % Then rotate in X
    angle_x_rad  = atan(focal_axis_new(2) / focal_axis_new(3)); 
    % Create another affine matrix for X
    vox_transf_X  = makehgtform('xrotate', angle_x_rad);

    % Combine transformations
    rotation_matrix = vox_transf_X * vox_transf_Y;

    % Adjust according to the scaling factor
    scale_matrix = makehgtform('scale', scale_factor);
    
    % Combine the scaling and rotation matrices
    rotation_and_scale_matrix = rotation_matrix*scale_matrix;
    
    % The rotated image size would be different from original;
    % here we find out the new dimensions
    TF = maketform('affine', rotation_and_scale_matrix');
    newdims = ceil(diff(findbounds(TF, [0 0 0; size(nii_image)])));

    % Now to keep the image in the field of view, it's better to rotate around
    % its center. By default, the transformations are done around the origin 
    % point (1, 1), so we first translate the image center to the origin point
    rotation_center = (size(nii_image)+1)/2;
    forward = makehgtform('translate', -rotation_center);

    % Then we translate it back to the center of the new dimensions
    backward = makehgtform('translate', (newdims+1)/2);

    % And now all the transformations are combined
    transformation_matrix = forward'*rotation_and_scale_matrix'*backward';

    % Then the image could be transformed
    TF = maketform('affine', transformation_matrix);

    rotated_img = tformarray(nii_image, TF, makeresampler('nearest', 'fill'), ...
        [1 2 3], [1 2 3], newdims, [], 0) ;

    % And the new positions for the transducer and the focus can be computed
    out_mat = round(tformfwd([trans_pos_grid; focus_pos_grid], TF));

    trans_pos_new = out_mat(1,:);
    focus_pos_new = out_mat(2,:);
    
    % Create plots of the original T1 image with the transducer in addition
    % to the rotated T1 image
    orig_with_transducer_img = plot_t1_with_transducer(nii_image, nii_header.PixelDimensions(1), trans_pos_grid, ...
            focus_pos_grid, parameters);
    rotated_with_transducer_img = plot_t1_with_transducer(rotated_img, nii_header.PixelDimensions(1)/scale_factor, trans_pos_new, focus_pos_new, parameters);

    % Present the two figures side-by-side
    montage_img = imtile({rotated_with_transducer_img, orig_with_transducer_img}, 'GridSize', [1 nan], 'BackgroundColor', [0.5 0.5 0.5]);
    
end