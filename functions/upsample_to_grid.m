function [upsampled_image, transformation_matrix, trans_pos_new, focus_pos_new] = ...
    upsample_to_grid(nii_image, previous_transformation_matrix, resample_factor, trans_pos_grid, focus_pos_grid)

% UPSAMPLE_TO_GRID Upsamples a 3D image to a higher resolution grid.
%
% This function resamples a 3D image (`nii_image`) to a higher resolution grid 
% using a specified resample factor. It also updates the transformation matrix 
% and computes the new positions for the transducer and focus points in the upsampled grid.
%
% Input:
%   nii_image                   - [Nx x Ny x Nz] matrix representing the 3D image to be upsampled.
%   previous_transformation_matrix - [4x4] affine transformation matrix used for the previous grid.
%   resample_factor             - Scalar specifying the resampling factor (e.g., 2 for doubling resolution).
%   trans_pos_grid              - [nx3] array specifying the transducer position in the original grid.
%   focus_pos_grid              - [nx3] array specifying the focus position in the original grid.
%
% Output:
%   upsampled_image             - [Mx x My x Mz] matrix representing the upsampled 3D image.
%   transformation_matrix       - [4x4] affine transformation matrix for the new grid.
%   trans_pos_new               - [nx3] array specifying the transducer position in the upsampled grid.
%   focus_pos_new               - [nx3] array specifying the focus position in the upsampled grid.

    % Update transformation matrix to include resampling
    transformation_matrix = diag([repmat(resample_factor, [1 3]), 1]) * previous_transformation_matrix;

    % Compute new dimensions of the upsampled image
    new_dims = ceil(diff(findbounds(maketform('affine', transformation_matrix), [0 0 0; size(nii_image)])));
    previous_dims = ceil(diff(findbounds(maketform('affine', previous_transformation_matrix), [0 0 0; size(nii_image)])));

    % Compute translation matrix to center the new image
    translation_matrix = makehgtform('translate', (new_dims - previous_dims) / 2);

    % Combine transformation and translation matrices
    TF = maketform('affine', transformation_matrix * translation_matrix');

    % Perform resampling using nearest-neighbor interpolation
    upsampled_image = tformarray(nii_image, TF, makeresampler('nearest', 'fill'), ...
        [1 2 3], [1 2 3], new_dims, [], 0);

    % Compute new positions for transducer and focus points in the upsampled grid
    nT = size(trans_pos_grid, 1);
    
    trans_pos_new = zeros(size(trans_pos_grid));
    focus_pos_new = zeros(size(focus_pos_grid));
    
    for i = 1:nT
        out_mat = round(tformfwd([trans_pos_grid(i, :); focus_pos_grid(i, :)], TF));        % [2n × 3]
    
        trans_pos_new(i, :) = out_mat(1, :);     % n×3
        focus_pos_new(i, :) = out_mat(2, :);     % n×3
    end

end