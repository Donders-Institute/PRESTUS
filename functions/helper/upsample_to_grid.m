function [upsampled_image, transformation_matrix, trans_pos_new, focus_pos_new] = ...
    upsample_to_grid(nii_image, previous_transformation_matrix, resample_factor, trans_pos_grid, focus_pos_grid)

% UPSAMPLE_TO_GRID  Resample a 3D image and update grid positions to a higher-resolution grid
%
% Scales the affine transform by resample_factor, computes new grid
% dimensions, resamples with nearest-neighbour interpolation, and
% transforms transducer/focus grid positions to the new coordinates.
%
% Use as:
%   [upsampled_image, transformation_matrix, trans_pos_new, focus_pos_new] = ...
%       upsample_to_grid(nii_image, previous_transformation_matrix, ...
%                        resample_factor, trans_pos_grid, focus_pos_grid)
%
% Input:
%   nii_image                      - [Nx x Ny x Nz] 3D image to resample
%   previous_transformation_matrix - [4x4] affine transform for the original grid
%   resample_factor                - resampling factor (e.g., 2 to double resolution)
%   trans_pos_grid                 - [Nx3] transducer positions in the original grid [voxels]
%   focus_pos_grid                 - [Nx3] focus positions in the original grid [voxels]
%
% Output:
%   upsampled_image       - [Mx x My x Mz] resampled 3D image
%   transformation_matrix - [4x4] affine transform for the new grid
%   trans_pos_new         - [Nx3] transducer positions in the upsampled grid [voxels]
%   focus_pos_new         - [Nx3] focus positions in the upsampled grid [voxels]
%
% See also: TFORMARRAY, MAKETFORM

arguments
    nii_image                      (:,:,:) {mustBeNumeric}
    previous_transformation_matrix (4,4)   {mustBeNumeric}
    resample_factor                (1,1)   {mustBeNumeric}
    trans_pos_grid                 (:,3)   {mustBeNumeric}
    focus_pos_grid                 (:,3)   {mustBeNumeric}
end

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