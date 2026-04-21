function [rotated_img, trans_pos_new, focus_pos_new, transformation_matrix, rotation_matrix, angle_x_rad, angle_y_rad, montage_img] = ...
    preproc_align_to_focal_axis(nii_image, nii_header, trans_pos_grid, focus_pos_grid, scale_factor, parameters)

% PREPROC_ALIGN_TO_FOCAL_AXIS  Rotate a 3D image so the focal axis aligns with the z-axis
%
% Rotates and optionally scales a 3D image such that the vector from
% trans_pos_grid to focus_pos_grid becomes aligned with the z-axis.
% Returns updated transducer/focus positions and the full affine transform.
%
% Use as:
%   [rotated_img, trans_pos_new, focus_pos_new, transformation_matrix, ...
%    rotation_matrix, angle_x_rad, angle_y_rad, montage_img] = ...
%       preproc_align_to_focal_axis(nii_image, nii_header, trans_pos_grid, ...
%           focus_pos_grid, scale_factor, parameters)
%
% Input:
%   nii_image      - [Nx x Ny x Nz] 3D image volume
%   nii_header     - NIfTI header struct (from niftiinfo)
%   trans_pos_grid - [1x3] transducer position in voxel coordinates
%   focus_pos_grid - [1x3] focus position in voxel coordinates
%   scale_factor   - image scaling factor (1 = no scaling)
%   parameters     - (1,1) simulation parameters struct
%
% Output:
%   rotated_img          - rotated and scaled image volume
%   trans_pos_new        - [1x3] transducer position after transformation
%   focus_pos_new        - [1x3] focus position after transformation
%   transformation_matrix - [4x4] combined affine transformation matrix
%   rotation_matrix      - [4x4] rotation-only matrix
%   angle_x_rad          - rotation angle around x-axis [rad]
%   angle_y_rad          - rotation angle around y-axis [rad]
%   montage_img          - comparison montage image
%
% Algorithm:
%   1. Focal vector d = focus_pos_grid - trans_pos_grid (homogeneous).
%   2. theta_y = atan2(d(1), d(3)); R_y(-theta_y) via makehgtform.
%   3. d' = R_y * d; theta_x = atan2(d'(2), d'(3)); R_x(theta_x).
%   4. R = R_x * R_y; S = scale; T_rs = R * S.
%   5. Bounds from tformfwd on corners; newdims = ceil(max-min).
%   6. Center c = (sz+1)/2; forward = translate(-c); backward = translate((newdims+1)/2).
%   7. T = forward' * T_rs' * backward'; rotated_img = tformarray(nii_image, T).
%   Assumes right-handed coordinates; angles ensure positive rotation.
%   Requires Image Processing Toolbox (makehgtform, tformarray).
%
% See also: PREPROC_HEAD, PREPROC_CROP_GRID

    arguments
        nii_image      (:,:,:) {mustBeNumeric}
        nii_header     (1,1) struct
        trans_pos_grid (1,3) double
        focus_pos_grid (1,3) double
        scale_factor   (1,1) double
        parameters     (1,1) struct
    end

    % Ensure inputs are on CPU (tformarray and other image processing functions don't support gpuArray)
    nii_image = gather(nii_image);
    trans_pos_grid = gather(trans_pos_grid);
    focus_pos_grid = gather(focus_pos_grid);

    %% Step 1: Compute focal axis
    % The focal axis is defined as a vector from `trans_pos_grid` to `focus_pos_grid`.
    % [Multi-transducer] the preprocessing will be based on the first transducer
    if parameters.transducer(1).align_to_focus 
        focal_axis = [focus_pos_grid - trans_pos_grid, 1]';
    else
        % Use the natural focus based on transducer curvature
        curv_radius = parameters.transducer(1).(parameters.transducer(1).type).curv_radius_mm;
        natural_focus = trans_pos_grid + [0, 0, curv_radius / parameters.grid_step_mm];
        focal_axis = [natural_focus - trans_pos_grid, 1]';
    end

    %% Step 2: Compute rotation angles
    % Rotate around the y-axis first to align focal axis with z-axis projection in x-z plane
    angle_y_rad  = atan(focal_axis(1) / focal_axis(3)); 
    if angle_y_rad < 0
        angle_y_rad = angle_y_rad + pi; % Ensure positive angles
    end
    % Create an affine matrix for y-rotation
    vox_transf_Y  = makehgtform('yrotate', -angle_y_rad); 

    % Compute rotated focal axis after y-rotation
    focal_axis_new = vox_transf_Y * focal_axis;

    % Rotate around x-axis to align remaining projection with z-axis
    angle_x_rad  = atan(focal_axis_new(2) / focal_axis_new(3)); 
    % Create an affine matrix for x-rotation
    vox_transf_X  = makehgtform('xrotate', angle_x_rad);

    %% Step 3: Combine rotation matrices
    % Combine sequential rotations into a single rotation matrix
    rotation_matrix = vox_transf_X * vox_transf_Y;

    %% Step 4: Apply scaling
    % Create a scaling matrix based on `scale_factor`
    scale_matrix = makehgtform('scale', scale_factor);

    % Combine scaling and rotation into a single transformation matrix
    rotation_and_scale_matrix = rotation_matrix * scale_matrix;

    %% Step 5: Compute new image dimensions
    % Determine new image dimensions after applying transformations.
    % We must ensure the new grid is large enough to contain the image AND the transducer/focus,
    % while keeping the image center at the center of the new grid.
    
    % Translate image center to origin before applying transformations
    rotation_center = (size(nii_image) + 1) / 2;
    forward = makehgtform('translate', -rotation_center);

    % Create the transform that centers the input, then rotates/scales
    % forward: translates input center to origin
    % rotation_and_scale_matrix: rotates/scales around origin
    T_centered = forward' * rotation_and_scale_matrix';
    TF_centered = maketform('affine', T_centered);
    
    % Get bounds of the transformed image corners
    img_bounds = findbounds(TF_centered, [0 0 0; size(nii_image)]);
    
    % Get transformed positions of transducer and focus
    tp_trans = tformfwd(trans_pos_grid, TF_centered);
    fp_trans = tformfwd(focus_pos_grid, TF_centered);
    
    % Combine all points to find the maximum extent from the center
    all_points = [img_bounds; tp_trans; fp_trans];
    max_abs = max(abs(all_points), [], 1);
    
    % Calculate symmetric dimensions to keep the image centered while fitting everything
    % Factor of 2 ensures coverage from -max to +max
    newdims = ceil(2 * max_abs) + 2; % +2 margin for safety

    %% Step 6: Center transformations around image center
    % (forward is already calculated above)

    % Translate back to center of new dimensions after transformations
    backward = makehgtform('translate', (newdims + 1) / 2);

    % Combine all transformations into a final affine transformation matrix
    transformation_matrix = forward' * rotation_and_scale_matrix' * backward';

    %% Step 7: Transform image
    % Apply affine transformation to rotate and scale the image
    % (nii_image already ensured to be on CPU)

    if numel(unique(nii_image))<20 % if the image is a mask use nearest neighbor
        parameters.interpolation = 'nearest';
    else
        parameters.interpolation = 'linear';
    end

    rotated_img = tformarray(nii_image, maketform('affine', transformation_matrix), ...
        makeresampler(parameters.interpolation, 'fill'), [1 2 3], [1 2 3], newdims, [], 0);

    %% Step 8: Update positions of transducer and focus
    % Transform transducer and focus positions using the affine transformation matrix
    out_mat = round(tformfwd([trans_pos_grid; focus_pos_grid], ...
        maketform('affine', transformation_matrix)));
    
    trans_pos_new = out_mat(1,:); % New transducer position
    focus_pos_new = out_mat(2,:); % New focus position

    %% Step 9: Generate visualization montage
    % Create plots of original and transformed images with overlaid transducer positions
    orig_with_transducer_img = plot_t1_with_transducer(nii_image, nii_header.PixelDimensions(1), ...
        trans_pos_grid, focus_pos_grid, parameters);
    
    rotated_with_transducer_img = plot_t1_with_transducer(rotated_img, ...
        nii_header.PixelDimensions(1) / scale_factor, trans_pos_new, focus_pos_new, parameters);

    % Combine original and transformed images into a side-by-side montage
    montage_img = imtile({orig_with_transducer_img, rotated_with_transducer_img}, ...
        'GridSize', [1 nan], 'BackgroundColor', [0.5 0.5 0.5]);
end
