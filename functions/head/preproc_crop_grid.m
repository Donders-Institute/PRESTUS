function [medium_masks, segmentation_crop, bone_crop, parameters, trans_pos_final, ...
         focus_pos_final, translation_matrix] = preproc_crop_grid( ...
         parameters, medium_masks, segmentation, bone_img, trans_pos_grid, focus_pos_grid)
% preproc_crop_grid Crops 3D head model grid to simulation-relevant region for k-Wave efficiency.
%
% Crops full-head model to CSF-guided + transducer-inclusive bounding box, excluding distant bone.
% Adds user-controlled symmetric padding BEFORE transducer_setup to prevent out-of-bounds.
% Optimizes grid dimensions for FFT performance. Computes skull edges and adjustments.
%
% Inputs:
%   parameters      - Struct with .pad_mm (symmetric padding mm), .grid_step_mm, .pml_size,
%                     .csf_mask_expansion_factor, .grid_max_expand, .debug, etc.
%   medium_masks    - 3D array: Current medium labels (non-CSF regions zeroed).
%   segmentation    - 3D array: Original segmentation/pseudoCT.
%   bone_img        - 3D array: Bone property image.
%   trans_pos_grid  - 1x3: Transducer position [x y z] grid indices.
%   focus_pos_grid  - 1x3: Focus position [x y z] grid indices.
%
% Outputs:
%   medium_masks          - Cropped 3D medium labels.
%   segmentation_crop     - Cropped original segmentation.
%   bone_crop             - Cropped bone image.
%   parameters            - Updated with .grid_dims = [Nx Ny Nz].
%   trans_pos_final       - 1x3: Adjusted transducer pos in cropped grid.
%   focus_pos_final       - 1x3: Adjusted focus pos in cropped grid.
%   translation_matrix    - 4x4: Homogeneous translation (original → cropped coordinates).
%

% === USER SYMMETRIC PADDING (BOTH SIDES) ===
if ~isfield(parameters, 'pad_mm') || isempty(parameters.pad_mm)
    parameters.pad_mm = 0;  % mm
end
pad_voxels = round(parameters.pad_mm / parameters.grid_step_mm);
pad_pre_post = [pad_voxels, pad_voxels, pad_voxels];  % Symmetric BOTH sides

% Track TOTAL pre-padding offset for translation matrix
total_pre_offset = zeros(1,3);

% Apply symmetric padding BEFORE transducer_setup
if any(pad_pre_post > 0)
    segmentation = padarray(segmentation, pad_pre_post, 0, 'both');
    medium_masks = padarray(medium_masks, pad_pre_post, 0, 'both');
    bone_img = padarray(bone_img, pad_pre_post, 0, 'both');
    total_pre_offset = total_pre_offset + pad_pre_post;  % User padding contribution
    % Positions shift by PRE-padding amount ('both' adds pre first)
    trans_pos_grid = trans_pos_grid + pad_pre_post;
    focus_pos_grid = focus_pos_grid + pad_pre_post;
end

% Expand CSF mask to crop the layered medium
if isfield(parameters, 'seg_labels') && any(strcmp(fieldnames(parameters.seg_labels), 'csf'))
    [medium_masks] = preproc_crop_eCSF(parameters, medium_masks, segmentation, trans_pos_grid);
else
    warning("CSF layer not specified or unknown ... will not use expanded CSF mask...");
end

% Set PML buffer as crop margin
crop_margin = parameters.pml_size + 1;

% Include transducer bowl geometry (safe with padding)
i_water = find(strcmp(fieldnames(parameters.medium), 'water'));
transducer_bowl = transducer_setup(parameters.transducer(1), trans_pos_grid, ...
    focus_pos_grid, size(segmentation), parameters.grid_step_mm);

% Compute crop bounds
orig_dims = size(medium_masks);
combinedmask = (medium_masks ~= i_water) | transducer_bowl;
[min_dims, max_dims, new_grid_dims] = get_crop_dims(double(combinedmask), crop_margin);
clear combinedmask;

% === CONDITIONAL PRE-PADDING if min_dims < 1 ===
if any(min_dims < 1)
    pad_amount = abs(min(min_dims, [1 1 1]));
    segmentation = padarray(segmentation, pad_amount, 0, 'pre');
    medium_masks = padarray(medium_masks, pad_amount, i_water, 'pre');
    bone_img = padarray(bone_img, pad_amount, 0, 'pre');
    total_pre_offset = total_pre_offset + pad_amount;  % Add conditional padding
    min_dims = max(min_dims, [1 1 1]);
    max_dims = max_dims + pad_amount;
    new_grid_dims = max_dims - min_dims + 1;
    trans_pos_grid = trans_pos_grid + pad_amount;
    focus_pos_grid = focus_pos_grid + pad_amount;
end

% Optimize dimensions for FFT
new_grid_dims(1) = find_min_factor(new_grid_dims(1), new_grid_dims(1) + parameters.grid_max_expand);
new_grid_dims(2) = find_min_factor(new_grid_dims(2), new_grid_dims(2) + parameters.grid_max_expand);
new_grid_dims(3) = find_min_factor(new_grid_dims(3), new_grid_dims(3) + parameters.grid_max_expand);
max_dims = min_dims + new_grid_dims - 1;

% === CONDITIONAL POST-PADDING if FFT expansion exceeds padding ===
if any(max_dims > size(medium_masks))
    pad_post_amount = max(0, max_dims - size(medium_masks));
    segmentation = padarray(segmentation, pad_post_amount, 0, 'post');
    medium_masks = padarray(medium_masks, pad_post_amount, i_water, 'post');
    bone_img = padarray(bone_img, pad_post_amount, 0, 'post');
    fprintf('Post-padding applied: [%d %d %d] voxels\n', pad_post_amount);
end

% Apply final crop
medium_masks = medium_masks(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
segmentation_crop = segmentation(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
bone_crop = bone_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));

% Update parameters
parameters.grid_dims = size(medium_masks);

% Final positions in cropped grid coordinates
trans_pos_final = trans_pos_grid - min_dims + 1;
focus_pos_final = focus_pos_grid - min_dims + 1;

% === COMPLETE TRANSLATION MATRIX ===
% Maps: original_world → final_cropped_world
% Offset = total_pre_padding + crop_offset
translation_matrix = makehgtform('translate', total_pre_offset + (1-min_dims));

% Display summary
fprintf('Grid: [%dx%dx%d] → [%dx%dx%d]\n', orig_dims, parameters.grid_dims);
fprintf('  User padding: %.1fmm = [%d %d %d] voxels\n', ...
        parameters.pad_mm, pad_pre_post);
fprintf('  Total pre-offset: [%d %d %d] voxels\n', total_pre_offset);
fprintf('  Crop bounds: min=[%d %d %d], max=[%d %d %d]\n', min_dims, max_dims);

end