function [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, parameters, trans_pos_final, ...
         focus_pos_final, translation_matrix] = preproc_crop_grid( ...
         parameters, medium_masks, segmentation, bone_mask_img, trans_pos_grid, focus_pos_grid, pseudoCT_img)
% PREPROC_CROP_GRID  Crop 3-D head model to the simulation-relevant region
%
% Crops the full-head model to a CSF-guided, transducer-inclusive bounding
% box, adds symmetric user padding, optimises dimensions for FFT, and
% applies conditional pre/post padding to prevent out-of-bounds access.
% The PML is always applied outside the returned grid (PMLInside=false).
% Symmetric padding (parameters.headmodel.head_pad_mm) is applied BEFORE
% transducer_setup to prevent out-of-bounds access. Grid dimensions are
% then adjusted to be FFT-efficient. Skull edges are computed to derive
% the crop bounding box.
%
% Use as:
%   [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, parameters, trans_pos_final, ...
%    focus_pos_final, translation_matrix] = preproc_crop_grid( ...
%       parameters, medium_masks, segmentation, bone_mask_img, trans_pos_grid, focus_pos_grid)
%   [medium_masks, segmentation_crop, bone_mask_crop, pseudoCT_crop, parameters, trans_pos_final, ...
%    focus_pos_final, translation_matrix] = preproc_crop_grid( ...
%       parameters, medium_masks, segmentation, bone_mask_img, trans_pos_grid, focus_pos_grid, pseudoCT_img)
%
% Input:
%   parameters     - (1,1) simulation configuration struct
%   medium_masks   - [Nx x Ny x Nz] medium label array
%   segmentation   - [Nx x Ny x Nz] tissue segmentation
%   bone_mask_img  - [Nx x Ny x Nz] binary skull mask
%   trans_pos_grid - [1x3] transducer position in voxel coordinates
%   focus_pos_grid - [1x3] focus position in voxel coordinates
%   pseudoCT_img   - [Nx x Ny x Nz] Hounsfield-unit skull image (optional, [] if unused)
%
% Output:
%   medium_masks       - cropped medium label array
%   segmentation_crop  - cropped tissue segmentation
%   bone_mask_crop     - cropped binary skull mask
%   pseudoCT_crop      - cropped pseudoCT image ([] when pseudoCT_img is not provided)
%   parameters         - updated struct with .grid.dims = [Nx Ny Nz]
%   trans_pos_final    - [1x3] transducer position in cropped grid
%   focus_pos_final    - [1x3] focus position in cropped grid
%   translation_matrix - [4x4] homogeneous translation matrix
%
% See also: PREPROC_CROP_ECSF, HEAD_SMOOTH_AND_CROP

if nargin < 7; pseudoCT_img = []; end

% === USER SYMMETRIC PADDING (BOTH SIDES) ===
if ~isfield(parameters.headmodel, 'head_pad_mm') || isempty(parameters.headmodel.head_pad_mm)
    parameters.headmodel.head_pad_mm = 0;  % mm
end
pad_voxels = round(parameters.headmodel.head_pad_mm / parameters.grid.resolution_mm);
pad_pre_post = [pad_voxels, pad_voxels, pad_voxels];  % Symmetric BOTH sides

% Track TOTAL pre-padding offset for translation matrix
total_pre_offset = zeros(1,3);

% Apply symmetric padding BEFORE transducer_setup
if any(pad_pre_post > 0)
    segmentation    = padarray(segmentation,    pad_pre_post, 0, 'both');
    medium_masks    = padarray(medium_masks,    pad_pre_post, 0, 'both');
    bone_mask_img   = padarray(bone_mask_img,   pad_pre_post, 0, 'both');
    if ~isempty(pseudoCT_img)
        pseudoCT_img = padarray(pseudoCT_img,   pad_pre_post, NaN, 'both');
    end
    total_pre_offset = total_pre_offset + pad_pre_post;  % User padding contribution
    % Positions shift by PRE-padding amount ('both' adds pre first)
    trans_pos_grid = trans_pos_grid + pad_pre_post;
    focus_pos_grid = focus_pos_grid + pad_pre_post;
end

% Expand CSF mask to crop the layered medium
[medium_masks] = preproc_crop_eCSF(parameters, medium_masks, segmentation, trans_pos_grid);

% Minimal fixed crop margin — PML is added outside the grid by k-Wave (PMLInside=false)
crop_margin = 1;

% Include transducer bowl geometry (safe with padding)
i_water = find(strcmp(fieldnames(parameters.medium_properties), 'water'));
transducer_bowl = transducer_setup(parameters.transducer(1), trans_pos_grid, ...
    focus_pos_grid, size(segmentation), parameters.grid.resolution_mm, parameters);

% Compute crop bounds
orig_dims = size(medium_masks);
combinedmask = (medium_masks ~= i_water) | transducer_bowl;
[min_dims, max_dims, new_grid_dims] = get_crop_dims(double(combinedmask), crop_margin);
clear combinedmask;

% === CONDITIONAL PRE-PADDING if min_dims < 1 ===
if any(min_dims < 1)
    pad_amount = abs(min(min_dims, [1 1 1]));
    segmentation  = padarray(segmentation,  pad_amount, 0,       'pre');
    medium_masks  = padarray(medium_masks,  pad_amount, i_water, 'pre');
    bone_mask_img = padarray(bone_mask_img, pad_amount, 0,       'pre');
    if ~isempty(pseudoCT_img)
        pseudoCT_img = padarray(pseudoCT_img, pad_amount, NaN, 'pre');
    end
    total_pre_offset = total_pre_offset + pad_amount;  % Add conditional padding
    min_dims = max(min_dims, [1 1 1]);
    max_dims = max_dims + pad_amount;
    new_grid_dims = max_dims - min_dims + 1;
    trans_pos_grid = trans_pos_grid + pad_amount;
    focus_pos_grid = focus_pos_grid + pad_amount;
end

% Optimize dimensions for FFT
new_grid_dims(1) = find_min_factor(new_grid_dims(1), new_grid_dims(1) + parameters.grid.max_expand);
new_grid_dims(2) = find_min_factor(new_grid_dims(2), new_grid_dims(2) + parameters.grid.max_expand);
new_grid_dims(3) = find_min_factor(new_grid_dims(3), new_grid_dims(3) + parameters.grid.max_expand);
max_dims = min_dims + new_grid_dims - 1;

% === CONDITIONAL POST-PADDING if FFT expansion exceeds padding ===
if any(max_dims > size(medium_masks))
    pad_post_amount = max(0, max_dims - size(medium_masks));
    segmentation  = padarray(segmentation,  pad_post_amount, 0,       'post');
    medium_masks  = padarray(medium_masks,  pad_post_amount, i_water, 'post');
    bone_mask_img = padarray(bone_mask_img, pad_post_amount, 0,       'post');
    if ~isempty(pseudoCT_img)
        pseudoCT_img = padarray(pseudoCT_img, pad_post_amount, NaN, 'post');
    end
    fprintf('Post-padding applied: [%d %d %d] voxels\n', pad_post_amount);
end

% Apply final crop
medium_masks    = medium_masks(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
segmentation_crop = segmentation(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
bone_mask_crop  = bone_mask_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
if ~isempty(pseudoCT_img)
    pseudoCT_crop = pseudoCT_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
else
    pseudoCT_crop = [];
end

% Update parameters
parameters.grid.dims = size(medium_masks);

% Final positions in cropped grid coordinates
trans_pos_final = trans_pos_grid - min_dims + 1;
focus_pos_final = focus_pos_grid - min_dims + 1;

% === COMPLETE TRANSLATION MATRIX ===
% Maps: original_world → final_cropped_world
% Offset = total_pre_padding + crop_offset
translation_matrix = makehgtform('translate', total_pre_offset + (1-min_dims));

% Display summary
fprintf('Grid: [%dx%dx%d] → [%dx%dx%d]\n', orig_dims, parameters.grid.dims);
fprintf('  User padding: %.1fmm = [%d %d %d] voxels\n', ...
        parameters.headmodel.head_pad_mm, pad_pre_post);
fprintf('  Total pre-offset: [%d %d %d] voxels\n', total_pre_offset);
fprintf('  Crop bounds: min=[%d %d %d], max=[%d %d %d]\n', min_dims, max_dims);

end