function [medium_masks, segmented_image_cropped, parameters, trans_pos_final, ...
         focus_pos_final, translation_matrix] = preproc_crop_grid( ...
         parameters, medium_masks, segmented_img, trans_pos_grid, focus_pos_grid)
% preproc_crop_grid Crops 3D head model grid to simulation-relevant region for k-Wave efficiency.
%
% Crops full-head model to CSF-guided + transducer-inclusive bounding box, excluding distant bone.
% Optimizes grid dimensions for FFT performance, pads as needed, computes skull edges and adjustments.
%
% Inputs:
%   parameters      - Struct with:
%                     .csf_mask_expansion_factor (mm), .pml_size, .grid_max_expand,
%                     .debug (1/0), .debug_dir, .subject_id, .simulation_medium, .results_filename_affix,
%                     .parameters.grid_step_mm (voxel size, used elsewhere).
%   medium_masks    - 3D array: Current medium labels (non-CSF regions zeroed).
%   csf_mask        - 3D binary: CSF mask (expanded for cropping guide).
%   segmented_img   - 3D array: Original segmentation/pseudoCT.
%   trans_pos_grid  - 1x3: Transducer position [x y z] grid indices.
%   focus_pos_grid  - 1x3: Focus position [x y z] grid indices.
%
% Outputs:
%   medium_masks          - Cropped 3D medium labels.
%   segmented_image_cropped - Cropped original segmentation.
%   skull_edge            - 3D binary: Skull edges (edge3 Canny).
%   parameters            - Updated with .grid_dims = [Nx Ny Nz].
%   trans_pos_final       - 1x3: Adjusted transducer pos in cropped grid.
%   focus_pos_final       - 1x3: Adjusted focus pos in cropped grid.
%   translation_matrix    - 4x4: Homogeneous translation (-min_dims).
%
% Dependencies:
%   Custom: get_crop_dims, transducer_setup, find_min_factor, makehgtform.
%   Image Processing: strel, imdilate, padarray, edge3, imshowpair.

    % Expand CSF mask to crop the layered medium
    if isfield(parameters, 'seg_labels') && any(strcmp(fieldnames(parameters.seg_labels), 'csf'))
        [medium_masks] = preproc_crop_eCSF(parameters, medium_masks, segmented_img, trans_pos_grid);
    else
        warning("CSF layer not specified or unknown ... will not use expanded CSF mask...");
    end

    % Set PML buffer as crop margin
    crop_margin = parameters.pml_size + 1;
    
    % Include transducer bowl geometry
    transducer_bowl = transducer_setup(parameters.transducer(1), trans_pos_grid, ...
        focus_pos_grid, size(segmented_img), parameters.grid_step_mm);

    % Compute crop bounds
    orig_dims = size(medium_masks);
    combinedmask = (medium_masks ~= 0) | transducer_bowl;  % Logical union of medium mask and transducer bowl
    [min_dims, max_dims, new_grid_dims] = get_crop_dims(double(combinedmask), crop_margin);
    clear combinedmask;

    % Pad if min_dims inside crop margin
    if any(min_dims < 1)
        pad_amount = abs(min(min_dims, [1 1 1]));
        segmented_img = padarray(segmented_img, pad_amount, 0, 'pre');
        medium_masks = padarray(medium_masks, pad_amount, 0, 'pre');
        min_dims = max(min_dims, [1 1 1]);
        max_dims = max_dims + pad_amount;
        new_grid_dims = max_dims - min_dims + 1;
    end

    % Optimize dimensions for FFT
    new_grid_dims(1) = find_min_factor(new_grid_dims(1), new_grid_dims(1) + parameters.grid_max_expand);
    new_grid_dims(2) = find_min_factor(new_grid_dims(2), new_grid_dims(2) + parameters.grid_max_expand);
    new_grid_dims(3) = find_min_factor(new_grid_dims(3), new_grid_dims(3) + parameters.grid_max_expand);
    max_dims = min_dims + new_grid_dims - 1;
    
    % Pad post if max_dims exceeds
    if any(max_dims > size(medium_masks))
        pad_amount = max(max_dims, size(segmented_img)) - size(segmented_img);
        segmented_img = padarray(segmented_img, pad_amount, 0, 'post');
        medium_masks = padarray(medium_masks, pad_amount, 0, 'post');
    end

    % Apply crop (create new grids)
    medium_masks = medium_masks(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));
    segmented_image_cropped = segmented_img(min_dims(1):max_dims(1), min_dims(2):max_dims(2), min_dims(3):max_dims(3));

    % Update parameters
    parameters.grid_dims = size(medium_masks);

    % Adjust positions
    trans_pos_final = trans_pos_grid - min_dims + 1;
    focus_pos_final = focus_pos_grid - min_dims + 1;
    
    % Translation matrix
    translation_matrix = makehgtform('translate', 1-min_dims);
    
    % Display adjusted parameters
    if orig_dims ~= new_grid_dims
        fprintf('Adjusted grid ... \n');
        fprintf('Grid dims: [%d %d %d] → [%d %d %d] \n', ...
                orig_dims, new_grid_dims);
    end

end