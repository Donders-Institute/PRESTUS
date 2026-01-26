function [medium_masks, skull_i] = skull_fill_holes(parameters, medium_masks, labels, focus_pos_grid, segmented_img)
% skull_fill_holes Expands skull mask and fills gaps to skin.
%
% This function processes a 3D medium mask used in ultrasound brain stimulation simulations (e.g., k-Wave).
% It ensures skull continuity by expanding the skull region using bone perimeter detection and filling 
% small gaps between skull and skin, while preserving trabecular bone and excluding eyes from bone labeling.
% Debug plots are generated if enabled, showing changes at the focus position slice.
%
% Inputs:
%   parameters      - Struct with simulation parameters:
%                     - smooth_window              (numeric) Window size for smoothing bone image.
%                     - smooth_threshold_skull     (numeric) Threshold for skull smoothing.
%                     - debug                      (1/0) Enable debug plots.
%                     - debug_dir                  (string) Directory for saving debug images.
%                     - subject_id                 (numeric) Subject ID for filename.
%                     - simulation_medium          (string) Medium type (e.g., 'water').
%                     - results_filename_affix     (string) Affix for output filenames.
%                     - seg_labels.eye             (optional int) Label index for eye tissue.
%   medium_masks    - 3D array of initial medium labels (updated in-place and returned).
%   labels          - Cell array of tissue label names (e.g., {'skin', 'skull_cortical'}).
%   focus_pos_grid  - 1x3 vector [x,y,z] indices of focus position for debug slice (y-slice used).
%   segmented_img   - 3D array of original segmentation labels.
%
% Outputs:
%   medium_masks    - Updated 3D medium mask with expanded skull and filled gaps.
%   skull_i         - Integer label value for skull tissue.
%
% Dependencies:
%   - smooth_img (custom function): Smooths bone image.
%   - Image Processing Toolbox: imerode, strel, imfill, bwlabeln, regionprops, label2rgb, imshowpair, montage.
%
% Example:
%   medium_masks = skull_fill_holes(params, medium_masks, labels, 2, bone_img, focus_pos, seg_img, trab_mask, 3);
%
% Notes:
%   - Assumes voxel-based 3D head model from neuroimaging (e.g., CT/MRI segmentation).
%   - Expansion uses morphological perimeter (3-voxel cube structuring element).
%   - Gap filling identifies largest non-skin-skull blob (likely CSF/air) and ignores it.
%   - Eyes default to water (label 0) to avoid erroneous bone assignment.

    if any(contains(labels, 'skull_cortical'))  
        % treat cortical bone as the base layer
        skull_i = find(ismember(labels, {'skull_cortical'; 'skull_trabecular'}));
    else
        skull_i = find(strcmp(labels,  'skull'));
    end

    % Retain pre-fill medium mask for later plotting
    medium_with_gaps = medium_masks;

    % Ensure continuous skull (fill small holes, connect thin regions)
    skull = ismember(medium_masks, skull_i);
    se = strel('sphere', 3);  % 3D ball for isotropic smoothing
    skull_continuous = imclose(skull, se);  % Dilate then erode
    skull_new = skull_continuous & ~skull;
    medium_masks(skull_new) = skull_i(1);  % Only add new voxels as cortical bone (if differentiated)

    % [DEBUG] Plot skull expansion at focus y-slice
    if parameters.debug == 1
        h = figure;
        montage({1-squeeze(skull(:,focus_pos_grid(2),:)), ...
            1-squeeze(skull_new(:,focus_pos_grid(2),:)), ...
            1-squeeze(skull_continuous(:,focus_pos_grid(2),:))}, 'Size', [1 3]);
        title('Skull (left), closed bone (center), continuous skull (right)');
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_skull_expansion%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png');
        close(h);
    end

    % Fill gaps between skull and skin (if skin label exists)
    if any(strcmp(labels, 'skin'))
        skin_i = find(strcmp(labels, 'skin'));
        skin = medium_masks == skin_i;
        skull = ismember(medium_masks, skull_i);
        skin_skull = skin | skull;
        % 3D hole-filling
        se = strel('sphere', 3);
        skin_skull_filled = imclose(skin_skull, se);  % Close small gaps + fill
        skin_skull_filled = imfill(skin_skull_filled, 'holes');  % Fill remaining

        % Preserve trabecular mask (if available)
        if any(contains(labels, 'skull_trabecular'))
            trabecular_i = find(strcmp(labels, 'skull_trabecular'));
            trabecular_mask = medium_masks(medium_masks==trabecular_i);
        end
        % update medium mask (with cortical skull in differentiated case)
        medium_masks((skin_skull_filled - skin_skull) > 0) = skull_i(1);
        % re-insert trabecular mask (if available)
        if any(contains(labels, 'skull_trabecular'))
            medium_masks(trabecular_mask ~= 0) = trabecular_i;
        end
    end % skull-skin

    % Label eye layer as water
    if isfield(parameters.seg_labels, 'eye')
        eye_i = parameters.seg_labels.eye;
        eye = segmented_img == eye_i;
        medium_masks(eye ~= 0) = 0; % Default to water
    end

    % [DEBUG] Plot before/after gap filling at focus y-slice
    if parameters.debug == 1
        medium_mask_updates = double(medium_masks-medium_with_gaps);
        medium_mask_updates(medium_mask_updates ~=0) = 1; % indicate where changes occured
        h = figure;
        imshowpair(label2rgb(squeeze(medium_mask_updates(:,focus_pos_grid(2),:)), 'parula'), ...
            label2rgb(squeeze(medium_masks(:,focus_pos_grid(2),:)), 'parula'), 'montage');
        title('Updated (left) and closed off (right) segmented images');
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_segmented_img_closing_gaps_changes%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png');
        close(h);
    end
end
