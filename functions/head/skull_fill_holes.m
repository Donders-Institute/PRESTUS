function [medium_masks, skull_i] = skull_fill_holes(parameters, medium_masks, focus_pos_grid, segmented_img)
% skull_fill_holes Expands skull mask and fills gaps to skin.
%
% This function processes a 3D medium mask used in ultrasound brain stimulation simulations (e.g., k-Wave).
% It ensures skull continuity by expanding the skull region using bone perimeter detection and filling 
% gaps, while preserving trabecular bone and excluding eyes from bone labeling.
% Debug plots are generated if enabled, showing changes at the focus position slice.
%
% Inputs:
%   parameters      - Struct with simulation parameters:
%                     - debug                      (1/0) Enable debug plots.
%                     - debug_dir                  (string) Directory for saving debug images.
%                     - subject_id                 (numeric) Subject ID for filename.
%                     - simulation.medium          (string) Medium type (e.g., 'water').
%                     - io.output_affix            (string) Affix for output filenames.
%                     - seg_labels.eye             (optional int) Label index for eye tissue.
%   medium_masks    - 3D array of initial medium labels (updated in-place and returned).
%   focus_pos_grid  - 1x3 vector [x,y,z] indices of focus position for debug slice (y-slice used).
%   segmented_img   - 3D array of original segmentation labels.
%
% Outputs:
%   medium_masks    - Updated 3D medium mask with expanded skull and filled gaps.
%   skull_i         - Integer label value for skull tissue.
%
% Dependencies:
%   - Image Processing Toolbox: imerode, strel, imfill, bwlabeln, regionprops, label2rgb, imshowpair, montage.
%
% Example:
%   medium_masks = skull_fill_holes(params, medium_masks, 2, bone_img, focus_pos, seg_img, trab_mask, 3);
%
% Notes:
%   - Assumes voxel-based 3D head model from neuroimaging (e.g., CT/MRI segmentation).
%   - Expansion uses morphological perimeter (3-voxel cube structuring element).
%   - Gap filling identifies largest non-skin-skull blob (likely CSF/air) and ignores it.
%   - Eyes default to water (label 0) to avoid erroneous bone assignment.

    labels_medium = fieldnames(parameters.medium_properties);
    labels_requested = fieldnames(parameters.layers);

    if parameters.pct.enabled == 0 && any(contains(labels_requested, 'skull_cortical'))  
        % treat cortical bone as the base layer
        skull_i = find(ismember(labels_medium, {'skull_cortical'; 'skull_trabecular'}));
    else
        skull_i = find(strcmp(labels_medium, 'skull'));
    end

    % Retain pre-fill medium mask for later plotting
    medium_with_gaps = medium_masks;

    % Ensure continuous skull (fill small holes, connect thin regions)
    skull = ismember(medium_masks, skull_i);
    if isfield(parameters.headmodel, 'skull_fill_method') && strcmp(parameters.headmodel.skull_fill_method, 'rubberwrap')
        % Local skull filling with rubber expansion
        skull_continuous = skull_rubber_wrap(parameters, skull, medium_masks, segmented_img);
    else
        % Isotropic dilation
        se = strel('sphere', 3);  % 3D ball for isotropic smoothing
        skull_continuous = imclose(skull, se);  % Dilate then erode
    end
    skull_new = skull_continuous & ~skull; % Identify voxels that were not part of the original skull mask
    medium_masks(skull_new) = skull_i(1);  % Only add new voxels as cortical bone (if differentiated)

    % [DEBUG] Plot skull expansion at focus y-slice
    if parameters.simulation.debug == 1
        h = figure;
        montage({1-squeeze(skull(:,focus_pos_grid(2),:)), ...
            1-squeeze(skull_new(:,focus_pos_grid(2),:)), ...
            1-squeeze(skull_continuous(:,focus_pos_grid(2),:))}, 'Size', [1 3]);
        title('Skull (left), added bone (center), continuous skull (right)');
        output_plot_filename = fullfile(parameters.io.debug_dir_preproc, ...
            sprintf('sub-%03d_%s_skull_expansion%s.png',...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
        saveas(h, output_plot_filename, 'png');
        close(h);
    end

    % Label eye layer as water
    seg_labels = charm_seg_labels();
    if isfield(seg_labels, 'eye')
        eye_i = seg_labels.eye;
        eye = segmented_img == eye_i;
        i_water = find(strcmp(fieldnames(parameters.medium_properties), 'water'));
        medium_masks(eye ~= 0) = i_water; % Default to water
    end

    % [DEBUG] Plot before/after gap filling at focus y-slice
    if parameters.simulation.debug == 1
        medium_mask_updates = double(medium_masks-medium_with_gaps);
        medium_mask_updates(medium_mask_updates ~=0) = 1; % indicate where changes occured
        h = figure;
        imshowpair(label2rgb(squeeze(medium_mask_updates(:,focus_pos_grid(2),:)), 'parula'), ...
            label2rgb(squeeze(medium_masks(:,focus_pos_grid(2),:)), 'parula'), 'montage');
        title('Updated (left) and closed off (right) segmented images');
        output_plot_filename = fullfile(parameters.io.debug_dir_preproc, ...
            sprintf('sub-%03d_%s_segmented_img_closing_gaps_changes%s.png',...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
        saveas(h, output_plot_filename, 'png');
        close(h);
    end
end
