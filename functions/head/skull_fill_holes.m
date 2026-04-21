function [medium_masks, skull_i] = skull_fill_holes(parameters, medium_masks, focus_pos_grid, segmented_img)
% SKULL_FILL_HOLES  Expand skull mask and fill gaps to maintain continuity
%
% Expands the skull region using morphological perimeter detection and
% fills gaps between bone and skin, while preserving trabecular bone and
% setting eye voxels to water. Debug slice figures are saved when enabled.
%
% Use as:
%   [medium_masks, skull_i] = skull_fill_holes(parameters, medium_masks, focus_pos_grid, segmented_img)
%
% Input:
%   parameters     - (1,1) simulation configuration struct
%   medium_masks   - [Nx x Ny x Nz] medium label array
%   focus_pos_grid - [1x3] focus position in voxel coordinates
%   segmented_img  - [Nx x Ny x Nz] SimNIBS tissue label volume
%
% Output:
%   medium_masks - updated medium label array with skull gaps filled
%   skull_i      - integer label index for skull tissue
%
% See also: SKULL_RUBBER_WRAP, PREPROC_MEDIUM_MASK

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
