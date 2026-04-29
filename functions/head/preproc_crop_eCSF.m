function [medium_masks] = preproc_crop_eCSF(parameters, medium_masks, segmented_img, trans_pos_grid)
% PREPROC_CROP_ECSF  Expand CSF mask to exclude distant bone from the simulation grid
%
% Dilates the CSF tissue label and sets all voxels outside the expanded
% CSF region to water, effectively removing distant bone from the medium
% masks before grid cropping.
%
% Use as:
%   medium_masks = preproc_crop_eCSF(parameters, medium_masks, segmented_img, trans_pos_grid)
%
% Input:
%   parameters     - (1,1) simulation configuration struct
%   medium_masks   - [Nx x Ny x Nz] medium label array
%   segmented_img  - [Nx x Ny x Nz] SimNIBS tissue label volume
%   trans_pos_grid - [1x3] transducer position in voxel coordinates
%
% Output:
%   medium_masks - updated medium label array with distant bone set to water
%
% See also: PREPROC_CROP_GRID, HEAD_SMOOTH_AND_CROP

    seg_labels = charm_seg_labels();
    if isfield(seg_labels, 'csf')
        % get CSF mask
        csf_mask = segmented_img==getidx(seg_labels,'csf');
        expansion_voxels = ceil(parameters.headmodel.csf_expansion / parameters.grid.resolution_mm);
        SE = strel('cube', expansion_voxels);
        csf_mask_expanded = imdilate(csf_mask, SE);
        % Mask out non-CSF-expanded regions (i.e., set to water layer)
        medium_masks(~csf_mask_expanded) = find(strcmp(fieldnames(parameters.medium_properties), 'water'));
        % [DEBUG] Visualize CSF expansion at transducer y-slice
        if parameters.simulation.debug == 1
            h = figure;
            imshowpair(squeeze(medium_masks(:, trans_pos_grid(2), :)), ...
                squeeze(csf_mask_expanded(:, trans_pos_grid(2), :)), 'falsecolor');
            title('CSF mask (pink) and the segmented image (white)');
            output_plot_filename = fullfile(parameters.io.dir_debug_preproc, ...
                sprintf('sub-%03d_%s_skull_csf_mask%s.png', ...
                parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
            saveas(h, output_plot_filename, 'png');
            close(h);
        end
    else
        warn("CSF layer not specified or unknown ... will not use expanded CSF mask...");
    end