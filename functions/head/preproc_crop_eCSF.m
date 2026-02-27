function [medium_masks] = preproc_crop_eCSF(parameters, medium_masks, segmented_img, trans_pos_grid)
% Expand CSF to mask areas away from the brain (e.g., exclude distant bone)

    if isfield(parameters, 'seg_labels') && any(strcmp(fieldnames(parameters.seg_labels), 'csf'))
        % get CSF mask
        csf_mask = segmented_img==getidx(parameters.seg_labels,'csf');
        expansion_voxels = ceil(parameters.csf_mask_expansion_factor / parameters.grid_step_mm);
        SE = strel('cube', expansion_voxels);
        csf_mask_expanded = imdilate(csf_mask, SE);
        % Mask out non-CSF-expanded regions (i.e., set to water layer)
        medium_masks(~csf_mask_expanded) = find(strcmp((fieldnames(parameters.medium), 'water'));
        % [DEBUG] Visualize CSF expansion at transducer y-slice
        if parameters.debug == 1
            h = figure;
            imshowpair(squeeze(medium_masks(:, trans_pos_grid(2), :)), ...
                squeeze(csf_mask_expanded(:, trans_pos_grid(2), :)), 'falsecolor');
            title('CSF mask (pink) and the segmented image (white)');
            output_plot_filename = fullfile(parameters.debug_dir, ...
                sprintf('sub-%03d_%s_skull_csf_mask%s.png', ...
                parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png');
            close(h);
        end
    else
        warning("CSF layer not specified or unknown ... will not use expanded CSF mask...");
    end