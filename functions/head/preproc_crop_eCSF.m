function [medium_masks] = preproc_crop_eCSF(parameters, medium_masks, segmented_img, trans_pos_grid, focus_pos_grid)
% PREPROC_CROP_ECSF  Expand CSF mask to exclude distant bone from the simulation grid
%
% Dilates the CSF tissue label and sets all voxels outside the expanded
% CSF region to water, effectively removing distant bone from the medium
% masks before grid cropping.
%
% Use as:
%   medium_masks = preproc_crop_eCSF(parameters, medium_masks, segmented_img, trans_pos_grid)
%   medium_masks = preproc_crop_eCSF(..., focus_pos_grid)
%
% Input:
%   parameters     - (1,1) simulation configuration struct
%   medium_masks   - [Nx x Ny x Nz] medium label array
%   segmented_img  - [Nx x Ny x Nz] SimNIBS tissue label volume
%   trans_pos_grid - [1x3] transducer position in voxel coordinates
%   focus_pos_grid - [1x3] focus position, where the head-model QC slices
%                    (SEG_QC) are taken (default: trans_pos_grid)
%
% Output:
%   medium_masks - updated medium label array with distant bone set to water
%
% See also: PREPROC_CROP_GRID, HEAD_SMOOTH_AND_CROP

    seg_labels = charm_seg_labels();
    if isfield(seg_labels, 'csf')
        % get CSF mask
        csf_mask = segmented_img==getidx(seg_labels,'csf');
        % Dilate by csf_expansion along each axis (cube of side 2r+1 voxels):
        % one moving maximum per axis, whose cost does not grow with r
        expansion_voxels = round(parameters.headmodel.csf_expansion / parameters.grid.resolution_mm);
        csf_mask_expanded = dilate_cube(csf_mask, expansion_voxels);
        before = [];
        if seg_qc('active'); before = medium_masks; end
        % Mask out non-CSF-expanded regions (i.e., set to water layer)
        medium_masks(~csf_mask_expanded) = find(strcmp(fieldnames(parameters.medium_properties), 'water'));
        if seg_qc('active')
            if nargin < 5 || isempty(focus_pos_grid); focus_pos_grid = trans_pos_grid; end
            seg_qc('add', 'csf_crop', before, medium_masks, focus_pos_grid, struct( ...
                'title', 'CSF-based crop', ...
                'description', sprintf(['Voxels farther than headmodel.csf_expansion (%g mm, %d voxels) ' ...
                'from CSF set to water, which removes distant bone and soft tissue.'], ...
                parameters.headmodel.csf_expansion, expansion_voxels), ...
                'grid', seg_qc('grid_note', parameters, 'placement'), 'label_space', 'medium'));
        end
        clear before
        % [DEBUG] Visualize CSF expansion at transducer y-slice
        if debug_on(parameters, 'preproc', 'figures')
            h = figure;
            imshowpair(squeeze(medium_masks(:, trans_pos_grid(2), :)), ...
                squeeze(csf_mask_expanded(:, trans_pos_grid(2), :)), 'falsecolor');
            title('CSF mask (pink) and the segmented image (white)');
            debug_save_figure(parameters, 'preproc', h, ...
                sprintf('sub-%03d_%s_skull_csf_mask%s', ...
                parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix), ...
                'Expanded CSF mask (pink) and medium mask (white), transducer y-slice');
        end
    else
        warn("CSF layer not specified or unknown ... will not use expanded CSF mask...");
    end
end

function M = dilate_cube(M, r)
% Dilation of a 3-D mask by a cube of side 2r+1 voxels (voxels outside the array are background)
    M = uint8(M);
    for d = 1:3
        M = movmax(M, [r r], d);
    end
    M = M > 0;
end
