function [medium_masks, skull_i] = skull_fill_holes(parameters, medium_masks, focus_pos_grid, segmented_img, qc_grid)
% SKULL_FILL_HOLES  Expand skull mask and fill gaps to maintain continuity
%
% Expands the skull region using morphological perimeter detection and
% fills gaps between bone and skin, while preserving trabecular bone and
% setting eye voxels to water. headmodel.skull_fill_method: 'off' keeps the
% skull mask as segmented (eye voxels still become water). Debug slice figures are saved when enabled.
%
% Use as:
%   [medium_masks, skull_i] = skull_fill_holes(parameters, medium_masks, focus_pos_grid, segmented_img)
%   [medium_masks, skull_i] = skull_fill_holes(..., qc_grid)
%
% Input:
%   parameters     - (1,1) simulation configuration struct
%   medium_masks   - [Nx x Ny x Nz] medium label array
%   focus_pos_grid - [1x3] focus position in voxel coordinates
%   segmented_img  - [Nx x Ny x Nz] SimNIBS tissue label volume
%   qc_grid        - grid description for the head-model QC records (SEG_QC);
%                    default: the placement grid
%
% Output:
%   medium_masks - updated medium label array with skull gaps filled
%   skull_i      - integer label index for skull tissue
%
% See also: SKULL_RUBBER_WRAP, PREPROC_MEDIUM_MASK

    if nargin < 5 || isempty(qc_grid)
        qc_grid = seg_qc('grid_note', parameters, 'placement');
    end

    labels_medium = fieldnames(parameters.medium_properties);
    labels_requested = fieldnames(parameters.layers);

    % Treat cortical bone as the base layer when it was requested and is
    % present in the mask. PREPROC_MEDIUM_MASK drops requested layers that the
    % segmentation lacks (e.g. bone labelled 4 instead of 7/8) without updating
    % parameters.layers, so the requested list alone would select labels that
    % do not occur and leave the skull unfilled.
    cortical_i = find(ismember(labels_medium, {'skull_cortical'; 'skull_trabecular'}));
    if parameters.pct.enabled == 0 && any(contains(labels_requested, 'skull_cortical')) && ...
            any(ismember(medium_masks(:), cortical_i))
        skull_i = cortical_i;
    else
        skull_i = find(strcmp(labels_medium, 'skull'));
    end

    % Retain pre-fill medium mask for later plotting
    medium_with_gaps = medium_masks;

    % Ensure continuous skull (fill small holes, connect thin regions)
    skull = ismember(medium_masks, skull_i);
    fill_method = 'rubberwrap';
    if isfield(parameters.headmodel, 'skull_fill_method') && ~isempty(parameters.headmodel.skull_fill_method)
        fill_method = char(parameters.headmodel.skull_fill_method);
    end
    if strcmp(fill_method, 'off')
        % No filling: the skull mask is used as segmented
        skull_continuous = skull;
    elseif strcmp(fill_method, 'rubberwrap')
        % Local skull filling with rubber expansion
        skull_continuous = skull_rubber_wrap(parameters, skull, medium_masks, segmented_img);
    else
        % Isotropic dilation
        [kx,ky,kz] = ndgrid(-3:3, -3:3, -3:3);
        se = (kx.^2 + ky.^2 + kz.^2) <= 9;     % sphere r=3
        skull_continuous = convn(logical(skull), se, 'same') > 0;            % dilate
        skull_continuous = convn(skull_continuous, se, 'same') == nnz(se);   % erode
    end
    skull_new = skull_continuous & ~skull; % Identify voxels that were not part of the original skull mask
    % Add new voxels as cortical bone when the skull is differentiated
    % (skull_i lists trabecular before cortical, in medium_properties order)
    fill_i = skull_i(1);
    if numel(skull_i) > 1
        fill_i = find(strcmp(labels_medium, 'skull_cortical'));
    end
    medium_masks(skull_new) = fill_i;

    if seg_qc('active')
        if strcmp(fill_method, 'off')
            fill_desc = 'Skull holes not filled (headmodel.skull_fill_method: off).';
        elseif strcmp(fill_method, 'rubberwrap')
            r_mm = 5;                                   % skull_rubber_wrap's default
            if isfield(parameters.headmodel, 'skull_wrap_radius_mm') && ~isempty(parameters.headmodel.skull_wrap_radius_mm)
                r_mm = parameters.headmodel.skull_wrap_radius_mm;
            end
            fill_desc = sprintf(['Skull holes filled by rubber wrapping (headmodel.skull_fill_method: rubberwrap, ' ...
                'skull_wrap_radius_mm %g mm): voxels between brain and skin inside the wrapped skull become bone.'], ...
                r_mm);
        else
            fill_desc = ['Skull holes filled by morphological closing (headmodel.skull_fill_method: imclose, ' ...
                'sphere of radius 3 voxels).'];
        end
        if numel(skull_i) > 1 && ~strcmp(fill_method, 'off')
            fill_desc = [fill_desc ' Added voxels are cortical bone.'];
        end
        seg_qc('add', 'skull_fill', medium_with_gaps, medium_masks, focus_pos_grid, struct( ...
            'title', 'Skull hole filling', 'description', fill_desc, ...
            'grid', qc_grid, 'label_space', 'medium'));
    end

    % [DEBUG] Plot skull expansion at focus y-slice
    if debug_on(parameters, 'preproc', 'figures')
        h = figure;
        montage({1-squeeze(skull(:,focus_pos_grid(2),:)), ...
            1-squeeze(skull_new(:,focus_pos_grid(2),:)), ...
            1-squeeze(skull_continuous(:,focus_pos_grid(2),:))}, 'Size', [1 3]);
        title('Skull (left), added bone (center), continuous skull (right)');
        debug_save_figure(parameters, 'preproc', h, ...
            sprintf('sub-%03d_%s_skull_expansion%s',...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix), ...
            'Skull (left), added bone (centre) and continuous skull (right), focus y-slice');
    end

    % Label eye layer as water
    seg_labels = charm_seg_labels();
    if isfield(seg_labels, 'eye')
        before_eyes = [];
        if seg_qc('active'); before_eyes = medium_masks; end
        eye_i = seg_labels.eye;
        eye = segmented_img == eye_i;
        i_water = find(strcmp(fieldnames(parameters.medium_properties), 'water'));
        medium_masks(eye ~= 0) = i_water; % Default to water
        if seg_qc('active')
            seg_qc('add', 'eyes', before_eyes, medium_masks, focus_pos_grid, struct( ...
                'title', 'Eyes to water', ...
                'description', sprintf('Voxels labelled eye by charm (label %d) set to water.', eye_i), ...
                'grid', qc_grid, 'label_space', 'medium'));
        end
        clear before_eyes
    end

    % [DEBUG] Plot before/after gap filling at focus y-slice
    if debug_on(parameters, 'preproc', 'figures')
        medium_mask_updates = double(medium_masks-medium_with_gaps);
        medium_mask_updates(medium_mask_updates ~=0) = 1; % indicate where changes occured
        h = figure;
        imshowpair(label2rgb(squeeze(medium_mask_updates(:,focus_pos_grid(2),:)), 'parula'), ...
            label2rgb(squeeze(medium_masks(:,focus_pos_grid(2),:)), 'parula'), 'montage');
        title('Updated (left) and closed off (right) segmented images');
        debug_save_figure(parameters, 'preproc', h, ...
            sprintf('sub-%03d_%s_segmented_img_closing_gaps_changes%s',...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix), ...
            'Changed voxels (left) and medium mask after hole filling and eyes to water (right), focus y-slice');
    end
end
