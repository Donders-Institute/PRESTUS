function medium_masks = preproc_medium_mask(segmented_img, parameters)
% PREPROC_MEDIUM_MASK  Convert SimNIBS segmentation to smoothed medium label indices
%
% Maps SimNIBS tissue labels to medium property indices defined in
% parameters.medium_properties. Applies layer-specific smoothing and
% handles cortical/trabecular skull differentiation when both are
% requested.
%
% Use as:
%   medium_masks = preproc_medium_mask(segmented_img, parameters)
%
% Input:
%   segmented_img - [Nx x Ny x Nz] or [Nx x Nz] tissue label volume
%   parameters    - (1,1) simulation configuration struct with .layers and .medium_properties
%
% Output:
%   medium_masks - [Nx x Ny x Nz] integer medium label array; voxels no
%                  smoothed layer claims are water
%
% See also: HEAD_SMOOTH_AND_CROP, PREPROC_HEAD

arguments
    segmented_img {mustBeNumericOrLogical}
    parameters    (1,1) struct
end

    % remove unavailable layers, unify skull for pCT
    [parameters] = check_layers(parameters, segmented_img);

    layer_labels = fieldnames(parameters.layers);
    medium_labels = fieldnames(parameters.medium_properties);
    medium_masks = zeros(size(segmented_img));
    
    % === SKULL CONFIGURATION DETECTION ===
    has_trabecular = ismember('skull_trabecular', layer_labels);
    has_cortical = ismember('skull_cortical', layer_labels);
    has_multiskull = has_trabecular && has_cortical;
    
    % Remove special skull layers from main loop
    loop_layer_labels = layer_labels;
    if has_multiskull
        loop_layer_labels = loop_layer_labels(~ismember(loop_layer_labels, ...
            {'skull', 'skull_cortical', 'skull_trabecular'}));
    end
    
    % === MAIN LOOP: Process non-special layers ===
    fwhm_mm  = parameters.headmodel.smooth_fwhm_mm;
    grid_mm  = parameters.grid.resolution_mm;
    sigma_vox = fwhm_mm ./ (grid_mm * 2.35482);
    fprintf('Smoothing medium masks: FWHM=[%s] mm → σ=[%s] vox\n', num2str(fwhm_mm), num2str(sigma_vox));
    for label_i = 1:length(loop_layer_labels)
        label_name = loop_layer_labels{label_i};
        sim_nibs_layers = parameters.layers.(label_name);
        layer_mask = ismember(segmented_img, sim_nibs_layers);
        
        % Tissue-specific smoothing
        if contains(label_name, 'skull')
            threshold = parameters.headmodel.smooth_threshold_skull;
        else
            threshold = parameters.headmodel.smooth_threshold_other;
        end
        
        layer_mask_smoothed = smooth_img(layer_mask, parameters.headmodel.smooth_fwhm_mm, parameters.grid.resolution_mm, ...
                                        threshold, parameters.headmodel.smooth_method);
        % assign tissue-specific medium ID
        medium_masks(layer_mask_smoothed ~= 0) = find(strcmp(medium_labels, label_name));
    end
    
    % === SPECIAL MULTI-LAYER SKULL HANDLING ===
    if has_multiskull  % Only when trabecular AND cortical present

        % Combined skull extent (skull + cortical + trabecular raw labels),
        % smoothed ONCE to get a single consistent outer/inner bone-shell
        % boundary. Smoothing the cortical/base mask and the trabecular mask
        % INDEPENDENTLY (each with its own Gaussian kernel + threshold) lets
        % their boundaries disagree at the smoothing edge: trabecular is a
        % physically-interior (diploe) layer and must never end up directly
        % adjacent to skin, but independent smoothing can let its boundary
        % bulge past where the separately-smoothed cortical boundary landed.
        % Constraining trabecular to be a subset of one combined shell fixes
        % this regardless of the smoothing kernel/threshold used.
        cortical_i   = find(strcmp(medium_labels, 'skull_cortical'));
        trabecular_i = find(strcmp(medium_labels, 'skull_trabecular'));

        combined_layers   = getidx(parameters.layers, {'skull', 'skull_cortical', 'skull_trabecular'});
        combined_mask     = ismember(segmented_img, combined_layers);
        combined_smoothed = smooth_img(combined_mask, parameters.headmodel.smooth_fwhm_mm, parameters.grid.resolution_mm, ...
                                        parameters.headmodel.smooth_threshold_skull, parameters.headmodel.smooth_method);

        % Step 1: the whole combined shell defaults to cortical.
        medium_masks(combined_smoothed ~= 0) = cortical_i;

        % Step 2: trabecular overlay, constrained to the combined shell so it
        % can never extend past the cortical envelope (e.g. touch skin).
        trabecular_mask = ismember(segmented_img, getidx(parameters.layers, 'skull_trabecular'));
        trabecular_mask_smoothed = smooth_img(trabecular_mask, parameters.headmodel.smooth_fwhm_mm, parameters.grid.resolution_mm, ...
                                             parameters.headmodel.smooth_threshold_skull, parameters.headmodel.smooth_method);
        trabecular_final = (trabecular_mask_smoothed ~= 0) & (combined_smoothed ~= 0);
        medium_masks(trabecular_final) = trabecular_i;
    end

    % Voxels no smoothed layer claims (e.g. where layers meet) are water,
    % as in MEDIUM_SETUP
    i_water = find(strcmp(medium_labels, 'water'), 1);
    if ~isempty(i_water)
        medium_masks(medium_masks == 0) = i_water;
    end

    medium_masks = single(medium_masks);

end