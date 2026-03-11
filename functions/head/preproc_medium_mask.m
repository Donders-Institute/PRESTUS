function medium_masks = preproc_medium_mask(segmented_img, parameters)
%MEDIUM_MASK_CREATE Convert SimNibs segmentation to smoothed label indices.
%
% INPUT:
%   segmented_img  - [Nx Ny Nz] SimNibs tissue layer_labels (1=wm,2=gm,3=csf,...)
%   parameters     - Struct with .layers, smoothing params
%
% OUTPUT:
%   medium_masks   - [Nx Ny Nz] Integer label indices (per layers order)
%
% LOGIC:
% Create a list of all requested layers. 
% Remove the water layer from the loop. 
% If skull_cortical is requested, remove skull and skull_trabecular from the loop. 
% In the skull_cortical loop proceed sequentially with whole skull smoothing and then insert trabecular layer.

    % remove unavailable layers, unify skull for pCT
    [parameters] = check_layers(parameters, segmented_img);

    layer_labels = fieldnames(parameters.layers);
    medium_labels = fieldnames(parameters.medium);
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
    for label_i = 1:length(loop_layer_labels)
        label_name = loop_layer_labels{label_i};
        sim_nibs_layers = parameters.layers.(label_name);
        layer_mask = ismember(segmented_img, sim_nibs_layers);
        
        % Tissue-specific smoothing
        if contains(label_name, 'skull')
            threshold = parameters.smooth_threshold_skull;
        else
            threshold = parameters.smooth_threshold_other;
        end
        
        layer_mask_smoothed = smooth_img(layer_mask, parameters.smooth_fwhm_mm, parameters.grid_step_mm, ...
                                        threshold, parameters.smooth_method);
        % assign tissue-specific medium ID
        medium_masks(layer_mask_smoothed ~= 0) = find(strcmp(medium_labels, label_name));
    end
    
    % === SPECIAL MULTI-LAYER SKULL HANDLING ===
    if has_multiskull  % Only when trabecular AND cortical present
        
        % Step 1: Assign base skull (skull + cortical) → cortical medium ID
        cortical_i = find(strcmp(medium_labels, 'skull_cortical'));
        skull_base_layers = getidx(parameters.layers, {'skull', 'skull_cortical'});
        layer_mask = ismember(segmented_img, skull_base_layers);
        layer_mask_smoothed = smooth_img(layer_mask, parameters.smooth_fwhm_mm, parameters.grid_step_mm, ...
                                        parameters.smooth_threshold_skull, parameters.smooth_method);
        medium_masks(layer_mask_smoothed ~= 0) = cortical_i;
        
        % Step 2: Overlay trabecular on top
        trabecular_i = find(strcmp(medium_labels, 'skull_trabecular'));
        trabecular_mask = ismember(segmented_img, getidx(parameters.layers, 'skull_trabecular'));
        trabecular_mask_smoothed = smooth_img(trabecular_mask, parameters.smooth_fwhm_mm, parameters.grid_step_mm, ...
                                             parameters.smooth_threshold_skull, parameters.smooth_method);
        medium_masks(trabecular_mask_smoothed ~= 0) = trabecular_i;
    end

    medium_masks = single(medium_masks);

end