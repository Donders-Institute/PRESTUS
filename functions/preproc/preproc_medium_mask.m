function medium_masks = preproc_medium_mask(segmented_img, parameters)
%MEDIUM_MASK_CREATE Convert SimNibs segmentation to smoothed label indices.
%
% INPUT:
%   segmented_img  - [Nx Ny Nz] SimNibs tissue labels (1=wm,2=gm,3=csf,...)
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

    labels = fieldnames(parameters.layers);
    medium_masks = zeros(size(segmented_img));
    
    % Skip water and handle multi-layer skull specially if present
    has_multiskull = any(strcmp(labels, 'skull_cortical'));
    loop_labels = labels(~ismember(labels, {'water'}));
    if has_multiskull
        loop_labels = loop_labels(~ismember(loop_labels, {'skull', 'skull_trabecular'}));
    end
    
    % Main loop: process non-special layers
    for label_i = 1:length(loop_labels)
        label_name = loop_labels{label_i};
        sim_nibs_layers = parameters.layers.(label_name);
        layer_mask = ismember(segmented_img, sim_nibs_layers);
        
        % Tissue-specific smoothing
        if contains(label_name, 'skull')
            threshold = parameters.smooth_threshold_skull;
        else
            threshold = parameters.smooth_threshold_other;
        end
        
        layer_mask_smoothed = smooth_img(layer_mask, parameters.smooth_window, ...
                                        threshold, parameters.smooth_method);
        medium_masks(layer_mask_smoothed ~= 0) = find(strcmp(labels, label_name));
    end
    
    % Special multi-layer skull handling (overrides base skull if present)
    if has_multiskull
        skull_i = find(strcmp(labels, 'skull_cortical'));
        
        % Combined cortical+trabecular skull base
        layer_mask = ismember(segmented_img, getidx(parameters.layers, 'skull'));
        layer_mask_smoothed = smooth_img(layer_mask, parameters.smooth_window, ...
                                        parameters.smooth_threshold_skull, parameters.smooth_method);
        medium_masks(layer_mask_smoothed ~= 0) = skull_i;
        
        % Overlay smoothed trabecular on top
        trabecular_i = find(strcmp(labels, 'skull_trabecular'));
        trabecular_mask = ismember(segmented_img, getidx(parameters.layers, 'skull_trabecular'));
        trabecular_mask_smoothed = smooth_img(trabecular_mask, parameters.smooth_window, ...
                                             parameters.smooth_threshold_skull, parameters.smooth_method);
        medium_masks(trabecular_mask_smoothed ~= 0) = trabecular_i;
    end

    medium_masks = single(medium_masks);
end