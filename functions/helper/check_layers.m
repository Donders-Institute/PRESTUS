function [parameters] = check_layers(parameters, segmentation)

% Function to check and remove unavailable segmentation layers.

    % Get all non-negative segmentation values
    seg_values = unique(segmentation(segmentation >= 0));

    % [WATER] force add layer incl. all segmentation (later remove ids specified in other layers)
    parameters.layers.water = seg_values;

    % [pCT] combine cortical and trabecular layers (if specified)
    if parameters.pct.enabled
        % Collect ALL skull layer indices
        skull_layers = [];
        if isfield(parameters.layers, 'skull')
            skull_layers = [skull_layers, parameters.layers.skull];
        end
        if isfield(parameters.layers, 'skull_cortical')
            skull_layers = [skull_layers, parameters.layers.skull_cortical];
        end
        if isfield(parameters.layers, 'skull_trabecular')
            skull_layers = [skull_layers, parameters.layers.skull_trabecular];
        end
        % Assign merged to skull
        parameters.layers.skull = unique(skull_layers);
        % Remove differentiated skull
        if isfield(parameters.layers, 'skull_cortical')
            parameters.layers = rmfield(parameters.layers, 'skull_cortical');
            warning('pCT requested with skull_cortical layer, unifying into a single skull layer ...');
        end
        if isfield(parameters.layers, 'skull_trabecular')
            parameters.layers = rmfield(parameters.layers, 'skull_trabecular');
            warning('pCT requested with skull_trabecular layer, unifying into a single skull layer ...');
        end
    end
    
    % Check all requested and available layers
    
    layers_requested = fieldnames(parameters.layers);
    layers_available = zeros(numel(layers_requested),1);
    
    for i_layer = 1:numel(layers_requested)
        if any(ismember(unique(segmentation), parameters.layers.(layers_requested{i_layer})))
            layers_available(i_layer) = 1;
        else
            warning(['[Layer] ', layers_requested{i_layer}, ' requested but not available in segmentation. Removing it...']);
        end
    end
    
    % fields to remove due to unavailable segmentations
    fields_to_remove = layers_requested(layers_available == 0);
    % remove them from the struct
    parameters.layers = rmfield(parameters.layers, fields_to_remove);

    % [WATER] remove segmentation ids attributed to non-water tissues
    layers_nonwater = layers_requested(layers_available & ~strcmp(layers_requested, 'water'));
    parameters.layers.water = parameters.layers.water(~ismember(parameters.layers.water, getidx(parameters.layers, layers_nonwater)));

end