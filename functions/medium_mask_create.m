function [medium_masks] = medium_mask_create(segmented_img, parameters, windowSize)
    arguments
        segmented_img (:,:,:) % tissue segmentation from SimNibs [alt: pseudoCT tissue mask]
        parameters struct
        windowSize % smoothing window size
    end
    
    % This function turns the original segmentations into medium masks such
    % that the setup_medium.m function can fill in the tissue-dependent parameters.
    % Note that tissue masks will assume the labels specified below.

    labels = fieldnames(parameters.layer_labels);

    % create an empty grid the size of the segmented image
    medium_masks = zeros(size(segmented_img));
    % add a smoothing threshold to bone and other non-water tissue.
    for label_i = 1:length(labels)
        if strcmp(labels{label_i}, 'water')
            continue
        end
        sim_nibs_layers = parameters.layer_labels.(labels{label_i});
        layer_mask = ismember(segmented_img, sim_nibs_layers);
        if contains(labels{label_i}, 'skull')
            smooth_threshold = parameters.skull_smooth_threshold;    
%             if strcmp(labels{label_i}, 'skull_trabecular') % two bone types are smoothed together later
%                 continue
%             end
        else
            smooth_threshold = parameters.other_smooth_threshold;
        end
        layer_mask_smoothed = smooth_img(layer_mask, windowSize, smooth_threshold);
        medium_masks(layer_mask_smoothed~=0) = label_i;
    end

end