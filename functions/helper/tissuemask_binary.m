function mask = tissuemask_binary(parameters, medium_masks)

    labels = fieldnames(parameters.layers);
    skull_i = find(strcmp(labels, 'skull'));
    cortical_i = find(strcmp(labels, 'skull_cortical'));
    trabecular_i = find(strcmp(labels, 'skull_trabecular'));
    all_skull_ids = [skull_i, cortical_i, trabecular_i];
    mask.skull = ismember(medium_masks,all_skull_ids);
    brain_i = find(strcmp(labels, 'brain'));
    mask.brain = ismember(medium_masks,brain_i);
    skin_i = find(strcmp(labels, 'skin'));
    mask.skin = ismember(medium_masks,skin_i);