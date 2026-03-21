function mask = tissuemask_binary(parameters, medium_masks)

    labels = fieldnames(parameters.medium_properties);
    
    skull_i = find(strcmp(labels, 'skull'));
    cortical_i = find(strcmp(labels, 'skull_cortical'));
    trabecular_i = find(strcmp(labels, 'skull_trabecular'));
    all_skull_ids = [skull_i, cortical_i, trabecular_i];
    mask.skull = ismember(medium_masks,all_skull_ids);
    mask.skull_cortical = ismember(medium_masks,cortical_i);
    mask.skull_trabecular = ismember(medium_masks,trabecular_i);

    brain_i = find(strcmp(labels, 'brain'));
    mask.brain = ismember(medium_masks,brain_i);

    skin_i = find(strcmp(labels, 'skin'));
    mask.skin = ismember(medium_masks,skin_i);