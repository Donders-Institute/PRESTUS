function mask = tissuemask_binary(parameters, medium_masks)
% TISSUEMASK_BINARY  Build per-tissue binary masks from a medium label map
%
% Derives skull, skull_cortical, skull_trabecular, brain, skin, and water
% binary masks by finding the label indices in parameters.medium_properties
% and applying ismember to medium_masks.
%
% Use as:
%   mask = tissuemask_binary(parameters, medium_masks)
%
% Input:
%   parameters   - PRESTUS config with medium_properties fieldnames used as tissue labels
%   medium_masks - medium label volume (values = medium index)
%
% Output:
%   mask - struct with binary fields: skull, skull_cortical, skull_trabecular,
%          brain, skin, water
%
% See also: ACOUSTIC_ANALYSIS, THERMAL_ANALYSIS, GETIDX

arguments
    parameters   (1,1) struct
    medium_masks {mustBeNumeric}
end

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