
%layers_skull = contains(fieldnames(parameters.layer_labels), 'skull')

% if (parameters.segmentation_software == 'charm')
%             skull_mask =  segmented_image_cropped==7 | segmented_image_cropped==8;
%         else
%             skull_mask =  segmented_image_cropped==4;
% end

%prova = getfield(parameters,'layer_labels',valu)


labels = fieldnames(parameters.layer_labels);

%  for label_i = 1:length(parameters.layer_labels)
%      if contains(fieldnames(parameters.layer_labels.(labels{label_i})), 'skull')
%          skull_layers{label_i} = parameters.layer_labels.(labels{label_i});
%      end
%  end
% 
%   sim_nibs_layers = parameters.layer_labels.(labels{label_i});
%         layer_mask = ismember(segmented_img, sim_nibs_layers);


  all_skull_ids = [];
        for label_i = find(contains(labels,  'skull'))'
            all_skull_ids = [all_skull_ids parameters.layer_labels.(labels{label_i})];
        end

  skull_mask = ismember(segmented_image_cropped,all_skull_ids);     