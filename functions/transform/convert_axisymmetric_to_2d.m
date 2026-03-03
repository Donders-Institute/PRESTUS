function [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
    convert_axisymmetric_to_2d(sensor_data, parameters, segmentation, medium_masks, kwave_medium, source, source_labels)
          
% Convert the data from symmetric about x=0 to a mirrored 2d setup.
sensor_data.p_final = cat(2, fliplr(sensor_data.p_final), sensor_data.p_final);
sensor_data.p_max_all = cat(2, fliplr(sensor_data.p_max_all), sensor_data.p_max_all);
segmentation = cat(2, fliplr(segmentation), segmentation);
medium_masks = cat(2, fliplr(medium_masks), medium_masks);
source.p_mask = cat(2, fliplr(source.p_mask), source.p_mask);
source_labels = cat(2, fliplr(source_labels), source_labels);
fields = fieldnames(kwave_medium);
for i = 1:numel(fields)
    if size(kwave_medium.(fields{i}),2)>1
        kwave_medium.(fields{i}) = cat(2, fliplr(kwave_medium.(fields{i})), kwave_medium.(fields{i}));
    end
end

% shift positions to diameter location
trans_pos_final = parameters.transducer(1).trans_pos;
focus_pos_final = parameters.transducer(1).focus_pos;

trans_pos_final(2) = trans_pos_final(2)+parameters.grid_dims(2);
focus_pos_final(2) = focus_pos_final(2)+parameters.grid_dims(2);

% convert radial dimension size into diameter
parameters.grid_dims(2) = parameters.grid_dims(2)*2;

% shift radial to x dim
parameters.grid_dims = fliplr(parameters.grid_dims);
trans_pos_final = fliplr(trans_pos_final);
focus_pos_final = fliplr(focus_pos_final);
segmentation = segmentation';
medium_masks = medium_masks';
sensor_data.p_max_all = sensor_data.p_max_all';
sensor_data.p_final = sensor_data.p_final';
source.p_mask = source.p_mask';
source_labels = source_labels';
fields = fieldnames(kwave_medium);
for i = 1:numel(fields)
    if size(kwave_medium.(fields{i}),2)>1
        kwave_medium.(fields{i}) = kwave_medium.(fields{i})';
    end
end

% Retain transducer and focus positions after all grid manipulations
parameters.transducer(1).trans_pos = trans_pos_final;
parameters.transducer(1).focus_pos = focus_pos_final;

% set up kgrid again for eventual heating sim
kgrid = kWaveGrid(parameters.grid_dims(1), parameters.grid_step_mm/1e3, ...
          parameters.grid_dims(2), parameters.grid_step_mm/1e3);

end