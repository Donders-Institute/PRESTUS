function [sensor_data, parameters, trans_pos_final, focus_pos_final, ...
    segmented_image_cropped, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
    convert_axisymmetric_to_3d(...
    sensor_data, parameters, trans_pos_final, focus_pos_final, ...
    segmented_image_cropped, medium_masks, kwave_medium, source, source_labels)

%% convert 2d axisymmetry images to 3d

sensor_data.p_final = radialExpand2DTo3D(sensor_data.p_final);
sensor_data.p_max_all = radialExpand2DTo3D(sensor_data.p_max_all);
segmented_image_cropped = radialExpand2DTo3D(segmented_image_cropped);
medium_masks = radialExpand2DTo3D(medium_masks);
source.p_mask = radialExpand2DTo3D(double(source.p_mask));
source_labels = radialExpand2DTo3D(double(source_labels));
fields = fieldnames(kwave_medium);
for i = 1:numel(fields)
    if size(kwave_medium.(fields{i}),2)>1
        kwave_medium.(fields{i}) = radialExpand2DTo3D(kwave_medium.(fields{i}));
    end
end

%% convert positions from 2d radial to 3d 

parameters.n_sim_dims = 3;

% shift positions to diameter location
trans_pos_final(2) = trans_pos_final(2)+parameters.grid_dims(2);
focus_pos_final(2) = focus_pos_final(2)+parameters.grid_dims(2);

% convert: radial x axial -> diameter x diameter x axial
parameters.grid_dims(2) = parameters.grid_dims(2)*2;
parameters.grid_dims = parameters.grid_dims';
parameters.grid_dims = [...
    parameters.grid_dims(2), ...
    parameters.grid_dims(2), ...
    parameters.grid_dims(1)];
parameters.default_grid_dims(2) = parameters.default_grid_dims(2)*2;
parameters.default_grid_dims = parameters.default_grid_dims';
parameters.default_grid_dims = [...
    parameters.default_grid_dims(2), ...
    parameters.default_grid_dims(2), ...
    parameters.default_grid_dims(1)];

% convert: transducer and focus positions
trans_pos_final = fliplr(trans_pos_final);
trans_pos_final = [trans_pos_final(1), trans_pos_final(1), trans_pos_final(2)];
focus_pos_final = fliplr(focus_pos_final);
focus_pos_final = [focus_pos_final(1), focus_pos_final(1), focus_pos_final(2)];

parameters.transducer.pos_grid = trans_pos_final;
parameters.focus_pos_grid = focus_pos_final;

%% set up 3d kgrid for follow-up heating simulations

kgrid = kWaveGrid(parameters.grid_dims(1), parameters.grid_step_m, ...
          parameters.grid_dims(2), parameters.grid_step_m, ...
          parameters.grid_dims(3), parameters.grid_step_m);

%% DEBUG: visually inspect

% % plot center of transducer
% figure; imagesc(squeeze(sensor_data.p_max_all(:,floor(parameters.grid_dims(1)/2),:)))
% 
% % plot each slice along the axial dimension
% figure;
% for index = 1:parameters.grid_dims(3)
%     imagesc(squeeze(sensor_data.p_max_all(:,:,index)))
%     pause(0.5)
% end

end