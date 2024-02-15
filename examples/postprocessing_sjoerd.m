
% start from tusim/code folder
cd /home/visual/andche/STAFF_SCI/andche_sandbox/TUS_sims/tusim/code/

parameters = load_parameters('sjoerd_config.yaml');
target_names = {'left_amygdala', 'right_amygdala'};
var_names = ["subject_id", "target", "medium", "layer", "max_pressure", "max_isppa",...
    "max_pos_x", "max_pos_y", "max_pos_z",   "max_press_pos_x", "max_press_pos_y", "max_press_pos_z", "max_isppa_at_expected_pos" ];
res_table = table('Size',[3*2*2*3 length(var_names)],'VariableTypes', ...
    [repmat("string", [1 4]) repmat("double", [1 2]) repmat("uint8", [1 6]) "double" ],'VariableNames', var_names); 
cur_row = 1;

% numbers
for subject_id = [1,5,6]
    for target_id = [1,2]
        results_filename_affix = sprintf('_target_%s', target_names{target_id});
        
        skull_obj = load(fullfile(parameters.data_path, sprintf('sub-%03d_after_cropping_and_smoothing%s.mat', ...
            subject_id, results_filename_affix)),'skull_mask', 'crop_translation_matrix','segmented_image_cropped', 'trans_pos_final', 'focus_pos_final');
        
        median_i = round(size(skull_obj.skull_mask,2)/2);
        imshowpair(squeeze(skull_obj.skull_mask(:,median_i,:)),squeeze(skull_obj.segmented_image_cropped(:,median_i,:)), 'montage')
        imshow(label2rgb(squeeze(skull_obj.segmented_image_cropped(:,median_i,:))))

        brain_only = skull_obj.segmented_image_cropped > 0 & skull_obj.segmented_image_cropped <4;
        skull_skin = skull_obj.segmented_image_cropped >= 4 & skull_obj.segmented_image_cropped <=5;

        for medium = ["water","water_and_skull"]
            res = load(fullfile(parameters.data_path, ...
                'sim_outputs', sprintf('sub-%03d_%s_results%s.mat', subject_id, medium, results_filename_affix)),...
                'sensor_data', 'parameters', 'kwave_medium');
            max_pressure_map = gather(res.sensor_data.p_max_all);
            Isppa_map = max_pressure_map.^2./(2*(res.kwave_medium.sound_speed.*res.kwave_medium.density)).*1e-4; 
            for layer = ["brain_only","skull_skin","all"]
                if strcmp(layer, 'all')
                    layer_mask = ones(size(Isppa_map));
                else
                    eval(strjoin(['layer_mask = ' layer ';']));
                end
                max_isppa_pos_layer = [];
                [max_Isppa_layer, max_x, max_y, max_z] = masked_max_3d(Isppa_map, layer_mask>0);
                [max_pressure_layer, max_px, max_py, max_pz] = masked_max_3d(max_pressure_map, layer_mask>0);
                max_pressure_at_expected_pos = Isppa_map(skull_obj.focus_pos_final(1), skull_obj.focus_pos_final(2), skull_obj.focus_pos_final(3));
                res_table(cur_row,:) = {subject_id target_names{target_id} medium layer max_pressure_layer*1e-6 ...
                    max_Isppa_layer max_x max_y max_z max_px max_py max_pz max_pressure_at_expected_pos };
                cur_row = cur_row+1;
            end

        end
    end
end
writetable(res_table, fullfile(parameters.data_path,'results_aggregated.csv'))

% images

for subject_id = [1,5,6]
    for target_id = [1,2]
        
        results_filename_affix = sprintf('_target_%s', target_names{target_id});
        
        load(fullfile(parameters.data_path, ...
                'sim_outputs', sprintf('sub-%03d_%s_results%s.mat', subject_id, 'water_and_skull', results_filename_affix)),...
                'parameters', 'sensor_data','kwave_medium');
            
        max_pressure_map = gather(sensor_data.p_max_all);
        Isppa_map = max_pressure_map.^2./(2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4; 

        t1_img = niftiread(fullfile(parameters.data_path, sprintf(parameters.t1_path_template, subject_id)));
        t1_hdr = niftiinfo(fullfile(parameters.data_path, sprintf(parameters.t1_path_template, subject_id)));
        rr_obj = load(fullfile(parameters.data_path, ...
            sprintf('sub-%03d_after_rotating_and_scaling%s.mat', subject_id, results_filename_affix)),'scale_rotate_recenter_matrix');
        
        skull_obj = load(fullfile(parameters.data_path, sprintf('sub-%03d_after_cropping_and_smoothing%s.mat',...
            subject_id, results_filename_affix)), 'crop_translation_matrix', 'trans_pos_final', 'focus_pos_final');

        final_transformation_matrix = rr_obj.scale_rotate_recenter_matrix*skull_obj.crop_translation_matrix';
        inv_final_transformation_matrix = maketform('affine', inv(final_transformation_matrix')');
        
        res_row = res_table(strcmp(res_table.target,target_names{target_id})&strcmp(res_table.subject_id, string(subject_id))&strcmp(res_table.medium, 'water_and_skull')&strcmp(res_table.layer,'brain_only'),:)
        max_isppa_pos = double(table2array(res_row(:,["max_pos_x","max_pos_y","max_pos_z"]))');
        
        orig_coordinates = [par_obj.parameters.transducer.pos_t1_grid;  par_obj.parameters.focus_pos_t1_grid];
        backtransf_coordinates = round(tformfwd([skull_obj.trans_pos_final skull_obj.focus_pos_final max_isppa_pos]', inv_final_transformation_matrix));
        
        isppa_map_backtransf = tformarray(Isppa_map, inv_final_transformation_matrix, ...
            makeresampler('linear', 'fill'), [1 2 3], [1 2 3], size(t1_img), [], 0) ;
        [~, source_labels] = transducer_setup(parameters.transducer, par_obj.parameters.transducer.pos_t1_grid, par_obj.parameters.focus_pos_t1_grid, ...
                                                            size(t1_img), t1_hdr.PixelDimensions(1));
        isppa_hdr = t1_hdr;
        
        isppa_hdr.Datatype = 'single';
        niftiwrite(isppa_map_backtransf, fullfile(parameters.data_path, sprintf('sub-%03d_final_isppa_orig_coord%s', subject_id, results_filename_affix)), isppa_hdr, 'Compressed', true)
        plot_isppa_over_image(isppa_map_backtransf, t1_img, source_labels, ...
            ones(size(t1_img)), {'y', backtransf_coordinates(3,2)}, backtransf_coordinates(1,:)', ...
            backtransf_coordinates(2,:)', backtransf_coordinates(3,:)')

        
   end
end
