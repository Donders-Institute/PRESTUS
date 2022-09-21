function pipeline_sjoerd(subject_id, target_id, parameters)

    % start from tusim/code folder
    cd /home/visual/andche/STAFF_SCI/andche_sandbox/TUS_sims/tusim/code/
    % add paths
    addpath('.');
    addpath('functions')
    addpath('toolboxes/kwave') % set your kwave path here
    addpath('toolboxes/Colormaps') % set your path to Colormaps files here
    addpath('toolboxes/export_fig') % set your path to export_fig files here
    addpath('toolboxes/yaml') % set your path to yaml files here
    
    amygdala_pos_allsubj = readmatrix(fullfile(parameters.data_path, 'amygdala_coordinates.csv'));

    filename_t1 = fullfile(parameters.data_path, sprintf(parameters.t1_path_template, subject_id));

    t1_header = niftiinfo(filename_t1);
    t1_image = niftiread(filename_t1);
    left_amygdala_ras = (amygdala_pos_allsubj(amygdala_pos_allsubj(:,1)==subject_id,2:4))';
    right_amygdala_ras  = (amygdala_pos_allsubj(amygdala_pos_allsubj(:,1)==subject_id,5:7))';
    
    left_amygdala_pos = ras_to_grid(left_amygdala_ras, t1_header);
    right_amygdala_pos = ras_to_grid(right_amygdala_ras, t1_header);

    subj_folder = fullfile(parameters.data_path, sprintf('simulations_Andrey/sub-%03d/', subject_id));
    trig_mark_files = dir(fullfile(subj_folder, '*.xml'));    

    reference_to_transducer_distance = -(parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm);
    
    left_trans_ras_pos = get_trans_pos_from_trigger_markers(fullfile(subj_folder, trig_mark_files(1).name), 10, ...
        reference_to_transducer_distance );
    left_trans_pos = ras_to_grid(left_trans_ras_pos, t1_header);

    right_trans_ras_pos = get_trans_pos_from_trigger_markers(fullfile(subj_folder, trig_mark_files(2).name), 10, ...
        reference_to_transducer_distance );
    right_trans_pos = ras_to_grid(right_trans_ras_pos, t1_header);

    imshowpair(plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), left_trans_pos, left_amygdala_pos, parameters), plot_t1_with_transducer(t1_image, t1_header.PixelDimensions(1), right_trans_pos, right_amygdala_pos, parameters),'montage');

    export_fig(fullfile(parameters.data_path, sprintf('sub-%03d_trans_positions.png', subject_id)), '-native')
    
    fprintf('Distance to left amygdala %03d to right amygdala %03d', norm(left_amygdala_ras-left_trans_ras_pos), norm(right_amygdala_ras-right_trans_ras_pos))
    transducers = [left_trans_pos right_trans_pos];
    targets = [left_amygdala_pos right_amygdala_pos];
    target_names = {'left_amygdala', 'right_amygdala'};
    
    %% Run the pipeline for each subject and target
        parameters.transducer.pos_t1_grid = transducers(:,target_id)';
        parameters.focus_pos_t1_grid = targets(:,target_id)';
        parameters.results_filename_affix = sprintf('_target_%s', target_names{target_id});
        
        if parameters.run_simulations_with_qsub
            % remember current folder to go back later
            current_dir = cd;
            cd(fullfile(parameters.data_path,'batch_job_logs'))
            parameters.interactive = 0;
            parameters.paths_to_add = path;
            parameters.overwrite_files = 'never';

            qsubfeval(@single_subject_pipeline, subject_id, parameters, 'timreq',  60*60*7,  'memreq',  20*(1024^3),  ...
                'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"', 'rerunable', 'yes');
            disp('Simulations job submitted to the cluster, stopping for now. Re-run the pipeline when the job finishes to do post-processing.')
            cd(current_dir);
        else
            [output_pressure_file_water_and_skull, parameters] = single_subject_pipeline(subject_id, parameters);
            % follow-up with water-only simulations and same locations
            new_parameters = parameters;
            new_parameters.simulation_medium = 'water';
            new_parameters.default_grid_dims = new_parameters.grid_dims;
            if subject_id == 5 && target_id == 2
                new_parameters.overwrite_files = 'always';
            end
            new_parameters.overwrite_simnibs = 0;
            
            output_pressure_file_water = single_subject_pipeline(subject_id, new_parameters);
        end
end

