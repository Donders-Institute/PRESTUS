function [output_pressure_file, parameters] = single_subject_pipeline(subject_id, parameters)
% add subject_id to parameters to pass arguments to functions more easily
parameters.subject_id = subject_id;

% if there are paths to be added, add them; this is mostly for batch runs
if isfield(parameters,'paths_to_add')
   path(path, parameters.paths_to_add)
end

% output directory

output_dir = fullfile(parameters.data_path, 'sim_outputs');
if ~exist(output_dir, 'file' )
    mkdir(output_dir)
end
output_pressure_file = fullfile(output_dir,sprintf('sub-%03d_%s_isppa%s.csv', subject_id, parameters.simulation_medium, parameters.results_filename_affix));

if strcmp(parameters.simulation_medium, 'water_and_skull')
    [skull_mask, segmented_image_cropped, skull_edge, trans_pos_final, focus_pos_final] = preprocess_brain(parameters, subject_id);
    if isempty(skull_mask)
        output_pressure_file = '';
        return;
    end
    parameters.grid_dims = size(skull_mask);
else
    % make sure that the grid dimensions are known
    assert(isfield(parameters, 'default_grid_dims'), 'The parameters structure should have the field grid_dims for the grid dimensions')
    parameters.grid_dims = parameters.default_grid_dims;
    skull_mask = [];
    segmented_image_cropped = zeros(parameters.grid_dims);
    if ~isfield(parameters.transducer, 'pos_grid') || ~isfield(parameters, 'focus_pos_grid')
        disp('Either grid or focus position is not set, positioning them arbitrarily based on the focal distance')
        % note that the focus position matters only for the orientation of
        % the transducer
        trans_pos_final = round([parameters.grid_dims(2:3)/2 10]);
        focus_pos_final = round([parameters.grid_dims(2:3)/2 10+(parameters.expected_focal_distance_mm)/parameters.grid_step_mm]);
    else
        trans_pos_final = parameters.transducer.pos_grid;
        focus_pos_final = parameters.focus_pos_grid;
    end
end

parameters.transducer.pos_grid = trans_pos_final;
parameters.focus_pos_grid = focus_pos_final;

%% SETUP MEDIUM
disp('Setting up kwave medium...')

kwave_medium = setup_medium(parameters, skull_mask);

%% SETUP SOURCE
disp('Setting up kwave source...')

max_sound_speed = max(kwave_medium.sound_speed(:));
[kgrid, source, sensor, source_labels] = setup_grid_source_sensor(parameters, max_sound_speed, trans_pos_final, focus_pos_final);

%% RUN ACOUSTIC SIMULATION
% =========================================================================
disp('Starting acoustic simulations...')

% pathname for the input and output files (used only for non-interactive
% computations)
parameters.kwave_input_filename  = fullfile(parameters.data_path, sprintf('sub-%03d_%s_input%s.h5', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
parameters.kwave_output_filename = fullfile(parameters.data_path, sprintf('sub-%03d_%s_output%s.h5', subject_id, parameters.simulation_medium, parameters.results_filename_affix));

kwave_input_args = struct('PMLInside', false, ...
    'PMLSize', 'auto', ...
    'PlotPML', false);

if strcmp(parameters.simulation_medium, 'water_and_skull')
   kwave_input_args.DisplayMask = skull_edge;
end

filename_sensor_data = fullfile(output_dir, sprintf('sub-%03d_%s_results%s.mat', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
if confirm_overwriting(filename_sensor_data, parameters)
    % double-check since segmentation takes a long time
    if parameters.interactive == 0 || confirmation_dlg('Running the simulations will take a long time, are you sure?', 'Yes', 'No')
%         
%         if isfield(parameters,'run_simulations_with_qsub') && parameters.run_simulations_with_qsub == 1
%             % remember current folder to later go back
%             current_dir = cd;
%             cd(fullfile(parameters.data_path,'batch_job_logs'))
%             parameters.interactive = 0;
%             parameters.paths_to_add = path;
%             qsubfeval(@run_simulations, kgrid, kwave_medium, source, sensor, kwave_input_args, parameters, 'timreq',  60*60*7,  'memreq',  20*(1024^3),  'options', '-l "nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=5.0"', 'rerunable', 'yes');
%             disp('Simulations job submitted to the cluster, stopping for now. Re-run the pipeline when the job finishes to do post-processing.')
%             cd(current_dir);
%             return;
%         else 
        sensor_data = run_simulations(kgrid, kwave_medium, source, sensor, kwave_input_args, parameters);
%         end
        save(filename_sensor_data, 'sensor_data', 'kgrid', 'kwave_medium', 'source', 'sensor', 'kwave_input_args', 'parameters' )
    end
else
    load(filename_sensor_data, 'sensor_data')
end

%% Process results
disp('Processing the results of acoustic simulations...')
data_max = gather(sensor_data.p_max_all); % gather is used since it could be a GPU array

comp_grid_size = size(sensor_data.p_max_all);

Isppa_map = data_max.^2./(2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4; 

max_Isppa = max(Isppa_map(:));
after_exit_plane_mask = ones(comp_grid_size);
bowl_depth_grid = round((parameters.transducer.curv_radius_mm-parameters.transducer.dist_to_plane_mm)/parameters.grid_step_mm);
if trans_pos_final(3) > comp_grid_size(3)/2
    after_exit_plane_mask(:,:,(trans_pos_final(3)-bowl_depth_grid):end) = 0;
else
    after_exit_plane_mask(:,:,1:(trans_pos_final(3)+bowl_depth_grid)) = 0;
end

[max_Isppa_after_exit_plane, Ix_eplane, Iy_eplane, Iz_eplane] = masked_max_3d(Isppa_map, after_exit_plane_mask);
max_isppa_eplane_pos = [Ix_eplane, Iy_eplane, Iz_eplane];
disp('Final transducer, expected focus, and max ISPPA positions')
[trans_pos_final, focus_pos_final, max_isppa_eplane_pos']
real_focal_distance = norm(max_isppa_eplane_pos-trans_pos_final)*parameters.grid_step_mm;

if strcmp(parameters.simulation_medium, 'water_and_skull')
    
    [max_Isppa_brain, Ix_brain, Iy_brain, Iz_brain] = masked_max_3d(Isppa_map, segmented_image_cropped>0 & segmented_image_cropped<3);
    [max_pressure_brain, Px_brain, Py_brain, Pz_brain] = masked_max_3d(data_max, segmented_image_cropped>0 & segmented_image_cropped<3);
    [max_Isppa_skull, Ix_skull, Iy_skull, Iz_skull] = masked_max_3d(Isppa_map, segmented_image_cropped==3);
    [max_pressure_skull, Px_skull, Py_skull, Pz_skull] = masked_max_3d(data_max, segmented_image_cropped==3);
    [max_Isppa_skin, Ix_skin, Iy_skin, Iz_skin] = masked_max_3d(Isppa_map, segmented_image_cropped==3);
    [max_pressure_skin, Px_skin, Py_skin, Pz_skin] = masked_max_3d(data_max, segmented_image_cropped==3);
    highlighted_pos = [Ix_brain, Iy_brain, Iz_brain];

    writematrix([subject_id max_Isppa max_Isppa_after_exit_plane real_focal_distance max_Isppa_brain max_Isppa_skull max_pressure_brain max_pressure_skull max_Isppa_skin max_pressure_skin], output_pressure_file);
else
    max_Isppa = max(Isppa_map(:));
    highlighted_pos = max_isppa_eplane_pos;

    writematrix([subject_id max_Isppa max_Isppa_after_exit_plane real_focal_distance], output_pressure_file);

end
output_plot = fullfile(output_dir,sprintf('sub-%03d_%s_isppa%s.png', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
plot_isppa_over_image(Isppa_map, segmented_image_cropped, source_labels, after_exit_plane_mask, {'y', trans_pos_final(2)}, trans_pos_final, focus_pos_final, highlighted_pos)

export_fig(output_plot, '-native')

if isfield(parameters, 'run_posthoc_water_sims') && parameters.run_posthoc_water_sims && strcmp(parameters.medium, 'water_and_skull')
    new_parameters = parameters;
    new_parameters.simulation_medium = 'water';
    new_parameters.default_grid_dims = new_parameters.grid_dims;
    single_subject_pipeline(subject_id, new_parameters);
end
end