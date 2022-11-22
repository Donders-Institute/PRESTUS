function [output_pressure_file, parameters] = single_subject_pipeline(subject_id, parameters)
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    %                       Single subject pipeline                     %
    %                                                                   %
    % This serves as the main pipeline for simulating the sonication    %
    % effects on the bone and neural tissue of individual subjects.     %
    % In practice, neither this pipeline nor any of the functions in    %
    % the functions folder need to be altered to run the simulations.   %
    %                                                                   %
    % Parameters are loaded in using a custom config file.              %
    % The file 'default_config' contains all possible to-be altered     %
    % parameters, but not all of them need to be used to succesfully    %
    % run the pipeline.                                                 %
    %                                                                   %
    % Some notes:                                                       %
    % - The pipeline is only able to simulate one transducer at a time, %
    % meaning that the pipeline has to be run once for each transducer  %
    % used in each subject.                                             %
    % - Matlab 2019b must be used since k-wave was no longer updated    %
    % for the release of Matlab 2020a onwards.                          %
    % - 'subject_id' must be a number.                                  %
    % - 'parameters' is a structure (see default_config for options)    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    fprintf('Starting processing for subject %i %s\n',subject_id, parameters.results_filename_affix)
    
    % Adds the paths to the 'functions' and 'toolboxes' folders
    currentLoc = fileparts(mfilename("fullpath"));
    functionsLoc = fullfile(currentLoc,'functions');
    toolboxesLoc = fullfile(currentLoc,'toolboxes');
    allPaths = regexp(path,pathsep,'Split');

    if ~any(ismember(functionsLoc,allPaths))
        addpath(functionsLoc);
    else
    end

    if ~any(ismember(toolboxesLoc,allPaths))
        addpath(genpath(toolboxesLoc));
    else
    end

    % If there are paths to be added, add them; this is mostly for batch runs
    if isfield(parameters,'paths_to_add') && ~isempty(parameters.paths_to_add)
        for nPaths = length(parameters.paths_to_add)
            addpath(parameters.paths_to_add(nPaths))
        end
    end

    % If the path and subpaths need to be added, use this instead
    if isfield(parameters,'subpaths_to_add') && ~isempty(parameters.subpaths_to_add)
        for nPaths = length(parameters.subpaths_to_add)
            addpath(genpath(parameters.subpaths_to_add(nPaths)))
        end
    end

    % Sanitize output file affix
    sanitized_affix = regexprep(parameters.results_filename_affix,'[^a-zA-Z0-9_]','_');
    if ~strcmp(sanitized_affix, parameters.results_filename_affix)
        fprintf('The original `results_filename_affix` was sanitized, "%s" will be used instead of "%s"\n', sanitized_affix, parameters.results_filename_affix)
        parameters.results_filename_affix = sanitized_affix;
    end

    % Create and use subfolders
    if isfield(parameters,'output_subfolder') && parameters.output_subfolder == 1
        output_dir = (sprintf('%ssub-%03d/', parameters.data_path, subject_id));
        if ~exist(output_dir, 'file' )
            mkdir(output_dir);
        end
    else
        output_dir = fullfile(parameters.data_path, 'sim_outputs/');
        if ~exist(output_dir, 'file')
            mkdir(output_dir)
        end
    end
    
    % save parameters to have a backlog
    parameters_file = fullfile(output_dir,sprintf('sub-%03d_parameters_%s.mat', subject_id, datestr(now,'dd_mm_yyyy_HHMMSS_FFF')));
    save(parameters_file, 'parameters')
    
    % Add subject_id and output_dir to parameters to pass arguments to functions more easily
    parameters.subject_id = subject_id;
    parameters.output_dir = output_dir;
    
    
    %% Start of simulations
    % Creates an output file to which output is written at a later stage
    output_pressure_file = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_isppa%s.csv', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
    % Tries an alternative method to calculate the expected focal distance
    % if none is entered into the config file
    if ~isfield(parameters, 'expected_focal_distance_mm')
        disp('Expected focal distance is not specified, trying to get it from positions on T1 grid')
        if ~isfield(parameters.transducer, 'pos_t1_grid') || ~isfield(parameters, 'focus_pos_t1_grid')
            error('Either the transducer position or the focus position on T1 grid are not specified, cannot compute the expected focal distance')
        end
        filename_t1 = fullfile(parameters.data_path, sprintf(parameters.t1_path_template, subject_id));
        if ~isfile(filename_t1)
            error('File does not exist: \r\n%s', filename_t1);
        end
        t1_info = niftiinfo(filename_t1);
        t1_grid_step_mm = t1_info.PixelDimensions(1);
        focal_distance_t1 = norm(parameters.focus_pos_t1_grid - parameters.transducer.pos_t1_grid);
        parameters.expected_focal_distance_mm = focal_distance_t1 * t1_grid_step_mm;
    end
    
    % Pre-processes the MRI data to segment the different forms of tissue
    % and visualise the position of the transducer with some help from SimNIBS.
    % For more documentation, see the 'preprocess_brain' function.
    if contains(parameters.simulation_medium, 'skull')|| strcmp(parameters.simulation_medium, 'layered')
        [medium_masks, segmented_image_cropped, skull_edge, trans_pos_final, focus_pos_final, t1_image_orig, t1_header, final_transformation_matrix, inv_final_transformation_matrix] = preprocess_brain(parameters, subject_id);
        if isempty(medium_masks)
            output_pressure_file = '';
            return;
        end
        parameters.grid_dims = size(medium_masks);
    else % In case simulations are not run in a skull of layered tissue,
        % alternative grid dimensions are set up
        assert(isfield(parameters, 'default_grid_dims'), 'The parameters structure should have the field grid_dims for the grid dimensions')
        parameters.grid_dims = parameters.default_grid_dims;
        if any(parameters.grid_dims==1)||length(parameters.grid_dims)==2
            parameters.grid_dims = squeeze(parameters.grid_dims);
            parameters.n_sim_dims = length(parameters.grid_dims);
            disp('One of the simulation grid dimensions is of length 1, assuming you want 2d simulations, dropping this dimension')
        end

        % Checks whether the transducer location and orientation are set,
        % and uses an arbitrary position if not
        medium_masks = [];
        segmented_image_cropped = zeros(parameters.grid_dims);
        if ~isfield(parameters.transducer, 'pos_grid') || ~isfield(parameters, 'focus_pos_grid')
            disp('Either grid or focus position is not set, positioning them arbitrarily based on the focal distance')
            % note that the focus position matters only for the orientation of the transducer
        end
        if ~isfield(parameters.transducer, 'pos_grid')
            trans_pos_final = round([parameters.grid_dims(1:(parameters.n_sim_dims-1))/2 parameters.pml_size+1]); % transducer positioned arbitrarily
        else
            trans_pos_final = parameters.transducer.pos_grid;
        end
        if ~isfield(parameters, 'focus_pos_grid')
            focus_pos_final = trans_pos_final;
            focus_pos_final(parameters.n_sim_dims) = round(focus_pos_final(parameters.n_sim_dims) + parameters.expected_focal_distance_mm/parameters.grid_step_mm);
        else
            focus_pos_final = parameters.focus_pos_grid;    
        end

    end

    % If dimension 1 is bigger than 2, the matrices of the focus and
    % transducer locations are transposed
    if size(trans_pos_final,1)>size(trans_pos_final, 2)
        focus_pos_final = focus_pos_final';
        trans_pos_final = trans_pos_final';
    end
    
    % If a PML layer is used to absorb waves reaching the edge of the grid,
    % this will check if there is enough room for a PML layer between the
    % transducer and the edge of the grid
    assert(min(abs([0,0,0;parameters.grid_dims]-trans_pos_final ),[],'all') > parameters.pml_size, 'The minimal distance between the transducer and the simulation grid boundary should be larger than the PML size. Adjust transducer position or the PML size')
    assert(min(abs([0,0,0;parameters.grid_dims]-focus_pos_final ),[],'all') > parameters.pml_size, 'The minimal distance between the focus position and the simulation grid boundary should be larger than the PML size. Adjust transducer position or the PML size')
    
    % Saves the new transducer and focus positions after all previous grid
    % manipulations
    parameters.transducer.pos_grid = trans_pos_final;
    parameters.focus_pos_grid = focus_pos_final;

    %% SETUP MEDIUM
    % For more documentation, see 'setup_medium'
    disp('Setting up kwave medium...')

    kwave_medium = setup_medium(parameters, medium_masks);

    %% SETUP SOURCE
    % For more documentation, see 'setup_grid_source_sensor'
    disp('Setting up kwave source...')

    if parameters.run_source_setup
        max_sound_speed = max(kwave_medium.sound_speed(:));
        [kgrid, source, sensor, source_labels] = setup_grid_source_sensor(parameters, max_sound_speed, trans_pos_final, focus_pos_final);
    end

    %% RUN ACOUSTIC SIMULATION
    % =========================================================================
    disp('Starting acoustic simulations...')

    % Pathname for the input and output files (used only for non-interactive computations)
    parameters.kwave_input_filename  = fullfile(parameters.data_path, sprintf('sub-%03d_%s_input%s.h5', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    parameters.kwave_output_filename = fullfile(parameters.data_path, sprintf('sub-%03d_%s_output%s.h5', subject_id, parameters.simulation_medium, parameters.results_filename_affix));

    % Defines the edge of the simulation as the edge of the PML layer (see line 148)
    kwave_input_args = struct('PMLInside', true, ...
        'PMLSize', parameters.pml_size, ...
        'PlotPML', true);

    if contains(parameters.simulation_medium, 'skull')|| strcmp(parameters.simulation_medium, 'layered')
       kwave_input_args.DisplayMask = skull_edge;
    end

    % Saves the sensory data in a .mat file with the simulation figures
    filename_sensor_data = fullfile(parameters.output_dir, sprintf('sub-%03d_%s_results%s.mat', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
    % Run the acoustic simulations
    % See 'run_simulations' for more documentation
    if parameters.run_acoustic_sims && confirm_overwriting(filename_sensor_data, parameters) && (parameters.interactive == 0 || confirmation_dlg('Running the simulations will take a long time, are you sure?', 'Yes', 'No'))
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
            save(filename_sensor_data, 'sensor_data', 'kgrid', 'kwave_medium', 'source', 'sensor', 'kwave_input_args', 'parameters' ,'-v7.3')
        
    else
        load(filename_sensor_data, 'sensor_data')
    end

    %% Process results
    disp('Processing the results of acoustic simulations...')

    % What is the highest pressure level for every gridpoint
    data_max = gather(sensor_data.p_max_all); % gather is used since it could be a GPU array

    % Calculates the Isppa for every gridpoint
    Isppa_map = data_max.^2./(2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4; 
    % Calculates the max Isppa
    max_Isppa = max(Isppa_map(:));

    % Creates the foundation for a mask before the exit plane to calculate max values outside of it
    comp_grid_size = size(sensor_data.p_max_all);
    after_exit_plane_mask = ones(comp_grid_size);

    bowl_depth_grid = round((parameters.transducer.curv_radius_mm-parameters.transducer.dist_to_plane_mm)/parameters.grid_step_mm);

    % Places the exit plane mask in the grid, adjusted to the amount of dimensions
    if parameters.n_sim_dims == 3
        if trans_pos_final(3) > comp_grid_size(3)/2
            after_exit_plane_mask(:,:,(trans_pos_final(parameters.n_sim_dims)-bowl_depth_grid):end) = 0;
        else
            after_exit_plane_mask(:,:,1:(trans_pos_final(parameters.n_sim_dims)+bowl_depth_grid)) = 0;
        end
    else
        if trans_pos_final(2) > comp_grid_size(2)/2
            after_exit_plane_mask(:,(trans_pos_final(parameters.n_sim_dims)-bowl_depth_grid):end) = 0;
        else
            after_exit_plane_mask(:,1:(trans_pos_final(parameters.n_sim_dims)+bowl_depth_grid)) = 0;
        end
    end

    % Calculates the X, Y and Z coordinates of the max. intensity
    [max_Isppa_after_exit_plane, Ix_eplane, Iy_eplane, Iz_eplane] = masked_max_3d(Isppa_map, after_exit_plane_mask);
    
    % Combines these coordinates into a point of max. intensity in the grid
    if parameters.n_sim_dims==3
        max_isppa_eplane_pos = [Ix_eplane, Iy_eplane, Iz_eplane];
    else 
        max_isppa_eplane_pos = [Ix_eplane, Iy_eplane];
    end
    disp('Final transducer, expected focus, and max ISPPA positions')

    % Calculates the average Isppa within a circle around the target
    [trans_pos_final', focus_pos_final', max_isppa_eplane_pos']
    real_focal_distance = norm(max_isppa_eplane_pos-trans_pos_final)*parameters.grid_step_mm;
    avg_radius = round(parameters.focus_area_radius/parameters.grid_step_mm); %grid
    avg_isppa_around_target = Isppa_map((focus_pos_final(1)-avg_radius):(focus_pos_final(1)+avg_radius),...
        (focus_pos_final(2)-avg_radius):(focus_pos_final(2)+avg_radius),...
        (focus_pos_final(3)-avg_radius):(focus_pos_final(3)+avg_radius));
    avg_isppa_around_target = mean(avg_isppa_around_target(:));
    
    % Reports the Isppa within the original stimulation target
    isppa_at_target = Isppa_map(focus_pos_final(1),focus_pos_final(2),focus_pos_final(3));
    
    % Overwrites the max Isppa by dividing it up into the max Isppa for
    % each layer in case a layered simulation_medium was selected
    if contains(parameters.simulation_medium, 'skull') || strcmp(parameters.simulation_medium, 'layered')
        [max_Isppa_brain, Ix_brain, Iy_brain, Iz_brain] = masked_max_3d(Isppa_map, segmented_image_cropped>0 & segmented_image_cropped<3);
        half_max = Isppa_map >= max_Isppa_brain/2 & segmented_image_cropped>0 & segmented_image_cropped<3;
        half_max_ISPPA_volume_brain = sum(half_max(:))*(parameters.grid_step_mm^3);

        [max_pressure_brain, Px_brain, Py_brain, Pz_brain] = masked_max_3d(data_max, segmented_image_cropped>0 & segmented_image_cropped<3);
        [max_Isppa_skull, Ix_skull, Iy_skull, Iz_skull] = masked_max_3d(Isppa_map, segmented_image_cropped==4);
        [max_pressure_skull, Px_skull, Py_skull, Pz_skull] = masked_max_3d(data_max, segmented_image_cropped==4);
        [max_Isppa_skin, Ix_skin, Iy_skin, Iz_skin] = masked_max_3d(Isppa_map, segmented_image_cropped==5);
        [max_pressure_skin, Px_skin, Py_skin, Pz_skin] = masked_max_3d(data_max, segmented_image_cropped==5);
        highlighted_pos = [Ix_brain, Iy_brain, Iz_brain];
        real_focal_distance = norm(highlighted_pos-trans_pos_final)*parameters.grid_step_mm;

        writetable(table(subject_id, max_Isppa, max_Isppa_after_exit_plane, real_focal_distance, max_Isppa_brain, max_Isppa_skull, max_pressure_brain, max_pressure_skull, max_Isppa_skin, max_pressure_skin, Ix_brain, Iy_brain, Iz_brain, trans_pos_final, focus_pos_final, isppa_at_target, avg_isppa_around_target, half_max_ISPPA_volume_brain), output_pressure_file);
    else % If no layered tissue was selected, the max Isppa is highlithed on the plane and written in a table.
        max_Isppa = max(Isppa_map(:)); %  Does this step need to be included? already done at line 225.
        highlighted_pos = max_isppa_eplane_pos;
        writetable(table(subject_id, max_Isppa, max_Isppa_after_exit_plane, real_focal_distance, trans_pos_final, focus_pos_final, isppa_at_target, avg_isppa_around_target), output_pressure_file);
    end

    % Plots the Isppa on the segmented figure
    output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_isppa%s.png', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
    if parameters.n_sim_dims==3
        plot_isppa_over_image(Isppa_map, segmented_image_cropped, source_labels, parameters, {'y', focus_pos_final(2)}, trans_pos_final, focus_pos_final, highlighted_pos);
    else
        plot_isppa_over_image_2d(Isppa_map, segmented_image_cropped, source_labels, parameters,  trans_pos_final, focus_pos_final, highlighted_pos);
    end
    export_fig(output_plot, '-native')
    close;
    
    % Plots the Isppa on the original T1 figure & in the MNI space for skull & layered mediums only
    if contains(parameters.simulation_medium, 'skull') || strcmp(parameters.simulation_medium, 'layered')
        backtransf_coordinates = round(tformfwd([trans_pos_final;  focus_pos_final; highlighted_pos], inv_final_transformation_matrix));

        isppa_orig_file = fullfile(parameters.output_dir, sprintf('sub-%03d_final_isppa_orig_coord%s.nii.gz',...
            subject_id, parameters.results_filename_affix));
        
        if confirm_overwriting(isppa_orig_file, parameters)
            % Transforms the Isppa map to the original T1 figure dimensions and orientation
            isppa_map_backtransf = tformarray(Isppa_map, inv_final_transformation_matrix, ...
                                    makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], size(t1_image_orig), [], 0) ;

            isppa_hdr = t1_header;
            isppa_hdr.Datatype = 'single';

            niftiwrite(isppa_map_backtransf, isppa_orig_file, isppa_hdr, 'Compressed', true)
        else 
            isppa_map_backtransf = niftiread(isppa_orig_file);
        end
        % Creates a visual overlay of the transducer
        [~, source_labels] = transducer_setup(parameters.transducer, backtransf_coordinates(1,:), backtransf_coordinates(2,:), ...
                                                    size(t1_image_orig), t1_header.PixelDimensions(1));

        % Plots the Isppa to fit the untransformed figure
        plot_isppa_over_image(isppa_map_backtransf, t1_image_orig, source_labels, ...
            parameters, {'y', backtransf_coordinates(2,2)}, backtransf_coordinates(1,:), ...
            backtransf_coordinates(2,:), backtransf_coordinates(3,:), 'show_rectangles', 0, 'grid_step', t1_header.PixelDimensions(1));
        output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_isppa_orig%s.png', subject_id, parameters.simulation_medium, parameters.results_filename_affix));

        export_fig(output_plot, '-native')
        close;
        
        headreco_folder = fullfile(parameters.data_path, sprintf('m2m_sub-%03d', subject_id));
        
        max_pressure_mni_file = fullfile(parameters.output_dir, sprintf('sub-%03d_final_pressure_MNI%s.nii.gz', subject_id, parameters.results_filename_affix));
        isppa_map_mni_file  = fullfile(parameters.output_dir, sprintf('sub-%03d_final_isppa_MNI%s.nii.gz', subject_id, parameters.results_filename_affix));
        segmented_image_mni_file = fullfile(parameters.output_dir, sprintf('sub-%03d_segmented_MNI%s.nii.gz', subject_id, parameters.results_filename_affix));

        isppa_map_mni = convert_final_to_MNI(Isppa_map, headreco_folder, inv_final_transformation_matrix, parameters, 'nifti_filename', isppa_map_mni_file);
        segmented_image_mni = convert_final_to_MNI(segmented_image_cropped, headreco_folder, inv_final_transformation_matrix,  parameters, 'nifti_filename', segmented_image_mni_file, 'nifti_data_type', 'uint8', 'BitsPerPixel', 8);
        max_pressure_mni = convert_final_to_MNI(data_max, headreco_folder, inv_final_transformation_matrix, parameters, 'nifti_filename', max_pressure_mni_file);
        
    end

    % Runs the heating simulation
    if isfield(parameters, 'run_heating_sims') && parameters.run_heating_sims 
        % Creates an output file to which output is written at a later stage
        filename_heating_data = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_heating_res%s.mat', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        
        % Clear sensor mask
        sensor.mask = zeros(size(sensor.mask));
        
        % Sets up parameters for the heating simulation and runs it
        if confirm_overwriting(filename_heating_data, parameters) && (parameters.interactive == 0 || confirmation_dlg('Running the thermal simulations will take a long time, are you sure?', 'Yes', 'No')) 
            
            % Set sensor along the focal axis 
            heating_window_dims = ones(2,3);
            for i = 1:2
                heating_window_dims(:,i) = [max(1, -parameters.thermal.sensor_xy_halfsize + parameters.transducer.pos_grid(i)), min(parameters.grid_dims(i), parameters.thermal.sensor_xy_halfsize + parameters.transducer.pos_grid(i))];
            end

            % Sets up a window to perform simulations in
            heating_window_dims(2,3) = parameters.grid_dims(3);
            sensor.mask(heating_window_dims(1,1):heating_window_dims(2,1), heating_window_dims(1,2):heating_window_dims(2,2), :) = 1;

            % For more documentation, see 'run_heating_simulations'
            [kwaveDiffusion, time_status_seq, maxT, focal_planeT]= run_heating_simulations(sensor_data, kgrid, kwave_medium, sensor, source, parameters, trans_pos_final);
            save(filename_heating_data, 'kwaveDiffusion','time_status_seq','heating_window_dims','sensor','maxT','focal_planeT','-v7.3');
        else 
            load(filename_heating_data);
        end

        % Sets up an empty medium mask if non is specified
        if isempty(medium_masks)
            medium_masks = zeros(parameters.grid_dims);
        end

        % Creates an output table for temperature readings
        output_table = readtable(output_pressure_file);
        output_table.maxT = gather(max(maxT, [], 'all'));

        % Overwrites the max temperature by dividing it up for each layer
        % in case a layered simulation_medium was selected
        if contains(parameters.simulation_medium, 'skull') || strcmp(parameters.simulation_medium, 'layered')
            output_table.maxT_brain = gather(masked_max_3d(maxT, segmented_image_cropped>0 & segmented_image_cropped<3));
            output_table.maxT_skull = gather(masked_max_3d(maxT, segmented_image_cropped==4));
            output_table.maxT_skin = gather(masked_max_3d(maxT, segmented_image_cropped==5));
        end
        writetable(output_table, output_pressure_file);

        % Creates a visual overlay of the transducer
        [~, source_labels] = transducer_setup(parameters.transducer, trans_pos_final, focus_pos_final, ...
                                                    size(segmented_image_cropped), t1_header.PixelDimensions(1));

        % Creates a line graph and a video of the heating effects
        plot_heating_sims(focal_planeT, time_status_seq, parameters, trans_pos_final, medium_masks);
        
        % Plots the maximum temperature in the segmented brain
        if max(maxT(:)) < 38
            temp_color_range = [37, 38];
        else
            temp_color_range = [37, max(maxT(:))];
        end
        maxT = gather(maxT);
        plot_isppa_over_image(maxT, segmented_image_cropped, source_labels, parameters, {'y', focus_pos_final(2)}, trans_pos_final, focus_pos_final, highlighted_pos, 'isppa_color_range', temp_color_range );
        output_plot = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_maxT%s.png', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        export_fig(output_plot, '-native')
        close;
        
        
        headreco_folder = fullfile(parameters.data_path, sprintf('m2m_sub-%03d', subject_id));
        heating_data_mni_file = fullfile(parameters.output_dir,sprintf('sub-%03d_heating_MNI%s.nii.gz', subject_id, parameters.results_filename_affix));
        convert_final_to_MNI(maxT, headreco_folder, inv_final_transformation_matrix, parameters, 'nifti_filename', heating_data_mni_file);

    end
    
    % Runs posthoc water simulations to check sonication parameters of the
    % transducer in free water
    if isfield(parameters, 'run_posthoc_water_sims') && parameters.run_posthoc_water_sims && ...
            (contains(parameters.simulation_medium, 'skull') || contains(parameters.simulation_medium, 'layered'))
        new_parameters = parameters;
        new_parameters.simulation_medium = 'water';
        new_parameters.run_heating_sims = 0;
        new_parameters.default_grid_dims = new_parameters.grid_dims;
        single_subject_pipeline(subject_id, new_parameters);
    end
    
end