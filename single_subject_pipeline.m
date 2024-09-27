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
    % - At least Matlab 2022b, SimNIBS 4.0 and k-Wave 1.4 must be used  %
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
        disp(['Adding ', functionsLoc]);
    else
    end

    if ~any(ismember(toolboxesLoc,allPaths))
        addpath(genpath(toolboxesLoc));
        disp(['Adding ', toolboxesLoc, 'and subfolders']);
    else
    end

    % If there are paths to be added, add them; this is mostly for batch runs
    if isfield(parameters,'paths_to_add') && ~isempty(parameters.paths_to_add)
        for nPaths = 1:length(parameters.paths_to_add)
            addpath(parameters.paths_to_add{nPaths})
            disp(['Adding ', parameters.paths_to_add{nPaths}]);
        end
    end

    % If the path and subpaths need to be added, use this instead
    if isfield(parameters,'subpaths_to_add') && ~isempty(parameters.subpaths_to_add)
        for nPaths = 1:length(parameters.subpaths_to_add)
            addpath(genpath(parameters.subpaths_to_add{nPaths}))
            disp(['Adding ', parameters.subpaths_to_add{nPaths}, 'and subfolders']);
        end
    end
    
    % test that kwave is added
    if ~exist('makeBowl','file')
        error('kwave not added');
    end

    % Make subfolder (if enabled) and check if directory exists
    if isfield(parameters,'subject_subfolder') && parameters.subject_subfolder == 1
        parameters.output_dir = fullfile(parameters.temp_output_dir, sprintf('sub-%03d', subject_id));
    else 
        parameters.output_dir = parameters.temp_output_dir;
    end

    % specify dedicated subfolder for debugging contents
    parameters.debug_dir = fullfile(parameters.output_dir, 'debug');

    if ~isfolder(parameters.output_dir)
        mkdir(parameters.output_dir);
    end
    if ~isfolder(parameters.debug_dir)
        mkdir(parameters.debug_dir);
    end
    if isfield(parameters,'seg_path') && ~isfolder(parameters.seg_path)
        mkdir(parameters.seg_path);
    end
    
    % Save parameters to have a backlog
    parameters_file = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_parameters%s_%s.mat', ...
        subject_id, parameters.results_filename_affix, ...
        datestr(now,'ddmmyy_HHMM')));
    save(parameters_file, 'parameters')
    
    % Add subject_id to parameters to pass arguments to functions more easily
    parameters.subject_id = subject_id;
    
    %% Extra settings needed for better usability
    warning('off','MATLAB:prnRenderer:opengl'); % suppress unneccessary warnings from export_fig when running without OpenGL
    
    % output GPU information (if requested)
    if strcmp(parameters.code_type, 'cuda') || strcmp(parameters.code_type, 'matlab_gpu')
        gpuDevice()
    end

    %% Start of simulations
    % Creates an output file to which output is written at a later stage
    output_pressure_file = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_output_table%s.csv', ...
        subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
    % Tries an alternative method to calculate the expected focal distance
    % if none is entered into the config file
    if ~isfield(parameters, 'expected_focal_distance_mm')
        disp('Expected focal distance is not specified, trying to get it from positions on T1 grid')
        if ~isfield(parameters.transducer, 'pos_t1_grid') || ~isfield(parameters, 'focus_pos_t1_grid')
            error('Either the transducer position or the focus position on T1 grid are not specified, cannot compute the expected focal distance')
        end
        filename_t1 = dir(fullfile(parameters.data_path, sprintf(parameters.t1_path_template, subject_id)));
        if isempty(filename_t1)
            error('File does not exist: \r\n%s', filename_t1);
        end
        % select first match in case multiple are found
        filename_t1 = fullfile(filename_t1(1).folder, filename_t1(1).name);
        t1_info = niftiinfo(filename_t1);
        t1_grid_step_mm = t1_info.PixelDimensions(1);
        focal_distance_t1 = norm(parameters.focus_pos_t1_grid - parameters.transducer.pos_t1_grid);
        parameters.expected_focal_distance_mm = focal_distance_t1 * t1_grid_step_mm; 
        clear filename_t1 t1_info focal_distance_t1 t1_grid_step_mm;
    end
    
    % if there is no specification of usepseudoCT, go for default of 0
    if ~isfield(parameters, 'usepseudoCT')
        parameters.usepseudoCT = 0;
    end

    % Pre-processes the MRI data to segment the different forms of tissue
    % and visualise the position of the transducer with some help from SimNIBS.
    % For more documentation, see the 'preprocess_brain' function.
    if contains(parameters.simulation_medium, 'skull')|| strcmp(parameters.simulation_medium, 'layered')
        [medium_masks, segmented_image_cropped, skull_edge, trans_pos_final, ...
        focus_pos_final, t1_image_orig, t1_header, final_transformation_matrix, ...
        inv_final_transformation_matrix] = preprocess_brain(parameters, subject_id);
        if isempty(medium_masks)
            output_pressure_file = '';
            return;
        end
        parameters.grid_dims = size(medium_masks);
    else % In case simulations are not run in a skull of layered tissue, alternative grid dimensions are set up
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

    if parameters.usepseudoCT == 1
        kwave_medium = setup_medium(parameters, medium_masks, segmented_image_cropped);
    else
        kwave_medium = setup_medium(parameters, medium_masks);
    end

    % split temp_0 from kwave_medium (due to kwave checks)
    temp_0 = kwave_medium.temp_0;
    kwave_medium = rmfield(kwave_medium, 'temp_0');
    
    % save medium images for debugging
    if (contains(parameters.simulation_medium, 'skull') || contains(parameters.simulation_medium, 'layered'))
                
        orig_hdr = t1_header; % header is based on original T1w (always present)
        orig_hdr.Datatype = 'single';

        filename_density = fullfile(parameters.debug_dir, sprintf('density'));
        density = single(tformarray(kwave_medium.density, inv_final_transformation_matrix, ...
            makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0)) ;
        if ~isfile(filename_density)
            niftiwrite(density, filename_density, orig_hdr, 'Compressed',true);
        end
        clear filename_density density;
        
        filename_sound_speed = fullfile(parameters.debug_dir, sprintf('sound_speed'));
        sound_speed = single(tformarray(kwave_medium.sound_speed, inv_final_transformation_matrix, ...
            makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0)) ;
        if ~isfile(filename_sound_speed)
            niftiwrite(sound_speed, filename_sound_speed, orig_hdr, 'Compressed',true);
        end
        clear filename_sound_speed sound_speed;
        
        filename_alpha_coeff = fullfile(parameters.debug_dir, sprintf('alpha_coeff'));
        alpha_coeff = single(tformarray(kwave_medium.alpha_coeff, inv_final_transformation_matrix, ...
            makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], orig_hdr.ImageSize, [], 0)) ;
        if ~isfile(filename_alpha_coeff)
            niftiwrite(alpha_coeff, filename_alpha_coeff, orig_hdr, 'Compressed',true);
        end
        clear filename_alpha_coeff alpha_coeff;
    end    

    %% SETUP SOURCE
    % For more documentation, see 'setup_grid_source_sensor'
    disp('Setting up kwave source...')

    if parameters.run_source_setup
        max_sound_speed = max(kwave_medium.sound_speed(:));
        [kgrid, source, sensor, source_labels] = setup_grid_source_sensor(parameters, max_sound_speed, trans_pos_final, focus_pos_final);
        % check stability
        % If estimated time step is smaller than the time step based on
        % default CFL, the estimated time step is used to redefine
        % transducer and sensor. Note: the estimated time step doesn't 
        % guarantee a stable simulation. If Nan numbers are acquired as
        % a result, you may want to try a time step smaller than the
        % estimated time step.	  
        disp('Check stability...')
        dt_stability_limit = checkStability(kgrid, kwave_medium);
        if ~isinf(dt_stability_limit) && kgrid.dt > dt_stability_limit
			disp('Adapt time step for simulation stability...')
            grid_time_step = dt_stability_limit*0.90; % use 90% of the limit (which are only an approximation in the heterogenous medium case: http://www.k-wave.org/documentation/checkStability.php)
            [kgrid, source, sensor, source_labels] = setup_grid_source_sensor(parameters, max_sound_speed, trans_pos_final, focus_pos_final, grid_time_step);
        end
    end

    %% RUN ACOUSTIC SIMULATION
    % =========================================================================
    disp('Starting acoustic simulations...')

    % Pathname for the input and output files (used only for non-interactive computations)
    parameters.kwave_input_filename  = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_%s_input%s.h5', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    parameters.kwave_output_filename = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_%s_output%s.h5', subject_id, parameters.simulation_medium, parameters.results_filename_affix));

    % Defines the edge of the simulation as the edge of the PML layer (see line 148)
    kwave_input_args = struct('PMLInside', true, ...
        'PMLSize', parameters.pml_size, ...
        'PlotPML', true);

    if contains(parameters.simulation_medium, 'skull')|| strcmp(parameters.simulation_medium, 'layered')
        kwave_input_args.DisplayMask = skull_edge;
    end

    % Looks up sensor data for use in simulations
    filename_sensor_data = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_%s_results%s.mat', ...
        subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
    % Run the acoustic simulations
    % See 'run_simulations' for more documentation
    if parameters.run_acoustic_sims && confirm_overwriting(filename_sensor_data, parameters) && (parameters.interactive == 0 || confirmation_dlg('Running the simulations will take a long time, are you sure?', 'Yes', 'No'))
        sensor_data = run_simulations(kgrid, kwave_medium, source, sensor, kwave_input_args, parameters);
        save(filename_sensor_data, 'sensor_data', 'kgrid', 'kwave_medium', 'source', 'sensor', 'kwave_input_args', 'parameters' ,'-v7.3')
    else
        disp('Skipping, the file already exists, loading it instead.')
        load(filename_sensor_data, 'sensor_data')
    end

    %% Process results
    disp('Processing the results of acoustic simulations...')

    % What is the highest pressure level for every gridpoint
    data_max = gather(sensor_data.p_max_all); % gather is used since it could be a GPU array
    max_pressure = max(data_max(:));

    % Calculates the Isppa for every gridpoint
    Isppa_map = data_max.^2./(2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4;
    % Calculates the max Isppa
    max_Isppa = max(Isppa_map(:));

    % Calculates the Mechanical Index for every gridpoint
    MI_map = (data_max/10^6)/sqrt((parameters.transducer.source_freq_hz/10^6));

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
    distance_target_real_maximum = norm(max_isppa_eplane_pos-focus_pos_final)*parameters.grid_step_mm;
    avg_radius = round(parameters.focus_area_radius/parameters.grid_step_mm); %grid
    avg_isppa_around_target = Isppa_map(...
        (focus_pos_final(1)-avg_radius):(focus_pos_final(1)+avg_radius),...
        (focus_pos_final(2)-avg_radius):(focus_pos_final(2)+avg_radius),...
        (focus_pos_final(3)-avg_radius):(focus_pos_final(3)+avg_radius));
    avg_isppa_around_target = mean(avg_isppa_around_target(:));
    
    % Reports the Isppa within the original stimulation target
    isppa_at_target = Isppa_map(focus_pos_final(1),focus_pos_final(2),focus_pos_final(3));
    
    % Creates a logical skull mask and register skull_ids
    labels = fieldnames(parameters.layer_labels);
    skull_i = find(strcmp(labels, 'skull_cortical'));
    trabecular_i = find(strcmp(labels, 'skull_trabecular'));
    all_skull_ids = [skull_i, trabecular_i];
    skull_mask = ismember(medium_masks,all_skull_ids);
    brain_i = find(strcmp(labels, 'brain'));
    brain_mask = ismember(medium_masks,brain_i);
    skin_i = find(strcmp(labels, 'skin'));
    skin_mask = ismember(medium_masks,skin_i);
    
    % Overwrites the max Isppa by dividing it up into the max Isppa for
    % each layer in case a layered simulation_medium was selected
    if contains(parameters.simulation_medium, 'skull') || strcmp(parameters.simulation_medium, 'layered')
        [max_Isppa_brain, Ix_brain, Iy_brain, Iz_brain] = masked_max_3d(Isppa_map, brain_mask);
        [min_Isppa_brain] = min(Isppa_map(brain_mask));
        half_max = Isppa_map >= max_Isppa_brain/2 & brain_mask;
        half_max_ISPPA_volume_brain = sum(half_max(:))*(parameters.grid_step_mm^3);
        [max_pressure_brain, Px_brain, Py_brain, Pz_brain] = masked_max_3d(data_max, brain_mask);
        [max_MI_brain, Px_brain, Py_brain, Pz_brain] = masked_max_3d(MI_map, brain_mask);
        
        [max_Isppa_skull, Ix_skull, Iy_skull, Iz_skull] = masked_max_3d(Isppa_map, skull_mask);
        [max_pressure_skull, Px_skull, Py_skull, Pz_skull] = masked_max_3d(data_max, skull_mask);
        [max_MI_skull, Px_skull, Py_skull, Pz_skull] = masked_max_3d(MI_map, skull_mask);
        
        [max_Isppa_skin, Ix_skin, Iy_skin, Iz_skin] = masked_max_3d(Isppa_map, skin_mask);
        [max_pressure_skin, Px_skin, Py_skin, Pz_skin] = masked_max_3d(data_max, skin_mask);
        [max_MI_skin, Px_skin, Py_skin, Pz_skin] = masked_max_3d(MI_map, skin_mask);

        highlighted_pos = [Ix_brain, Iy_brain, Iz_brain];
        real_focal_distance = norm(highlighted_pos-trans_pos_final)*parameters.grid_step_mm;

        writetable(table(subject_id, max_Isppa, max_Isppa_after_exit_plane, real_focal_distance, max_Isppa_skin, max_Isppa_skull, max_Isppa_brain, max_pressure_skin, max_pressure_skull, max_pressure_brain, max_MI_skin, max_MI_skull, max_MI_brain, Ix_brain, Iy_brain, Iz_brain, trans_pos_final, focus_pos_final, isppa_at_target, avg_isppa_around_target, half_max_ISPPA_volume_brain), output_pressure_file);
    else % If no layered tissue was selected, the max Isppa is highlighted on the plane and written in a table.
        max_Isppa = max(Isppa_map(:)); %  Does this step need to be included? already done at line 225.
        highlighted_pos = max_isppa_eplane_pos;
        writetable(table(subject_id, max_Isppa, max_Isppa_after_exit_plane, max_pressure, real_focal_distance, trans_pos_final, focus_pos_final, isppa_at_target, avg_isppa_around_target), output_pressure_file);
    end

    % Plots the Isppa on the segmented image
    
    if parameters.n_sim_dims==3
        %options.isppa_color_range = [0.5, max_Isppa_brain];
        [~,~,~,~,~,~,~,h]=plot_isppa_over_image(...
            Isppa_map, segmented_image_cropped, source_labels, parameters, ...
            {'y', focus_pos_final(2)}, trans_pos_final, focus_pos_final, highlighted_pos);
    else
        [~,~,h]=plot_isppa_over_image_2d(Isppa_map, segmented_image_cropped, source_labels, parameters,  trans_pos_final, focus_pos_final, highlighted_pos);
    end
    output_plot = fullfile(parameters.output_dir,...
        sprintf('sub-%03d_%s_isppa%s.png', ...
        subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    saveas(h, output_plot, 'png')
    close(h);

    %% RUN HEATING SIMULATIONS
    % =========================================================================
    if isfield(parameters, 'run_heating_sims') && parameters.run_heating_sims 
        disp('Starting heating simulations...')
        % Creates an output file to which output is written at a later stage
        filename_heating_data = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_heating_res%s.mat', subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        
        % Clear sensor mask
        sensor.mask = zeros(size(sensor.mask));
        
        % Sets up parameters for the heating simulation and runs it
        if confirm_overwriting(filename_heating_data, parameters) && (parameters.interactive == 0 || confirmation_dlg('Running the thermal simulations will take a long time, are you sure?', 'Yes', 'No')) 
            
            % Set sensor along the focal axis 
            heating_window_dims = ones(2,3);
            for i = 1:2
                heating_window_dims(:,i) = [...
                    max(1, -parameters.thermal.sensor_xy_halfsize + parameters.transducer.pos_grid(i)), ...
                    min(parameters.grid_dims(i), parameters.thermal.sensor_xy_halfsize + parameters.transducer.pos_grid(i))];
            end

            % Sets up a window to perform simulations in
            heating_window_dims(2,3) = parameters.grid_dims(3);
            sensor.mask(heating_window_dims(1,1):heating_window_dims(2,1), heating_window_dims(1,2):heating_window_dims(2,2), :) = 1;

            % add starting temperature
            kwave_medium.temp_0 = temp_0;

            % For more documentation, see 'run_heating_simulations'
            [kwaveDiffusion, time_status_seq, maxT, focal_planeT, maxCEM43, focal_planeCEM43]= ...
                run_heating_simulations(sensor_data, kgrid, kwave_medium, sensor, source, parameters, trans_pos_final);
            
            % apply gather in case variables are GPU arrays
            maxT = gather(maxT);
            focal_planeT = gather(focal_planeT);
            maxCEM43 = gather(maxCEM43);
            focal_planeCEM43 = gather(focal_planeCEM43);

            save(filename_heating_data, 'kwaveDiffusion','time_status_seq',...
                'heating_window_dims','sensor','maxT','focal_planeT','maxCEM43','focal_planeCEM43','-v7.3');

        else 
            disp('Skipping, the file already exists, loading it instead.')
            load(filename_heating_data);
        end

        % variable may have been saved under a different name in the past (backcomp)
        if exist('CEM43')
            focal_planeCEM43=CEM43;
        end

        % Sets up an empty medium mask if non is specified
        if isempty(medium_masks)
            medium_masks = zeros(parameters.grid_dims);
        end

        % Creates an output table for temperature readings
        output_table = readtable(output_pressure_file);
        
        % convert GPU data 
        maxT = gather(maxT);
        CEM43 = gather(maxCEM43);

        output_table.maxT = max(maxT, [], 'all');
        output_table.maxCEM43 = max(CEM43, [], 'all');
        % Overwrites the max temperature by dividing it up for each layer
        % in case a layered simulation_medium was selected
        if contains(parameters.simulation_medium, 'skull') || strcmp(parameters.simulation_medium, 'layered')
            output_table.maxT_brain = masked_max_3d(maxT, brain_mask);
            output_table.maxT_skull = masked_max_3d(maxT, skull_mask); 
            output_table.maxT_skin = masked_max_3d(maxT, skin_mask);
            output_table.riseT_brain = masked_max_3d(maxT, brain_mask)-parameters.thermal.temp_0.brain;
            output_table.riseT_skull = masked_max_3d(maxT, skull_mask)-parameters.thermal.temp_0.skull; 
            output_table.riseT_skin = masked_max_3d(maxT, skin_mask)-parameters.thermal.temp_0.skin;
            output_table.CEM43_brain = masked_max_3d(CEM43, brain_mask);
            output_table.CEM43_skull = masked_max_3d(CEM43, skull_mask); 
            output_table.CEM43_skin = masked_max_3d(CEM43, skin_mask);
        end
        writetable(output_table, output_pressure_file);

        % Creates a visual overlay of the transducer
        [~, source_labels] = transducer_setup(parameters.transducer, trans_pos_final, focus_pos_final, ...
                                                    size(segmented_image_cropped), t1_header.PixelDimensions(1));
        % Creates a line graph and a video of the heating effects
        plot_heating_sims(focal_planeT, time_status_seq, parameters, trans_pos_final, medium_masks, focal_planeCEM43);
                
        % Plots the maximum temperature in the segmented brain
        if output_table.maxT < 38
            temp_color_range = [37, 38];
        else
            temp_color_range = [37, output_table.maxT];
        end

        % plot ISPPA superimposed on segmentation
        [~,~,~,~,~,~,~,h]=plot_isppa_over_image(...
            maxT, segmented_image_cropped, source_labels, parameters, ...
            {'y', focus_pos_final(2)}, trans_pos_final, focus_pos_final, ...
            highlighted_pos, 'isppa_color_range', temp_color_range );
        
        output_plot_filename = fullfile(parameters.output_dir,...
            sprintf('sub-%03d_%s_maxT%s.png',...
            subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);
    end

    %% Plots the data on the original T1 image and in MNI space
    if contains(parameters.simulation_medium, 'skull') || strcmp(parameters.simulation_medium, 'layered')
        backtransf_coordinates = round(tformfwd([trans_pos_final;  focus_pos_final; highlighted_pos], inv_final_transformation_matrix));

        data_types = ["isppa","MI","pressure","medium_masks"];
        if  isfield(parameters, 'run_heating_sims') && parameters.run_heating_sims 
            data_types  = [data_types, "heating", "heatrise", "CEM43"];
        end
        for data_type = data_types
            orig_file = fullfile(parameters.output_dir, sprintf('sub-%03d_final_%s_orig_coord%s',...
                subject_id, data_type, parameters.results_filename_affix));
            mni_file  = fullfile(parameters.output_dir, sprintf('sub-%03d_final_%s_MNI%s.nii.gz',...
                subject_id, data_type, parameters.results_filename_affix));

            if strcmp(data_type, "isppa")
                data = single(Isppa_map);
            elseif strcmp(data_type, "MI")
                data = single(MI_map);
            elseif strcmp(data_type, "pressure")
                data = single(data_max);
            elseif strcmp(data_type, "medium_masks")
                data = medium_masks;
            elseif strcmp(data_type, "heating")
                data = single(maxT);
            elseif strcmp(data_type, "heatrise")
                data = single(maxT-temp_0);
            elseif strcmp(data_type, "CEM43")
                data = single(CEM43);
            end
            orig_file_with_ext = strcat(orig_file, '.nii.gz');
    
            if confirm_overwriting(orig_file_with_ext, parameters)
                % Transforms the data to original T1 image dimensions and orientation
                orig_hdr = t1_header;

                if strcmp(data_type, "medium_masks")
                    data_backtransf = tformarray(uint8(data), inv_final_transformation_matrix, ...
                                            makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], size(t1_image_orig), [], 0) ;

                    orig_hdr.Datatype = 'uint8';
                    orig_hdr.BitsPerPixel = 8;
                else
                    data_backtransf = tformarray(data, inv_final_transformation_matrix, ...
                                            makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], size(t1_image_orig), [], 0) ;
                    
                    orig_hdr = t1_header;
                    orig_hdr.Datatype = 'single';
                end
                niftiwrite(data_backtransf, orig_file, orig_hdr, 'Compressed', true)
            else 
                data_backtransf = niftiread(orig_file_with_ext);
            end
            if strcmp(data_type, "isppa")
                % Creates a visual overlay of the transducer
                [~, source_labels] = transducer_setup(parameters.transducer, backtransf_coordinates(1,:), backtransf_coordinates(2,:), ...
                                                            size(t1_image_orig), t1_header.PixelDimensions(1));
                % Plots the Isppa over the untransformed image
                [~,~,~,~,~,~,~,h]=plot_isppa_over_image(data_backtransf, t1_image_orig, source_labels, ...
                    parameters, {'y', backtransf_coordinates(2,2)}, backtransf_coordinates(1,:), ...
                    backtransf_coordinates(2,:), backtransf_coordinates(3,:), ...
                    'show_rectangles', 0, 'grid_step', t1_header.PixelDimensions(1), ...
                    'isppa_color_range', [0 max_Isppa_brain], ...
                    'isppa_threshold_low', 0, 'isppa_threshold_high', max_Isppa_brain, ...
                    'rotation', 0); % rotation = 90 not implemented for transducer overlay
        
                output_plot_filename = fullfile(parameters.output_dir,...
                    sprintf('sub-%03d_%s_isppa_t1%s.png', ...
                    subject_id, parameters.simulation_medium, ...
                    parameters.results_filename_affix));
                saveas(h, output_plot_filename, 'png')
                close(h);
            end
            
            m2m_folder= fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', subject_id));
                
            if ~confirm_overwriting(mni_file, parameters)
                continue
            end
            if strcmp(parameters.segmentation_software, 'headreco')
                
                if strcmp(data_type, "medium_masks")
                   convert_final_to_MNI_matlab(data, m2m_folder, inv_final_transformation_matrix, parameters, 'nifti_filename', mni_file,  'nifti_data_type', 'uint8', 'BitsPerPixel', 8);
                else
                   convert_final_to_MNI_matlab(data, m2m_folder, inv_final_transformation_matrix, parameters, 'nifti_filename', mni_file);
                end
            elseif strcmp(parameters.segmentation_software, 'charm')
                convert_final_to_MNI_simnibs(orig_file_with_ext , m2m_folder, mni_file, parameters, 'interpolation_order', 0);
            end
        
        end
        
        % Since charm does not transform the T1 into MNI space, one is manually created here
        if strcmp(parameters.segmentation_software, 'charm')
            path_to_input_img = fullfile(m2m_folder,'T1.nii.gz');
            path_to_output_img = fullfile(m2m_folder,'toMNI','T1_to_MNI_post-hoc.nii.gz');

            if ~exist(path_to_output_img,'file')
                convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters)
            end
        end
    end
    
    %% Runs posthoc water simulations
    % To check sonication parameters of the transducer in free water
    if isfield(parameters, 'run_posthoc_water_sims') && parameters.run_posthoc_water_sims && ...
            (contains(parameters.simulation_medium, 'skull') || contains(parameters.simulation_medium, 'layered'))
        new_parameters = parameters;
        new_parameters.simulation_medium = 'water';
        new_parameters.run_heating_sims = 0;
        new_parameters.default_grid_dims = new_parameters.grid_dims;
        % restore subject-specific path to original path if done earlier in this function
        if isfield(new_parameters,'subject_subfolder') && new_parameters.subject_subfolder == 1
            new_parameters.output_dir = fileparts(new_parameters.output_dir);
        end
        single_subject_pipeline(subject_id, new_parameters);
    end

    disp('Pipeline finished successfully');
end
