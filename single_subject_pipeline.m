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
        parameters.output_dir = fullfile(parameters.sim_path, sprintf('sub-%03d', subject_id));
    else 
        parameters.output_dir = parameters.sim_path;
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
    
    % deactivate pseudoCT by default
    if ~isfield(parameters, 'usepseudoCT')
        parameters.usepseudoCT = 0;
    end

    % Pre-processes MRI data to segment the different forms of tissue
    % and visualise the position of the transducer with some help from SimNIBS.
    % For more documentation, see the 'preprocess_brain' function.
    if contains(parameters.simulation_medium, 'skull') || strcmp(parameters.simulation_medium, 'layered')
        [medium_masks, segmented_image_cropped, skull_edge, trans_pos_final, ...
        focus_pos_final, t1_image_orig, t1_header, final_transformation_matrix, ...
        inv_final_transformation_matrix] = preprocess_brain(parameters, subject_id);
        if isempty(medium_masks)
            output_pressure_file = '';
            return;
        end
        parameters.grid_dims = size(medium_masks);
    else % In case simulations are not run in a skull of layered tissue, alternative grid dimensions are set up
        assert(isfield(parameters, 'default_grid_dims'), ...
            'The parameters structure should have the field grid_dims for the grid dimensions')
        parameters.grid_dims = parameters.default_grid_dims;
        if any(parameters.grid_dims==1)||length(parameters.grid_dims)==2
            parameters.grid_dims = squeeze(parameters.grid_dims);
            parameters.n_sim_dims = length(parameters.grid_dims);
            disp('One of the simulation grid dimensions is of length 1. Assuming you want 2D simulations ... dropping this dimension')
        end

        % Checks whether the transducer location and orientation are set,
        % and uses an arbitrary position if not
        medium_masks = [];
        if strcmp(parameters.simulation_medium, 'phantom')
            segmentation_folder = fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', subject_id));
            filename_segmented = fullfile(segmentation_folder, 'final_tissues.nii.gz');
            segmented_img = niftiread(filename_segmented);
            if size(segmented_img) == parameters.grid_dims
                disp('Check passed: phantom dimensions fit requested grid...');
            elseif size(segmented_img') == parameters.grid_dims
                segmented_img = segmented_img';
                disp('Check passed: phantom dimensions fit requested grid...');
            else
                warning('WARNING: phantom dimensions DO NOT fit requested grid...')
            end
            % create medium mask according to indices in parameters.layer_labels (see smooth_and_crop.m)
            [medium_masks] = medium_mask_create(segmented_img, parameters, 1);
            clear segmented_img;
        end
        segmented_image_cropped = zeros(parameters.grid_dims);
        if ~isfield(parameters.transducer, 'pos_grid') || ~isfield(parameters, 'focus_pos_grid')
            disp('Either grid or focus position is not set, positioning them arbitrarily based on the focal distance')
            % note that the focus position matters only for the orientation of the transducer
        end
        % set transducer position in grid
        if ~isfield(parameters.transducer, 'pos_grid')
            % transducer positioned arbitrarily (2D only)
            % y: first position beyond pml layer
            % x: halfway
            trans_pos_final = round(...
                [parameters.grid_dims(1:(parameters.n_sim_dims-1))/2, ...
                parameters.pml_size+1]);
        else
            trans_pos_final = parameters.transducer.pos_grid;
            % Adjust if the positions are transposed
            if size(trans_pos_final,1)>size(trans_pos_final, 2)
                warning('Specified transducer position appears transposed...adjusting');
                trans_pos_final = trans_pos_final';
            end
        end
        % set focus position in grid
        if ~isfield(parameters, 'focus_pos_grid')
            % no focus point specified
            % position focus at expected distance from transducer
            % index dimension depends on 2D/3D
            % this already accounts for PML size
            focus_pos_final = trans_pos_final;
            n_dim = numel(focus_pos_final);
            focus_pos_final(n_dim) = ...
                round(focus_pos_final(n_dim) + ...
                parameters.expected_focal_distance_mm/parameters.grid_step_mm);
            clear n_dim;
        else
            focus_pos_final = parameters.focus_pos_grid;
            % Adjust if the positions are transposed (2D only)
            if numel(focus_pos_final) == 2 && size(focus_pos_final,1)>size(focus_pos_final, 2)
                warning('Specified focus position appears transposed...adjusting');
                focus_pos_final = focus_pos_final';
            end
        end
    end
    
    % If a PML layer is used to absorb waves reaching the edge of the grid,
    % this will check if there is enough room for a PML layer between the
    % transducer and the edge of the grid
    assert(min(abs([repmat(0, 1, numel(parameters.grid_dims));parameters.grid_dims]-...
        trans_pos_final ),[],'all') > parameters.pml_size, ...
        'The minimal distance between the transducer and the simulation grid boundary should be larger than the PML size. Adjust transducer position or the PML size')
    assert(min(abs([repmat(0, 1, numel(parameters.grid_dims));parameters.grid_dims]-...
        focus_pos_final ),[],'all') > parameters.pml_size, ...
        'The minimal distance between the focus position and the simulation grid boundary should be larger than the PML size. Adjust transducer position or the PML size')
       
    % adapt grid dimensions to axisymmetry if requested
    % grid should be specified as [axial, radial x 2]
    if numel(focus_pos_final) == 2 && ...
            isfield(parameters, 'axisymmetric') && parameters.axisymmetric == 1
        % ensure that radial(y) dim is shorter than axial (x) dim
        if parameters.grid_dims(2) > parameters.grid_dims(1)
            parameters.grid_dims = fliplr(parameters.grid_dims);
            parameters.default_grid_dims = fliplr(parameters.default_grid_dims);
            trans_pos_final = fliplr(trans_pos_final);
            focus_pos_final = fliplr(focus_pos_final);
            segmented_image_cropped = segmented_image_cropped';
            medium_masks = medium_masks';
        end
        % halve the grid along the radial axis
        % see http://www.k-wave.org/documentation/kspaceFirstOrderAS.php
        Ny_half = floor(parameters.grid_dims(2)/2);
        parameters.grid_dims(2) = Ny_half;
        parameters.default_grid_dims(2) = Ny_half;
        segmented_image_cropped = segmented_image_cropped(:,Ny_half+1:end);
        medium_masks = medium_masks(:,Ny_half+1:end);
        % set transducer and focus position to the radial midline
        trans_pos_final(2) = 1; 
        focus_pos_final(2) = 1; 
    end

    % Retain transducer and focus positions after all grid manipulations
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
    
    % save images of assigned medium properties
    if (contains(parameters.simulation_medium, 'skull') || ...
            contains(parameters.simulation_medium, 'layered'))
        medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, 'sound_speed')
        medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, 'density')
        medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, 'alpha_coeff')
        medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, 'alpha_power')
        medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, 'thermal_conductivity')
        medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, 'specific_heat')
        medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, 'perfusion_coeff')
        medium_properties_nifti(parameters, kwave_medium, inv_final_transformation_matrix, t1_header, 'absorption_fraction')
    end

    % split temp_0 & absorption_fraction from kwave_medium (due to kwave checks)
    temp_0 = kwave_medium.temp_0;
    kwave_medium = rmfield(kwave_medium, 'temp_0');
    absorption_fraction = kwave_medium.absorption_fraction;
    kwave_medium = rmfield(kwave_medium, 'absorption_fraction');

    %% SETUP SOURCE
    % =========================================================================
    % For more documentation, see 'setup_grid_source_sensor'

    if parameters.run_source_setup
        disp('Setting up kwave source...')

        max_sound_speed = max(kwave_medium.sound_speed(:));
        [kgrid, source, sensor, source_labels] = ...
            setup_grid_source_sensor(...
            parameters, ...
            max_sound_speed, ...
            trans_pos_final, ...
            focus_pos_final);
        
        % check stability
        % If estimated time step is smaller than the time step based on
        % default CFL, the estimated time step is used to redefine
        % transducer and sensor. Note: the estimated time step doesn't 
        % guarantee a stable simulation. If NaN numbers result, a smaller
        % time step than estimated may be optimal.	  
        disp('Check stability...')
        dt_stability_limit = checkStability(kgrid, kwave_medium);
        if ~isinf(dt_stability_limit) && kgrid.dt > dt_stability_limit
			disp('Adapt time step for simulation stability...')
            grid_time_step = dt_stability_limit*0.90; % use 90% of the limit (which are only an approximation in the heterogenous medium case: http://www.k-wave.org/documentation/checkStability.php)
            [kgrid, source, sensor, source_labels] = setup_grid_source_sensor(parameters, max_sound_speed, trans_pos_final, focus_pos_final, grid_time_step);
        end
    end

    %% ACOUSTIC SIMULATION
    % =========================================================================
    % See 'run_simulations' for more documentation

    filename_sensor_data = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_%s_results%s.mat', ...
        subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
    parameters.acoustics_available = 0;
    if parameters.run_acoustic_sims && ...
            confirm_overwriting(filename_sensor_data, parameters) && ...
            (parameters.interactive == 0 || ...
            confirmation_dlg('Running the simulations will take a long time, are you sure?', 'Yes', 'No'))

        disp('Starting acoustic simulations...')

        % Pathname for the input and output files (used only for non-interactive computations)
        parameters.kwave_input_filename  = fullfile(parameters.output_dir, ...
            sprintf('sub-%03d_%s_input%s.h5', subject_id, ...
            parameters.simulation_medium, parameters.results_filename_affix));
        parameters.kwave_output_filename = fullfile(parameters.output_dir, ...
            sprintf('sub-%03d_%s_output%s.h5', subject_id, ...
            parameters.simulation_medium, parameters.results_filename_affix));

        % Defines the edge of the simulation as the edge of the PML layer
%         if strcmp(parameters.simulation_medium, 'phantom')
%             plminside = false; % phantom does not yet include pml
%             if isfield(parameters, 'axisymmetric') && parameters.axisymmetric == 1
%                 trans_pos_final(1) = trans_pos_final(1)-parameters.pml_size;
%                 focus_pos_final(1) = focus_pos_final(1)-parameters.pml_size;
%             else
%                 trans_pos_final(2) = trans_pos_final(2)-parameters.pml_size;
%                 focus_pos_final(2) = focus_pos_final(2)-parameters.pml_size;
%             end
%         else
            plminside = true;
%         end
        kwave_input_args = struct('PMLInside', plminside, ...
            'PMLSize', parameters.pml_size, ...
            'PlotPML', true);

        if contains(parameters.simulation_medium, 'skull')|| ...
                strcmp(parameters.simulation_medium, 'layered')
            kwave_input_args.DisplayMask = skull_edge;
        end

        if parameters.run_source_setup==0
            error('Source setup not requested. Not able to proceed with acoustic simulation.')
        end

        sensor_data = run_simulations(kgrid, kwave_medium, source, sensor, ...
            kwave_input_args, parameters);

        % re-add fields for future loading
        kwave_medium.temp_0 = temp_0; clear temp_0;
        kwave_medium.absorption_fraction = absorption_fraction; clear absorption_fraction;

        % if using axisymmetric settings, rearrange axes
        if numel(focus_pos_final) == 2 && ...
                isfield(parameters, 'axisymmetric') && parameters.axisymmetric == 1
            [sensor_data, parameters, trans_pos_final, focus_pos_final, ...
                segmented_image_cropped, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
                convert_axisymmetric_to_2d(...
                sensor_data, parameters, trans_pos_final, focus_pos_final, ...
                segmented_image_cropped, medium_masks, kwave_medium, source, source_labels);
        end

        if isfield(parameters, 'savemat') && parameters.savemat==0
            disp("Not saving acoustic output matrices ...")
        else
            save(filename_sensor_data, 'sensor_data', 'kgrid', ...
                'kwave_medium', 'source', 'sensor', 'kwave_input_args', ...
                'parameters' ,'-v7.3')
        end
        parameters.acoustics_available = 1;
    elseif exist(filename_sensor_data, 'file')
        disp('Skipping, the file already exists, loading it instead.')
        load(filename_sensor_data, 'sensor_data')
        parameters.acoustics_available = 1;
    else
        parameters.acoustics_available = 0;
    end

    %% PROCESS ACOUSTIC RESULTS
    % =========================================================================

    if parameters.acoustics_available == 1
        disp('Processing the results of acoustic simulations...')
    
        % What is the highest pressure level for every gridpoint
        acoustic_pressure = gather(sensor_data.p_max_all); % gather is used since it could be a GPU array
        max_pressure = max(acoustic_pressure(:));
    
        % Calculates the Isppa for every gridpoint
        acoustic_isppa = acoustic_pressure.^2./...
            (2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4;
        % Calculates the max Isppa
        max_Isppa = max(acoustic_isppa(:));
    
        % Calculates the Mechanical Index for every gridpoint
        acoustic_MI = (acoustic_pressure/10^6)/...
            sqrt((parameters.transducer.source_freq_hz/10^6));
    
        % Creates the foundation for a mask before the exit plane to calculate max values outside of it
        comp_grid_size = size(sensor_data.p_max_all);
        after_exit_plane_mask = ones(comp_grid_size);
        bowl_depth_grid = round((parameters.transducer.curv_radius_mm-...
            parameters.transducer.dist_to_plane_mm)/parameters.grid_step_mm);
        % Places the exit plane mask in the grid, adjusted to the amount of dimensions
        if parameters.n_sim_dims == 3
            if trans_pos_final(3) > comp_grid_size(3)/2
                after_exit_plane_mask(:,:,(trans_pos_final(parameters.n_sim_dims)-...
                    bowl_depth_grid):end) = 0;
            else
                after_exit_plane_mask(:,:,1:(trans_pos_final(parameters.n_sim_dims)+...
                    bowl_depth_grid)) = 0;
            end
        else
            if trans_pos_final(2) > comp_grid_size(2)/2
                after_exit_plane_mask(:,(trans_pos_final(parameters.n_sim_dims)-...
                    bowl_depth_grid):end) = 0;
            else
                after_exit_plane_mask(:,1:(trans_pos_final(parameters.n_sim_dims)+...
                    bowl_depth_grid)) = 0;
            end
        end
    
        % Calculates the X, Y and Z coordinates of the max. intensity
        [max_Isppa_after_exit_plane, Ix_eplane, Iy_eplane, Iz_eplane] = ...
            masked_max_3d(acoustic_isppa, after_exit_plane_mask);
        
        % Combines these coordinates into a point of max. intensity in the grid
        if parameters.n_sim_dims==3
            max_isppa_eplane_pos = [Ix_eplane, Iy_eplane, Iz_eplane];
        else 
            max_isppa_eplane_pos = [Ix_eplane, Iy_eplane];
        end
        disp('Final transducer, expected focus, and max ISPPA positions')
    
        % Calculates the average Isppa within a circle around the target
%         [trans_pos_final', focus_pos_final', max_isppa_eplane_pos']
        real_focal_distance = norm(max_isppa_eplane_pos-trans_pos_final)*parameters.grid_step_mm; % [mm]
        distance_target_real_maximum = norm(max_isppa_eplane_pos-focus_pos_final)*parameters.grid_step_mm; % [mm]
        % convert the radius from mm to voxels
        avg_radius = round(parameters.focus_area_radius/parameters.grid_step_mm); % [voxel]
        idx = arrayfun(@(d) max(1,focus_pos_final(d)-avg_radius):...
            min(size(acoustic_isppa,d),focus_pos_final(d)+avg_radius), ...
               1:ndims(acoustic_isppa), 'UniformOutput', false);
        avg_isppa_around_target = acoustic_isppa(idx{:});
        avg_isppa_around_target = mean(avg_isppa_around_target(:));
        
        % Reports the Isppa within the original stimulation target
        idx = num2cell(focus_pos_final);
        isppa_at_target = acoustic_isppa(idx{:});

        % Creates a logical skull mask and register skull_ids
        labels = fieldnames(parameters.layer_labels);
        skull_i = find(strcmp(labels, 'skull'));
        cortical_i = find(strcmp(labels, 'skull_cortical'));
        trabecular_i = find(strcmp(labels, 'skull_trabecular'));
        all_skull_ids = [skull_i, cortical_i, trabecular_i];
        mask_skull = ismember(medium_masks,all_skull_ids);
        brain_i = find(strcmp(labels, 'brain'));
        mask_brain = ismember(medium_masks,brain_i);
        skin_i = find(strcmp(labels, 'skin'));
        mask_skin = ismember(medium_masks,skin_i);
        
        % Layer-specific outcomes (in case a layered simulation)
        if contains(parameters.simulation_medium, 'skull') || ...
                strcmp(parameters.simulation_medium, 'layered') || ...
                strcmp(parameters.simulation_medium, 'phantom')

            % calculate max. isppa and location across full space
            [~, Ix, Iy, Iz] = masked_max_3d(acoustic_isppa, ones(size(medium_masks)));
            
            % extract indices in brain medium
            [max_Isppa_brain, Ix_brain, Iy_brain, Iz_brain] = ...
                masked_max_3d(acoustic_isppa, mask_brain);
            [min_Isppa_brain] = min(acoustic_isppa(mask_brain));
            half_max = acoustic_isppa >= max_Isppa_brain/2 & mask_brain;
            half_max_ISPPA_volume_brain = sum(half_max(:))*(parameters.grid_step_mm^3);
            [max_pressure_brain] = masked_max_3d(acoustic_pressure, mask_brain);
            [max_MI_brain] = masked_max_3d(acoustic_MI, mask_brain);
            
            % extract indices in skull medium
            [max_Isppa_skull] = masked_max_3d(acoustic_isppa, mask_skull);
            [max_pressure_skull] = masked_max_3d(acoustic_pressure, mask_skull);
            [max_MI_skull] = masked_max_3d(acoustic_MI, mask_skull);
            
            % extract indices in skin medium
            [max_Isppa_skin] = masked_max_3d(acoustic_isppa, mask_skin);
            [max_pressure_skin] = masked_max_3d(acoustic_pressure, mask_skin);
            [max_MI_skin] = masked_max_3d(acoustic_MI, mask_skin);
    
            % calculate real focal distance based on max isppa in whole medium
            highlighted_pos = [Ix, Iy, Iz];
            highlighted_pos = highlighted_pos(1:numel(trans_pos_final));
            real_focal_distance = norm(highlighted_pos-trans_pos_final)*parameters.grid_step_mm;
    
            writetable(table(subject_id, max_Isppa, max_Isppa_after_exit_plane, real_focal_distance, max_Isppa_skin, max_Isppa_skull, max_Isppa_brain, max_pressure_skin, max_pressure_skull, max_pressure_brain, max_MI_skin, max_MI_skull, max_MI_brain, Ix_brain, Iy_brain, Iz_brain, trans_pos_final, focus_pos_final, isppa_at_target, avg_isppa_around_target, half_max_ISPPA_volume_brain), output_pressure_file);
        else % If no layered tissue was selected, the max Isppa is highlighted on the plane and written in a table.
            max_Isppa = max(acoustic_isppa(:)); %  Does this step need to be included? already done at line 225.
            highlighted_pos = max_isppa_eplane_pos;
            writetable(table(subject_id, max_Isppa, max_Isppa_after_exit_plane, max_pressure, real_focal_distance, trans_pos_final, focus_pos_final, isppa_at_target, avg_isppa_around_target), output_pressure_file);
        end
    
        % Plot intensity on the segmented image
        
        if parameters.n_sim_dims==3
            %options.isppa_color_range = [0.5, max_Isppa_brain];
            [~,~,~,~,~,~,~,h]=plot_isppa_over_image(...
                acoustic_isppa, ...
                segmented_image_cropped, ...
                source_labels, ...
                parameters, ...
                {'y', focus_pos_final(2)}, ...
                trans_pos_final, ...
                focus_pos_final, ...
                highlighted_pos);
        else
            h = plot_isppa_over_image_2d(...
                acoustic_isppa, ...
                segmented_image_cropped, ...
                source_labels, ...
                after_exit_plane_mask, ...
                trans_pos_final, ...
                focus_pos_final, ...
                highlighted_pos);
        end
        output_plot = fullfile(parameters.output_dir,...
            sprintf('sub-%03d_%s_isppa%s.png', ...
            subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        set(h, 'InvertHardcopy', 'off'); % keep original colours
        saveas(h, output_plot, 'png')
        close(h);
    else
        disp('No acoustic simulation results available. Skipping analysis...')
    end

    %% HEATING SIMULATIONS
    % =========================================================================

    parameters.heating_available = 0;
    if isfield(parameters, 'run_heating_sims') && parameters.run_heating_sims && parameters.acoustics_available == 1
        disp('Starting heating simulations...')
        % Creates an output file to which output is written at a later stage
        filename_heating_data = fullfile(parameters.output_dir,sprintf('sub-%03d_%s_heating_res%s.mat', ...
            subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        
        % Clear sensor mask
        sensor.mask = zeros(size(sensor.mask));
        
        % Sets up parameters for the heating simulation and runs it
        if confirm_overwriting(filename_heating_data, parameters) && (parameters.interactive == 0 || ...
            confirmation_dlg('Running the thermal simulations will take a long time, are you sure?', 'Yes', 'No')) 
            
            % Set sensor along the focal axis 
            heating_window_dims = ones(2,3);
            for i = 1:2
                heating_window_dims(:,i) = [...
                    max(1, -parameters.thermal.sensor_xy_halfsize + parameters.transducer.pos_grid(i)), ...
                    min(parameters.grid_dims(i), parameters.thermal.sensor_xy_halfsize + parameters.transducer.pos_grid(i))];
            end

            % Sets up a window to perform simulations in
            if length(size(parameters.grid_dims))>2
                heating_window_dims(2,3) = parameters.grid_dims(3);
            end
            sensor.mask(heating_window_dims(1,1):heating_window_dims(2,1), heating_window_dims(1,2):heating_window_dims(2,2), :) = 1;

            % if k-plan pseudoCT setup is used, density and sound speed in bone are fixed for heating sims
            % https://dispatch.k-plan.io/static/docs/simulation-pipeline.html
            if parameters.usepseudoCT ==1 && strcmp(parameters.pseudoCT_variant, 'k-plan')
                kwave_medium.density(medium_masks==mask_skull) = 1850;
                kwave_medium.sound_speed(medium_masks==mask_skull) = ...
                    1.33*kwave_medium.density(medium_masks==mask_skull)+167; 
            end

            % convert attenuation into absorption
            kwave_medium.alpha_coeff = kwave_medium.alpha_coeff .* kwave_medium.absorption_fraction;

            % For more documentation, see 'run_heating_simulations'
            [kwaveDiffusion, ...
                time_status_seq, ...
                heating_maxT, ...
                heating_focal_planeT, ...
                heating_CEM43, ...
                heating_focal_planeCEM43] = ...
                run_heating_simulations(sensor_data, ...
                kgrid, ...
                kwave_medium, ...
                sensor, ...
                source, ...
                parameters, ...
                trans_pos_final);
            
            % apply gather in case variables are GPU arrays
            heating_maxT = gather(heating_maxT);
            heating_focal_planeT = gather(heating_focal_planeT);
            heating_CEM43 = gather(heating_CEM43);
            heating_focal_planeCEM43 = gather(heating_focal_planeCEM43);

            if isfield(parameters, 'savemat') && parameters.savemat==1
                disp("Not saving heating output matrices ...")
            else
                save(filename_heating_data, 'kwaveDiffusion','time_status_seq',...
                    'heating_window_dims','sensor','heating_maxT','heating_focal_planeT','heating_CEM43','heating_focal_planeCEM43','-v7.3');
            end
            parameters.heating_available = 1;
        elseif exist(filename_heating_data, 'file')
            disp('Skipping, the file already exists, loading it instead.')
            load(filename_heating_data);
            parameters.heating_available = 1;
        else 
            warning('Heating simulations requested, but no acoustic results available. Other misspecification is possible.')
            parameters.heating_available = 0;
        end
    end

    %% PROCESS HEATING RESULTS
    % ================================================================

    if parameters.heating_available == 1
        disp('Processing the results of heating simulations...')

        % Sets up an empty medium mask if none is specified
        if isempty(medium_masks)
            medium_masks = zeros(parameters.grid_dims);
        end

        % Creates an output table for temperature readings
        output_table = readtable(output_pressure_file);

        output_table.maxT = max(heating_maxT, [], 'all');
        output_table.maxCEM43 = max(heating_CEM43, [], 'all');
        % Overwrites the max temperature by dividing it up for each layer
        % in case a layered simulation_medium was selected
        if contains(parameters.simulation_medium, 'skull') || ...
                strcmp(parameters.simulation_medium, 'layered') || ...
                strcmp(parameters.simulation_medium, 'phantom')
            output_table.maxT_brain = masked_max_3d(heating_maxT, mask_brain);
            output_table.maxT_skull = masked_max_3d(heating_maxT, mask_skull); 
            output_table.maxT_skin = masked_max_3d(heating_maxT, mask_skin);
            output_table.riseT_brain = masked_max_3d(heating_maxT, mask_brain)-parameters.thermal.temp_0.brain;
            output_table.riseT_skull = masked_max_3d(heating_maxT, mask_skull)-parameters.thermal.temp_0.skull; 
            output_table.riseT_skin = masked_max_3d(heating_maxT, mask_skin)-parameters.thermal.temp_0.skin;
            output_table.CEM43_brain = masked_max_3d(heating_CEM43, mask_brain);
            output_table.CEM43_skull = masked_max_3d(heating_CEM43, mask_skull); 
            output_table.CEM43_skin = masked_max_3d(heating_CEM43, mask_skin);
        end
        writetable(output_table, output_pressure_file);

        % Creates a visual overlay of the transducer (if 3D T1 image is available)
        if exist('t1_header')
            grid_step_mm = t1_header.PixelDimensions(1);
            [~, source_labels] = transducer_setup(...
                parameters.transducer, ...
                trans_pos_final, ...
                focus_pos_final, ...
                size(segmented_image_cropped), ...
                grid_step_mm);
        else
            source_labels = zeros(size(segmented_image_cropped));
        end

        % Creates a line graph and a video of the heating effects
        plot_heating_sims(...
            heating_focal_planeT, ...
            time_status_seq, ...
            parameters, ...
            trans_pos_final, ...
            medium_masks, ...
            heating_focal_planeCEM43);
                
        % Plots the maximum temperature in the segmented brain
        if output_table.maxT < 38
            temp_color_range = [37, 38];
        else
            temp_color_range = [37, output_table.maxT];
        end

        clear output_table;

        % plot heating superimposed on segmentation
        if ndims(heating_maxT) == 3
            [~,~,~,~,~,~,~,h]=plot_isppa_over_image(...
                    heating_maxT, ...
                    segmented_image_cropped, ...
                    source_labels, ...
                    parameters, ...
                    {'y', focus_pos_final(2)}, ...
                    trans_pos_final, ...
                    focus_pos_final, ...
                    highlighted_pos, ...
                    'overlay_color_range', temp_color_range);
        elseif ndims(heating_maxT) == 2
            [h]=plot_isppa_over_image_2d(...
                    heating_maxT, ...
                    medium_masks, ...
                    source_labels, ...
                    after_exit_plane_mask, ...
                    trans_pos_final, ...
                    focus_pos_final, ...
                    highlighted_pos, ...
                    'overlay_color_range', temp_color_range, ...
                    'bg_bw_range', [0, numel(fieldnames(parameters.layer_labels))]);
        end
        output_plot_filename = fullfile(parameters.output_dir,...
            sprintf('sub-%03d_%s_maxT%s.png',...
            subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);
    else
        disp('No heating simulation results available. Skipping analysis...')
    end

    %% CREATE NIFTI IMAGES
    % ============================================================================
    % plot various metrics on both the subject-space T1 image and in MNI space

    if contains(parameters.simulation_medium, 'skull') || ...
            strcmp(parameters.simulation_medium, 'layered') || ...
            strcmp(parameters.simulation_medium, 'phantom')

        data_types = "medium_masks";
        if parameters.acoustics_available == 1 
            data_types  = [data_types, "isppa","MI","pressure"];
        end
        if isfield(parameters, 'run_heating_sims') && parameters.run_heating_sims 
            data_types  = [data_types, "heating", "heatrise", "CEM43"];
        end
        for data_type = data_types
            orig_file = fullfile(parameters.output_dir, sprintf('sub-%03d_final_%s_orig_coord%s',...
                subject_id, data_type, parameters.results_filename_affix));
            mni_file  = fullfile(parameters.output_dir, sprintf('sub-%03d_final_%s_MNI%s.nii.gz',...
                subject_id, data_type, parameters.results_filename_affix));

            if strcmp(data_type, "isppa")
                data = single(acoustic_isppa);
            elseif strcmp(data_type, "MI")
                data = single(acoustic_MI);
            elseif strcmp(data_type, "pressure")
                data = single(acoustic_pressure);
            elseif strcmp(data_type, "medium_masks")
                data = medium_masks;
            elseif strcmp(data_type, "heating")
                data = single(heating_maxT);
            elseif strcmp(data_type, "heatrise")
                data = single(heating_maxT-kwave_medium.temp_0);
            elseif strcmp(data_type, "CEM43")
                data = single(heating_CEM43);
            end
            orig_file_with_ext = strcat(orig_file, '.nii.gz');
    
            if confirm_overwriting(orig_file_with_ext, parameters)
                if ~strcmp(parameters.simulation_medium, 'phantom')
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
                        
                        orig_hdr.Datatype = 'single';
                    end
                    niftiwrite(data_backtransf, orig_file, orig_hdr, 'Compressed', true);
                else
                    niftiwrite(data, orig_file, 'Compressed', true);
                end
            else 
                data_backtransf = niftiread(orig_file_with_ext);
            end

            if strcmp(data_type, "isppa") && ~strcmp(parameters.simulation_medium, 'phantom')
                % Creates a visual overlay of the transducer
                backtransf_coordinates = round(tformfwd([trans_pos_final;  focus_pos_final; highlighted_pos], inv_final_transformation_matrix));
                [~, source_labels] = transducer_setup(parameters.transducer, backtransf_coordinates(1,:), backtransf_coordinates(2,:), ...
                                                            size(t1_image_orig), t1_header.PixelDimensions(1));
                % Plots the Isppa over the untransformed image
                backtransf_coordinates = round(tformfwd([trans_pos_final;  focus_pos_final; highlighted_pos], inv_final_transformation_matrix));
                [~,~,~,~,~,~,~,h]=plot_isppa_over_image(...
                    data_backtransf, ...
                    t1_image_orig, ...
                    source_labels, ...
                    parameters, ...
                    {'y', backtransf_coordinates(2,2)}, ...
                    backtransf_coordinates(1,:), ...
                    backtransf_coordinates(2,:), ...
                    backtransf_coordinates(3,:), ...
                    'show_rectangles', 0, ...
                    'grid_step', t1_header.PixelDimensions(1), ...
                    'overlay_color_range', [0 max_Isppa_brain], ...
                    'overlay_threshold_low', 0, ...
                    'overlay_threshold_high', max_Isppa_brain, ...
                    'rotation', 0); % rotation = 90 not implemented for transducer overlay
        
                output_plot_filename = fullfile(parameters.output_dir,...
                    sprintf('sub-%03d_%s_isppa_t1%s.png', ...
                    subject_id, parameters.simulation_medium, ...
                    parameters.results_filename_affix));
                saveas(h, output_plot_filename, 'png')
                close(h);
            end
            
            m2m_folder= fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', subject_id));
            
            % transform outputs to MNI space (using SimNibs or applying transformation matrix)
            if ~confirm_overwriting(mni_file, parameters) || strcmp(parameters.simulation_medium, 'phantom')
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

        % cleanup to reduce RAM load
        clear data acoustic_* heating_* mask_* temp_0
        
        % Since charm does not transform the T1 into MNI space, one is manually created here
        if strcmp(parameters.segmentation_software, 'charm') && ~strcmp(parameters.simulation_medium, 'phantom')
            path_to_input_img = fullfile(m2m_folder,'T1.nii.gz');
            path_to_output_img = fullfile(m2m_folder,'toMNI','T1_to_MNI_post-hoc.nii.gz');

            if ~exist(path_to_output_img,'file')
                convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters)
            end
        end
    end
    
    %% Post-hoc acoustic simulations in water
    % ============================================================================
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
