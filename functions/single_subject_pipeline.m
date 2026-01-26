function [parameters] = single_subject_pipeline(subject_id, parameters, options)
    arguments
        subject_id 
        parameters struct
        options.adopted_heatmap (:,:,:) = []
        options.sequential_configs struct = struct()
    end

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
    % - Matlab 2022b+, SimNIBS 4.0, and k-Wave 1.4 have been tested     %
    % - 'subject_id' must be a number.                                  %
    % - 'parameters' is a structure (see default_config for options)    %
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % ====================================================================
    %% PATH & LOG setup
    % ====================================================================

    currentLoc = fileparts(mfilename("fullpath"));
    % add functions here to detect path setup function
    addpath(genpath(fullfile(currentLoc, '..', 'functions')));
    [parameters] = path_log_setup(parameters, get_prestus_path, subject_id);

    fprintf('Starting processing for subject %i %s\n',...
        parameters.subject_id, parameters.results_filename_affix)
    
    % ====================================================================
    %% SEGMENT planning image (structural MRI) with SimNIBS
    % ====================================================================
    % if segmentation is not yet available

    fprintf('========================================\n');
    fprintf('SEGMENTATION \n');
    fprintf('========================================\n\n');

    if contains(parameters.simulation_medium, {'skull'; 'layered'})
        log_timer('start','segmentation', parameters.seg_path);
        preproc_segmentation(parameters)
        log_timer('stop','segmentation');
    else
        disp('No head segmentation necessary...')
    end

    % EXIT is no simulation is requested (= segmentation only)
    if ~any([parameters.run_acoustic_sims, ...
            parameters.run_heating_sims, ...
            parameters.run_posthoc_water_sims])
        disp(newline)
        disp('No simulation requested...')
        exit;
    end

    % ====================================================================
    %% GRID: PREPROCESS structural MRI & POSITION transducer + target
    % ====================================================================
    % reorient image, determine transducer & target position in image
    % For more documentation, see the 'preproc_head' function.

    fprintf('========================================\n');
    fprintf('GRID SETUP & HEAD PREPROC \n');
    fprintf('========================================\n\n');

    log_timer('start','preproc', parameters.output_dir);

    % Focal distance calculation (if not specified)
    parameters = focal_distance_calculation(parameters);

    % Set up grid by preprocessing the planning image or reading in phantom
    [parameters, medium_masks, segmentation, skull_edge, planimg] = ...
        grid_tissue_setup(parameters);

    % Position the transducer(s) in the grid
    [parameters] = grid_transducer_location(parameters, planimg);

    % Adapt grid to axisymmetry (if requested)
    [parameters, segmentation, medium_masks] = ...
        grid_axisymmetry(parameters, segmentation, medium_masks);

    % Extract variables for quick access
    trans_pos = parameters.transducer(1).trans_pos;
    focus_pos = parameters.transducer(1).focus_pos;

    log_timer('stop','preproc');
    
    % ====================================================================
    %% SETUP MEDIUM
    % ====================================================================
    % For more documentation, see 'medium_setup'
    
    fprintf('========================================\n');
    fprintf('MEDIUM PROPERTY MAPPING \n');
    fprintf('========================================\n\n');

    log_timer('start','medium', parameters.output_dir);

    if parameters.usepseudoCT == 1
        kwave_medium = medium_setup(parameters, medium_masks, segmentation);
    else
        kwave_medium = medium_setup(parameters, medium_masks);
    end
    
    % save images of assigned medium properties
    if contains(parameters.simulation_medium, {'skull'; 'layered'}) && parameters.debug == 1
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'sound_speed')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'density')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_coeff')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_power')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'thermal_conductivity')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'specific_heat')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'perfusion_coeff')
        medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'absorption_fraction')
    end

    % split temp_0 & absorption_fraction from kwave_medium (to pass internal kwave checks)
    if isfield(parameters, 'adopted_heatmap') && parameters.adopted_heatmap == 1
        heatmap_image = niftiread(parameters.adopted_heatmap);
        medium_plus.temp_0 = double(tformarray(heatmap_image, maketform("affine", planimg.transf), ...
            makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], size(medium_masks), [], 0));
    else
        medium_plus.temp_0 = kwave_medium.temp_0;
    end
    kwave_medium = rmfield(kwave_medium, 'temp_0');
    medium_plus.absorption_fraction = kwave_medium.absorption_fraction;
    kwave_medium = rmfield(kwave_medium, 'absorption_fraction');

    log_timer('stop','medium');

    % ====================================================================
    %% SETUP SOURCE
    % ====================================================================
    % For more documentation, see 'source_sensor_setup'.
    % Precomputed source can be loaded (if available).

    fprintf('========================================\n');
    fprintf('K-WAVE SOURCE SETUP \n');
    fprintf('========================================\n\n');

    if parameters.run_source_setup
        log_timer('start','source', parameters.output_dir);

        max_sound_speed = max(kwave_medium.sound_speed(:));
        [kgrid, source, sensor, source_labels] = ...
            source_sensor_setup(...
            parameters, ...
            max_sound_speed, ...
            trans_pos, ...
            focus_pos);
        
        % Check stability & adjust source time step if necessary
        % If estimated time step is smaller than the time step based on default CFL, the estimated time step is used to redefine transducer and sensor. 
        % Note: The estimated time step does not guarantee a stable simulation. If NaN numbers result, a smaller time step than estimated may be optimal.	  
        disp('Check stability...')
        dt_stability_limit = checkStability(kgrid, kwave_medium);
        if ~isinf(dt_stability_limit) && kgrid.dt > dt_stability_limit
			disp('Adapt time step for simulation stability...')
            % Use 90% of the theoretical limit (which is only an approximation in the heterogenous medium case: http://www.k-wave.org/documentation/checkStability.php)
            grid_time_step = dt_stability_limit*0.90;
            [kgrid, source, sensor, source_labels] = source_sensor_setup(parameters, max_sound_speed, trans_pos, focus_pos, grid_time_step);
        end
        log_timer('stop', 'source');
    else
        disp('No source setup requested... no simulations will be performed.')
    end

    % ====================================================================
    %% ACOUSTIC SIMULATION
    % ====================================================================
    % See 'acoustic_simulation' for more documentation

    fprintf('========================================\n');
    fprintf('ACOUSTIC SIMULATION \n');
    fprintf('========================================\n\n');

    log_timer('start','acoustic', parameters.output_dir);

    filename_sensor_data = fullfile(parameters.output_dir, ...
        sprintf('sub-%03d_%s_results%s.mat', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    
    parameters.acoustics_available = 0;
    if parameters.run_acoustic_sims && ...
            confirm_overwriting(filename_sensor_data, parameters) && ...
            (parameters.interactive == 0 || ...
            confirmation_dlg('Running the simulations will take a long time, are you sure?', 'Yes', 'No'))

        [sensor_data] = ...
            acoustic_wrapper(...
            parameters, ...
            kgrid, ...
            kwave_medium, ...
            source, ...
            sensor, ...
            skull_edge, ...
            filename_sensor_data);

        parameters.acoustics_available = 1;
    elseif exist(filename_sensor_data, 'file')
        disp('Skipping acoustic simulation, loading existing output file.')
        load(filename_sensor_data);
        parameters.acoustics_available = 1;
    else
        parameters.acoustics_available = 0;
    end

    log_timer('stop', 'acoustic');

    % =========================================================================
    %% ACOUSTIC ANALYSIS
    % =========================================================================

    fprintf('========================================\n');
    fprintf('ACOUSTIC ANALYSIS \n');
    fprintf('========================================\n\n');

    if parameters.acoustics_available == 1
        log_timer('start','acoustic_analysis', parameters.output_dir);

        % Pass thermally relevant (but kwave-irregular) medium fields
        % these need to be loaded here to potentially transform them from axisymmetry below
        kwave_medium.temp_0 = medium_plus.temp_0;
        kwave_medium.absorption_fraction = medium_plus.absorption_fraction; 
        clear medium_plus;

        % convert media and results to 2D/3D (if axisymmetry was used)
        if numel(focus_pos) == 2 && isfield(parameters, 'axisymmetric') && parameters.axisymmetric == 1
            [sensor_data, parameters, trans_pos, focus_pos, segmentation, ...
                medium_masks, kwave_medium, kgrid, source, source_labels] = ...
                acoustic_convert_axisymmetry(...
                parameters, sensor_data, trans_pos, focus_pos, segmentation, ...
                medium_masks, kwave_medium, source, source_labels);
        end

        % perform acoustic analysis
        [results_acoustic, acoustic_isppa, acoustic_MI, acoustic_pressure, highlighted_pos] = ...
            acoustic_analysis(parameters, kwave_medium, medium_masks, ...
                        sensor_data, trans_pos, focus_pos, segmentation, ...
                        source_labels);
        log_timer('stop', 'acoustic_analysis');
    else
        disp('No acoustic simulation results available. Skipping analysis...')
    end

    % =========================================================================
    %% THERMAL SIMULATIONS
    % =========================================================================

    fprintf('========================================\n');
    fprintf('THERMAL SIMULATIONS \n');
    fprintf('========================================\n\n');

    log_timer('start','thermal', parameters.output_dir);

    parameters.heating_available = 0;
    if isfield(parameters, 'run_heating_sims') && parameters.run_heating_sims && parameters.acoustics_available == 1
        
        disp('Starting thermal simulations...')
        % Name of thermal simulation output file
        filename_heating_data = fullfile(parameters.output_dir,...
            sprintf('sub-%03d_%s_heating_res%s.mat', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        
        % Check whether thermal results axist and - if so - should be overwritten
        if confirm_overwriting(filename_heating_data, parameters) && (parameters.interactive == 0 || ...
            confirmation_dlg('Running the thermal simulations will take a long time, are you sure?', 'Yes', 'No')) 

            % Set up and run thermal simulation
            [kwaveDiffusion, ...
                time_status_seq, ...
                results_heating.maxT, ...
                results_heating.focal_planeT, ...
                results_heating.CEM43, ...
                results_heating.focal_planeCEM43] = ...
                thermal_simulation(...
                parameters, ...
                sensor_data, ...
                kgrid, ...
                kwave_medium, ...
                sensor, ...
                source, ...
                trans_pos, ...
                planimg.transf, ...
                medium_masks);

            if isfield(parameters, 'savemat') && parameters.savemat==0
                disp("Not saving thermal simulation output matrices ...")
            else
                save(filename_heating_data, ...
                    'kwaveDiffusion',...
                    'time_status_seq',...
                    'sensor',...
                    'results_heating',...
                    '-v7.3');
            end
            parameters.heating_available = 1;
        elseif exist(filename_heating_data, 'file')
            disp('Skipping thermal simulation, loading existing output file.')
            load(filename_heating_data);
            parameters.heating_available = 1;
        else 
            warning('Heating simulations requested, but no acoustic results available. Other misspecification is possible.')
            parameters.heating_available = 0;
        end
    else
        parameters.heating_available = 0;
        results_heating = [];
    end
    log_timer('stop','thermal');

    % ================================================================
    %% THERMAL ANALYSIS
    % ================================================================

    fprintf('========================================\n');
    fprintf('THERMAL ANALYSIS \n');
    fprintf('========================================\n\n');

    if parameters.heating_available == 1
        log_timer('start','thermal_analysis', parameters.output_dir);
        thermal_analysis(parameters, results_heating, time_status_seq, ...
            medium_masks, trans_pos, focus_pos, highlighted_pos, segmentation);
        log_timer('stop','thermal_analysis');
    else
        disp('No heating simulation results available. Skipping thermal analysis...')
    end

    % ================================================================
    %% CREATE NIFTI IMAGES
    % ================================================================
    % plot various metrics on both the subject-space T1 image & MNI space
    
    fprintf('========================================\n');
    fprintf('NIFTI IMAGES \n');
    fprintf('========================================\n\n');
    
    log_timer('start','nifti', parameters.output_dir);

    simulation_nifti(parameters, planimg, results_acoustic, ...
                            acoustic_isppa, acoustic_MI, acoustic_pressure, ...
                            medium_masks, results_heating, kwave_medium, trans_pos, ...
                            focus_pos, highlighted_pos)

    log_timer('stop','nifti');

    % cleanup to reduce RAM load
    clear acoustic_* heating_*

    % ====================================================================
    %% END OF THIS SIMULATIN
    % ====================================================================

    fprintf('========================================\n');
    fprintf('END \n');
    fprintf('========================================\n\n');

    % capture time, RAM, & GB load of pipeline
    log_timer('stop','single_subject_pipeline')

    % indicate success
    disp('Pipeline finished successfully');

    % end logging
    diary('off')

    % ====================================================================
    %% POST-HOC ACOUSTIC WATER SIMULATION
    % ====================================================================
    % To check sonication parameters of the transducer in free water

    if isfield(parameters, 'run_posthoc_water_sims') && parameters.run_posthoc_water_sims && ...
            contains(parameters.simulation_medium, {'layered'})

        fprintf('POST-HOC ACOUSTIC WATER SIMULATION \n');

        if numel(parameters.transducer) > 1
            warning(['Post-hoc water simulations are not implemented for multiple transducers. ' ...
                     'Post-hoc water simulation will be run only for the first specified transducer. ' ...
                     'Consider running separate configs to test individual transducers.']);
        end
        water_parameters = parameters;
        water_parameters.simulation_medium = 'water';
        water_parameters.run_heating_sims = 0;
        water_parameters.run_posthoc_water_sims = 0;
        % run with the same grid dimension as the real simulation (overkill?)
        water_parameters.default_grid_dims = water_parameters.grid_dims;
        % restore subject-specific path to original path if done earlier in this function
        if isfield(water_parameters,'subject_subfolder') && water_parameters.subject_subfolder == 1
            water_parameters.output_dir = fileparts(water_parameters.output_dir);
        end
        % inherit submit medium from main pipeline
        switch parameters.hpc_submit_medium
            case 'slurm'
                single_subject_pipeline_with_slurm(subject_id, water_parameters, 0);
            case 'matlab'
                single_subject_pipeline(parameters.subject_id, water_parameters);
        end
        clear water_parameters;
    end

    % ====================================================================
    %% FOLLOW-UP SIMULATION WITH IDENTICAL MEDIUM
    % ====================================================================

    if ~isempty(fieldnames(options.sequential_configs))
        fprintf('========================================\n');
        fprintf('FOLLOW-UP SIMULATION WITH IDENTICAL MEDIUM \n');
        fprintf('========================================\n\n');
        % Select the config next in line
        sequential_configs = options.sequential_configs;
        fields = fieldnames(sequential_configs);
        numbers = cellfun(@(x) sscanf(x, 'config_%d'), fields);
        [~, minIdx] = min(numbers);
        lowestField = fields{minIdx};
        sequential_parameters = sequential_configs.(lowestField);
        sequential_configs = rmfield(sequential_configs, lowestField);
        % restore subject-specific path to original path if done earlier in this function
        sequential_parameters.adopted_heatmap = fullfile(parameters.output_dir, sprintf('sub-%03d_final_%s_orig_coord%s',...
                parameters.subject_id, 'heating', parameters.results_filename_affix));
        sequential_parameters.adopted_cumulative_heat = fullfile(parameters.output_dir, sprintf('sub-%03d_final_%s_orig_coord%s',...
                parameters.subject_id, 'CEM43', parameters.results_filename_affix));
        adopted_heatmap = ones(2,2,2);
        fprintf('Running subsequent heating simulation on %s\n', lowestField);
        % Limitation: currently only support SLURM
        if ~isempty(fieldnames(sequential_configs))
            single_subject_pipeline_with_slurm(parameters.subject_id, sequential_parameters, false, '24:00:00', 64, 'adopted_heatmap', adopted_heatmap, 'sequential_configs', sequential_configs);
        else
            single_subject_pipeline_with_slurm(parameters.subject_id, sequential_parameters, false, '24:00:00', 64, 'adopted_heatmap', adopted_heatmap);
        end
    end

end