function [parameters] = prestus_pipeline(parameters, options)
% PRESTUS_PIPELINE  End-to-end transcranial focused ultrasound simulation pipeline
%
% Executes the full PRESTUS workflow for a single subject and transducer
% configuration, from head model construction through acoustic and thermal
% simulation to output report generation.
%
% Pipeline stages (each independently toggled via parameters.modules.*):
%   1. Segmentation      SimNIBS charm segmentation of structural MRI
%   2. Head processing   Reorientation, rescaling, smoothing, and cropping
%                        of the head model to the simulation grid
%   3. Source setup      Construction of the transducer source array
%                        (kWaveArray or analytical; annular or matrix)
%   4. Acoustic sim      k-Wave pressure field simulation
%   5. Acoustic analysis ISPPA, ISPTA, MI, and focal metrics
%   6. Thermal sim       Bio-heat equation solved over the sonication protocol
%   7. Thermal analysis  Peak temperature, CEM43, and safety metrics
%   8. NIfTI export      Pressure and temperature maps written as NIfTI volumes
%   9. Free-water sim    Repeat acoustic simulation in homogeneous water medium
%  10. Report            Self-contained HTML summary report
%
% Use as:
%   parameters = prestus_pipeline(parameters)
%   parameters = prestus_pipeline(parameters, options)
%
%   Typically invoked via prestus_pipeline_start, which handles platform
%   selection (MATLAB / SLURM / qsub) and uncertainty mode dispatch.
%
% Input:
%   parameters - struct built from default_config.yaml and a subject-specific
%                override config; must include subject_id, path.sim, and all
%                transducer/grid settings; see configs/default_config.yaml for
%                the full list of available fields
%   options    - (optional) struct; currently supports:
%                  .sequential_configs  cell array of additional YAML config
%                                       paths merged in sequence before running
%
% Output:
%   parameters - updated struct with fields added by path_log_setup and
%                intermediate processing steps (e.g. io.output_dir,
%                io.kwave_source_filename, transducer positions)
%
% Dependencies: SimNIBS 4, k-Wave >= 1.4.1, MATLAB >= 2023b
%
% See also: PRESTUS_PIPELINE_START, LOAD_PARAMETERS, UNCERTAINTY_PIPELINE

    arguments
        parameters struct
        options struct = struct()
    end
    
    % ====================================================================
    %% PATH & LOG setup
    % ====================================================================

    currentLoc = fileparts(mfilename("fullpath"));
    safe_addpath(fullfile(currentLoc, '..', 'functions'));
    [parameters] = path_log_setup(parameters, get_prestus_path);

    % ====================================================================
    %% TELEMETRY (opt-in, anonymous)
    % ====================================================================
    telemetry_setup();                    % no-op after first run
    track_usage('run_start', parameters);
    telemetry_t0 = tic;
    try

    % ====================================================================
    %% UNCERTAINTY MODE (MATLAB platform only)
    % ====================================================================
    if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'uncertainty') && ...
            parameters.simulation.uncertainty
        uncertainty_pipeline(parameters, options);
        return;
    end

    % ====================================================================
    %% MULTI-ISPPA MODE (MATLAB platform only)
    % ====================================================================
    if is_multi_isppa_mode(parameters)
        multi_isppa_pipeline(parameters, options);
        return;
    end

    fprintf('Starting processing for subject %i %s\n',...
        parameters.subject_id, parameters.io.output_affix)
    
    % ====================================================================
    %% SEGMENT planning image (structural MRI) with SimNIBS
    % ====================================================================

    fprintf('========================================\n');
    fprintf('SEGMENTATION \n');
    fprintf('========================================\n\n');
    log_timer('start','segmentation', parameters.path.seg);

    if contains(parameters.simulation.medium, {'layered'})
        parameters = preproc_segmentation(parameters);
    else
        disp('No head segmentation necessary...')
    end
    log_timer('stop','segmentation');

    if isfield(parameters.modules, 'segmentation_only') && parameters.modules.segmentation_only
        fprintf('Only segmentation requested: finishing.\n');
        log_timer('stop','prestus_pipeline')
        diary('off')
        return;
    end

    % ====================================================================
    %% TRANSDUCER PLACEMENT
    % ====================================================================

    fprintf('========================================\n');
    fprintf('TRANSDUCER PLACEMENT \n');
    fprintf('========================================\n\n');
    log_timer('start','placement', parameters.io.output_dir);

    if ~isfield(parameters.modules, 'run_transducer_placement') || ...
            parameters.modules.run_transducer_placement == 1
        parameters = preproc_transducer_placement(parameters);
    else
        disp('Transducer placement stage skipped.');
    end
    log_timer('stop','placement');

    % ====================================================================
    %% GRID: PREPROCESS structural MRI & POSITION transducer + target
    % ====================================================================

    fprintf('========================================\n');
    fprintf('GRID SETUP & HEAD PREPROC \n');
    fprintf('========================================\n\n');
    log_timer('start','preproc', parameters.io.output_dir);

    if isfield(parameters.io, 'preproc_affix')
        preproc_affix = parameters.io.preproc_affix;
    else
        preproc_affix = parameters.io.output_affix;
    end
    filename_grid_cache = fullfile(parameters.io.cache_dir, ...
        sprintf('sub-%03d_%s%s_grid_cache.mat', ...
        parameters.subject_id, parameters.simulation.medium, preproc_affix));

    if ~isfield(parameters.modules, 'run_grid_setup') || parameters.modules.run_grid_setup == 1
        if ~confirm_overwriting(filename_grid_cache, parameters)
            disp('Loading grid cache — skipping head preprocessing.')
            load(filename_grid_cache);
            parameters.grid       = grid_cache_info.parameters.grid;
            parameters.transducer = grid_cache_info.parameters.transducer;
            clear grid_cache_info
            % Re-derive computed transducer fields absent in older cache versions.
            parameters = load_transducer_parameters(parameters);
            trans_pos = parameters.transducer(1).trans_pos;
            focus_pos = parameters.transducer(1).focus_pos;
        else
            parameters = focal_distance_calculation(parameters);

            [parameters, medium_masks, segmentation, bone, planimg] = ...
                grid_tissue_setup(parameters);

            [parameters] = grid_transducer_location(parameters, planimg);

            [parameters, segmentation, bone, medium_masks] = ...
                grid_axisymmetry(parameters, segmentation, bone, medium_masks);

            trans_pos = parameters.transducer(1).trans_pos;
            focus_pos = parameters.transducer(1).focus_pos;

            if should_save_output(parameters.io, 'save_grid_cache')
                grid_cache_info.parameters = parameters;
                save(filename_grid_cache, ...
                    'planimg', 'medium_masks', 'segmentation', 'bone', ...
                    'trans_pos', 'focus_pos', 'grid_cache_info', '-v7.3');
                clear grid_cache_info
            end
        end
    else
        disp('No grid setup requested...no simulations will be performed.')
    end
    log_timer('stop','preproc');

    % ====================================================================
    %% SETUP MEDIUM
    % ====================================================================

    fprintf('========================================\n');
    fprintf('MEDIUM PROPERTY MAPPING \n');
    fprintf('========================================\n\n');
    log_timer('start','medium', parameters.io.output_dir);

    filename_medium_cache = fullfile(parameters.io.cache_dir, ...
        sprintf('sub-%03d_%s%s_medium_cache.mat', ...
        parameters.subject_id, parameters.simulation.medium, preproc_affix));

    if ~isfield(parameters.modules, 'run_medium_setup') || parameters.modules.run_medium_setup == 1
        if ~confirm_overwriting(filename_medium_cache, parameters)
            disp('Loading medium cache — skipping medium property mapping.')
            load(filename_medium_cache);
        else
            if parameters.pct.enabled == 1
                [kwave_medium, medium_plus] = medium_setup(parameters, medium_masks, planimg, bone);
            else
                [kwave_medium, medium_plus] = medium_setup(parameters, medium_masks, planimg);
            end

            if should_save_output(parameters.io, 'save_medium_cache')
                save(filename_medium_cache, 'kwave_medium', 'medium_plus', '-v7.3');
            end
        end
    else
        disp('No medium mapping requested...no simulations will be performed.')
    end
    log_timer('stop','medium');

    % ====================================================================
    %% SETUP SOURCE
    % ====================================================================

    fprintf('========================================\n');
    fprintf('K-WAVE SOURCE SETUP \n');
    fprintf('========================================\n\n');
    log_timer('start','source', parameters.io.output_dir);

    if ~isfield(parameters.modules, 'run_source_setup') || parameters.modules.run_source_setup==1
        max_sound_speed = max(kwave_medium.sound_speed(:));
        min_sound_speed = min(kwave_medium.sound_speed(:));
        [kgrid, source, sensor, source_labels] = ...
            source_sensor_setup(...
            parameters, ...
            max_sound_speed, ...
            trans_pos, ...
            focus_pos, ...
            [], ...
            min_sound_speed);

        if isfield(parameters.grid, 'source_limit_fraction') && parameters.grid.source_limit_fraction ~=0
            disp('Check stability...')
            dt_stability_limit = checkStability(kgrid, kwave_medium);
            fprintf('Stability limit estimate for time step: %.1d.\n', dt_stability_limit);
            if ~isinf(dt_stability_limit) && kgrid.dt > dt_stability_limit
			    disp('Adapt time step for simulation stability...')
                clear source sensor source_labels
                % checkStability is only an approximation for heterogeneous media
                grid_time_step = dt_stability_limit*parameters.grid.source_limit_fraction;
                [kgrid, source, sensor, source_labels] = source_sensor_setup(parameters, max_sound_speed, trans_pos, focus_pos, grid_time_step, min_sound_speed);
            end
        end
    else
        disp('No source setup requested...no simulations will be performed.')
    end
    log_timer('stop', 'source');

    % ====================================================================
    %% FREE-WATER BASELINE (provenance for post-hoc ISPPA scaling)
    % ====================================================================

    acoustic_provenance = struct();
    run_baseline = (isfield(parameters.modules, 'run_water_baseline') && ...
                    parameters.modules.run_water_baseline == 1) || ...
                   (isfield(parameters, 'calibration') && ...
                    isfield(parameters.calibration, 'target_isppa_wcm2') && ...
                    ~isempty(parameters.calibration.target_isppa_wcm2));
    if run_baseline && ...
            (~isfield(parameters.modules, 'run_source_setup') || parameters.modules.run_source_setup == 1)
        log_timer('start', 'freefield_baseline', parameters.io.output_dir);
        acoustic_provenance.freefield_isppa_wcm2 = ...
            water_baseline(parameters, kgrid, source, sensor);
        log_timer('stop', 'freefield_baseline');
    end

    % ====================================================================
    %% ACOUSTIC SIMULATION
    % ====================================================================

    fprintf('========================================\n');
    fprintf('ACOUSTIC SIMULATION \n');
    fprintf('========================================\n\n');
    log_timer('start','acoustic', parameters.io.output_dir);

    % Thermal jobs use a per-target output_affix; acoustic_cache_affix points back to the base run.
    if isfield(parameters.io, 'acoustic_cache_affix')
        acoustic_file_affix = parameters.io.acoustic_cache_affix;
    else
        acoustic_file_affix = parameters.io.output_affix;
    end
    filename_sensor_data = fullfile(parameters.io.cache_dir, ...
        sprintf('sub-%03d_%s%s_results.mat',...
        parameters.subject_id, parameters.simulation.medium, acoustic_file_affix));

    parameters.state.acoustics_available = 0;
    if isfield(parameters.modules, 'run_acoustic_sims') && parameters.modules.run_acoustic_sims &&...
        confirm_overwriting(filename_sensor_data, parameters) && ...
        (parameters.simulation.interactive == 0 || ...
        confirmation_dlg('Running the simulations will take a long time, are you sure?', 'Yes', 'No'))

        [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
            acoustic_wrapper(...
            parameters, ...
            kgrid, ...
            kwave_medium, ...
            source, ...
            sensor, ...
            medium_masks, ...
            filename_sensor_data, ...
            segmentation, ...
            source_labels, ...
            acoustic_provenance);

        parameters.state.acoustics_available = 1;

    elseif exist(filename_sensor_data, 'file')
        disp('Skipping acoustic simulation, loading existing output file.')
        load(filename_sensor_data);
        parameters.state.acoustics_available = 1;
        sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters);
    else
        disp('No acoustic simulation available or requested ... skipping analysis')
        parameters.state.acoustics_available = 0;
        parameters.modules.run_acoustic_analysis = 0;
    end
    log_timer('stop', 'acoustic');

    % =========================================================================
    %% ACOUSTIC ANALYSIS
    % =========================================================================

    fprintf('========================================\n');
    fprintf('ACOUSTIC ANALYSIS \n');
    fprintf('========================================\n\n');
    log_timer('start','acoustic_analysis', parameters.io.output_dir);

    if (~isfield(parameters.modules, 'run_acoustic_analysis') || parameters.modules.run_acoustic_analysis)
        [results_acoustic, acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos] = ...
            acoustic_analysis(parameters, kwave_medium, medium_masks, sensor_data, segmentation, source_labels);
    else
        disp('No acoustic simulation results available (or requested). Skipping analysis...')
        results_acoustic = [];
        acoustic_Ipa = [];
        acoustic_MI = [];
        acoustic_pressure = [];
        highlighted_pos = [];
    end
    log_timer('stop', 'acoustic_analysis');

    clear source_labels

    % =========================================================================
    %% THERMAL SIMULATIONS
    % =========================================================================

    fprintf('========================================\n');
    fprintf('THERMAL SIMULATIONS \n');
    fprintf('========================================\n\n');
    log_timer('start','thermal', parameters.io.output_dir);

    parameters.state.heating_available = 0;
    if isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims && parameters.state.acoustics_available == 1

        disp('Starting thermal simulations...')
        filename_heating_data = fullfile(parameters.io.cache_dir,...
            sprintf('sub-%03d_%s%s_heating_res.mat',...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
        
        if confirm_overwriting(filename_heating_data, parameters) && (parameters.simulation.interactive == 0 || ...
            confirmation_dlg('Running the thermal simulations will take a long time, are you sure?', 'Yes', 'No'))

            kwave_medium.temp_0              = medium_plus.temp_0;
            kwave_medium.absorption_fraction = medium_plus.absorption_fraction;
            clear medium_plus;
            [kwaveDiffusion, ...
                time_status_seq, ...
                results_heating.maxT, ...
                results_heating.focal_planeT, ...
                results_heating.heating_endT, ...
                results_heating.CEM43, ...
                results_heating.focal_planeCEM43, ...
                results_heating.CEM43_end, ...
                results_heating.timeseries, ...
                results_heating.CEM43_iso, ...
                results_heating.focal_planeCEM43_iso, ...
                results_heating.CEM43_iso_end] = ...
                thermal_simulation(...
                parameters, ...
                sensor_data, ...
                kgrid, ...
                kwave_medium, ...
                sensor, ...
                source, ...
                planimg.transf, ...
                medium_masks);

            if ~should_save_output(parameters.io, 'save_thermal_matrices')
                disp('Not saving thermal simulation output matrices ...')
            else
                save(filename_heating_data, ...
                    'kwaveDiffusion',...
                    'time_status_seq',...
                    'sensor',...
                    'results_heating',...
                    'kwave_medium', ...
                    '-v7.3');
            end

            clear sensor_data source sensor kgrid kwaveDiffusion

            parameters.state.heating_available = 1;
        elseif exist(filename_heating_data, 'file')
            disp('Skipping thermal simulation, loading existing output file.')
            load(filename_heating_data);
            parameters.state.heating_available = 1;
        else
            warning('Heating simulations requested, but no acoustic results available. Other misspecification is possible.')
            parameters.state.heating_available = 0;
            parameters.modules.run_thermal_analysis = 0;
            results_heating = struct();
        end
    else
        parameters.state.heating_available = 0;
        results_heating = struct();
    end
    log_timer('stop','thermal');

    % ================================================================
    %% THERMAL ANALYSIS
    % ================================================================

    fprintf('========================================\n');
    fprintf('THERMAL ANALYSIS \n');
    fprintf('========================================\n\n');
    log_timer('start','thermal_analysis', parameters.io.output_dir);

    if parameters.state.heating_available == 1 && ...
            (~isfield(parameters.modules, 'run_thermal_analysis') || parameters.modules.run_thermal_analysis)
        thermal_analysis(parameters, results_heating, time_status_seq, ...
            medium_masks, highlighted_pos, segmentation);
    else
        disp('No heating simulation results available (or requested). Skipping thermal analysis...')
    end
    log_timer('stop','thermal_analysis');

    clear time_status_seq segmentation bone

    % ================================================================
    %% CREATE NIFTI IMAGES
    % ================================================================
    
    fprintf('========================================\n');
    fprintf('NIFTI IMAGES \n');
    fprintf('========================================\n\n');
    log_timer('start','nifti', parameters.io.output_dir);

    if ~isfield(parameters.modules, 'run_nifti_creation') || parameters.modules.run_nifti_creation==1
        simulation_nifti(parameters, planimg, results_acoustic, ...
                                acoustic_Ipa, acoustic_MI, acoustic_pressure, ...
                                medium_masks, results_heating, kwave_medium, highlighted_pos)
    else
        disp('No nifti creation requested...')
    end

    log_timer('stop','nifti');

    clear acoustic_* results_heating medium_masks kwave_medium planimg

    % ====================================================================
    %% END OF THIS SIMULATION
    % ====================================================================

    fprintf('========================================\n');
    fprintf('END \n');
    fprintf('========================================\n\n');

    log_timer('stop','prestus_pipeline')
    disp('Pipeline finished successfully');

    if isfield(parameters.modules, 'generate_report') && parameters.modules.generate_report
        generate_simulation_report(parameters);
    end

    if isfield(parameters.modules, 'multi_isppa_report') && parameters.modules.multi_isppa_report
        if isfield(parameters, 'multi_isppa')
            generate_multi_isppa_report(parameters, parameters.multi_isppa);
        else
            warning('multi_isppa_report flag set but parameters.multi_isppa is missing — skipping.');
        end
    end

    if isfield(parameters.modules, 'uncertainty_report') && parameters.modules.uncertainty_report
        if isfield(parameters, 'uncertainty') && isfield(parameters.uncertainty, 'affixes')
            generate_uncertainty_report(parameters, parameters.uncertainty.affixes);
            if ~should_save_output(parameters.io, 'save_matrices')
                cleanup_uncertainty_intermediates(parameters, parameters.uncertainty.affixes);
            end
        else
            generate_uncertainty_report(parameters);
        end
    end

    diary('off')

    % ====================================================================
    %% POST-HOC ACOUSTIC WATER SIMULATION
    % ====================================================================

    if isfield(parameters.modules, 'run_posthoc_water_sims') && parameters.modules.run_posthoc_water_sims && ...
            contains(parameters.simulation.medium, {'layered', 'phantom'})

        fprintf('POST-HOC ACOUSTIC WATER SIMULATION \n');

        if numel(parameters.transducer) > 1
            warning(['Post-hoc water simulations are not implemented for multiple transducers. ' ...
                     'Post-hoc water simulation will be run only for the first specified transducer. ' ...
                     'Consider running separate configs to test individual transducers.']);
        end
        water_parameters = parameters;
        water_parameters.simulation.medium = 'water';
        water_parameters.modules.run_heating_sims = 0;
        water_parameters.modules.run_posthoc_water_sims = 0;
        water_parameters.simulation.debug = 0;
        water_parameters.grid.default_dims = water_parameters.grid.dims;
        water_parameters.hpc.timelimit = '05:00:00';
        water_parameters.hpc.memorylimit = 40;
        water_parameters.hpc.wait_for_job = false;
        prestus_pipeline_start(water_parameters);
        clear water_parameters;
    end

    % ====================================================================
    %% FOLLOW-UP SIMULATION WITH IDENTICAL MEDIUM
    % ====================================================================

    if any(strcmp(fieldnames(options), 'sequential_configs'))
        fprintf('========================================\n');
        fprintf('FOLLOW-UP SIMULATION WITH IDENTICAL MEDIUM \n');
        fprintf('========================================\n\n');
        sequential_configs = options.sequential_configs;
        fields = fieldnames(sequential_configs);
        numbers = cellfun(@(x) sscanf(x, 'config_%d'), fields);
        [~, minIdx] = min(numbers);
        lowestField = fields{minIdx};
        sequential_parameters = sequential_configs.(lowestField);
        sequential_configs = rmfield(sequential_configs, lowestField);
        sequential_parameters.io.adopted_heatmap = fullfile(parameters.io.output_dir, sprintf('sub-%03d_final_%s_orig_coord%s',...
                parameters.subject_id, 'heating_end', parameters.io.output_affix));
        sequential_parameters.io.adopted_cem43 = fullfile(parameters.io.output_dir, sprintf('sub-%03d_final_%s_orig_coord%s',...
                parameters.subject_id, 'CEM43_end', parameters.io.output_affix));
        fprintf('Running subsequent heating simulation on %s\n', lowestField);
        if ~isempty(fieldnames(sequential_configs))
            options.sequential_configs = sequential_configs;
        else
            options = rmfield(options, 'sequential_configs');
        end
        prestus_pipeline_start(sequential_parameters, options)
    end

    track_usage('run_end', parameters, struct('duration_s', toc(telemetry_t0), 'status', 'success'));

    catch me_
        track_usage('run_error', parameters, struct( ...
            'duration_s', toc(telemetry_t0), ...
            'status',     'error', ...
            'error_id',   me_.identifier));
        rethrow(me_);
    end

end
