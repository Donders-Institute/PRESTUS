function [parameters] = prestus_pipeline(parameters, options)
% PRESTUS_PIPELINE  End-to-end transcranial focused ultrasound simulation pipeline
%
% Executes the full PRESTUS workflow for a single subject and transducer
% configuration, from head model construction through acoustic and thermal
% simulation to output report generation.
%
% Pipeline stages (each independently toggled via parameters.modules.*):
%   1. Segmentation      SimNIBS charm segmentation of structural MRI
%   2. Transducer placement  Neuronavigation-guided entry/target point setup
%   3. Grid & head model Reorientation, rescaling, smoothing, and cropping
%                        of the segmented head to the simulation grid
%   4. Medium setup      Assignment of acoustic/thermal tissue properties
%                        (density, sound speed, attenuation, heat capacity …)
%   5. Source setup      Construction of the transducer source array
%                        (kWaveArray or analytical; annular or matrix)
%   6. Acoustic sim      k-Wave pressure field simulation
%   7. Acoustic analysis ISPPA, ISPTA, MI, and focal metrics; NIfTI export
%   8. Thermal sim       Bio-heat equation solved over the sonication protocol
%   9. Thermal analysis  Peak temperature, CEM43, and safety metrics; NIfTI export
%  10. Report            Self-contained HTML summary report
%  11. Post-hoc water sim  Optional repeat of acoustic sim in homogeneous water
%                          (for free-field reference without skull)
%  12. Sequential follow-up  Optional chained simulation using end-temperature
%                             of this run as the starting state for the next
%
% Cache files written at each stage allow individual stages to be skipped on
% re-runs.  Stage-level toggles live in parameters.modules.*; see
% config/config_default.yaml for all available flags.
%
% Use as:
%   parameters = prestus_pipeline(parameters)
%   parameters = prestus_pipeline(parameters, options)
%
%   Typically invoked via prestus_pipeline_start, which handles platform
%   selection (MATLAB / SLURM / qsub) and special-mode dispatch.
%
% Input:
%   parameters - struct built from config_default.yaml and a subject-specific
%                override config; must include subject_id, path.sim, and all
%                transducer/grid settings; see config/config_default.yaml for
%                the full list of available fields
%   options    - (optional) struct; currently supports:
%                  .sequential_configs  struct of follow-up configs to run in
%                                       sequence after this simulation, each
%                                       starting from the current run's
%                                       end-temperature and CEM43 state
%
% Output:
%   parameters - updated struct with fields added by path_log_setup and
%                intermediate processing steps (e.g. io.dir_output,
%                io.filename_kwave_source, transducer positions)
%
% Dependencies: SimNIBS 4, k-Wave >= 1.4.1, MATLAB >= 2023b
%
% See also: PRESTUS_PIPELINE_START, LOAD_PARAMETERS, UNCERTAINTY_PIPELINE,
%           SEQUENTIAL_PIPELINE

    arguments
        parameters struct
        options struct = struct()
    end

    % ====================================================================
    %% PATH & LOG SETUP
    % ====================================================================

    currentLoc = fileparts(mfilename("fullpath"));
    safe_addpath(fullfile(currentLoc, '..', 'functions'));
    [parameters] = path_log_setup(parameters, get_prestus_path);

    % ====================================================================
    %% TELEMETRY  (opt-in, fully anonymous)
    % ====================================================================

    telemetry_setup();   % reads/writes consent flag on first call, no-op after
    telemetry_t0  = tic;
    telemetry_rid = generate_run_id();
    try

    % ====================================================================
    %% SPECIAL-MODE INTERCEPTS
    %
    % Several orchestration modes reuse prestus_pipeline as their entry
    % point but redirect execution to a dedicated pipeline function before
    % any simulation work begins.  Each returns immediately after dispatch
    % so the rest of this function is not reached.
    % ====================================================================

    % Uncertainty mode: runs default / liberal / conservative tissue-property
    % variants and combines them into a single uncertainty HTML report.
    if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'uncertainty') && ...
            parameters.simulation.uncertainty
        uncertainty_pipeline(parameters, options);
        return;
    end

    % Multi-ISPPA mode: sweeps a list of target intensities and compares
    % focal metrics across the sweep without re-running the full pipeline
    % for each point.
    if is_multi_isppa_mode(parameters)
        multi_isppa_pipeline(parameters, options);
        return;
    end

    % Async-combine stage: a lightweight HPC job that merges intensity
    % volumes from independently simulated transducers.  No k-Wave is
    % invoked — the job only calls combine_async_intensity and exits.
    if isfield(parameters, 'modules') && isfield(parameters.modules, 'combine_async') && ...
            parameters.modules.combine_async
        targets = [];
        if isfield(parameters.async_combine, 'targets_wcm2')
            targets = parameters.async_combine.targets_wcm2;
        end
        combine_async_intensity( ...
            parameters.async_combine.files_in, ...
            parameters.async_combine.file_out, ...
            targets);
        return;
    end

    fprintf('Starting processing for subject %i %s\n', ...
        parameters.subject_id, parameters.io.output_affix)

    % ====================================================================
    %% STAGE 1 — SEGMENTATION
    %
    % Runs SimNIBS charm to segment the structural MRI into tissue classes
    % (skin, skull, CSF, brain, …).  Only needed for layered simulations;
    % water / phantom runs skip this stage entirely.
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
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
    %% STAGE 2 — TRANSDUCER PLACEMENT
    %
    % Resolves entry and target coordinates from neuronavigation data
    % (e.g. Localite) or from coordinates provided directly in the config.
    % Writes a placement QC figure to the output folder.
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
    fprintf('========================================\n\n');
    log_timer('start','placement', parameters.io.dir_output);

    if ~isfield(parameters.modules, 'run_transducer_placement') || ...
            parameters.modules.run_transducer_placement == 1
        parameters = preproc_transducer_placement(parameters);
    else
        disp('Transducer placement stage skipped.');
    end
    log_timer('stop','placement');

    % ====================================================================
    %% STAGE 3 — GRID SETUP & HEAD PREPROCESSING
    %
    % Reorients and crops the segmented head, places the transducer and
    % target within the k-Wave grid, and optionally applies axisymmetric
    % reduction for 2-D simulations.  Results are cached so that
    % re-runs with the same config reuse the grid without recomputing.
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
    fprintf('========================================\n\n');
    log_timer('start','preproc', parameters.io.dir_output);

    % sequential_pipeline sets preproc_affix to the base run's output_affix so
    % follow-up simulations reuse the existing grid without recomputing it.
    if isfield(parameters.io, 'preproc_affix')
        preproc_affix = parameters.io.preproc_affix;
    else
        preproc_affix = parameters.io.output_affix;
    end
    filename_grid_cache = fullfile(parameters.io.dir_cache, ...
        sprintf('sub-%03d_%s_grid_cache%s.mat', ...
        parameters.subject_id, parameters.simulation.medium, preproc_affix));

    if ~isfield(parameters.modules, 'run_grid_setup') || parameters.modules.run_grid_setup == 1
        if ~confirm_overwriting(filename_grid_cache, parameters)
            disp('Loading grid cache - skipping head preprocessing.')
            load(filename_grid_cache);
            parameters.grid       = grid_cache_info.parameters.grid;
            parameters.transducer = grid_cache_info.parameters.transducer;
            clear grid_cache_info
            % Cache was written by an older PRESTUS version that did not store
            % all derived transducer fields; recompute them from the loaded config.
            parameters = load_transducer_parameters(parameters);
            trans_pos = parameters.transducer(1).trans_pos;
            focus_pos = parameters.transducer(1).focus_pos;
        else
            parameters = load_transducer_parameters(parameters);

            [parameters, medium_masks, segmentation, bone_mask, pseudoCT, planimg] = ...
                grid_tissue_setup(parameters);

            [parameters] = grid_transducer_location(parameters, planimg);

            [parameters, segmentation, bone_mask, pseudoCT, medium_masks] = ...
                grid_axisymmetry(parameters, segmentation, bone_mask, pseudoCT, medium_masks);

            trans_pos = parameters.transducer(1).trans_pos;
            focus_pos = parameters.transducer(1).focus_pos;

            if should_save_output(parameters.io, 'save_grid_cache')
                grid_cache_info.parameters = parameters;
                save(filename_grid_cache, ...
                    'planimg', 'medium_masks', 'segmentation', 'bone_mask', 'pseudoCT', ...
                    'trans_pos', 'focus_pos', 'grid_cache_info', '-v7.3');
                clear grid_cache_info
            end
        end
    else
        disp('No grid setup requested...no simulations will be performed.')
    end
    log_timer('stop','preproc');

    % ====================================================================
    %% STAGE 4 — MEDIUM PROPERTY MAPPING
    %
    % Maps tissue-class labels from the segmentation onto acoustic and
    % thermal material properties (density, sound speed, attenuation,
    % heat capacity, perfusion rate, …) using the lookup tables in the
    % config.  Produces kwave_medium and medium_plus structs consumed by
    % the acoustic and thermal stages respectively.
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
    fprintf('========================================\n\n');
    log_timer('start','medium', parameters.io.dir_output);

    filename_medium_cache = fullfile(parameters.io.dir_cache, ...
        sprintf('sub-%03d_%s_medium_cache%s.mat', ...
        parameters.subject_id, parameters.simulation.medium, preproc_affix));

    if ~isfield(parameters.modules, 'run_medium_setup') || parameters.modules.run_medium_setup == 1
        if ~confirm_overwriting(filename_medium_cache, parameters)
            disp('Loading medium cache - skipping medium property mapping.')
            load(filename_medium_cache);
        else
            [kwave_medium, medium_plus] = medium_setup(parameters, medium_masks, planimg, pseudoCT);

            if should_save_output(parameters.io, 'save_medium_cache')
                save(filename_medium_cache, 'kwave_medium', 'medium_plus', '-v7.3');
            end
        end
    else
        disp('No medium mapping requested...no simulations will be performed.')
    end
    log_timer('stop','medium');

    log_timer('start', 'nifti_medium', parameters.io.dir_output);
    nifti_medium(parameters, planimg, medium_masks, kwave_medium, pseudoCT);
    log_timer('stop', 'nifti_medium');

    % ====================================================================
    %% STAGE 5 — SOURCE & SENSOR SETUP
    %
    % Constructs the k-Wave source matrix from the transducer geometry and
    % focus/steering settings.  Optionally runs a CFL stability check and
    % reduces the time step if the default kgrid.dt would produce an
    % unstable simulation (relevant for high-attenuation skull models).
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
    fprintf('========================================\n\n');
    log_timer('start','source', parameters.io.dir_output);

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

        if isfield(parameters.grid, 'source_limit_fraction') && parameters.grid.source_limit_fraction ~= 0
            disp('Checking CFL stability condition...')
            dt_stability_limit = checkStability(kgrid, kwave_medium);
            fprintf('Stability limit estimate for time step: %.1d.\n', dt_stability_limit);
            if ~isinf(dt_stability_limit) && kgrid.dt > dt_stability_limit
                disp('Adapting time step for simulation stability...')
                clear source sensor source_labels
                % checkStability is an approximation for heterogeneous media;
                % source_limit_fraction adds a safety margin below the estimate.
                grid_time_step = dt_stability_limit * parameters.grid.source_limit_fraction;
                [kgrid, source, sensor, source_labels] = source_sensor_setup( ...
                    parameters, max_sound_speed, trans_pos, focus_pos, grid_time_step, min_sound_speed);
            end
        end
    else
        disp('No source setup requested...no simulations will be performed.')
    end
    log_timer('stop', 'source');

    % ====================================================================
    %% FREE-WATER BASELINE
    %
    % Runs a brief in-water simulation with the same source to measure the
    % unobstructed free-field ISPPA.  This is used as a provenance record
    % and, when per-transducer target intensities are specified, to derive
    % the per-element amplitude scaling that maps the simulated output onto
    % the requested in-situ intensity.
    % ====================================================================

    acoustic_provenance = struct();
    has_per_transducer_target = isfield(parameters, 'transducer') && ...
        any(arrayfun(@(t) isfield(t, 'target_isppa_wcm2') && ...
                          ~isempty(t.target_isppa_wcm2) && ...
                          any(isfinite(t.target_isppa_wcm2)), parameters.transducer));
    run_baseline = (isfield(parameters.modules, 'run_water_baseline') && ...
                    parameters.modules.run_water_baseline == 1) || ...
                   has_per_transducer_target;
    if run_baseline && ...
            (~isfield(parameters.modules, 'run_source_setup') || parameters.modules.run_source_setup == 1)
        log_timer('start', 'freefield_baseline', parameters.io.dir_output);
        acoustic_provenance.freefield_isppa_wcm2 = ...
            water_baseline(parameters, kgrid, source, sensor);
        log_timer('stop', 'freefield_baseline');
    end

    % ====================================================================
    %% STAGE 6 — ACOUSTIC SIMULATION
    %
    % Runs k-Wave to propagate the pressure field through the medium.
    % If a cached result file already exists (and overwriting is not
    % requested), it is loaded directly.  The acoustic_cache_affix field
    % lets sequential follow-up runs point back to the base run's result
    % rather than re-running the simulation.
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
    fprintf('========================================\n\n');
    log_timer('start','acoustic', parameters.io.dir_output);

    % Thermal jobs use a per-target output_affix; acoustic_cache_affix points back to the base run.
    if isfield(parameters.io, 'acoustic_cache_affix')
        acoustic_file_affix = parameters.io.acoustic_cache_affix;
    else
        acoustic_file_affix = parameters.io.output_affix;
    end
    filename_sensor_data = fullfile(parameters.io.dir_cache, ...
        sprintf('sub-%03d_%s_results%s.mat',...
        parameters.subject_id, parameters.simulation.medium, acoustic_file_affix));

    parameters.state.acoustics_available = 0;
    if isfield(parameters.modules, 'run_acoustic_sims') && parameters.modules.run_acoustic_sims && ...
        confirm_overwriting(filename_sensor_data, parameters) && ...
        (parameters.simulation.interactive == 0 || ...
        confirmation_dlg('Running the simulations will take a long time, are you sure?', 'Yes', 'No'))

        [sensor_data, parameters, segmentation, medium_masks, kwave_medium, kgrid, source, source_labels] = ...
            acoustic_wrapper( ...
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

        % acoustic_wrapper saves the sensor data after expanding the
        % axisymmetric result to 2-D / 3-D, but the current parameters
        % struct still holds the pre-expansion axisymmetric grid dims and
        % transducer positions (because we bypassed acoustic_wrapper).
        % Detect the mismatch from the loaded array size and re-apply the
        % same [Nz x Nr] → [2*Nr x Nz] coordinate transformation so that
        % all downstream code sees consistent positions.
        if numel(parameters.grid.dims) == 2 && ...
                isfield(parameters.grid, 'axisymmetric') && parameters.grid.axisymmetric
            expand_to_3d = strcmp(parameters.simulation.medium, 'phantom') || ...
                           (isfield(parameters.modules, 'run_heating_sims') && ...
                            parameters.modules.run_heating_sims == 1);
            if ~expand_to_3d
                Nr_orig = parameters.grid.dims(2);
                Nz_orig = parameters.grid.dims(1);
                if size(sensor_data.p_max_all, 1) == Nr_orig * 2
                    for ti_load = 1:numel(parameters.transducer)
                        tp = parameters.transducer(ti_load).trans_pos;
                        fp = parameters.transducer(ti_load).focus_pos;
                        tp(2) = tp(2) + Nr_orig;
                        fp(2) = fp(2) + Nr_orig;
                        parameters.transducer(ti_load).trans_pos = fliplr(tp);
                        parameters.transducer(ti_load).focus_pos = fliplr(fp);
                    end
                    parameters.grid.dims = [Nr_orig * 2, Nz_orig];
                end
            end
        end
        sensor_data = apply_isppa_scaling(sensor_data, acoustic_provenance, parameters);
    else
        disp('No acoustic simulation available or requested ... skipping analysis')
        parameters.state.acoustics_available = 0;
        parameters.modules.run_acoustic_analysis = 0;
    end
    log_timer('stop', 'acoustic');

    % ====================================================================
    %% STAGE 7 — ACOUSTIC ANALYSIS & NIFTI EXPORT
    %
    % Derives ISPPA, ISPTA, mechanical index, focal distance, and skull
    % normalised ratio from the pressure field.  Writes the intensity and
    % pressure maps to NIfTI (T1w space and optionally MNI space).
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
    fprintf('========================================\n\n');
    log_timer('start','acoustic_analysis', parameters.io.dir_output);

    if (~isfield(parameters.modules, 'run_acoustic_analysis') || parameters.modules.run_acoustic_analysis)
        [results_acoustic, acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos] = ...
            acoustic_analysis(parameters, kwave_medium, medium_masks, sensor_data, segmentation, source_labels);
    else
        disp('No acoustic simulation results available (or requested). Skipping analysis...')
        results_acoustic = [];
        acoustic_Ipa     = [];
        acoustic_MI      = [];
        acoustic_pressure = [];
        highlighted_pos  = [];
    end
    log_timer('stop', 'acoustic_analysis');

    if parameters.state.acoustics_available
        log_timer('start', 'nifti_acoustic', parameters.io.dir_output);
        nifti_acoustic(parameters, planimg, results_acoustic, ...
            acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos);
        log_timer('stop', 'nifti_acoustic');
    end
    % acoustic pressure maps no longer needed after NIfTI export;
    % highlighted_pos still required by thermal_analysis below
    clear results_acoustic acoustic_Ipa acoustic_MI acoustic_pressure

    clear source_labels

    % ====================================================================
    %% STAGE 8 — THERMAL SIMULATION
    %
    % Solves the Pennes bio-heat equation using k-Wave's thermal diffusion
    % solver over the full sonication sequence (on/off pulses, ISI, …).
    % Heat deposition is proportional to the acoustic absorption at each
    % voxel derived from the pressure field above.
    %
    % In sequential mode (follow-up run), the temperature and CEM43 maps
    % from the previous run are used as initial conditions so that
    % cumulative heating is captured correctly across stimulation bouts.
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
    fprintf('========================================\n\n');
    log_timer('start','thermal', parameters.io.dir_output);

    parameters.state.heating_available = 0;
    if isfield(parameters.modules, 'run_heating_sims') && parameters.modules.run_heating_sims && parameters.state.acoustics_available == 1

        disp('Starting thermal simulations...')
        % sequential_pipeline gives each follow-up run a unique thermal_cache_affix
        % (e.g. <base_affix>_seq1) so the "already done" cache check does not
        % silently load the base run's thermal result.
        if isfield(parameters.io, 'thermal_cache_affix')
            thermal_file_affix = parameters.io.thermal_cache_affix;
        else
            thermal_file_affix = parameters.io.output_affix;
        end
        filename_heating_data = fullfile(parameters.io.dir_cache,...
            sprintf('sub-%03d_%s_heating_res%s.mat',...
            parameters.subject_id, parameters.simulation.medium, thermal_file_affix));
        
        if confirm_overwriting(filename_heating_data, parameters) && (parameters.simulation.interactive == 0 || ...
            confirmation_dlg('Running the thermal simulations will take a long time, are you sure?', 'Yes', 'No'))

            kwave_medium.temp_0              = medium_plus.temp_0;
            kwave_medium.absorption_fraction = medium_plus.absorption_fraction;
            clear medium_plus;

            [kwaveDiffusion, time_status_seq, results_heating] = ...
                thermal_simulation(...
                parameters, sensor_data, kgrid, kwave_medium, sensor, source, planimg.transf, medium_masks);

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
            warning('prestus_pipeline:thermalNoAcoustics', ...
                ['Heating simulations requested but no acoustic results are available. ' ...
                 'Check that run_acoustic_sims is enabled and the acoustic cache file exists.']);
            parameters.state.heating_available = 0;
            parameters.modules.run_thermal_analysis = 0;
            results_heating = struct();
        end
    else
        parameters.state.heating_available = 0;
        results_heating = struct();
    end
    log_timer('stop','thermal');

    % ====================================================================
    %% STAGE 9 — THERMAL ANALYSIS & NIFTI EXPORT
    %
    % Extracts safety metrics (peak temperature rise, CEM43, ISO CEM43)
    % per tissue layer and compares them against regulatory limits.
    % Writes heating and CEM43 maps as NIfTI volumes in T1w (and
    % optionally MNI) space.
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
    fprintf('========================================\n\n');
    log_timer('start','thermal_analysis', parameters.io.dir_output);

    if parameters.state.heating_available == 1 && ...
            (~isfield(parameters.modules, 'run_thermal_analysis') || parameters.modules.run_thermal_analysis)
        thermal_analysis(parameters, results_heating, time_status_seq, ...
            medium_masks, highlighted_pos, segmentation);
    else
        disp('No heating simulation results available (or requested). Skipping thermal analysis...')
    end
    log_timer('stop','thermal_analysis');
    % highlighted_pos and segmentation consumed by thermal_analysis above
    clear highlighted_pos segmentation

    if parameters.state.heating_available
        log_timer('start', 'nifti_thermal', parameters.io.dir_output);
        nifti_thermal(parameters, planimg, results_heating, kwave_medium);
        log_timer('stop', 'nifti_thermal');
    end

    clear time_status_seq segmentation bone_mask pseudoCT acoustic_* results_heating medium_masks kwave_medium planimg

    % ====================================================================
    %% STAGE 10 — REPORTS
    %
    % Generates one or more HTML summary reports depending on which report
    % flags are set.  Reports are written after all simulation outputs are
    % on disk so they can embed figures and metrics from this run.
    % ====================================================================

    fprintf('========================================\n');
    fprintf('\n');
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
            warning('prestus_pipeline:missingMultiIsppa', ...
                'multi_isppa_report flag set but parameters.multi_isppa is missing - skipping.');
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
    %% STAGE 11 — POST-HOC WATER SIMULATION  (optional)
    %
    % Re-runs the acoustic simulation in a homogeneous water medium using
    % the same transducer and grid geometry.  Provides a skull-free
    % reference for normalisation (e.g. skull transmission ratio).
    % Only applicable to layered / phantom media; water simulations skip
    % this stage automatically.
    % ====================================================================

    if isfield(parameters.modules, 'run_posthoc_water_sims') && parameters.modules.run_posthoc_water_sims && ...
            contains(parameters.simulation.medium, {'layered', 'phantom'})

        fprintf('========================================\n');
    fprintf('\n');
    fprintf('========================================\n\n');

        if numel(parameters.transducer) > 1
            warning('prestus_pipeline:posthocMultiTransducer', ...
                ['Post-hoc water simulations are not implemented for multiple transducers. ' ...
                 'Only the first transducer will be used. ' ...
                 'Consider running separate single-transducer configs instead.']);
        end
        water_parameters = parameters;
        water_parameters.simulation.medium = 'water';
        water_parameters.modules.run_heating_sims = 0;
        water_parameters.modules.run_posthoc_water_sims = 0;
        water_parameters.simulation.debug = 0;
        water_parameters.grid.default_dims = water_parameters.grid.dims;
        water_parameters.io.log_file = '';  % new log file so water run does not append to the layered diary
        water_parameters.hpc.timelimit = '05:00:00';
        water_parameters.hpc.memorylimit = 40;
        water_parameters.hpc.wait_for_job = false;
        prestus_pipeline_start(water_parameters);
        clear water_parameters;
    end

    % ====================================================================
    %% STAGE 12 — SEQUENTIAL FOLLOW-UP SIMULATION  (optional)
    %
    % When options.sequential_configs is present, dispatches the next
    % simulation in the sequence with the end-temperature and CEM43 maps
    % of this run wired up as initial conditions.  The entire chain
    % recurses through sequential_pipeline until all configs are consumed,
    % then a multi-run summary report is generated.
    % ====================================================================

    if any(strcmp(fieldnames(options), 'sequential_configs'))
        fprintf('========================================\n');
        fprintf('FOLLOW-UP SIMULATION WITH IDENTICAL MEDIUM\n');
        fprintf('========================================\n\n');
        sequential_pipeline(parameters, options);
    end

    % ====================================================================
    %% CLEANUP
    % Clear simulation workspace variables that are no longer needed after
    % all pipeline stages and dispatches have completed.
    % ====================================================================

    clear kgrid source sensor sensor_data trans_pos focus_pos ...
          acoustic_provenance results_acoustic ...
          acoustic_Ipa acoustic_MI acoustic_pressure highlighted_pos

    track_usage('run_end', parameters, struct('duration_s', toc(telemetry_t0), 'status', 'success', 'run_id', telemetry_rid));

    catch me_
        track_usage('run_error', parameters, struct( ...
            'duration_s', toc(telemetry_t0), ...
            'status',     'error', ...
            'run_id',     telemetry_rid, ...
            'error_id',   me_.identifier));
        rethrow(me_);
    end

end
