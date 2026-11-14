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
%   3. Grid & head model grid_stage — reorientation, rescaling, smoothing,
%                        and cropping of the segmented head to the simulation grid
%   4. Medium setup      medium_stage — assignment of acoustic/thermal tissue
%                        properties (density, sound speed, attenuation, …)
%   5. Source setup      Construction of the transducer source array
%                        (kWaveArray or analytical; annular or matrix)
%   6. Acoustic sim      acoustic_stage — k-Wave / jwave / BabelBrain
%   7. Acoustic analysis ISPPA, ISPTA, MI, focal metrics; NIfTI export
%   8. Thermal sim       thermal_stage — bio-heat equation over sonication protocol
%   9. Thermal analysis  Peak temperature, CEM43, safety metrics; NIfTI export
%  10. Report            Self-contained HTML summary report
%  11. Post-hoc water sim  Optional repeat of acoustic sim in homogeneous water
%  12. Sequential follow-up  Optional chained simulation (end-temperature handoff)
%
% Cache files written at each stage allow individual stages to be skipped on
% re-runs.  Stage-level toggles live in parameters.modules.*; see
% config/config_default.yaml for all available flags.
%
% Use as:
%   parameters = prestus_pipeline(parameters)
%   parameters = prestus_pipeline(parameters, options)
%
% Input:
%   parameters - struct built from config_default.yaml and subject-specific
%                override config; must include subject_id, path.sim, and all
%                transducer/grid settings
%   options    - (optional) struct; supports:
%                  .sequential_configs  struct of follow-up configs for chained runs
%
% Output:
%   parameters - updated struct with pipeline outputs (io.dir_output, etc.)
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

    % Stop before any computation if a tool a requested stage needs is missing
    % (SimNIBS binaries, k-Wave binary, GPU, JVM).
    check_tool_availability(parameters);

    % ====================================================================
    %% TELEMETRY  (opt-in, fully anonymous)
    % ====================================================================

    telemetry_setup();
    telemetry_t0  = tic;
    telemetry_rid = generate_run_id();
    run_stats('reset');
    try

    % ====================================================================
    %% SPECIAL-MODE INTERCEPTS
    % ====================================================================

    if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'uncertainty') && ...
            parameters.simulation.uncertainty
        uncertainty_pipeline(parameters, options);
        return;
    end

    if is_multi_isppa_mode(parameters)
        multi_isppa_pipeline(parameters, options);
        return;
    end

    if is_pose_workflow_mode(parameters)
        pose_workflow_pipeline(parameters, options);
        return;
    end

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
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
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
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
    log_timer('start','placement', parameters.io.dir_output);

    if ~isfield(parameters.modules, 'run_transducer_placement') || ...
            parameters.modules.run_transducer_placement == 1
        parameters = preproc_transducer_placement(parameters);
    else
        % Without the placement stage RAS/MNI mm would be read as voxels
        if ~strcmp(manual_coordinate_space(parameters), 'voxel')
            error('PRESTUS:placement:unconvertedPositions', ...
                ['placement.manual.coordinate_space = ''%s'' requires ' ...
                 'modules.run_transducer_placement = 1 (positions are converted there).'], ...
                manual_coordinate_space(parameters));
        end
        disp('Transducer placement stage skipped.');
    end
    log_timer('stop','placement');

    % Planned positions in T1 voxels, kept for the targeting analysis: the
    % grid stage replaces transducer(i).trans_pos/focus_pos with grid indices.
    % For a non-RAS+ T1, GRID_TISSUE_SETUP re-records them in the reoriented grid.
    parameters.state.planned_trans_t1_vox = planned_positions(parameters, 'trans_pos');
    parameters.state.planned_focus_t1_vox = planned_positions(parameters, 'focus_pos');

    % ====================================================================
    %% STAGE 3 — GRID SETUP & HEAD PREPROCESSING
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
    log_timer('start','preproc', parameters.io.dir_output);

    [parameters, planimg, medium_masks, segmentation, bone_mask, pseudoCT, trans_pos, focus_pos] = ...
        grid_stage(parameters);

    log_timer('stop','preproc');

    % ====================================================================
    %% STAGE 4 — MEDIUM PROPERTY MAPPING
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
    log_timer('start','medium', parameters.io.dir_output);

    [kwave_medium, medium_plus, planimg, medium_masks] = ...
        medium_stage(parameters, medium_masks, planimg, pseudoCT);

    log_timer('stop','medium');

    log_timer('start', 'nifti_medium', parameters.io.dir_output);
    % property maps include the effective absorption maps kept in medium_plus
    medium_maps = kwave_medium;
    for f = intersect(fieldnames(medium_plus), {'absorption_fraction', 'absorption_np'})'
        medium_maps.(f{1}) = medium_plus.(f{1});
    end
    nifti_medium(parameters, planimg, medium_masks, medium_maps, pseudoCT);
    clear medium_maps
    log_timer('stop', 'nifti_medium');
    clear bone_mask

    % ====================================================================
    %% ACOUSTIC FOV CROP  (between medium and source setup)
    %
    % Resolves 'auto' FOV settings, applies jwave PML expansion if needed,
    % and crops medium/masks/segmentation to the beam-axis FOV sub-grid.
    % fov_state carries the full-grid backups for thermal restore and
    % nifti_acoustic back-projection.
    % ====================================================================

    [parameters, kwave_medium, medium_masks, segmentation, trans_pos, focus_pos, fov_state] = ...
        acoustic_fov_crop(parameters, kwave_medium, medium_masks, segmentation, trans_pos, focus_pos);

    % Record the FOV-crop voxel offset on planimg so any downstream stage that
    % exports the (possibly FOV-cropped) grid directly to RAS space — e.g.
    % babelbrain_wrapper's write_babelbrain_input — can build a correct affine
    % via grid_to_ras_affine without needing fov_state threaded separately.
    if fov_state.active
        planimg.fov_offset_ac = fov_state.offset;
    else
        planimg.fov_offset_ac = [1 1 1];
    end

    % ====================================================================
    %% STAGE 5 — SOURCE & SENSOR SETUP
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
    log_timer('start','source', parameters.io.dir_output);

    use_babelbrain_acoustic = isfield(parameters, 'acoustic_backend') && ...
        strcmpi(parameters.acoustic_backend, 'babelbrain');
    coherent_drive = [];

    if (~isfield(parameters.modules, 'run_source_setup') || parameters.modules.run_source_setup == 1) && ...
            ~use_babelbrain_acoustic
        max_sound_speed = max(kwave_medium.sound_speed(:));
        min_sound_speed = min(kwave_medium.sound_speed(:));
        % Coherent multi-transducer: set each transducer's drive to its own
        % ISPPA target (free-water baseline of that transducer alone)
        if numel(parameters.transducer) > 1 && ~is_async_mode(parameters)
            simulate_baseline = @(p) single_source_water_baseline(p, max_sound_speed, ...
                min_sound_speed, trans_pos, focus_pos);
            [parameters, coherent_drive] = coherent_drive_calibration(parameters, simulate_baseline);
            clear simulate_baseline
        end
        [kgrid, source, sensor, source_labels, src_weights] = source_sensor_setup( ...
            parameters, max_sound_speed, trans_pos, focus_pos, [], min_sound_speed);
        if isfield(parameters.grid, 'source_limit_fraction') && parameters.grid.source_limit_fraction ~= 0
            disp('Checking CFL stability condition...')
            dt_stability_limit = checkStability(kgrid, kwave_medium);
            fprintf('Stability limit estimate for time step: %.1d.\n', dt_stability_limit);
            if ~isinf(dt_stability_limit) && kgrid.dt > dt_stability_limit
                disp('Adapting time step for simulation stability...')
                clear source sensor source_labels
                grid_time_step = dt_stability_limit * parameters.grid.source_limit_fraction;
                [kgrid, source, sensor, source_labels, src_weights] = source_sensor_setup( ...
                    parameters, max_sound_speed, trans_pos, focus_pos, grid_time_step, min_sound_speed);
            end
        end
    elseif use_babelbrain_acoustic
        disp('[pipeline] BabelBrain acoustic backend: skipping source setup.')
        kgrid = []; source = []; sensor = []; source_labels = []; src_weights = [];
    else
        disp('No source setup requested...no simulations will be performed.')
        kgrid = []; source = []; sensor = []; source_labels = []; src_weights = [];
    end
    log_timer('stop', 'source');

    % ====================================================================
    %% ACOUSTIC PROVENANCE INIT & FREE-WATER BASELINE
    % ====================================================================

    acoustic_provenance = struct();
    if fov_state.active
        acoustic_provenance.fov_offset_ac = fov_state.offset;
        acoustic_provenance.full_ac_dims  = fov_state.full_dims;
    end
    tr1 = parameters.transducer(1);
    if isfield(tr1, 'tx_calibration') && isfield(tr1.tx_calibration, 'isppa_wcm2') && ...
            isscalar(tr1.tx_calibration.isppa_wcm2) && isfinite(tr1.tx_calibration.isppa_wcm2)
        acoustic_provenance.tx_calibration = tr1.tx_calibration;
    end
    clear tr1
    if ~isempty(coherent_drive)
        acoustic_provenance.coherent_drive = coherent_drive;
    end

    tr_only = isfield(parameters, 'time_reversal') && parameters.time_reversal.enabled && ...
        isfield(parameters.time_reversal, 'only') && parameters.time_reversal.only;
    if ~tr_only
        acoustic_provenance = acoustic_water_baseline( ...
            parameters, acoustic_provenance, kgrid, source, sensor, source_labels);
    end

    % ====================================================================
    %% TIME REVERSAL  (optional: refocusing phases, aberration metrics, pose search)
    % ====================================================================
    % After the free-water baseline, so that ISPPA targeting keeps the drive
    % amplitudes of the geometric drive when time-reversal phases are applied.

    if isfield(parameters, 'time_reversal') && parameters.time_reversal.enabled
        if isempty(kgrid)
            warning('prestus_pipeline:timeReversalNoGrid', ...
                ['time_reversal needs the k-Wave source setup; skipped (BabelBrain has its own ' ...
                 'refocusing: babelbrain.do_refocusing).']);
        else
            log_timer('start', 'time_reversal', parameters.io.dir_output);
            % The forward source (its signals: GBs at fine resolution) is not
            % needed by the reverse run; the stage rebuilds it afterwards
            source = struct();
            [parameters, source, source_labels] = time_reversal_stage(parameters, kgrid, kwave_medium, ...
                medium_masks, sensor, src_weights, trans_pos, focus_pos, source, source_labels, planimg, fov_state);
            log_timer('stop', 'time_reversal');
        end
    end
    clear src_weights

    if tr_only
        fprintf('time_reversal.only: no acoustic or thermal simulation. Finishing.\n');
        log_timer('stop','prestus_pipeline')
        diary('off')
        run_stats('write', parameters, toc(telemetry_t0));
        return;
    end

    % ====================================================================
    %% STAGE 6 — ACOUSTIC SIMULATION
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
    log_timer('start','acoustic', parameters.io.dir_output);

    [sensor_data, parameters, kgrid, kwave_medium, medium_masks, segmentation, ...
     source, source_labels, acoustic_provenance] = ...
        acoustic_stage(parameters, kgrid, kwave_medium, source, sensor, ...
                       medium_masks, segmentation, source_labels, acoustic_provenance, ...
                       planimg, trans_pos, focus_pos, pseudoCT);
    clear pseudoCT

    log_timer('stop', 'acoustic');

    if parameters.state.acoustic_dry_run
        fprintf('Acoustic dry run: solver inputs written, no field simulated. Finishing.\n');
        if isfield(parameters.modules, 'store_acoustic_results') && parameters.modules.store_acoustic_results
            parameters.benchmark_sensor_data = sensor_data;
        end
        log_timer('stop','prestus_pipeline')
        diary('off')
        return;
    end

    % ====================================================================
    %% STAGE 7 — ACOUSTIC ANALYSIS & NIFTI EXPORT
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
    log_timer('start','acoustic_analysis', parameters.io.dir_output);

    if ~isfield(parameters.modules, 'run_acoustic_analysis') || parameters.modules.run_acoustic_analysis
        [results_acoustic, acoustic_Itp, acoustic_Ipa, acoustic_Ita, acoustic_MI, acoustic_pressure, highlighted_pos] = ...
            acoustic_analysis(parameters, kwave_medium, medium_masks, sensor_data, ...
                              segmentation, source_labels, [], acoustic_provenance);
    else
        disp('No acoustic simulation results available (or requested). Skipping analysis...')
        results_acoustic  = [];
        acoustic_Itp      = [];
        acoustic_Ipa      = [];
        acoustic_Ita      = [];
        acoustic_MI       = [];
        acoustic_pressure = [];
        highlighted_pos   = [];
    end
    log_timer('stop', 'acoustic_analysis');

    if parameters.state.acoustics_available
        log_timer('start', 'nifti_acoustic', parameters.io.dir_output);
        if fov_state.active
            % Back-project FOV-sized analysis volumes onto the full grid before
            % NIfTI export so outputs are in full simulation-grid space.
            params_nii           = parameters;
            params_nii.grid.dims = fov_state.full_dims;
            params_nii = fov_restore_transducers(params_nii, fov_state);
            % FOV box in full-grid voxels, outlined on the T1 overlay.
            params_nii.grid.acoustic_fov_offset_vox = fov_state.offset;
            params_nii.grid.acoustic_fov_dims_vox   = parameters.grid.dims;
            for fmf = {'focal_mask_3dB', 'focal_mask_6dB'}
                fn = fmf{1};
                if isfield(results_acoustic, fn) && ~isempty(results_acoustic.(fn))
                    results_acoustic.(fn) = fov_backproject_volume( ...
                        uint8(results_acoustic.(fn)), fov_state.offset, fov_state.full_dims) > 0;
                end
            end
            acoustic_Ita_full = [];
            if ~isempty(acoustic_Ita)
                acoustic_Ita_full = fov_backproject_volume(acoustic_Ita, fov_state.offset, fov_state.full_dims);
            end
            nifti_acoustic(params_nii, planimg, results_acoustic, ...
                fov_backproject_volume(acoustic_Itp,      fov_state.offset, fov_state.full_dims), ...
                fov_backproject_volume(acoustic_Ipa,      fov_state.offset, fov_state.full_dims), ...
                acoustic_Ita_full, ...
                fov_backproject_volume(acoustic_MI,       fov_state.offset, fov_state.full_dims), ...
                fov_backproject_volume(acoustic_pressure, fov_state.offset, fov_state.full_dims), ...
                fov_to_grid_pos(highlighted_pos, fov_state.offset));
            clear params_nii acoustic_Ita_full
        else
            nifti_acoustic(parameters, planimg, results_acoustic, ...
                acoustic_Itp, acoustic_Ipa, acoustic_Ita, acoustic_MI, acoustic_pressure, highlighted_pos);
        end
        log_timer('stop', 'nifti_acoustic');

        % Sidecar of the pressure map, so a later run can start from it
        % (path.acoustic_input) instead of simulating
        write_acoustic_sidecar(parameters, acoustic_provenance, sensor_data);

        % Targeting accuracy from the focal-mask NIfTIs just written
        if (~isfield(parameters.modules, 'run_targeting_analysis') || parameters.modules.run_targeting_analysis) && ...
                (~isfield(parameters.modules, 'run_nifti_creation') || parameters.modules.run_nifti_creation) && ...
                contains(parameters.simulation.medium, {'layered', 'phantom'})
            try
                targeting_analysis(parameters);
            catch ME
                warning('prestus_pipeline:targeting', 'Targeting analysis failed: %s', ME.message);
            end
        end
    end

    if isfield(parameters, 'modules') && isfield(parameters.modules, 'store_acoustic_results') && ...
            parameters.modules.store_acoustic_results
        if exist('results_acoustic', 'var'), parameters.results_acoustic      = results_acoustic; end
        if exist('sensor_data',      'var'), parameters.benchmark_sensor_data = sensor_data;      end
    end
    clear results_acoustic acoustic_Itp acoustic_Ipa acoustic_Ita acoustic_MI acoustic_pressure

    % Optional BabelBrain input export (io.export_babelbrain or acoustic_backend='babelbrain').
    do_bb_export = (isfield(parameters.io, 'export_babelbrain') && parameters.io.export_babelbrain) || ...
                   (isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'babelbrain'));
    if do_bb_export
        if ~isfield(parameters.io, 'dir_babelbrain') || isempty(parameters.io.dir_babelbrain)
            parameters.io.dir_babelbrain = fullfile(parameters.io.dir_output, 'babelbrain');
        end
        export_babelbrain_inputs(parameters, medium_masks, focus_pos, planimg, struct());
    end

    clear source_labels

    % ====================================================================
    %% STAGE 8 — THERMAL SIMULATION
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
    log_timer('start','thermal', parameters.io.dir_output);

    [results_heating, time_status_seq, kwave_medium, medium_masks, segmentation, ...
     planimg, highlighted_pos, parameters] = ...
        thermal_stage(parameters, sensor_data, kgrid, kwave_medium, medium_plus, ...
                      medium_masks, segmentation, planimg, source, sensor, ...
                      trans_pos, focus_pos, acoustic_provenance, fov_state, highlighted_pos);

    log_timer('stop','thermal');
    clear sensor_data source sensor kgrid medium_plus acoustic_provenance

    % ====================================================================
    %% STAGE 9 — THERMAL ANALYSIS & NIFTI EXPORT
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
    log_timer('start','thermal_analysis', parameters.io.dir_output);

    if parameters.state.heating_available == 1 && ...
            (~isfield(parameters.modules, 'run_thermal_analysis') || parameters.modules.run_thermal_analysis)
        thermal_analysis(parameters, results_heating, time_status_seq, ...
            medium_masks, highlighted_pos, segmentation);
    else
        disp('No heating simulation results available (or requested). Skipping thermal analysis...')
    end
    log_timer('stop','thermal_analysis');
    clear highlighted_pos segmentation

    if parameters.state.heating_available
        log_timer('start', 'nifti_thermal', parameters.io.dir_output);
        nifti_thermal(parameters, planimg, results_heating, kwave_medium);
        log_timer('stop', 'nifti_thermal');
    end

    clear time_status_seq results_heating medium_masks kwave_medium planimg fov_state

    % ====================================================================
    %% STAGE 10 — REPORTS
    % ====================================================================

    fprintf('========================================\n\n========================================\n\n');
    log_timer('stop','prestus_pipeline')
    disp('Pipeline finished successfully');

    if isfield(parameters.modules, 'generate_report') && parameters.modules.generate_report
        % Final parameter state as seen by the report, so the report can be
        % re-rendered later without re-running (see preview_report_styles).
        try
            save(fullfile(parameters.io.dir_cache, sprintf('sub-%03d_%s_report_parameters%s.mat', ...
                parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix)), 'parameters');
        catch ME
            warning('prestus_pipeline:reportParams', 'Could not save report parameters: %s', ME.message);
        end
        generate_simulation_report(parameters);
    end

    if isfield(parameters.modules, 'generate_group_report') && parameters.modules.generate_group_report
        try
            group_subjects = discover_group_subjects(parameters);
            if ~isempty(group_subjects)
                generate_group_report(parameters, group_subjects);
            end
        catch ME
            warning('prestus_pipeline:groupReport', ...
                'Group report generation failed: %s', ME.message);
        end
    end

    if isfield(parameters.modules, 'multi_isppa_report') && parameters.modules.multi_isppa_report
        if isfield(parameters, 'multi_isppa')
            generate_multi_isppa_report(parameters, parameters.multi_isppa);
        else
            warning('prestus_pipeline:missingMultiIsppa', ...
                'multi_isppa_report flag set but parameters.multi_isppa is missing — skipping.');
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

    % Record this run's size and resources before post-hoc water or
    % follow-up runs (which may run inline) start records of their own.
    run_stats('write', parameters, toc(telemetry_t0));

    % ====================================================================
    %% STAGE 11 — POST-HOC WATER SIMULATION  (optional)
    % ====================================================================

    if isfield(parameters.modules, 'run_posthoc_water_sims') && ...
            parameters.modules.run_posthoc_water_sims && ...
            contains(parameters.simulation.medium, {'layered', 'phantom'})

        fprintf('========================================\n\n========================================\n\n');
        transducer_N_water = numel(parameters.transducer);
        for it_water = 1:transducer_N_water
            if transducer_N_water > 1
                fprintf('[pipeline] Post-hoc water sim — transducer %d / %d\n', it_water, transducer_N_water);
            end
            water_parameters = posthoc_water_parameters(parameters, it_water);
            prestus_pipeline_start(water_parameters);
        end
        clear water_parameters
    end

    % ====================================================================
    %% STAGE 12 — SEQUENTIAL FOLLOW-UP SIMULATION  (optional)
    % ====================================================================

    if any(strcmp(fieldnames(options), 'sequential_configs'))
        fprintf('========================================\n');
        fprintf('FOLLOW-UP SIMULATION WITH IDENTICAL MEDIUM\n');
        fprintf('========================================\n\n');
        sequential_pipeline(parameters, options);
    end

    track_usage('run_end', parameters, ...
        struct('duration_s', toc(telemetry_t0), 'status', 'success', 'run_id', telemetry_rid), options);

    catch me_
        diary('off')   % do not keep logging into this run's log after it failed
        track_usage('run_error', parameters, struct( ...
            'duration_s', toc(telemetry_t0), ...
            'status',     'error', ...
            'run_id',     telemetry_rid, ...
            'error_id',   me_.identifier), options);
        rethrow(me_);
    end

end

function pos = planned_positions(parameters, field)
% [N×3] transducer(i).<field> for all transducers; [] when any is unset.
    pos = [];
    if ~isfield(parameters, 'transducer') || ~isfield(parameters.transducer, field)
        return
    end
    for ti = 1:numel(parameters.transducer)
        v = parameters.transducer(ti).(field);
        if ~isnumeric(v) || numel(v) ~= 3
            pos = [];
            return
        end
        pos(ti, :) = double(v(:)'); %#ok<AGROW>
    end
end

function isppa = single_source_water_baseline(parameters, max_sound_speed, min_sound_speed, trans_pos, focus_pos)
% Free-water ISPPA of the source in parameters (COHERENT_DRIVE_CALIBRATION
% passes all transducers but one at zero amplitude): same grid and geometry as
% the run, built and simulated in water like the run's own baseline. The
% source gets its own affix, is always rebuilt and never saved, so it can
% neither come from nor replace the run's cached source (kwave_source*.mat).
    parameters.io.output_affix         = [parameters.io.output_affix, '_drivecal'];
    parameters.io.overwrite_files      = 'always';
    parameters.io.save_source_matrices = 0;
    [kgrid, source, sensor, source_labels] = source_sensor_setup( ...
        parameters, max_sound_speed, trans_pos, focus_pos, [], min_sound_speed);
    isppa = water_baseline(parameters, kgrid, source, sensor, source_labels);
end
