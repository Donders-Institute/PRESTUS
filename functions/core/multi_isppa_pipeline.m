function jobs = multi_isppa_pipeline(parameters, options)
% MULTI_ISPPA_PIPELINE  Run acoustic simulation once, then thermal analysis at
%                       multiple free-water target intensities in parallel
%
% When transducer(1).target_isppa_wcm2 contains more than one value, the
% pipeline automatically splits into:
%
%   Stage 1  Acoustic simulation + free-water baseline measurement
%            Cached with acoustic_provenance.freefield_isppa_wcm2 recorded
%            alongside the sensor_data. No thermal simulation.
%
%   Stages 2…N  One thermal analysis job per target_isppa_wcm2 value.
%               Each job loads the cached acoustic result, scales
%               sensor_data.p_max_all by sqrt(target/baseline), runs the
%               thermal simulation and analysis, and generates an HTML report.
%               Output files carry an affix _desc-isppaNNN (NNN = integer W/cm²).
%
%   Stage N+1   (optional) Summary report aggregating all ISPPA variants.
%               Active when options.generate_summary_report = true (default).
%
% Platform behaviour:
%   matlab     Stages run sequentially in the current MATLAB session.
%   slurm/qsub Stage 1 is submitted; stages 2…N carry an afterok dependency
%              on stage 1; stage N+1 depends on all thermal jobs. With
%              thermal_backend 'babelbrain' each thermal job also waits for
%              the previous one (BabelBrain's thermal output file is named
%              after the shared acoustic file, not the target).
%              All jobs are submitted immediately without blocking.
%
% Overwrite semantics (see RESOLVE_CHAINED_OVERWRITE):
%   io.overwrite_files = 'never'   Stage 1 resumes from its acoustic cache;
%                                  thermal/summary stages are skipped when
%                                  their outputs exist.
%   io.overwrite_files = 'always'  Stage 1 resumes; thermal/summary stages
%                                  always re-run and overwrite their outputs.
%   io.overwrite_files = 'force'   every stage re-runs.
%   Thermal/summary sub-jobs always load the Stage-1 head/grid/medium caches
%   (io.reuse_upstream_caches).
%
% Use as:
%   multi_isppa_pipeline(parameters)
%   multi_isppa_pipeline(parameters, options)
%
% Input:
%   parameters - base PRESTUS config; transducer(1).target_isppa_wcm2 must be
%                a numeric vector with at least two elements [W/cm²];
%                do NOT set io.output_affix — this function manages it
%   options    - pipeline knobs (optional, default: struct()):
%                .generate_summary_report  - true (default) / false
%                .acoustic_timelimit       - HPC wall time for stage 1
%                                            (default '06:00:00')
%                .thermal_timelimit        - HPC wall time per thermal job
%                                            (default '02:00:00')
%                .summary_timelimit        - HPC wall time for summary stage
%                                            (default '00:30:00')
%                .acoustic_memorylimit     - RAM in GB for stage 1
%                .thermal_memorylimit      - RAM in GB per thermal job
%                .summary_memorylimit      - RAM in GB for summary stage
%                .acoustic_partition       - HPC partition for stage 1
%                .summary_partition        - HPC partition for summary stage
%                .depend_job_ids           - scheduler job IDs the first submitted
%                                            stage waits for (default []); set by a
%                                            caller that nests this run, e.g.
%                                            UNCERTAINTY_PIPELINE
%                .job_name_prefix          - inserted after 'PRESTUS-' in job names
%                                            (default '')
%
% Output:
%   jobs       - struct: .acoustic, .thermal, .summary (scheduler job IDs; [] for
%                skipped stages and in MATLAB), .all (every submitted ID), and
%                .output_affixes (1×M cell, io.output_affix of each target's outputs)
%
% See also: PRESTUS_PIPELINE_START, PRESTUS_PIPELINE,
%           MEASURE_FREEFIELD_BASELINE, APPLY_ISPPA_SCALING

arguments
    parameters (1,1) struct
    options    (1,1) struct = struct()
end

% =========================================================================
%% Validate
% =========================================================================

% Coherent multi-transducer runs set each transducer's drive from its own target
% before the simulation (COHERENT_DRIVE_CALIBRATION): one simulation per sweep point.
if numel(parameters.transducer) > 1
    error('prestus:isppa:coherentSweep', ...
        ['ISPPA sweeps are not available for coherent multi-transducer runs: each transducer''s ' ...
         'drive is set from its own target before the simulation, so every sweep point needs its ' ...
         'own simulation. Run one config per point, or use simulation.transducer_coupling: async.']);
end

% Targets are read from transducer(1).target_isppa_wcm2 (canonical location).
% calibration.target_isppa_wcm2 is the deprecated fallback handled by load_parameters.
if ~isfield(parameters, 'transducer') || isempty(parameters.transducer) || ...
        ~isfield(parameters.transducer(1), 'target_isppa_wcm2') || ...
        numel(parameters.transducer(1).target_isppa_wcm2) < 2 || ...
        ~any(isfinite(parameters.transducer(1).target_isppa_wcm2))
    error(['multi_isppa_pipeline: transducer(1).target_isppa_wcm2 must contain at least ' ...
           'two finite values. Set it in the transducer config or YAML.']);
end

targets = parameters.transducer(1).target_isppa_wcm2(:)';

% =========================================================================
%% Resolve options
% =========================================================================

if ~isfield(options, 'generate_summary_report'); options.generate_summary_report = true;   end
if ~isfield(options, 'acoustic_timelimit');      options.acoustic_timelimit      = '06:00:00'; end
if ~isfield(options, 'thermal_timelimit');       options.thermal_timelimit       = '02:00:00'; end
if ~isfield(options, 'summary_timelimit');       options.summary_timelimit       = '00:30:00'; end
if ~isfield(options, 'acoustic_memorylimit');    options.acoustic_memorylimit    = []; end
if ~isfield(options, 'thermal_memorylimit');     options.thermal_memorylimit     = []; end
if ~isfield(options, 'summary_memorylimit');     options.summary_memorylimit     = []; end
if ~isfield(options, 'acoustic_partition');      options.acoustic_partition      = []; end
if ~isfield(options, 'summary_partition');       options.summary_partition       = []; end
if ~isfield(options, 'depend_job_ids');          options.depend_job_ids          = []; end
if ~isfield(options, 'job_name_prefix');         options.job_name_prefix         = ''; end

% =========================================================================
%% Resolve platform
% =========================================================================

if ~isfield(parameters, 'platform')
    parameters.platform = 'auto';
end
platform = parameters.platform;
if strcmp(platform, 'auto')
    platform = hpc_detect_system();
end

% =========================================================================
%% Print header
% =========================================================================

fprintf('\n');
fprintf('╔══════════════════════════════════════╗\n');
fprintf('║   PRESTUS Multi-ISPPA Pipeline       ║\n');
fprintf('╠══════════════════════════════════════╣\n');
fprintf('║  Subject : sub-%03d                  ║\n', parameters.subject_id);
fprintf('║  Medium  : %-26s ║\n',                     parameters.simulation.medium);
fprintf('║  Platform: %-26s ║\n',                     upper(platform));
fprintf('║  Targets : %-26s ║\n', ...
    strjoin(arrayfun(@(x) sprintf('%.0f', x), targets, 'UniformOutput', false), ', '));
fprintf('╚══════════════════════════════════════╝\n\n');

% =========================================================================
%% Build parameter structs
% =========================================================================

base_affix = '';
if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix') && ...
        ~isempty(parameters.io.output_affix)
    base_affix = parameters.io.output_affix;
end

ow = resolve_chained_overwrite(parameters);

p_acoustic = make_acoustic_params(parameters, base_affix, ow);

p_thermal  = cell(1, numel(targets));
affixes    = cell(1, numel(targets));
for ti = 1:numel(targets)
    affixes{ti}   = sprintf('%s_desc-isppa%03d', base_affix, round(targets(ti)));
    p_thermal{ti} = make_thermal_params(parameters, targets(ti), affixes{ti}, base_affix, ow);
end

p_summary = make_summary_params(parameters, targets, affixes, base_affix, ow);

jobs = struct('acoustic', [], 'thermal', [], 'summary', [], 'all', [], ...
              'output_affixes', {affixes});

% =========================================================================
%% Execute
% =========================================================================

switch platform

    % ---------------------------------------------------------------------
    case 'matlab'
    % Stages run sequentially in the current MATLAB session.
    % ---------------------------------------------------------------------

        run_stage(1, 'Acoustic simulation + free-water baseline', ...
                  @() prestus_pipeline(p_acoustic));

        for ti = 1:numel(targets)
            run_stage(ti + 1, sprintf('Thermal — %.0f W/cm²', targets(ti)), ...
                      @() prestus_pipeline(p_thermal{ti}));
        end

        if options.generate_summary_report
            assert_thermal_outputs_exist(parameters, affixes);
            run_stage(numel(targets) + 2, 'Multi-ISPPA summary report', ...
                      @() prestus_pipeline(p_summary));
        end

    % ---------------------------------------------------------------------
    case {'slurm', 'qsub'}
    % Jobs are submitted immediately with scheduler afterok dependencies.
    % ---------------------------------------------------------------------

        subj       = sprintf('sub-%03d', parameters.subject_id);
        output_dir = get_output_dir(parameters);
        subject_id = parameters.subject_id;
        medium     = parameters.simulation.medium;

        % ── Stage 1: Acoustic ────────────────────────────────────────────
        acoustic_sentinel = fullfile(output_dir, 'cache', ...
            sprintf('sub-%03d_%s_results%s.mat', subject_id, medium, base_affix));
        p_acoustic.hpc.timelimit = options.acoustic_timelimit;
        p_acoustic.hpc.job_name  = sprintf('PRESTUS-%smi1-acoustic_%s', options.job_name_prefix, subj);
        p_acoustic               = hpc_add_dependency(p_acoustic, options.depend_job_ids);
        p_acoustic               = apply_memorylimit(p_acoustic, options.acoustic_memorylimit);
        p_acoustic               = apply_partition(p_acoustic, options.acoustic_partition);

        if isfile(acoustic_sentinel) && ~ow.rerun_intermediates
            fprintf('[Stage 1] Acoustic cache exists — skipping submission.\n');
            job_id_acoustic = [];
        else
            fprintf('[Stage 1] Submitting acoustic job...\n');
            job_id_acoustic = prestus_pipeline_start(p_acoustic);
            fprintf('          → job %s\n', num2str(job_id_acoustic));
        end
        fprintf('\n');

        % ── Stages 2…N: Thermal per target ───────────────────────────────
        job_ids_thermal  = zeros(1, numel(targets));
        thermal_submitted = false(1, numel(targets));
        bb_thermal = isfield(parameters, 'thermal_backend') && ...
            strcmpi(parameters.thermal_backend, 'babelbrain');

        for ti = 1:numel(targets)
            affix    = affixes{ti};
            csv_file = fullfile(output_dir, ...
                sprintf('sub-%03d_%s%s.csv', subject_id, medium, affix));
            stage_label = sprintf('Stage %d', ti + 1);
            fprintf('[%s] Submitting thermal job (%.0f W/cm²)  ', stage_label, targets(ti));

            if isfile(csv_file) && ~ow.rerun_outputs
                fprintf('— output exists, skipping.\n');
            else
                pt = p_thermal{ti};
                pt.hpc.timelimit = options.thermal_timelimit;
                pt.hpc.job_name  = sprintf('PRESTUS-%smi%d-thermal_isppa%03d_%s', ...
                                           options.job_name_prefix, ti + 1, round(targets(ti)), subj);
                pt = apply_memorylimit(pt, options.thermal_memorylimit);
                if ~isempty(job_id_acoustic)
                    pt = hpc_add_dependency(pt, job_id_acoustic);
                else
                    pt = hpc_add_dependency(pt, options.depend_job_ids);
                end
                % BabelBrain thermal writes <acoustic stem>-ThermalField_AllCombinations.h5
                % next to the shared Stage-1 acoustic file and reads it back, whatever
                % the target: parallel jobs would overwrite each other's result.
                if bb_thermal && any(thermal_submitted)
                    pt = hpc_add_dependency(pt, job_ids_thermal(find(thermal_submitted, 1, 'last')));
                end
                job_ids_thermal(ti)   = prestus_pipeline_start(pt);
                thermal_submitted(ti) = true;
                fprintf('→ job %s\n', num2str(job_ids_thermal(ti)));
            end
        end
        fprintf('\n');

        % ── Stage N+1: Summary ───────────────────────────────────────────
        if options.generate_summary_report
            summary_file = fullfile(output_dir, ...
                sprintf('sub-%03d_%s_desc-multiIsppa_report%s.html', subject_id, medium, base_affix));
            fprintf('[Stage %d] Submitting summary report  ', numel(targets) + 2);
            if isfile(summary_file) && ~ow.rerun_outputs
                fprintf('— report exists, skipping.\n');
            else
                p_summary.hpc.timelimit = options.summary_timelimit;
                p_summary.hpc.job_name  = sprintf('PRESTUS-%smi-summary_%s', options.job_name_prefix, subj);
                p_summary               = apply_memorylimit(p_summary, options.summary_memorylimit);
                p_summary               = apply_partition(p_summary, options.summary_partition);
                pending_ids = job_ids_thermal(thermal_submitted);
                if ~isempty(pending_ids)
                    p_summary = hpc_add_dependency(p_summary, pending_ids);
                else
                    p_summary = hpc_add_dependency(p_summary, options.depend_job_ids);
                end
                job_id_summary = prestus_pipeline_start(p_summary);
                fprintf('→ job %s\n', num2str(job_id_summary));
            end
        end
        fprintf('\n');

        % ── Summary ──────────────────────────────────────────────────────
        fprintf('All jobs submitted. Monitor: squeue --user $USER -o "%%.18i %%.9P %%.50j %%.8u %%.2t %%.10M %%.6D %%R"\n');
        all_ids = [job_id_acoustic, job_ids_thermal(thermal_submitted)];
        if options.generate_summary_report && exist('job_id_summary', 'var')
            all_ids = [all_ids, job_id_summary];
        end
        all_ids = all_ids(all_ids > 0);
        jobs.acoustic = job_id_acoustic;
        jobs.thermal  = job_ids_thermal(thermal_submitted);
        if exist('job_id_summary', 'var'); jobs.summary = job_id_summary; end
        jobs.all      = all_ids;
        if ~isempty(all_ids)
            fprintf('Submitted job IDs : %s\n', num2str(all_ids));
            fprintf('To cancel         : scancel %s\n', num2str(all_ids));
        end

    % ---------------------------------------------------------------------
    otherwise
        error('multi_isppa_pipeline:unknownPlatform', ...
            'Unknown platform ''%s''. Expected ''matlab'', ''slurm'', or ''qsub''.', platform);
end

end % multi_isppa_pipeline

% =========================================================================
%% Stage runner (MATLAB platform only)
% =========================================================================

function run_stage(n, label, fn)
    fprintf('────────────────────────────────────────\n');
    fprintf('[Stage %d] %s\n', n, label);
    fprintf('────────────────────────────────────────\n');
    fn();
    fprintf('\n');
end

% =========================================================================
%% Parameter builders
% =========================================================================

function p = clear_multi_isppa_targets(p)
% Replace each transducer's target vector with NaN so that:
%   (a) is_multi_isppa_mode returns false → no recursion on compute nodes
%   (b) the field is present so apply_isppa_scaling can detect it as unset
    if isfield(p, 'transducer')
        for ti = 1:numel(p.transducer)
            p.transducer(ti).target_isppa_wcm2 = NaN;
        end
    end
end

function p = make_acoustic_params(base, base_affix, ow)
% Stage 1: run the full acoustic simulation and free-water baseline.
% Thermal simulation is disabled — each thermal job runs separately.
    p = clear_multi_isppa_targets(base);
    p.modules.run_source_setup       = 1;
    p.modules.run_acoustic_sims      = 1;
    p.modules.run_water_baseline = 1;
    p.modules.run_heating_sims       = 0;
    p.modules.run_thermal_analysis   = 0;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 0;
    p.io.output_affix                = base_affix;
    p.io.overwrite_files             = ow.intermediates;
    p.io.overwrite_simnibs           = 0;
    p.io.save_acoustic_matrices      = 1;   % cache must persist for thermal jobs
    p.io.save_thermal_matrices       = 0;
    p.simulation.interactive         = 0;
end

function p = make_thermal_params(base, target_isppa, affix, base_affix, ow)
% Stages 2…N: load cached acoustic results, apply pressure scaling to
% target_isppa, run thermal simulation and analysis, generate HTML report.
    p = clear_multi_isppa_targets(base);
    p.transducer(1).target_isppa_wcm2 = target_isppa;  % scalar for this job
    p.io.acoustic_cache_affix        = base_affix;     % points at Stage-1 cache file
    p.modules.run_source_setup       = 0;   % kgrid/source/sensor come from acoustic cache
    p.modules.run_acoustic_sims      = 0;   % load from cache
    p.modules.run_water_baseline     = 0;
    p.modules.run_heating_sims       = 1;
    run_thermal = ~isfield(base.modules, 'run_thermal_analysis') || ...
                  base.modules.run_thermal_analysis;
    p.modules.run_thermal_analysis   = run_thermal;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 1;
    p.io.output_affix                = affix;
    % Point at the stage-1 preprocessing cache: the acoustic stage's own affix, or
    % the caches it inherited (an uncertainty variant reuses the shared stage-1 head/grid).
    p.io.preproc_affix               = base_affix;
    if isfield(base, 'io') && isfield(base.io, 'preproc_affix')
        p.io.preproc_affix           = base.io.preproc_affix;
    end
    p.io.overwrite_files             = ow.outputs;
    p.io.reuse_upstream_caches       = {'head', 'grid', 'medium'};
    p.io.overwrite_simnibs           = 0;
    p.io.save_acoustic_matrices      = 0;
    p.io.save_thermal_matrices       = isfield(base.modules, 'run_heating_sims') && ...
                                       base.modules.run_heating_sims;
    p.simulation.interactive         = 0;
end

function p = make_summary_params(base, targets, affixes, base_affix, ow)
% Stage N+1: generate a summary report comparing all ISPPA variants.
    p = clear_multi_isppa_targets(base);
    p.modules.run_source_setup       = 0;
    p.modules.run_acoustic_sims      = 0;
    p.modules.run_water_baseline   = 0;
    p.modules.run_heating_sims       = 0;
    p.modules.run_thermal_analysis   = 0;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 0;
    p.modules.multi_isppa_report     = 1;
    p.io.output_affix                = base_affix;
    p.io.overwrite_files             = ow.outputs;
    p.io.reuse_upstream_caches       = {'head', 'grid', 'medium'};
    p.simulation.interactive         = 0;
    p.multi_isppa.targets_wcm2       = targets;
    p.multi_isppa.affixes            = affixes;
end

% =========================================================================
%% Pre-report validation (MATLAB platform)
% =========================================================================

function assert_thermal_outputs_exist(parameters, affixes)
    % Filename must match the thermal stage's actual output (see the csv_file
    % built in the slurm/qsub branch above, and resolve_io_dirs.m's
    % io.filename_table derivation) — not '_output_table', which no stage
    % here ever writes.
    subject_id = parameters.subject_id;
    medium     = parameters.simulation.medium;
    output_dir = get_output_dir(parameters);
    for ti = 1:numel(affixes)
        csv = fullfile(output_dir, ...
            sprintf('sub-%03d_%s%s.csv', subject_id, medium, affixes{ti}));
        if ~isfile(csv)
            warn('multi_isppa_pipeline:missingOutput', ...
                'Thermal output not found for affix ''%s'' — summary will be incomplete.\n  Expected: %s', ...
                affixes{ti}, csv);
        end
    end
end

% =========================================================================
%% HPC helpers
% =========================================================================

function p = apply_memorylimit(p, limit)
    if ~isempty(limit)
        p.hpc.memorylimit = limit;
    end
end

function p = apply_partition(p, partition)
    if ~isempty(partition)
        p.hpc.partition = partition;
    elseif isfield(p.hpc, 'partition')
        p.hpc = rmfield(p.hpc, 'partition');
    end
end
