function multi_isppa_pipeline(parameters, options)
% MULTI_ISPPA_PIPELINE  Run acoustic simulation once, then thermal analysis at
%                       multiple free-water target intensities in parallel
%
% When calibration.target_isppa_wcm2 contains more than one value, the
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
%               Output files carry an affix _isppaNNN (NNN = integer W/cm²).
%
%   Stage N+1   (optional) Summary report aggregating all ISPPA variants.
%               Active when options.generate_summary_report = true (default).
%
% Platform behaviour:
%   matlab     Stages run sequentially in the current MATLAB session.
%   slurm/qsub Stage 1 is submitted; stages 2…N carry an afterok dependency
%              on stage 1; stage N+1 depends on all thermal jobs.
%              All jobs are submitted immediately without blocking.
%
% Use as:
%   multi_isppa_pipeline(parameters)
%   multi_isppa_pipeline(parameters, options)
%
% Input:
%   parameters - base PRESTUS config; calibration.target_isppa_wcm2 must be
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

if ~isfield(parameters, 'calibration') || ...
        ~isfield(parameters.calibration, 'target_isppa_wcm2') || ...
        numel(parameters.calibration.target_isppa_wcm2) < 2
    error('multi_isppa_pipeline: calibration.target_isppa_wcm2 must contain at least two values.');
end

targets = parameters.calibration.target_isppa_wcm2(:)';

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

p_acoustic = make_acoustic_params(parameters, base_affix);

p_thermal  = cell(1, numel(targets));
affixes    = cell(1, numel(targets));
for ti = 1:numel(targets)
    affixes{ti}   = sprintf('%s_isppa%03d', base_affix, round(targets(ti)));
    p_thermal{ti} = make_thermal_params(parameters, targets(ti), affixes{ti}, base_affix);
end

p_summary = make_summary_params(parameters, targets, affixes, base_affix);

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
        p_acoustic.hpc.job_name  = sprintf('PRESTUS-mi1-acoustic_%s', subj);
        p_acoustic               = apply_memorylimit(p_acoustic, options.acoustic_memorylimit);
        p_acoustic               = apply_partition(p_acoustic, options.acoustic_partition);

        if isfile(acoustic_sentinel)
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

        for ti = 1:numel(targets)
            affix    = affixes{ti};
            csv_file = fullfile(output_dir, ...
                sprintf('sub-%03d_%s_output_table%s.csv', subject_id, medium, affix));
            stage_label = sprintf('Stage %d', ti + 1);
            fprintf('[%s] Submitting thermal job (%.0f W/cm²)  ', stage_label, targets(ti));

            if isfile(csv_file)
                fprintf('— output exists, skipping.\n');
            else
                pt = p_thermal{ti};
                pt.hpc.timelimit = options.thermal_timelimit;
                pt.hpc.job_name  = sprintf('PRESTUS-mi%d-thermal_isppa%03d_%s', ...
                                           ti + 1, round(targets(ti)), subj);
                pt = apply_memorylimit(pt, options.thermal_memorylimit);
                if ~isempty(job_id_acoustic)
                    pt.hpc.depend_job_id = job_id_acoustic;
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
                sprintf('sub-%03d_%s_multi_isppa_report%s.html', subject_id, medium, base_affix));
            fprintf('[Stage %d] Submitting summary report  ', numel(targets) + 2);
            if isfile(summary_file)
                fprintf('— report exists, skipping.\n');
            else
                p_summary.hpc.timelimit = options.summary_timelimit;
                p_summary.hpc.job_name  = sprintf('PRESTUS-mi-summary_%s', subj);
                p_summary               = apply_memorylimit(p_summary, options.summary_memorylimit);
                p_summary               = apply_partition(p_summary, options.summary_partition);
                pending_ids = job_ids_thermal(thermal_submitted);
                if ~isempty(pending_ids)
                    p_summary.hpc.depend_job_ids = pending_ids;
                end
                job_id_summary = prestus_pipeline_start(p_summary);
                fprintf('→ job %s\n', num2str(job_id_summary));
            end
        end
        fprintf('\n');

        % ── Summary ──────────────────────────────────────────────────────
        fprintf('All jobs submitted. Monitor: squeue -u $USER\n');
        all_ids = [job_id_acoustic, job_ids_thermal(thermal_submitted)];
        if options.generate_summary_report && exist('job_id_summary', 'var')
            all_ids = [all_ids, job_id_summary];
        end
        all_ids = all_ids(all_ids > 0);
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
% Replace the vector of targets with a scalar NaN so that:
%   (a) is_multi_isppa_mode returns false → no recursion on compute nodes
%   (b) the target field is present so apply_isppa_scaling can check it
    if isfield(p, 'calibration') && isfield(p.calibration, 'target_isppa_wcm2')
        p.calibration.target_isppa_wcm2 = NaN;
    end
end

function p = make_acoustic_params(base, base_affix)
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
    p.io.overwrite_files             = 'never';
    p.io.overwrite_simnibs           = 0;
    p.io.save_acoustic_matrices      = 1;   % cache must persist for thermal jobs
    p.io.save_thermal_matrices       = 0;
    p.simulation.interactive         = 0;
end

function p = make_thermal_params(base, target_isppa, affix, base_affix)
% Stages 2…N: load cached acoustic results, apply pressure scaling to
% target_isppa, run thermal simulation and analysis, generate HTML report.
    p = clear_multi_isppa_targets(base);
    p.calibration.target_isppa_wcm2  = target_isppa;   % scalar for this job
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
    p.io.preproc_affix               = '';   % point at stage-1 preprocessing cache
    p.io.overwrite_files             = 'never';
    p.io.overwrite_simnibs           = 0;
    p.io.save_acoustic_matrices      = 0;
    p.io.save_thermal_matrices       = isfield(base.modules, 'run_heating_sims') && ...
                                       base.modules.run_heating_sims;
    p.simulation.interactive         = 0;
end

function p = make_summary_params(base, targets, affixes, base_affix)
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
    p.io.overwrite_files             = 'never';
    p.simulation.interactive         = 0;
    p.multi_isppa.targets_wcm2       = targets;
    p.multi_isppa.affixes            = affixes;
end

% =========================================================================
%% Pre-report validation (MATLAB platform)
% =========================================================================

function assert_thermal_outputs_exist(parameters, affixes)
    subject_id = parameters.subject_id;
    medium     = parameters.simulation.medium;
    output_dir = get_output_dir(parameters);
    for ti = 1:numel(affixes)
        csv = fullfile(output_dir, ...
            sprintf('sub-%03d_%s_output_table%s.csv', subject_id, medium, affixes{ti}));
        if ~isfile(csv)
            warning('multi_isppa_pipeline:missingOutput', ...
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
