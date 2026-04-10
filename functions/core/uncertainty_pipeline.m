function uncertainty_pipeline(parameters, options)
% UNCERTAINTY_PIPELINE  Run the full PRESTUS uncertainty quantification workflow.
%
% Rather than a single simulation with one set of tissue properties, this
% pipeline runs three variants — default, liberal, conservative — and
% combines them into a unified uncertainty HTML report.
%
% Pipeline stages
% ---------------
%   Stage 1  Preprocessing & source setup          (always serial)
%   Stage 2  Default simulation      ─┐
%   Stage 3  Liberal simulation       ├─ parallel on HPC, sequential in MATLAB
%   Stage 4  Conservative simulation ─┘
%   Stage 5  Uncertainty report generation         (after stages 2–4)
%
% Variant definitions
% -------------------
%   Default      Best-estimate medium properties (reference)
%   Liberal      Low skull impedance / attenuation → higher intracranial
%                intensity; represents the worst-case safety scenario
%   Conservative High skull impedance / attenuation → lower intracranial
%                intensity; represents the lower plausible bound
%
% Usage
% -----
%   uncertainty_pipeline(parameters)
%   uncertainty_pipeline(parameters, options)
%
%   Typically called indirectly via prestus_pipeline_start when
%   parameters.simulation.uncertainty = true.
%
% Inputs
% ------
%   parameters   Base PRESTUS parameters struct. Must have subject_id,
%                simulation.medium, path.sim, and all transducer/grid
%                settings configured.
%                Do NOT set io.output_affix — this function manages it.
%
%   options      Optional struct:
%     .affixes             Struct with fields .default / .liberal / .conservative
%                          controlling the filename suffix for each variant.
%                          Defaults: '' / '_liberal' / '_conservative'
%     .liberal_config      Path to liberal medium YAML override file.
%                          Default: <prestus_root>/configs/uncertainty/
%                                   config_medium_liberal.yaml
%     .conservative_config Path to conservative medium YAML override file.
%                          Default: <prestus_root>/configs/uncertainty/
%                                   config_medium_conservative.yaml
%     .stage1_timelimit    HPC wall time for stage 1    (default '12:00:00')
%     .sim_timelimit       HPC wall time for stages 2–4 (default '06:00:00')
%     .report_timelimit    HPC wall time for stage 5    (default '00:30:00')
%     .stage1_memorylimit  RAM in GB for stage 1        (default: 16 GB)
%     .stage1_partition   Partition for stage 1         (default: scheduler default)
%     .report_partition   Partition for stage 5         (default: scheduler default)
%                          Set either to 'gpu' to avoid a busy batch queue.
%     .sim_memorylimit     RAM in GB for stages 2–4     (default: inherit)
%     .report_memorylimit  RAM in GB for stage 5        (default: inherit)
%
% Platform behaviour
% ------------------
%   matlab       Stages run sequentially in the current MATLAB session.
%   slurm/qsub   All five jobs are submitted immediately. Stages 2–4 carry
%                an afterok dependency on stage 1; stage 5 carries afterok
%                dependencies on all three simulation jobs. No MATLAB session
%                blocking is required.
%
% Note: on HPC, prestus_pipeline_start must return the submitted job id as
% its first output argument.

arguments
    parameters   struct
    options      struct = struct()
end

% =========================================================================
%% Resolve options — fill in defaults for any unset fields
% =========================================================================

uncertainty_configs = fullfile(get_prestus_path(), 'configs', 'uncertainty');

if ~isfield(options, 'affixes')
    % If the caller has already set io.output_affix on the base parameters,
    % use it as the default variant's affix so existing naming is preserved.
    % Liberal and conservative affixes are always appended on top of that base.
    base_affix = '';
    if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix') && ...
            ~isempty(parameters.io.output_affix)
        base_affix = parameters.io.output_affix;
    end
    options.affixes = struct( ...
        'default',      base_affix, ...
        'liberal',      [base_affix '_liberal'], ...
        'conservative', [base_affix '_conservative']);
end
if ~isfield(options, 'liberal_config')
    options.liberal_config      = fullfile(uncertainty_configs, 'config_medium_liberal.yaml');
end
if ~isfield(options, 'conservative_config')
    options.conservative_config = fullfile(uncertainty_configs, 'config_medium_conservative.yaml');
end
if ~isfield(options, 'stage1_timelimit');    options.stage1_timelimit    = '12:00:00'; end
if ~isfield(options, 'sim_timelimit');       options.sim_timelimit       = '06:00:00'; end
if ~isfield(options, 'report_timelimit');    options.report_timelimit    = '00:30:00'; end
% Memory limits per stage (GB).
% Simulation and report stages inherit hpc.memorylimit from the base config.
if ~isfield(options, 'stage1_memorylimit'); options.stage1_memorylimit  = []; end  % resolved below after platform detection
if ~isfield(options, 'stage1_partition');   options.stage1_partition   = []; end  % [] = scheduler default; set e.g. 'gpu' to avoid a long batch queue
if ~isfield(options, 'report_partition');   options.report_partition   = []; end  % [] = scheduler default; set e.g. 'gpu' to avoid a long batch queue
if ~isfield(options, 'sim_memorylimit');    options.sim_memorylimit     = []; end
if ~isfield(options, 'report_memorylimit'); options.report_memorylimit  = []; end

% =========================================================================
%% Parameters overridden by this pipeline
% =========================================================================
% The builders below forcibly set the following fields regardless of what
% is in the base parameters.  This is intentional — they are pipeline-
% managed settings that must not be left to per-subject configs:
%
%   io.output_affix     Variant-specific suffix (managed via options.affixes).
%                       The base value is absorbed into options.affixes.default
%                       above so it is not silently discarded.
%
%   io.preproc_affix    Set to '' for simulation variants so that head-
%                       preprocessing file lookups (reoriented/scaled data,
%                       cropped/smoothed skull) point at the stage-1 cache
%                       rather than looking for per-variant files that do not
%                       exist. preproc_head.m honours this field.
%                       NOTE: kwave source files are NOT shared — each variant
%                       computes its own source using io.output_affix so that
%                       the time axis matches the variant's medium sound speed.
%
%   io.overwrite_files  Hardcoded to 'never' so a restarted pipeline resumes
%                       from where it left off rather than re-running finished
%                       variants from scratch.
%
%   io.overwrite_simnibs  Hardcoded to 0 for the same reason.
%
%   io.save_grid_cache         Stage 1: 1 (write cache for stages 2–4)
%                              Stages 2–4: 0 (already written)
%   io.save_source_matrices    Stage 1: 1 (write source for stages 2–4)
%                              Stages 2–4: 0 (already written)
%   io.save_acoustic_matrices  All stages: 0 (large; not needed after analysis)
%   io.save_thermal_matrices   Stages 2–4: inherited from base modules.run_heating_sims
%                              Stage 1: 0 (no thermal sim)
%
%   simulation.interactive  Hardcoded to 0; HPC jobs cannot show dialogs.
%
%   modules.*           All module flags are set explicitly by each builder —
%                       this is the whole point of the pipeline, so user-set
%                       module flags in the base parameters are ignored here.
%
%   hpc.timelimit       Set per stage from options.*_timelimit.
%
% Parameters that ARE inherited unchanged from the base:
%   transducer, grid, path, hpc (except timelimit/job_name/depend*),
%   simulation.medium, simulation.code_type, pct, startup, subject_id.

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

fprintf('\n');
fprintf('╔══════════════════════════════════════╗\n');
fprintf('║   PRESTUS Uncertainty Pipeline       ║\n');
fprintf('╠══════════════════════════════════════╣\n');
fprintf('║  Subject : sub-%03d                  ║\n', parameters.subject_id);
fprintf('║  Medium  : %-26s ║\n',                     parameters.simulation.medium);
fprintf('║  Platform: %-26s ║\n',                     upper(platform));
fprintf('╚══════════════════════════════════════╝\n\n');

% =========================================================================
%% Build one parameter struct per stage
%
%   Each builder calls clear_uncertainty_flag so that when a struct is
%   serialised to a temp .mat file for HPC submission, the flag is absent
%   and prestus_pipeline runs normally on the compute node (no recursion).
% =========================================================================

p_stage1       = make_stage1_params(parameters);

p_default      = make_sim_params(parameters, ...
                     options.affixes.default, []);
p_liberal      = make_sim_params(parameters, ...
                     options.affixes.liberal, options.liberal_config);
p_conservative = make_sim_params(parameters, ...
                     options.affixes.conservative, options.conservative_config);

% Pre-assign deterministic log file paths for all five stages so that
% generate_uncertainty_report can locate each stage's timing data.
% Using a single pipeline-run timestamp keeps file names stable across
% MATLAB and HPC (where the paths must be encoded before job submission).
run_ts     = string(datetime('now'), 'yyMMdd_HHmm');
output_dir = get_output_dir(parameters);
subj       = sprintf('sub-%03d', parameters.subject_id);
medium     = parameters.simulation.medium;

log_files.stage1      = fullfile(output_dir, sprintf('%s_%s_stage1_%s.txt',       subj, medium, run_ts));
log_files.default     = fullfile(output_dir, sprintf('%s_%s%s_%s.txt',            subj, medium, options.affixes.default,      run_ts));
log_files.liberal     = fullfile(output_dir, sprintf('%s_%s%s_%s.txt',            subj, medium, options.affixes.liberal,      run_ts));
log_files.conservative= fullfile(output_dir, sprintf('%s_%s%s_%s.txt',            subj, medium, options.affixes.conservative, run_ts));
log_files.report      = fullfile(output_dir, sprintf('%s_%s_stage5_report_%s.txt',subj, medium, run_ts));

p_stage1.io.log_file       = log_files.stage1;
p_default.io.log_file      = log_files.default;
p_liberal.io.log_file      = log_files.liberal;
p_conservative.io.log_file = log_files.conservative;

p_report       = make_report_params(parameters, options.affixes, log_files);

% =========================================================================
%% Execute
% =========================================================================

switch platform

    % ---------------------------------------------------------------------
    case 'matlab'
    % All five stages run sequentially in the current MATLAB session.
    % ---------------------------------------------------------------------

        run_stage(1, 'Preprocessing & source setup', @() prestus_pipeline(p_stage1));
        run_stage(2, 'Default simulation',           @() prestus_pipeline(p_default));
        run_stage(3, 'Liberal simulation',           @() prestus_pipeline(p_liberal));
        run_stage(4, 'Conservative simulation',      @() prestus_pipeline(p_conservative));

        % Warn about any missing variant outputs before attempting report
        assert_variant_outputs_exist(parameters, options.affixes);
        run_stage(5, 'Uncertainty report',           @() prestus_pipeline(p_report));

        % Clean up intermediate files when save_matrices = 0
        if ~should_save_output(parameters.io, 'save_matrices')
            cleanup_uncertainty_intermediates(parameters, options.affixes);
        end

    % ---------------------------------------------------------------------
    case {'slurm', 'qsub'}
    % Jobs are submitted immediately with scheduler afterok dependencies.
    % No MATLAB session blocking is required.
    % ---------------------------------------------------------------------

        subj       = sprintf('sub-%03d', parameters.subject_id);
        output_dir = get_output_dir(parameters);
        subject_id = parameters.subject_id;
        medium     = parameters.simulation.medium;

        % ── Stage 1 ──────────────────────────────────────────────────────
        % Skip if the head preprocessing cache already exists.
        % Sentinel: cache/sub-NNN_<medium>_after_cropping_and_smoothing.mat
        % (written by preproc_head with preproc_affix = '').
        % When skipped, stages 2–4 are submitted without an afterok
        % dependency so they can start immediately.
        preproc_sentinel = fullfile(output_dir, 'cache', ...
            sprintf('sub-%03d_%s_after_cropping_and_smoothing.mat', subject_id, medium));
        p_stage1.hpc.timelimit = options.stage1_timelimit;
        p_stage1.hpc.job_name  = sprintf('PRESTUS-u1-preproc_%s', subj);
        p_stage1               = apply_memorylimit(p_stage1, options.stage1_memorylimit);
        p_stage1               = strip_gpu_requirements(p_stage1, options.stage1_partition);
        if isfile(preproc_sentinel)
            fprintf('[Stage 1] Preprocessing cache exists — skipping submission.\n');
            job_id_stage1 = [];
        else
            fprintf('[Stage 1] Submitting preprocessing job...\n');
            job_id_stage1 = prestus_pipeline_start(p_stage1);
            fprintf('          → job %s  (mem %iG, time %s)\n', ...
                num2str(job_id_stage1), p_stage1.hpc.memorylimit, p_stage1.hpc.timelimit);
        end
        fprintf('\n');

        % ── Stages 2–4 ───────────────────────────────────────────────────
        % Each variant is skipped if its output CSV already exists.
        % If not skipped, it depends on stage 1 (when stage 1 was submitted).
        sim_stages = { ...
            'Stage 2', 'default',      p_default,      options.affixes.default;      ...
            'Stage 3', 'liberal',      p_liberal,      options.affixes.liberal;       ...
            'Stage 4', 'conservative', p_conservative, options.affixes.conservative   ...
        };
        job_names = { ...
            sprintf('PRESTUS-u2-sim-default_%s',      subj); ...
            sprintf('PRESTUS-u3-sim-liberal_%s',       subj); ...
            sprintf('PRESTUS-u4-sim-conservative_%s',  subj)  ...
        };
        job_ids_sim = zeros(1, 3);
        sim_submitted = false(1, 3);
        for si = 1:3
            stage_label = sim_stages{si, 1};
            variant     = sim_stages{si, 2};
            p_sim       = sim_stages{si, 3};
            affix       = sim_stages{si, 4};
            csv_file    = fullfile(output_dir, ...
                sprintf('sub-%03d_%s_output_table%s.csv', subject_id, medium, affix));
            fprintf('[%s] Submitting %s simulation    ', stage_label, variant);
            if isfile(csv_file)
                fprintf('— output exists, skipping.\n');
            else
                p_sim.hpc.timelimit = options.sim_timelimit;
                p_sim.hpc.job_name  = job_names{si};
                p_sim               = apply_memorylimit(p_sim, options.sim_memorylimit);
                if ~isempty(job_id_stage1)
                    p_sim.hpc.depend_job_id = job_id_stage1;
                end
                job_ids_sim(si)     = prestus_pipeline_start(p_sim);
                sim_submitted(si)   = true;
                fprintf('→ job %s  (mem %iG)\n', ...
                    num2str(job_ids_sim(si)), p_sim.hpc.memorylimit);
            end
        end
        fprintf('\n');
        job_id_default      = job_ids_sim(1);
        job_id_liberal      = job_ids_sim(2);
        job_id_conservative = job_ids_sim(3);

        % ── Stage 5 ──────────────────────────────────────────────────────
        % Skip if the report already exists.
        % Depends on whichever simulation jobs were actually submitted.
        report_file = fullfile(output_dir, ...
            sprintf('sub-%03d_%s_uncertainty_report.html', subject_id, medium));
        fprintf('[Stage 5] Submitting uncertainty report      ');
        if isfile(report_file)
            fprintf('— report exists, skipping.\n');
            job_id_report = [];
        else
            p_report.hpc.timelimit  = options.report_timelimit;
            p_report.hpc.job_name   = sprintf('PRESTUS-u5-report_%s', subj);
            p_report                = apply_memorylimit(p_report, options.report_memorylimit);
            p_report                = strip_gpu_requirements(p_report, options.report_partition);
            pending_sim_ids = job_ids_sim(sim_submitted);
            if ~isempty(pending_sim_ids)
                p_report.hpc.depend_job_ids = pending_sim_ids;
            end
            job_id_report = prestus_pipeline_start(p_report);
            fprintf('→ job %s  (mem %iG)\n', num2str(job_id_report), p_report.hpc.memorylimit);
        end
        fprintf('\n');

        % ── Summary ──────────────────────────────────────────────────────
        fprintf('All jobs submitted. Monitor: squeue -u $USER\n');
        submitted_ids = [job_id_stage1, job_ids_sim(sim_submitted), job_id_report];
        submitted_ids = submitted_ids(submitted_ids > 0);
        if ~isempty(submitted_ids)
            fprintf('Submitted job IDs : %s\n', num2str(submitted_ids));
            fprintf('To cancel         : scancel %s\n', num2str(submitted_ids));
        end
        fprintf('Logs              : %s/%s/batch_job_logs/\n', parameters.path.sim, subj);

    % ---------------------------------------------------------------------
    otherwise
        error('uncertainty_pipeline:unknownPlatform', ...
            'Unknown platform ''%s''. Expected ''matlab'', ''slurm'', or ''qsub''.', platform);
end

end % uncertainty_pipeline

% =========================================================================
%% Stage runner (MATLAB platform only)
% =========================================================================

function run_stage(n, label, fn)
% Print a consistent stage header and run fn().
    fprintf('────────────────────────────────────────\n');
    fprintf('[Stage %d] %s\n', n, label);
    fprintf('────────────────────────────────────────\n');
    fn();
    fprintf('\n');
end

% =========================================================================
%% Parameter builders
% =========================================================================

function p = clear_uncertainty_flag(p)
% Prevent variant parameter structs from re-triggering uncertainty mode.
% When prestus_pipeline_start serialises p to a .mat file for HPC, this
% ensures the compute node runs the normal pipeline, not another uncertainty
% pipeline (which would cause infinite job submission).
    if isfield(p, 'simulation') && isfield(p.simulation, 'uncertainty')
        p.simulation.uncertainty = false;
    end
end

function p = make_stage1_params(base)
% Stage 1: preprocessing and source setup only — no simulation.
%
% Produces the shared head-preprocessing cache used by all three simulation
% variants (skull, tissue masks, etc.). Each variant computes its own kwave
% source matrix (via io.output_affix) because the time axis depends on the
% variant's medium sound speed.
%
% GPU is not required: no k-Wave simulation runs in this stage.
    p = clear_uncertainty_flag(base);
    p.modules.run_source_setup         = 1;
    p.modules.run_acoustic_sims        = 0;
    p.modules.run_heating_sims         = 0;
    p.modules.run_thermal_analysis     = 0;
    p.modules.run_posthoc_water_sims   = 0;
    p.modules.generate_report          = 0;
    p.io.output_affix                  = '';
    p.io.overwrite_files               = 'never';
    p.io.overwrite_simnibs             = 0;
    p.io.save_source_matrices          = 0;   % source is recomputed per variant; no need to cache here
    p.io.save_acoustic_matrices        = 0;   % no acoustic sim in stage 1
    p.io.save_thermal_matrices         = 0;   % no thermal sim in stage 1
    p.simulation.interactive           = 0;
    p                                  = strip_gpu_requirements(p);
end

function p = make_sim_params(base, affix, medium_config)
% Stages 2–4: full acoustic + thermal simulation.
%
% If medium_config is a valid YAML path, its properties are merged on top
% of the base parameters (only the overridden fields change; all transducer,
% grid, and path settings are inherited).
%
% run_source_setup = 1 is intentional: each variant checks for cached source
% files from stage 1 and skips recomputation (overwrite_files = 'never').
    p = clear_uncertainty_flag(base);

    if ~isempty(medium_config)
        if ~(ischar(medium_config) || isstring(medium_config)) || ~isfile(medium_config)
            error('uncertainty_pipeline:configNotFound', ...
                'Medium config file not found: %s', medium_config);
        end
        extra = yaml.loadFile(char(medium_config), 'ConvertToArray', true);
        p     = MergeStruct(p, extra);
    end

    p.modules.run_source_setup         = 1;
    p.modules.run_acoustic_sims        = 1;
    % Thermal steps: inherit from base config rather than forcing on.
    % run_heating_sims = 0 in the base will produce acoustic-only uncertainty
    % variants, which is valid when thermal estimation is not requested.
    run_thermal = isfield(base.modules, 'run_heating_sims') && base.modules.run_heating_sims;
    p.modules.run_heating_sims         = run_thermal;
    p.modules.run_thermal_analysis     = run_thermal && ...
                                         isfield(base.modules, 'run_thermal_analysis') && ...
                                         base.modules.run_thermal_analysis;
    p.modules.run_posthoc_water_sims   = 0;
    p.modules.generate_report          = 1;
    p.io.output_affix                  = affix;
    % Both the grid cache and kwave source were written by stage 1 with
    % affix ''. Point lookups at those files rather than recomputing.
    p.io.preproc_affix                 = '';
    p.io.overwrite_files               = 'never';
    p.io.overwrite_simnibs             = 0;
    p.io.save_source_matrices          = 0;   % source is recomputed at the start of each variant; no benefit to caching
    p.io.save_acoustic_matrices        = 0;   % large; not needed after analysis
    p.io.save_thermal_matrices         = run_thermal;  % only needed when thermal ran
    p.simulation.interactive           = 0;
end

function p = make_report_params(base, affixes, log_files)
% Stage 5: generate the uncertainty report only — no simulation.
%
% The affixes struct is passed through parameters.uncertainty.affixes so
% that generate_uncertainty_report can locate the three variant output files.
% log_files is a struct with fields stage1, default, liberal, conservative,
% report — pre-assigned paths that allow the report to parse per-stage timing.
    p = clear_uncertainty_flag(base);
    p.modules.run_source_setup       = 0;
    p.modules.run_acoustic_sims      = 0;
    p.modules.run_heating_sims       = 0;
    p.modules.run_thermal_analysis   = 0;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 0;
    p.modules.uncertainty_report     = 1;
    p.io.output_affix                = '';
    p.io.log_file                    = log_files.report;
    p.io.overwrite_files             = 'never';
    p.simulation.interactive         = 0;
    p.uncertainty.affixes            = affixes;
    p.uncertainty.log_files          = log_files;
end

% =========================================================================
%% Pre-report validation
% =========================================================================

function p = apply_memorylimit(p, limit)
% Override hpc.memorylimit for this stage if an explicit limit is provided.
% When limit is empty, the value already in parameters.hpc.memorylimit is
% kept unchanged (i.e. inherited from the base config).
    if ~isempty(limit)
        p.hpc.memorylimit = limit;
    end
end

function p = strip_gpu_requirements(p, partition_override)
% Clear GPU compute requirements (gres, n_gpu) and force CPU code path.
% The partition is left unchanged unless partition_override is provided,
% allowing CPU-only stages to still be submitted to the GPU partition
% (e.g., to avoid a long batch queue).
    if nargin < 2; partition_override = []; end
    gpu_fields = {'gpu', 'n_gpu'};
    for f = gpu_fields
        if isfield(p.hpc, f{1})
            p.hpc = rmfield(p.hpc, f{1});
        end
    end
    if ~isempty(partition_override)
        p.hpc.partition = partition_override;
    else
        % Remove partition so scheduler uses its default (typically batch)
        if isfield(p.hpc, 'partition')
            p.hpc = rmfield(p.hpc, 'partition');
        end
    end
    % Force CPU code path so write_slurm_script does not re-add GPU directives
    p.simulation.code_type = 'matlab_cpu';
end

function assert_variant_outputs_exist(parameters, affixes)
% Warn (not error) for each simulation variant whose output CSV is missing.
% Called before stage 5 on the MATLAB platform so the user knows in advance
% that the report will be incomplete rather than discovering it after the
% fact inside generate_uncertainty_report.
    subject_id = parameters.subject_id;
    medium     = parameters.simulation.medium;
    output_dir = get_output_dir(parameters);

    variants = { ...
        'Conservative', affixes.conservative; ...
        'Default',      affixes.default;      ...
        'Liberal',      affixes.liberal        ...
    };

    for v = 1:size(variants, 1)
        label  = variants{v, 1};
        affix  = variants{v, 2};
        csv    = fullfile(output_dir, ...
            sprintf('sub-%03d_%s_output_table%s.csv', subject_id, medium, affix));
        if ~isfile(csv)
            warning('uncertainty_pipeline:missingOutput', ...
                '%s variant output not found — report will be incomplete.\n  Expected: %s', ...
                label, csv);
        end
    end
end
