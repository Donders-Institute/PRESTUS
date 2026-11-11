function jobs = uncertainty_pipeline(parameters, options)
% UNCERTAINTY_PIPELINE  Run the full PRESTUS uncertainty quantification workflow
%
% Rather than a single simulation with one set of tissue properties, this
% pipeline runs five variants (UNCERTAINTY_VARIANTS: default, Temp-H, Temp-L,
% Trans-H, Trans-L) and combines them into a unified uncertainty HTML report.
%
% Pipeline stages:
%   Stage 1  Preprocessing & source setup          (always serial)
%   Stage 2  Default simulation      ─┐
%   Stage 3  Temp-H simulation        │
%   Stage 4  Temp-L simulation        ├─ parallel on HPC, sequential in MATLAB
%   Stage 5  Trans-H simulation       │  (Trans-H/L: acoustic only)
%   Stage 6  Trans-L simulation      ─┘
%   Stage 7  Uncertainty report generation         (after stages 2–6)
%
% Nested chained modes: when the base config is itself a multi-ISPPA sweep
% (vector transducer.target_isppa_wcm2) or an async multi-transducer run
% (simulation.transducer_coupling: async), each variant runs that pipeline
% (MULTI_ISPPA_PIPELINE / ASYNC_TRANSDUCER_PIPELINE) with the variant's medium,
% and Stage 7 writes one uncertainty report per sweep point. Uncertainty is
% always the outer mode. A multi-ISPPA variant reuses the Stage-1 head/grid
% caches; an async variant builds its own (async Stage 0), so Stage 1 is not
% run. Inner-pipeline options (e.g. .thermal_timelimit) pass through options.
% Sequential follow-up runs are not supported in nested mode.
%
% Variant definitions (ITRUSST tissue-property bounds, config/uncertainty/):
%   Default      PRESTUS default medium properties; with pct.enabled, the
%                skull is mapped from the pseudo-CT
%   Temp-H/L     Per property, the bound that maximises / minimises brain
%                heating; acoustic + thermal
%   Trans-H/L    Per acoustic property, the bound that maximises / minimises
%                brain intensity; acoustic only
%   The pairs differ because brain sound speed, density and attenuation and
%   skull sound speed raise heating and intensity in opposite directions.
%   With pct.enabled, the four bound variants share the pseudo-CT run's head
%   geometry but give the skull layer the homogeneous values of their file.
%
% Platform behaviour:
%   matlab       Stages run sequentially in the current MATLAB session.
%   slurm/qsub   All seven jobs are submitted immediately. Stages 2–6 carry
%                an afterok dependency on stage 1; stage 7 carries afterok
%                dependencies on all five simulation jobs. In nested mode
%                the variants' inner job graphs are submitted from here, and
%                each per-point report waits for every job they submitted.
%                prestus_pipeline_start must return the submitted job id as
%                its first output argument.
%
% Use as:
%   uncertainty_pipeline(parameters)
%   uncertainty_pipeline(parameters, options)
%
% Input:
%   parameters - base PRESTUS config; must have subject_id, simulation.medium,
%                path.sim, and transducer/grid settings;
%                do NOT set io.output_affix — this function manages it
%   options    - pipeline knobs (optional, default: struct()):
%                .affixes             — struct with one field per variant name
%                                       (default, temp_high, temp_low, trans_high,
%                                       trans_low): filename suffixes per variant
%                .<name>_config       — path to the medium YAML of variant <name>
%                                       (e.g. .temp_high_config; default:
%                                       simulation.uncertainty_<name>, else the
%                                       shipped file; see UNCERTAINTY_MEDIUM_CONFIG)
%                .stage1_timelimit    — HPC wall time for stage 1 (default '12:00:00')
%                .sim_timelimit       — HPC wall time for stages 2-6 (default '06:00:00')
%                .report_timelimit    — HPC wall time for stage 7 (default '00:30:00')
%                .stage1_memorylimit  — RAM in GB for stage 1 (default: scheduler default)
%                .stage1_partition    — partition for stage 1 (default: scheduler default)
%                .report_partition    — partition for the report stage (default: scheduler default)
%                .sim_memorylimit     — RAM in GB for stages 2-6 (default: inherit)
%                .report_memorylimit  — RAM in GB for stage 7 (default: inherit)
%                .sequential_configs  — struct of follow-up run parameter structs
%                                       (config_2, config_3, …) for sequential
%                                       thermal continuation with matched uncertainty
%                                       variants (the thermal outputs of the default
%                                       and Temp variants are passed to the same
%                                       variant of the next run; Trans variants
%                                       re-run acoustic only)
%
% Output:
%   jobs       - struct with .all: every scheduler job ID submitted (stage 1,
%                variants or their inner jobs, reports, sequential runs);
%                [] when run in MATLAB. Lets a caller make later jobs depend
%                on this run.
%
% See also: PRESTUS_PIPELINE_START, PRESTUS_PIPELINE, GENERATE_UNCERTAINTY_REPORT,
%           UNCERTAINTY_VARIANTS

arguments
    parameters   struct
    options      struct = struct()
end

% =========================================================================
%% Resolve options — fill in defaults for any unset fields
% =========================================================================

variants = uncertainty_variants();
reject_legacy_options(options);
if ~isfield(options, 'affixes')
    % If the caller has already set io.output_affix on the base parameters,
    % use it as the default variant's affix so existing naming is preserved.
    % The other variants' descriptors are always appended on top of that base.
    base_affix = '';
    if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix') && ...
            ~isempty(parameters.io.output_affix)
        base_affix = parameters.io.output_affix;
    end
    options.affixes = struct();
    for v = 1:numel(variants)
        options.affixes.(variants(v).name) = [base_affix variants(v).desc];
    end
end
% Variant medium files: options.* > simulation.uncertainty_<variant> > shipped file.
for v = 2:numel(variants)
    key = [variants(v).name '_config'];
    if ~isfield(options, key)
        options.(key) = uncertainty_medium_config(parameters, variants(v).name);
    end
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
%   io.overwrite_files  Mapped per stage by RESOLVE_CHAINED_OVERWRITE:
%                       'never'  stage 1 resumes from its preprocessing cache;
%                                variants/reports are skipped when their
%                                outputs exist (a restarted pipeline resumes).
%                       'always' stage 1 resumes; variants/reports always
%                                re-run and overwrite their outputs.
%                       'force'  every stage re-runs.
%                       Variants always load the stage-1 head/grid caches
%                       (io.reuse_upstream_caches).
%
%   io.overwrite_simnibs  Hardcoded to 0 for the same reason.
%
%   io.save_grid_cache         Stage 1: 1 (write cache for stages 2–6)
%                              Stages 2–6: 0 (already written)
%   io.save_source_matrices    Stage 1: 1 (write source for stages 2–6)
%                              Stages 2–6: 0 (already written)
%   io.save_acoustic_matrices  All stages: 0 (large; not needed after analysis)
%   io.save_thermal_matrices   Default and Temp variants: inherited from base
%                              modules.run_heating_sims; Trans variants and
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

% Nested chained mode: '' | 'multi_isppa' | 'async' (see header).
nested = nested_mode(parameters);
if strcmp(nested, 'multi_isppa') && numel(parameters.transducer) > 1
    % Fail before any job is submitted (MULTI_ISPPA_PIPELINE would reject it later).
    error('prestus:isppa:coherentSweep', ...
        ['ISPPA sweeps are not available for coherent multi-transducer runs, with or ' ...
         'without uncertainty. Run one config per point, or use ' ...
         'simulation.transducer_coupling: async.']);
end
if ~isempty(nested)
    fprintf('Nested mode: each variant runs a %s pipeline.\n\n', strrep(nested, '_', '-'));
end

% =========================================================================
%% Build one parameter struct per stage
%
%   Each builder calls clear_uncertainty_flag so that when a struct is
%   serialised to a temp .mat file for HPC submission, the flag is absent
%   and prestus_pipeline runs normally on the compute node (no recursion).
% =========================================================================

ow = resolve_chained_overwrite(parameters);
jobs = struct('all', []);

p_stage1   = make_stage1_params(parameters, ow);

n_var      = numel(variants);
p_variants = cell(1, n_var);
for v = 1:n_var
    p_variants{v} = make_sim_params(parameters, options.affixes.(variants(v).name), ...
                        get_variant_medium_config(variants(v).name, options), ...
                        variants(v).heating, ow);
end

% Pre-assign deterministic log file paths for all seven stages so that
% generate_uncertainty_report can locate each stage's timing data.
% Using a single pipeline-run timestamp keeps file names stable across
% MATLAB and HPC (where the paths must be encoded before job submission).
run_ts     = string(datetime('now'), 'yyMMdd_HHmm');
output_dir = get_output_dir(parameters);
subj       = sprintf('sub-%03d', parameters.subject_id);
medium     = parameters.simulation.medium;

logs_dir = fullfile(output_dir, 'log');
if ~isfolder(logs_dir); mkdir(logs_dir); end
log_files.stage1 = fullfile(logs_dir, sprintf('%s_%s_stage1_%s.txt', subj, medium, run_ts));
for v = 1:n_var
    log_files.(variants(v).name) = fullfile(logs_dir, sprintf('%s_%s%s_%s.txt', ...
        subj, medium, options.affixes.(variants(v).name), run_ts));
end
log_files.report = fullfile(logs_dir, sprintf('%s_%s_desc-report_%s.txt', subj, medium, run_ts));

p_stage1.io.log_file = log_files.stage1;
if isempty(nested)
    % Nested variants are job graphs: each inner job writes its own log.
    for v = 1:n_var
        p_variants{v}.io.log_file = log_files.(variants(v).name);
    end
end

p_report       = make_report_params(parameters, options.affixes, log_files, ow, options);

% =========================================================================
%% Sequential config pre-processing
% =========================================================================

has_sequential = isfield(options, 'sequential_configs') && ...
                 isstruct(options.sequential_configs) && ...
                 ~isempty(fieldnames(options.sequential_configs));
if has_sequential && ~isempty(nested)
    error('uncertainty_pipeline:nestedSequential', ...
        ['Sequential follow-up runs are not supported when the uncertainty variants ' ...
         'run a %s pipeline. Run the sequential chain with a scalar target and ' ...
         'coherent coupling, or drop options.sequential_configs.'], strrep(nested, '_', '-'));
end

% run_index_start: the run number assigned to the FIRST sequential follow-up.
% Run 1 is the base uncertainty run, so sequential runs start at 2.
run_index_start = 2;

% =========================================================================
%% Execute
% =========================================================================

switch platform

    % ---------------------------------------------------------------------
    case 'matlab'
    % All seven stages run sequentially in the current MATLAB session.
    % ---------------------------------------------------------------------

        if ~isempty(nested)
            run_nested_matlab(parameters, options, nested, platform, ...
                p_variants, p_stage1, log_files, ow);
            return
        end

        run_stage(1, 'Preprocessing & source setup', @() prestus_pipeline(p_stage1));
        for v = 1:n_var
            run_stage(v + 1, variant_stage_label(variants(v)), ...
                      @() prestus_pipeline(p_variants{v}));
        end

        % Warn about any missing variant outputs before attempting report
        assert_variant_outputs_exist(parameters, options.affixes);
        run_stage(n_var + 2, 'Uncertainty report', @() prestus_pipeline(p_report));

        % Clean up intermediate files when save_matrices = 0
        if ~should_save_output(parameters.io, 'save_matrices')
            cleanup_uncertainty_intermediates(parameters, options.affixes);
        end

        % Sequential runs: each follow-up run repeats the uncertainty
        % variants; the default and Temp variants inherit the thermal
        % end-state from the matched variant of the prior run.
        if has_sequential
            seq_run_params = run_sequential_matlab(parameters, options, run_index_start);
            % Generate combined sequential uncertainty report.
            % Prepend run-1 params so the base run appears as the first entry.
            % Its bands span the Temp pair (the Trans variants have no heating).
            all_runs = struct();
            for v = 1:n_var
                all_runs.(variants(v).name) = [p_variants(v), seq_run_params.(variants(v).name)];
            end
            try
                generate_sequential_report( ...
                    all_runs.default, {}, ...
                    struct('high_params_list', {all_runs.temp_high}, ...
                           'low_params_list',  {all_runs.temp_low}, ...
                           'high_label', 'Temp-H', 'low_label', 'Temp-L'));
            catch ME_seq_rep
                warning('uncertainty_pipeline:seqReport', ...
                    'Sequential uncertainty report failed: %s', ME_seq_rep.message);
            end
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

        % Run prefix: 'r1-' when sequential configs are present, '' otherwise.
        % This keeps job names unambiguous when multiple runs are chained
        % (r1-u2-default, r2-u2-default, …); single-run submissions are
        % named u1-preproc, u2-sim-default, …
        has_sequential = isfield(options, 'sequential_configs') && ...
                         ~isempty(fieldnames(options.sequential_configs));
        if has_sequential
            run_prefix = 'r1-';
        else
            run_prefix = '';
        end

        if ~isempty(nested)
            jobs.all = submit_nested_hpc(parameters, options, nested, platform, ...
                p_variants, p_stage1, log_files, ow);
            return
        end

        % ── Stage 1 ──────────────────────────────────────────────────────
        % Skip if the head preprocessing cache already exists.
        % Sentinel: cache/sub-NNN_<medium>_cropped_smoothed.mat
        % (written by preproc_head with preproc_affix = '').
        % When skipped, stages 2–6 are submitted without an afterok
        % dependency so they can start immediately.
        preproc_sentinel = fullfile(output_dir, 'cache', ...
            sprintf('sub-%03d_%s_cropped_smoothed.mat', subject_id, medium));
        p_stage1.hpc.timelimit = options.stage1_timelimit;
        p_stage1.hpc.job_name  = sprintf('PRESTUS-%su1-preproc_%s', run_prefix, subj);
        p_stage1               = apply_memorylimit(p_stage1, options.stage1_memorylimit);
        p_stage1               = strip_gpu_requirements(p_stage1, options.stage1_partition);
        if isfile(preproc_sentinel) && ~ow.rerun_intermediates
            fprintf('[Stage 1] Preprocessing cache exists — skipping submission.\n');
            job_id_stage1 = [];
        else
            fprintf('[Stage 1] Submitting preprocessing job...\n');
            job_id_stage1 = prestus_pipeline_start(p_stage1);
            fprintf('          → job %s  (mem %iG, time %s)\n', ...
                num2str(job_id_stage1), p_stage1.hpc.memorylimit, p_stage1.hpc.timelimit);
        end
        fprintf('\n');

        % ── Stages 2–6 ───────────────────────────────────────────────────
        % Each variant is skipped if its output CSV already exists.
        % If not skipped, it depends on stage 1 (when stage 1 was submitted).
        [job_ids_sim, sim_submitted] = submit_variant_sims(p_variants, options.affixes, ...
            job_id_stage1, run_prefix, '[Stage %d]', parameters, subj, options, ow);

        % ── Stage 7 ──────────────────────────────────────────────────────
        % Skip if the report already exists.
        % Depends on whichever simulation jobs were actually submitted.
        report_file = fullfile(output_dir, ...
            sprintf('sub-%03d_%s_desc-uncertainty_report%s.html', subject_id, medium, options.affixes.default));
        job_id_report = submit_report(p_report, report_file, job_ids_sim(sim_submitted), ...
            sprintf('[Stage %d]', n_var + 2), ...
            sprintf('PRESTUS-%su%d-report_%s', run_prefix, n_var + 2, subj), options, ow);

        % ── Sequential runs (HPC) ────────────────────────────────────────
        % Each sequential run submits its variants with afterok on the
        % prior run's report job, then submits its own report job.
        all_submitted_ids = [job_id_stage1, job_ids_sim(sim_submitted), job_id_report];
        if has_sequential
            [seq_submitted_ids, ~] = submit_sequential_hpc( ...
                parameters, options, job_id_report, subj, run_index_start);
            all_submitted_ids = [all_submitted_ids, seq_submitted_ids];
        end

        % ── Summary ──────────────────────────────────────────────────────
        fprintf('All jobs submitted. Monitor: squeue --user $USER -o "%%.18i %%.9P %%.50j %%.8u %%.2t %%.10M %%.6D %%R"\n');
        submitted_ids = all_submitted_ids(all_submitted_ids > 0);
        jobs.all = submitted_ids;
        if ~isempty(submitted_ids)
            fprintf('Submitted job IDs : %s\n', num2str(submitted_ids));
            fprintf('To cancel         : scancel %s\n', num2str(submitted_ids));
        end
        fprintf('HPC logs          : %s/%s/log_hpc/\n', parameters.path.sim, subj);

    % ---------------------------------------------------------------------
    otherwise
        error('uncertainty_pipeline:unknownPlatform', ...
            'Unknown platform ''%s''. Expected ''matlab'', ''slurm'', or ''qsub''.', platform);
end

end % uncertainty_pipeline

% =========================================================================
%% Sequential helpers
% =========================================================================

function seq_affixes = make_seq_affixes(seq_number, base_affixes)
% Derive variant affixes for sequential run <seq_number>: each variant's
% affix with the seq suffix appended.
    seq_suffix = sprintf('_seq%d', seq_number);
    seq_affixes = struct();
    for f = fieldnames(base_affixes)'
        seq_affixes.(f{1}) = [base_affixes.(f{1}) seq_suffix];
    end
end

function p = wire_thermal_handoff(p, nii_dir, prior_affix, subject_id, medium)
% Set adopted_heatmap / adopted_cem43 / adopted_cem43_iso from prior variant.
    p.io.adopted_heatmap    = fullfile(nii_dir, ...
        sprintf('sub-%03d_%s_T1w%s_heating_end.nii.gz',  subject_id, medium, prior_affix));
    p.io.adopted_cem43      = fullfile(nii_dir, ...
        sprintf('sub-%03d_%s_T1w%s_CEM43_end.nii.gz',    subject_id, medium, prior_affix));
    p.io.adopted_cem43_iso  = fullfile(nii_dir, ...
        sprintf('sub-%03d_%s_T1w%s_CEM43_iso_end.nii.gz',subject_id, medium, prior_affix));
end

function [sorted_fields, sorted_numbers] = sort_seq_configs(seq_configs)
% Return sequential config field names and their numeric indices in sorted order.
    fields   = fieldnames(seq_configs);
    numbers  = cellfun(@(x) sscanf(x, 'config_%d'), fields);
    [sorted_numbers, idx] = sort(numbers);
    sorted_fields = fields(idx);
end

function medium_config = get_variant_medium_config(variant_name, options)
% Return the YAML override path for a given variant name ([] for default).
    if strcmp(variant_name, 'default')
        medium_config = [];
    else
        medium_config = options.([variant_name '_config']);
    end
end

function label = variant_stage_label(variant)
    label = sprintf('%s simulation', variant.label);
    if ~variant.heating
        label = [label ' (acoustic only)'];
    end
end

function reject_legacy_options(options)
% The liberal / conservative variants of v0.6 were replaced by the Temp and
% Trans pairs; their options do not map onto either pair.
    legacy = {'liberal_config', 'conservative_config'};
    used = legacy(isfield(options, legacy));
    if isfield(options, 'affixes') && isstruct(options.affixes) && ...
            any(isfield(options.affixes, {'liberal', 'conservative'}))
        used{end+1} = 'affixes.liberal / affixes.conservative';
    end
    if ~isempty(used)
        error('uncertainty_pipeline:legacyVariant', ...
            ['options.%s belong to the removed liberal/conservative variants. Use ' ...
             'options.temp_high_config, .temp_low_config, .trans_high_config, ' ...
             '.trans_low_config and options.affixes.<variant> (see UNCERTAINTY_VARIANTS ' ...
             'and doc_uncertainty.md).'], strjoin(used, ', options.'));
    end
end

% -------------------------------------------------------------------------

function seq_run_params = run_sequential_matlab(base_parameters, options, run_index_start)
% Execute sequential uncertainty runs on the MATLAB platform.
% Returns a struct with one field per variant name, each a cell array of
% parameter structs (one per sequential run), for passing to
% generate_sequential_report.

    subject_id = base_parameters.subject_id;
    medium     = base_parameters.simulation.medium;
    nii_dir    = fullfile(get_output_dir(base_parameters), 'nii');
    run_ts     = string(datetime('now'), 'yyMMdd_HHmm');
    logs_dir   = fullfile(get_output_dir(base_parameters), 'log');

    [sorted_fields, sorted_numbers] = sort_seq_configs(options.sequential_configs);
    ow = resolve_chained_overwrite(base_parameters);

    variants       = uncertainty_variants();
    n_var          = numel(variants);
    seq_run_params = struct();
    for vi = 1:n_var
        seq_run_params.(variants(vi).name) = {};
    end
    prior_affixes  = options.affixes;

    for si = 1:numel(sorted_fields)
        run_idx    = run_index_start + si - 1;
        seq_p_base = options.sequential_configs.(sorted_fields{si});
        seq_num    = sorted_numbers(si);
        seq_aff    = make_seq_affixes(seq_num, options.affixes);

        p_variants = cell(1, n_var);

        for vi = 1:n_var
            vname       = variants(vi).name;
            cur_affix   = seq_aff.(vname);
            prior_affix = prior_affixes.(vname);
            med_cfg     = get_variant_medium_config(vname, options);

            p_v = make_sim_params(seq_p_base, cur_affix, med_cfg, variants(vi).heating, ow);
            if variants(vi).heating
                p_v = wire_thermal_handoff(p_v, nii_dir, prior_affix, subject_id, medium);
            end

            % Assign log file
            p_v.io.log_file = fullfile(logs_dir, ...
                sprintf('sub-%03d_%s%s_r%d_%s.txt', subject_id, medium, cur_affix, run_idx, run_ts));

            p_variants{vi} = p_v;
        end

        % Build report params for this sequential run
        seq_log_files.stage1 = '';
        for vi = 1:n_var
            seq_log_files.(variants(vi).name) = p_variants{vi}.io.log_file;
        end
        seq_log_files.report = fullfile(logs_dir, ...
            sprintf('sub-%03d_%s_desc-report_r%d_%s.txt', subject_id, medium, run_idx, run_ts));
        p_rep = make_report_params(seq_p_base, seq_aff, seq_log_files, ow, options);

        % Run
        for vi = 1:n_var
            run_stage(vi + 1, sprintf('%s (r%d)', variant_stage_label(variants(vi)), run_idx), ...
                @() prestus_pipeline(p_variants{vi}));
        end
        assert_variant_outputs_exist(seq_p_base, seq_aff);
        run_stage(n_var + 2, sprintf('Uncertainty report (r%d)', run_idx), ...
            @() prestus_pipeline(p_rep));

        for vi = 1:n_var
            seq_run_params.(variants(vi).name){end+1} = p_variants{vi};
        end

        prior_affixes = seq_aff;
    end
end

% -------------------------------------------------------------------------

function [all_ids, prior_job_id_report] = submit_sequential_hpc( ...
        base_parameters, options, initial_report_job_id, subj, run_index_start)
% Submit sequential uncertainty runs to the HPC scheduler.
% Each run's variants depend on the prior run's report job.
% Returns all submitted job IDs and the final report job ID.

    subject_id = base_parameters.subject_id;
    medium     = base_parameters.simulation.medium;
    nii_dir    = fullfile(get_output_dir(base_parameters), 'nii');
    run_ts     = string(datetime('now'), 'yyMMdd_HHmm');
    logs_dir   = fullfile(get_output_dir(base_parameters), 'log');
    output_dir = get_output_dir(base_parameters);

    [sorted_fields, sorted_numbers] = sort_seq_configs(options.sequential_configs);
    ow = resolve_chained_overwrite(base_parameters);

    all_ids            = [];
    prior_job_id_report = initial_report_job_id;
    prior_affixes      = options.affixes;

    variants = uncertainty_variants();
    n_var    = numel(variants);

    for si = 1:numel(sorted_fields)
        run_idx    = run_index_start + si - 1;
        run_prefix = sprintf('r%d-', run_idx);
        seq_p_base = options.sequential_configs.(sorted_fields{si});
        seq_num    = sorted_numbers(si);
        seq_aff    = make_seq_affixes(seq_num, options.affixes);

        seq_log_files = struct('stage1', '');
        p_seq = cell(1, n_var);
        for vi = 1:n_var
            vname       = variants(vi).name;
            cur_affix   = seq_aff.(vname);
            med_cfg     = get_variant_medium_config(vname, options);

            p_v = make_sim_params(seq_p_base, cur_affix, med_cfg, variants(vi).heating, ow);
            if variants(vi).heating
                p_v = wire_thermal_handoff(p_v, nii_dir, prior_affixes.(vname), subject_id, medium);
            end
            p_v.io.log_file = fullfile(logs_dir, ...
                sprintf('sub-%03d_%s%s_r%d_%s.txt', subject_id, medium, cur_affix, run_idx, run_ts));
            seq_log_files.(vname) = p_v.io.log_file;
            p_seq{vi} = p_v;
        end
        [job_ids_sim, sim_submitted] = submit_variant_sims(p_seq, seq_aff, prior_job_id_report, ...
            run_prefix, sprintf('[r%d-u%%d]', run_idx), base_parameters, subj, options, ow);

        % Report job for this sequential run
        seq_log_files.report = fullfile(logs_dir, ...
            sprintf('sub-%03d_%s_desc-report_r%d_%s.txt', subject_id, medium, run_idx, run_ts));

        p_rep = make_report_params(seq_p_base, seq_aff, seq_log_files, ow, options);
        report_file = fullfile(output_dir, ...
            sprintf('sub-%03d_%s_desc-uncertainty_report%s.html', subject_id, medium, seq_aff.default));

        prior_job_id_report = submit_report(p_rep, report_file, job_ids_sim(sim_submitted), ...
            sprintf('[r%d-u%d]', run_idx, n_var + 2), ...
            sprintf('PRESTUS-%su%d-report_%s', run_prefix, n_var + 2, subj), options, ow);

        all_ids = [all_ids, job_ids_sim(sim_submitted), prior_job_id_report];
        prior_affixes = seq_aff;
    end
end

% =========================================================================
%% HPC submission helpers
% =========================================================================

function [job_ids, submitted] = submit_variant_sims(p_variants, affixes, depend_id, ...
        run_prefix, label_fmt, base, subj, options, ow)
% Submit one simulation job per uncertainty variant. A variant whose output
% CSV exists is skipped unless the run re-runs outputs. label_fmt holds one
% %d for the stage number (e.g. '[Stage %d]'). Returns the job IDs (0 where
% skipped) and which variants were submitted.
    variants   = uncertainty_variants();
    n_var      = numel(variants);
    output_dir = get_output_dir(base);
    job_ids    = zeros(1, n_var);
    submitted  = false(1, n_var);
    for vi = 1:n_var
        name     = variants(vi).name;
        p        = p_variants{vi};
        csv_file = fullfile(output_dir, sprintf('sub-%03d_%s%s.csv', ...
            base.subject_id, base.simulation.medium, affixes.(name)));
        fprintf('%s Submitting %s simulation    ', sprintf(label_fmt, vi + 1), name);
        if isfile(csv_file) && ~ow.rerun_outputs
            fprintf('— output exists, skipping.\n');
            continue
        end
        p.hpc.timelimit = options.sim_timelimit;
        p.hpc.job_name  = sprintf('PRESTUS-%su%d-sim-%s_%s', run_prefix, vi + 1, name, subj);
        p               = apply_memorylimit(p, options.sim_memorylimit);
        if ~isempty(depend_id)
            p.hpc.depend_job_id = depend_id;
        end
        job_ids(vi)   = prestus_pipeline_start(p);
        submitted(vi) = true;
        fprintf('→ job %s  (mem %iG)\n', num2str(job_ids(vi)), p.hpc.memorylimit);
    end
    fprintf('\n');
end

function job_id = submit_report(p, report_file, depend_ids, label, job_name, options, ow)
% Submit the uncertainty report job after the given simulation jobs; []
% when the report exists and the run does not re-run outputs.
    fprintf('%s Submitting uncertainty report      ', label);
    job_id = [];
    if isfile(report_file) && ~ow.rerun_outputs
        fprintf('— report exists, skipping.\n\n');
        return
    end
    p.hpc.timelimit = options.report_timelimit;
    p.hpc.job_name  = job_name;
    p               = apply_memorylimit(p, options.report_memorylimit);
    p               = strip_gpu_requirements(p, options.report_partition);
    if ~isempty(depend_ids)
        p.hpc.depend_job_ids = depend_ids;
    end
    job_id = prestus_pipeline_start(p);
    fprintf('→ job %s  (mem %iG)\n\n', num2str(job_id), p.hpc.memorylimit);
end

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

function mode = nested_mode(p)
% Chained mode each uncertainty variant runs: 'async' | 'multi_isppa' | ''.
% Async absorbs a vector target as its own sweep, so it takes precedence.
    if is_async_mode(p)
        mode = 'async';
    elseif is_multi_isppa_mode(p)
        mode = 'multi_isppa';
    else
        mode = '';
    end
end

function p = clear_nested_modes(p)
% Stage 1 and the report run once for all variants (and all sweep points):
% clear sweep targets and async coupling so neither re-enters a chained
% pipeline. Neither stage uses the drive level (CACHE_SIGNATURE excludes it).
    if isfield(p, 'transducer')
        for ti = 1:numel(p.transducer)
            if isfield(p.transducer(ti), 'target_isppa_wcm2')
                p.transducer(ti).target_isppa_wcm2 = NaN;
            end
        end
    end
    if isfield(p, 'simulation') && isfield(p.simulation, 'transducer_coupling')
        p.simulation.transducer_coupling = 'coherent';
    end
end

function run = run_variant(p, nested, platform, inner_options)
% Run one variant's inner chained pipeline; returns its job IDs and the
% output affix of each sweep point (run.output_affixes).
    p.platform = platform;
    switch nested
        case 'multi_isppa'
            run = multi_isppa_pipeline(p, inner_options);
        case 'async'
            run = async_transducer_pipeline(p, inner_options);
    end
    if ~isfield(run, 'all')
        run.all = [run.preproc(:)', run.acoustic(:)', run.combine(:)', ...
                   run.thermal(:)', run.water(:)'];
    end
    run.all = run.all(run.all > 0);
end

function inner = inner_options(options, prefix, depend_ids)
% Pass the caller's options through to the inner pipeline (it reads the
% fields it knows, e.g. .thermal_timelimit), plus nesting knobs.
    inner = options;
    inner.job_name_prefix = prefix;
    inner.depend_job_ids  = depend_ids;
end

function reports = make_point_reports(parameters, options, runs, log_files, ow)
% One uncertainty report per sweep point j, comparing the variants'
% outputs at that point. Returns a struct array with .params, .suffix, .file.
    n_points = numel(runs{1}.output_affixes);
    if any(cellfun(@(r) numel(r.output_affixes), runs) ~= n_points)
        error('uncertainty_pipeline:pointMismatch', ...
            'The uncertainty variants produced different numbers of sweep points.');
    end
    names = {uncertainty_variants().name};
    reports = struct('params', {}, 'suffix', {}, 'file', {});
    for j = 1:n_points
        % Each point's affix is <variant affix><point suffix>; the suffix is shared.
        suffix = point_suffix(runs{1}.output_affixes{j}, options.affixes.(names{1}));
        affixes = struct();
        for v = 1:numel(names)
            if ~strcmp(point_suffix(runs{v}.output_affixes{j}, options.affixes.(names{v})), suffix)
                error('uncertainty_pipeline:pointMismatch', ...
                    'Sweep point %d has different output affixes across variants.', j);
            end
            affixes.(names{v}) = runs{v}.output_affixes{j};
        end
        lf = log_files;
        [d, n, e] = fileparts(log_files.report);
        lf.report = fullfile(d, [n suffix e]);
        p = make_report_params(parameters, affixes, lf, ow, options);
        reports(j).params = p;
        reports(j).suffix = suffix;
        reports(j).file   = fullfile(get_output_dir(parameters), sprintf( ...
            'sub-%03d_%s_desc-uncertainty_report%s.html', ...
            parameters.subject_id, parameters.simulation.medium, p.uncertainty.report_suffix));
    end
end

function suffix = point_suffix(full_affix, variant_affix)
    if ~startsWith(full_affix, variant_affix)
        error('uncertainty_pipeline:pointMismatch', ...
            'Output affix ''%s'' does not start with variant affix ''%s''.', full_affix, variant_affix);
    end
    suffix = full_affix(numel(variant_affix)+1:end);
end

function run_nested_matlab(parameters, options, nested, platform, p_variants, p_stage1, log_files, ow)
% MATLAB platform, nested mode: Stage 1 (multi-ISPPA only), each variant's
% inner pipeline, then one report per sweep point.
    if strcmp(nested, 'multi_isppa')
        run_stage(1, 'Preprocessing & source setup', @() prestus_pipeline(p_stage1));
    end
    variants = uncertainty_variants();
    runs = cell(1, numel(variants));
    for v = 1:numel(variants)
        fprintf('────────────────────────────────────────\n');
        fprintf('[Stage %d] %s (%s)\n', v + 1, variant_stage_label(variants(v)), strrep(nested, '_', '-'));
        fprintf('────────────────────────────────────────\n');
        runs{v} = run_variant(p_variants{v}, nested, platform, inner_options(options, '', []));
    end
    reports = make_point_reports(parameters, options, runs, log_files, ow);
    for j = 1:numel(reports)
        assert_variant_outputs_exist(parameters, reports(j).params.uncertainty.affixes);
        run_stage(numel(variants) + 2, sprintf('Uncertainty report%s', reports(j).suffix), ...
                  @() prestus_pipeline(reports(j).params));
    end
end

function submitted_ids = submit_nested_hpc(parameters, options, nested, platform, p_variants, p_stage1, log_files, ow)
% HPC, nested mode: submit Stage 1 (multi-ISPPA only) and each variant's inner
% job graph from here, so every report job can wait for the jobs that write its inputs.
    subj = sprintf('sub-%03d', parameters.subject_id);
    job_id_stage1 = [];
    if strcmp(nested, 'multi_isppa')
        preproc_sentinel = fullfile(get_output_dir(parameters), 'cache', ...
            sprintf('sub-%03d_%s_cropped_smoothed.mat', parameters.subject_id, parameters.simulation.medium));
        if isfile(preproc_sentinel) && ~ow.rerun_intermediates
            fprintf('[Stage 1] Preprocessing cache exists — skipping submission.\n');
        else
            p_stage1.hpc.timelimit = options.stage1_timelimit;
            p_stage1.hpc.job_name  = sprintf('PRESTUS-u1-preproc_%s', subj);
            p_stage1               = apply_memorylimit(p_stage1, options.stage1_memorylimit);
            p_stage1               = strip_gpu_requirements(p_stage1, options.stage1_partition);
            job_id_stage1          = prestus_pipeline_start(p_stage1);
            fprintf('[Stage 1] Preprocessing → job %s\n', num2str(job_id_stage1));
        end
    else
        fprintf('[Stage 1] Not needed: each async variant builds its own grid and medium.\n');
    end
    fprintf('\n');

    variants = {uncertainty_variants().name};
    report_stage = numel(variants) + 2;
    runs = cell(1, numel(variants));
    for v = 1:numel(variants)
        fprintf('[Stage %d] %s variant (%s)\n', v + 1, variants{v}, strrep(nested, '_', '-'));
        prefix  = sprintf('u%d-%s-', v + 1, variants{v});
        runs{v} = run_variant(p_variants{v}, nested, platform, ...
                              inner_options(options, prefix, job_id_stage1));
        fprintf('\n');
    end
    variant_ids = [job_id_stage1, cellfun(@(r) r.all, runs, 'UniformOutput', false)];
    variant_ids = [variant_ids{:}];

    reports = make_point_reports(parameters, options, runs, log_files, ow);
    report_ids = [];
    for j = 1:numel(reports)
        fprintf('[Stage %d] Uncertainty report%s  ', report_stage, reports(j).suffix);
        if isfile(reports(j).file) && ~ow.rerun_outputs
            fprintf('— report exists, skipping.\n');
            continue
        end
        p = reports(j).params;
        p.hpc.timelimit = options.report_timelimit;
        p.hpc.job_name  = sprintf('PRESTUS-u%d-report%s_%s', report_stage, reports(j).suffix, subj);
        p               = apply_memorylimit(p, options.report_memorylimit);
        p               = strip_gpu_requirements(p, options.report_partition);
        p               = hpc_add_dependency(p, variant_ids);
        report_ids(end+1) = prestus_pipeline_start(p); %#ok<AGROW>
        fprintf('→ job %s\n', num2str(report_ids(end)));
    end
    fprintf('\n');

    submitted_ids = [variant_ids, report_ids];
    submitted_ids = submitted_ids(submitted_ids > 0);
    fprintf('All jobs submitted. Monitor: squeue --user $USER -o "%%.18i %%.9P %%.50j %%.8u %%.2t %%.10M %%.6D %%R"\n');
    if ~isempty(submitted_ids)
        fprintf('Submitted job IDs : %s\n', num2str(submitted_ids));
        fprintf('To cancel         : scancel %s\n', num2str(submitted_ids));
    end
end

function p = make_stage1_params(base, ow)
% Stage 1: preprocessing and source setup only — no simulation.
%
% Produces the shared head-preprocessing cache used by all simulation
% variants (skull, tissue masks, etc.). Each variant computes its own kwave
% source matrix (via io.output_affix) because the time axis depends on the
% variant's medium sound speed.
%
% GPU is not required: no k-Wave simulation runs in this stage.
    p = clear_uncertainty_flag(base);
    p = clear_nested_modes(p);
    p.modules.run_source_setup         = 1;
    p.modules.run_acoustic_sims        = 0;
    p.modules.run_heating_sims         = 0;
    p.modules.run_thermal_analysis     = 0;
    p.modules.run_posthoc_water_sims   = 0;
    p.modules.generate_report          = 0;
    p.io.output_affix                  = '';
    p.io.overwrite_files               = ow.intermediates;
    p.io.overwrite_simnibs             = 0;
    p.io.save_source_matrices          = 0;   % source is recomputed per variant; no need to cache here
    p.io.save_acoustic_matrices        = 0;   % no acoustic sim in stage 1
    p.io.save_thermal_matrices         = 0;   % no thermal sim in stage 1
    p.simulation.interactive           = 0;
    p                                  = strip_gpu_requirements(p);
end

function p = make_sim_params(base, affix, medium_config, heating, ow)
% Stages 2–6: acoustic simulation, plus thermal when heating is true and the
% base config runs it (heating = false for the acoustic-only Trans variants).
%
% If medium_config is a valid YAML path, its properties are merged on top
% of the base parameters (only the overridden fields change; all transducer,
% grid, and path settings are inherited).
%
% run_source_setup = 1 is intentional: each variant checks for cached source
% files from stage 1 and skips recomputation unless overwrite_files = 'always'.
    p = clear_uncertainty_flag(base);

    if ~isempty(medium_config)
        if ~(ischar(medium_config) || isstring(medium_config)) || ~isfile(medium_config)
            error('uncertainty_pipeline:configNotFound', ...
                'Medium config file not found: %s', medium_config);
        end
        extra = yaml.loadFile(char(medium_config), 'ConvertToArray', true);
        p     = mergestruct(p, extra);
        p     = medium_absorption_layer(p, extra, char(medium_config));
    end

    p.modules.run_source_setup         = 1;
    p.modules.run_acoustic_sims        = 1;
    % Thermal steps: inherit from base config rather than forcing on.
    % run_heating_sims = 0 in the base will produce acoustic-only uncertainty
    % variants, which is valid when thermal estimation is not requested.
    run_thermal = heating && isfield(base.modules, 'run_heating_sims') && base.modules.run_heating_sims;
    % Nested multi-ISPPA / async pipelines run their own thermal jobs; this
    % flag keeps them acoustic only for a non-heating variant.
    p.uncertainty.acoustic_only        = ~heating;
    % With a pseudo-CT, only the default variant maps the skull from it; the
    % bound variants keep its head geometry with a homogeneous skull layer
    % (PCT_SKULL_MAPPED).
    p.uncertainty.homogeneous_skull    = ~isempty(medium_config);
    p.modules.run_heating_sims         = run_thermal;
    p.modules.run_thermal_analysis     = run_thermal && ...
                                         isfield(base.modules, 'run_thermal_analysis') && ...
                                         base.modules.run_thermal_analysis;
    p.modules.run_posthoc_water_sims   = 0;
    % Free water does not depend on the medium variant, but each variant runs
    % its own acoustic simulation. 'auto' runs the baseline only when a target
    % ISPPA needs one (no calibration, or target differs from it), so a target
    % is never reported unscaled; see NEEDS_WATER_BASELINE.
    p.modules.run_water_baseline       = 'auto';
    p.modules.generate_report          = 1;
    p.io.output_affix                  = affix;
    % Both the grid cache and kwave source were written by stage 1 with
    % affix ''. Point lookups at those files rather than recomputing.
    p.io.preproc_affix                 = '';
    % The medium is variant-specific: keep each variant's medium cache apart.
    p.io.medium_affix                  = affix;
    p.io.overwrite_files               = ow.outputs;
    % Medium is variant-specific, so only the stage-1 head/grid caches are reused.
    p.io.reuse_upstream_caches         = {'head', 'grid'};
    p.io.overwrite_simnibs             = 0;
    p.io.save_source_matrices          = 0;   % source is recomputed at the start of each variant; no benefit to caching
    p.io.save_acoustic_matrices        = 0;   % large; not needed after analysis
    p.io.save_thermal_matrices         = run_thermal;  % only needed when thermal ran
    p.simulation.interactive           = 0;
end

function p = make_report_params(base, affixes, log_files, ow, options)
% Stage 7: generate the uncertainty report only — no simulation.
%
% The affixes struct is passed through parameters.uncertainty.affixes so
% that generate_uncertainty_report can locate the variant output files.
% log_files is a struct with fields stage1, report and one per variant name
% — pre-assigned paths that allow the report to parse per-stage timing.
    p = clear_uncertainty_flag(base);
    if ~isempty(nested_mode(base))
        % The report reads the variants' outputs only; an async run has no
        % shared Stage-1 grid to load, so do not rebuild one.
        p.modules.run_grid_setup     = 0;
        p.modules.run_medium_setup   = 0;
    end
    p = clear_nested_modes(p);
    p.modules.run_source_setup       = 0;
    p.modules.run_acoustic_sims      = 0;
    p.modules.run_heating_sims       = 0;
    p.modules.run_thermal_analysis   = 0;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.run_water_baseline     = 0;
    p.modules.generate_report        = 0;
    p.modules.uncertainty_report     = 1;
    p.io.output_affix                = '';
    p.io.log_file                    = log_files.report;
    p.io.overwrite_files             = ow.outputs;
    p.io.reuse_upstream_caches       = {'head', 'grid', 'medium'};
    p.simulation.interactive         = 0;
    p.uncertainty.affixes            = affixes;
    p.uncertainty.log_files          = log_files;
    % Medium file of each bound variant, for the pseudo-CT comparison
    p.uncertainty.medium_files       = struct();
    for v = uncertainty_variants()
        if ~strcmp(v.name, 'default')
            p.uncertainty.medium_files.(v.name) = get_variant_medium_config(v.name, options);
        end
    end
    % The report is named after the default variant's affix (base affix and,
    % in nested mode, the sweep point), so runs with different affixes coexist.
    p.uncertainty.report_suffix      = affixes.default;
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
% Called before the report stage on the MATLAB platform so the user knows in
% advance that the report will be incomplete rather than discovering it after
% the fact inside generate_uncertainty_report.
    subject_id = parameters.subject_id;
    medium     = parameters.simulation.medium;
    output_dir = get_output_dir(parameters);

    variants = uncertainty_variants();
    for v = 1:numel(variants)
        label  = variants(v).label;
        affix  = affixes.(variants(v).name);
        csv    = fullfile(output_dir, ...
            sprintf('sub-%03d_%s%s.csv', subject_id, medium, affix));
        if ~isfile(csv)
            warn('uncertainty_pipeline:missingOutput', ...
                '%s variant output not found — report will be incomplete.\n  Expected: %s', ...
                label, csv);
        end
    end
end
