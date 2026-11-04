function jobs = async_transducer_pipeline(parameters, options)
% ASYNC_TRANSDUCER_PIPELINE  Simulate N asynchronously firing transducers
%                            with independent acoustic fields and a joint
%                            thermal model.
%
% Asynchronously firing transducers (e.g., alternating pulses, non-overlapping
% duty cycles) do not produce coherent acoustic interference. Each transducer
% is first scaled to its own target free-water ISPPA, then heat deposition
% is the incoherent sum of the scaled intensity fields:
%
%   p_i_scaled = p_i · sqrt(target_i / baseline_i)
%   Q_total   ∝ Σ  p_i_scaled² / (2·ρ·c)
%
% Multi-ISPPA sweep: when any transducer carries a vector target
% (transducer(i).target_isppa_wcm2 with more than one value), one
% combine+thermal pair is run per sweep point.  All transducers must carry
% either a scalar target (held fixed across the sweep) or a vector of the
% same length as the sweep.
%
% Activated automatically by PRESTUS_PIPELINE_START when
% parameters.simulation.transducer_coupling = 'async' and
% numel(parameters.transducer) > 1.
%
% Pipeline stages:
%
%   Stage 1…N      One acoustic simulation per transducer (parallel on HPC).
%                  Affixes: _trans-01, _trans-02, …, _trans-NN
%
%   Stage N+1…N+M  One COMBINE_ASYNC_INTENSITY call per sweep point j.
%                  Passes the j-th target for each transducer; writes a
%                  combined cache with affix _desc-asyncCombinedSjj (or
%                  _desc-asyncCombined when M=1).
%
%   Stage N+M+1…   One joint (thermal + intensity) job per sweep point,
%   N+2M           reading the combined cache.  Affixes: _desc-asyncJointSjj
%                  (or _desc-asyncJoint when M=1).
%
%   Stage N+2M+1…  (optional) One free-water acoustic simulation per
%   N+2M+N         transducer, when modules.run_posthoc_water_sims == 1.
%                  Uses the same driving parameters as Stage 1…N but with
%                  simulation.medium = 'water'.  Affixes match Stage 1…N.
%
% Overwrite semantics (see RESOLVE_CHAINED_OVERWRITE):
%   io.overwrite_files = 'never'   intermediate stages (preproc, acoustic,
%                                  combine) resume from existing caches;
%                                  output stages (joint, water validation)
%                                  are skipped when their CSV exists.
%   io.overwrite_files = 'always'  intermediates resume; output stages always
%                                  re-run and overwrite CSV/NIfTI/report.
%   io.overwrite_files = 'force'   every stage re-runs.
%   Acoustic and joint sub-jobs always load the shared head/grid/medium
%   caches built by Stage 0 (io.reuse_upstream_caches), so parallel jobs
%   never rebuild them.
%
% Use as:
%   async_transducer_pipeline(parameters)
%   async_transducer_pipeline(parameters, options)
%   jobs = async_transducer_pipeline(...)
%
% Input:
%   parameters - PRESTUS config with parameters.transducer as a struct array
%                (N ≥ 2) and simulation.transducer_coupling = 'async'.
%                Each transducer may carry a scalar or vector
%                target_isppa_wcm2.  io.output_affix, if set, prefixes every
%                stage's affix (e.g. _desc-run1_trans-01, _desc-run1_desc-asyncJoint).
%                simulation.async_grid_mode (optional): 'shared' (default) —
%                one grid sized to fit every transducer and oriented by
%                grid.mode (default 'focal_axis': mean beam axis); or
%                'per_transducer' — each transducer gets its own optimally
%                aligned, minimally cropped grid, resampled onto
%                transducer(1)'s grid at combine time. See doc_multitransducer.md.
%   options    - (optional) struct with HPC knobs:
%                .acoustic_timelimit   wall time per acoustic job ('06:00:00')
%                .combine_timelimit    wall time per combine job  ('00:30:00')
%                .thermal_timelimit    wall time per thermal job  ('04:00:00')
%                .acoustic_memorylimit / combine_memorylimit / thermal_memorylimit
%                .acoustic_partition  / combine_partition  / thermal_partition
%                .depend_job_ids       scheduler job IDs every submitted job also
%                                      waits for (default []); set by a caller
%                                      that nests this run, e.g. UNCERTAINTY_PIPELINE
%                .job_name_prefix      inserted after 'PRESTUS-' in job names ('')
%
% Output:
%   jobs       - scheduler job IDs submitted per stage (fields preproc,
%                acoustic, combine, thermal, water; [] for stages that were
%                skipped because their output exists, and all [] when run in
%                MATLAB). Lets a caller make later jobs depend on this run.
%                .output_affixes (1×M cell) is the io.output_affix of each sweep
%                point's joint (thermal) outputs.
%
% See also: PRESTUS_PIPELINE_START, COMBINE_ASYNC_INTENSITY,
%           IS_ASYNC_MODE, MULTI_ISPPA_PIPELINE

arguments
    parameters (1,1) struct
    options    (1,1) struct = struct()
end

% =========================================================================
%% Validate
% =========================================================================

N = numel(parameters.transducer);
if N < 2
    error('async_transducer_pipeline:singleTransducer', ...
        'async mode requires at least 2 transducers; only %d found.', N);
end
% BabelBrain backends support one transducer per run (errors here, before any job)
if (isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'babelbrain')) || ...
        (isfield(parameters, 'thermal_backend') && strcmpi(parameters.thermal_backend, 'babelbrain'))
    babelbrain_check_geometry(parameters);
end

% =========================================================================
%% Resolve options
% =========================================================================

if ~isfield(options, 'preproc_timelimit');    options.preproc_timelimit    = '01:00:00'; end
if ~isfield(options, 'acoustic_timelimit');   options.acoustic_timelimit   = '06:00:00'; end
if ~isfield(options, 'combine_timelimit');    options.combine_timelimit    = '00:30:00'; end
if ~isfield(options, 'thermal_timelimit');    options.thermal_timelimit    = '04:00:00'; end
if ~isfield(options, 'water_timelimit');      options.water_timelimit      = '06:00:00'; end
if ~isfield(options, 'preproc_memorylimit');  options.preproc_memorylimit  = [];         end
if ~isfield(options, 'acoustic_memorylimit'); options.acoustic_memorylimit = [];         end
if ~isfield(options, 'combine_memorylimit');  options.combine_memorylimit  = [];         end
if ~isfield(options, 'thermal_memorylimit');  options.thermal_memorylimit  = [];         end
if ~isfield(options, 'water_memorylimit');    options.water_memorylimit    = [];         end
if ~isfield(options, 'preproc_partition');    options.preproc_partition    = [];         end
if ~isfield(options, 'acoustic_partition');   options.acoustic_partition   = [];         end
if ~isfield(options, 'combine_partition');    options.combine_partition    = [];         end
if ~isfield(options, 'thermal_partition');    options.thermal_partition    = [];         end
if ~isfield(options, 'water_partition');      options.water_partition      = [];         end
if ~isfield(options, 'depend_job_ids');       options.depend_job_ids       = [];         end
if ~isfield(options, 'job_name_prefix');      options.job_name_prefix      = '';         end

% =========================================================================
%% Resolve multi-ISPPA sweep points
% =========================================================================

% Collect target vectors per transducer; scalars broadcast across the sweep.
target_vecs = cell(1, N);
n_points    = 1;
for ti = 1:N
    if isfield(parameters.transducer(ti), 'target_isppa_wcm2') && ...
            ~isempty(parameters.transducer(ti).target_isppa_wcm2) && ...
            any(isfinite(parameters.transducer(ti).target_isppa_wcm2))
        target_vecs{ti} = parameters.transducer(ti).target_isppa_wcm2(:)';
        n_points = max(n_points, numel(target_vecs{ti}));
    else
        target_vecs{ti} = NaN;  % no scaling for this transducer
    end
end

% Validate lengths: each must be scalar (1) or equal to n_points.
for ti = 1:N
    if numel(target_vecs{ti}) > 1 && numel(target_vecs{ti}) ~= n_points
        error('async_transducer_pipeline:sweepLengthMismatch', ...
            ['transducer(%d).target_isppa_wcm2 has %d values but the sweep has ' ...
             '%d points. All transducers must carry a scalar target or a vector ' ...
             'of the same length.'], ti, numel(target_vecs{ti}), n_points);
    end
end

% Build targets matrix [N × n_points]; scalars broadcast across sweep points.
targets_matrix = NaN(N, n_points);
for ti = 1:N
    targets_matrix(ti, :) = target_vecs{ti};
end

% =========================================================================
%% Resolve platform
% =========================================================================

platform = parameters.platform;
if strcmp(platform, 'auto')
    platform = hpc_detect_system();
end

% =========================================================================
%% Build affixes
% =========================================================================

base_affix = '';
if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix') && ...
        ~isempty(parameters.io.output_affix)
    base_affix = parameters.io.output_affix;
end

affix_per_t = cell(1, N);
for ti = 1:N
    affix_per_t{ti} = [base_affix, trans_suffix(ti, N)];
end

% Grid geometry mode: 'shared' (default) builds one grid that fits every
% transducer (oriented by grid.mode, union bounding box); 'per_transducer' gives
% each transducer its own optimally-aligned, minimally-cropped grid (matching
% solo/sequential mode) and resamples onto transducer(1)'s grid at combine
% time. See doc_multitransducer.md.
async_grid_mode = 'shared';
if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'async_grid_mode') && ...
        ~isempty(parameters.simulation.async_grid_mode)
    async_grid_mode = parameters.simulation.async_grid_mode;
end
if ~ismember(async_grid_mode, {'shared', 'per_transducer'})
    error('async_transducer_pipeline:unknownGridMode', ...
        'simulation.async_grid_mode must be ''shared'' or ''per_transducer'', got ''%s''.', ...
        async_grid_mode);
end
per_transducer_grid = strcmp(async_grid_mode, 'per_transducer');

% Preproc affix(es): shared mode builds one grid sized to fit every
% transducer; per_transducer mode builds N independent grids, each aligned
% and cropped to just its own transducer.
if per_transducer_grid
    affix_preproc = cell(1, N);
    for ti = 1:N
        affix_preproc{ti} = [base_affix, sprintf('_desc-asyncPreproc-t%02d', ti)];
    end
else
    affix_preproc = [base_affix, '_desc-asyncPreproc'];
end

% Combined and thermal affixes — include sweep index when M > 1.
affix_combined = cell(1, n_points);
affix_thermal  = cell(1, n_points);
for j = 1:n_points
    if n_points == 1
        affix_combined{j} = [base_affix, '_desc-asyncCombined'];
        affix_thermal{j}  = [base_affix, '_desc-asyncJoint'];
    else
        affix_combined{j} = sprintf('%s_desc-asyncCombinedS%02d', base_affix, j);
        affix_thermal{j}  = sprintf('%s_desc-asyncJointS%02d',    base_affix, j);
    end
end

% =========================================================================
%% Build parameter structs
% =========================================================================

ow = resolve_chained_overwrite(parameters);

p_acoustic = cell(1, N);
for ti = 1:N
    p_acoustic{ti} = make_acoustic_params(parameters, ti, affix_per_t{ti}, ...
        preproc_affix_for(affix_preproc, ti), ow);
end

% Thermal always keys off transducer(1)'s grid — matches the reference frame
% combine_async_intensity uses for both 'shared' and 'per_transducer' modes.
p_thermal = cell(1, n_points);
for j = 1:n_points
    p_thermal{j} = make_thermal_params(parameters, affix_combined{j}, affix_thermal{j}, ...
        preproc_affix_for(affix_preproc, 1), ow);
end

% =========================================================================
%% Resolve file paths
% =========================================================================

subject_id = parameters.subject_id;
medium     = parameters.simulation.medium;
output_dir = get_output_dir(parameters);
cache_dir  = fullfile(output_dir, 'cache');

if per_transducer_grid
    files_preproc = cell(1, N);
    for ti = 1:N
        files_preproc{ti} = fullfile(cache_dir, ...
            sprintf('sub-%03d_%s_grid_cache%s.mat', subject_id, medium, affix_preproc{ti}));
    end
else
    files_preproc = {fullfile(cache_dir, ...
        sprintf('sub-%03d_%s_grid_cache%s.mat', subject_id, medium, affix_preproc))};
end

files_acoustic = cell(1, N);
for ti = 1:N
    files_acoustic{ti} = fullfile(cache_dir, ...
        sprintf('sub-%03d_%s_results%s.mat', subject_id, medium, affix_per_t{ti}));
end
files_combined = cell(1, n_points);
files_thermal  = cell(1, n_points);
for j = 1:n_points
    files_combined{j} = fullfile(cache_dir, ...
        sprintf('sub-%03d_%s_results%s.mat', subject_id, medium, affix_combined{j}));
    files_thermal{j}  = fullfile(output_dir, ...
        sprintf('sub-%03d_%s%s.csv', subject_id, medium, affix_thermal{j}));
end

run_water = isfield(parameters, 'modules') && ...
            isfield(parameters.modules, 'run_posthoc_water_sims') && ...
            parameters.modules.run_posthoc_water_sims == 1;

files_water = cell(1, N);
for ti = 1:N
    files_water{ti} = fullfile(output_dir, ...
        sprintf('sub-%03d_water%s.csv', subject_id, affix_per_t{ti}));
end

% =========================================================================
%% Print header
% =========================================================================

fprintf('\n');
fprintf('╔══════════════════════════════════════════╗\n');
fprintf('║   PRESTUS Async Transducer Pipeline      ║\n');
fprintf('╠══════════════════════════════════════════╣\n');
fprintf('║  Subject      : sub-%03d                ║\n', subject_id);
fprintf('║  Medium       : %-26s ║\n', medium);
fprintf('║  Platform     : %-26s ║\n', upper(platform));
fprintf('║  Transducers  : %-26d ║\n', N);
fprintf('║  Sweep points : %-26d ║\n', n_points);
for ti = 1:N
    targets_str = strjoin(arrayfun(@(x) sprintf('%.1f', x), targets_matrix(ti,:), ...
        'UniformOutput', false), ', ');
    fprintf('║  t%02d target   : %-26s ║\n', ti, [targets_str, ' W/cm²']);
end
fprintf('║  Water valid. : %-26s ║\n', mat2str(run_water));
fprintf('╚══════════════════════════════════════════╝\n\n');

% =========================================================================
%% Execute
% =========================================================================

jobs = struct('preproc', [], 'acoustic', [], 'combine', [], 'thermal', [], 'water', [], ...
              'output_affixes', {affix_thermal});

switch platform

    % ---------------------------------------------------------------------
    case 'matlab'
    % ---------------------------------------------------------------------

        if per_transducer_grid
            for ti = 1:N
                run_stage(ti, sprintf('Preproc — transducer %d/%d (own grid)', ti, N), ...
                          @() prestus_pipeline(make_preproc_params(parameters, affix_preproc{ti}, ti, ow)));
            end
        else
            run_stage(0, 'Preproc — shared grid for all transducers', ...
                      @() prestus_pipeline(make_preproc_params(parameters, affix_preproc, [], ow)));
        end

        for ti = 1:N
            run_stage(ti, sprintf('Acoustic — transducer %d/%d', ti, N), ...
                      @() prestus_pipeline(p_acoustic{ti}));
        end

        for j = 1:n_points
            label = sprintf('Combine (%d transducers)', N);
            if n_points > 1; label = [label, sprintf(' — sweep %d/%d', j, n_points)]; end %#ok<AGROW>
            run_stage(N + j, label, ...
                @() combine_async_intensity(files_acoustic, files_combined{j}, ...
                                            targets_matrix(:, j)'));
        end

        for j = 1:n_points
            label = 'Thermal on combined field';
            if n_points > 1; label = [label, sprintf(' — sweep %d/%d', j, n_points)]; end %#ok<AGROW>
            run_stage(N + n_points + j, label, ...
                      @() prestus_pipeline(p_thermal{j}));
        end

        if run_water
            for ti = 1:N
                run_stage(N + 2*n_points + ti, ...
                    sprintf('Water validation — transducer %d/%d', ti, N), ...
                    @() prestus_pipeline(make_water_params(parameters, ti, affix_per_t{ti}, ow)));
            end
        end

    % ---------------------------------------------------------------------
    case {'slurm', 'qsub'}
    % ---------------------------------------------------------------------

        subj = sprintf('sub-%03d', subject_id);

        % ── Stage 0: Preproc grid(s) ──────────────────────────────────────
        job_ids_preproc   = zeros(1, N);
        preproc_submitted = false(1, N);

        if per_transducer_grid
            for ti = 1:N
                fprintf('[Stage 0.%d] Preproc — transducer %d/%d (own grid)\n', ti, ti, N);
                if isfile(files_preproc{ti}) && ~ow.rerun_intermediates
                    fprintf('             Cache exists — skipping.\n');
                else
                    p_prep = make_preproc_params(parameters, affix_preproc{ti}, ti, ow);
                    p_prep.hpc.timelimit = options.preproc_timelimit;
                    p_prep.hpc.job_name  = sprintf('PRESTUS-%sat%02d-pp_%s', options.job_name_prefix, ti, subj);
                    p_prep = apply_memorylimit(p_prep, options.preproc_memorylimit);
                    p_prep = apply_partition(p_prep, options.preproc_partition);
                    p_prep = hpc_add_dependency(p_prep, options.depend_job_ids);
                    job_ids_preproc(ti)   = prestus_pipeline_start(p_prep);
                    preproc_submitted(ti) = true;
                    fprintf('             → job %s\n', num2str(job_ids_preproc(ti)));
                end
            end
        else
            fprintf('[Stage 0] Preproc — shared grid\n');
            if isfile(files_preproc{1}) && ~ow.rerun_intermediates
                fprintf('           Cache exists — skipping.\n');
            else
                p_prep = make_preproc_params(parameters, affix_preproc, [], ow);
                p_prep.hpc.timelimit = options.preproc_timelimit;
                p_prep.hpc.job_name  = sprintf('PRESTUS-%sat-pp_%s', options.job_name_prefix, subj);
                p_prep = apply_memorylimit(p_prep, options.preproc_memorylimit);
                p_prep = apply_partition(p_prep, options.preproc_partition);
                p_prep = hpc_add_dependency(p_prep, options.depend_job_ids);
                job_id_shared        = prestus_pipeline_start(p_prep);
                job_ids_preproc(:)   = job_id_shared;
                preproc_submitted(:) = true;
                fprintf('           → job %s\n', num2str(job_id_shared));
            end
        end
        fprintf('\n');

        % ── Stages 1…N: Acoustic (parallel) ──────────────────────────────
        job_ids_acoustic   = zeros(1, N);
        acoustic_submitted = false(1, N);

        for ti = 1:N
            fprintf('[Stage %d/%d] Acoustic — transducer %d\n', ti, N, ti);
            if isfile(files_acoustic{ti}) && ~ow.rerun_intermediates
                fprintf('             Cache exists — skipping.\n');
            else
                pa = p_acoustic{ti};
                pa.hpc.timelimit = options.acoustic_timelimit;
                pa.hpc.job_name  = sprintf('PRESTUS-%sat%02d-ac_%s', options.job_name_prefix, ti, subj);
                pa = apply_memorylimit(pa, options.acoustic_memorylimit);
                pa = apply_partition(pa, options.acoustic_partition);
                if preproc_submitted(ti)
                    pa = hpc_add_dependency(pa, job_ids_preproc(ti));
                end
                pa = hpc_add_dependency(pa, options.depend_job_ids);
                job_ids_acoustic(ti)   = prestus_pipeline_start(pa);
                acoustic_submitted(ti) = true;
                fprintf('             → job %s\n', num2str(job_ids_acoustic(ti)));
            end
        end
        fprintf('\n');

        pending_acoustic = job_ids_acoustic(acoustic_submitted);
        jobs.preproc  = unique(job_ids_preproc(preproc_submitted));
        jobs.acoustic = pending_acoustic;

        % ── Stages N+1…N+M: Combine per sweep point ──────────────────────
        job_ids_combine   = zeros(1, n_points);
        combine_submitted = false(1, n_points);

        for j = 1:n_points
            stage_label = sprintf('Stage %d', N + j);
            if n_points > 1
                fprintf('[%s] Combine — sweep %d/%d\n', stage_label, j, n_points);
            else
                fprintf('[%s] Combine\n', stage_label);
            end
            if isfile(files_combined{j}) && ~ow.rerun_intermediates
                fprintf('           Cache exists — skipping.\n');
            else
                p_comb = make_combine_params(parameters, files_acoustic, ...
                                             files_combined{j}, targets_matrix(:,j)');
                p_comb.hpc.timelimit = options.combine_timelimit;
                p_comb.hpc.job_name  = sprintf('PRESTUS-%sat-cb%02d_%s', options.job_name_prefix, j, subj);
                p_comb = apply_memorylimit(p_comb, options.combine_memorylimit);
                p_comb = apply_partition(p_comb, options.combine_partition);
                if ~isempty(pending_acoustic)
                    p_comb = hpc_add_dependency(p_comb, pending_acoustic);
                end
                p_comb = hpc_add_dependency(p_comb, options.depend_job_ids);
                job_ids_combine(j)   = prestus_pipeline_start(p_comb);
                combine_submitted(j) = true;
                fprintf('           → job %s\n', num2str(job_ids_combine(j)));
            end
        end
        fprintf('\n');
        jobs.combine = job_ids_combine(combine_submitted);

        % ── Stages N+M+1…N+2M: Thermal per sweep point ───────────────────
        for j = 1:n_points
            stage_label = sprintf('Stage %d', N + n_points + j);
            if n_points > 1
                fprintf('[%s] Thermal — sweep %d/%d\n', stage_label, j, n_points);
            else
                fprintf('[%s] Thermal\n', stage_label);
            end
            if isfile(files_thermal{j}) && ~ow.rerun_outputs
                fprintf('           Output exists — skipping.\n');
            else
                pt = p_thermal{j};
                pt.hpc.timelimit = options.thermal_timelimit;
                pt.hpc.job_name  = sprintf('PRESTUS-%sat-th%02d_%s', options.job_name_prefix, j, subj);
                pt = apply_memorylimit(pt, options.thermal_memorylimit);
                pt = apply_partition(pt, options.thermal_partition);
                if combine_submitted(j)
                    pt = hpc_add_dependency(pt, job_ids_combine(j));
                end
                pt = hpc_add_dependency(pt, options.depend_job_ids);
                job_id_th = prestus_pipeline_start(pt);
                jobs.thermal(end+1) = job_id_th;
                fprintf('           → job %s\n', num2str(job_id_th));
            end
        end
        fprintf('\n');

        % ── Stages N+2M+1…N+2M+N: Water validation (optional) ───────────
        if run_water
            % Water jobs depend on the acoustic jobs (need source/preproc cache).
            pending_for_water = job_ids_acoustic(acoustic_submitted);

            for ti = 1:N
                stage_label = sprintf('Stage %d', N + 2*n_points + ti);
                fprintf('[%s] Water validation — transducer %d/%d\n', stage_label, ti, N);
                if isfile(files_water{ti}) && ~ow.rerun_outputs
                    fprintf('           Output exists — skipping.\n');
                else
                    pw = make_water_params(parameters, ti, affix_per_t{ti}, ow);
                    pw.hpc.timelimit = options.water_timelimit;
                    pw.hpc.job_name  = sprintf('PRESTUS-%sat%02d-wv_%s', options.job_name_prefix, ti, subj);
                    pw = apply_memorylimit(pw, options.water_memorylimit);
                    pw = apply_partition(pw, options.water_partition);
                    if ~isempty(pending_for_water)
                        pw = hpc_add_dependency(pw, pending_for_water);
                    end
                    pw = hpc_add_dependency(pw, options.depend_job_ids);
                    job_id_wv = prestus_pipeline_start(pw);
                    jobs.water(end+1) = job_id_wv;
                    fprintf('           → job %s\n', num2str(job_id_wv));
                end
            end
            fprintf('\n');
        end

        fprintf('All jobs submitted. Monitor: squeue --user $USER -o "%%.18i %%.9P %%.50j %%.8u %%.2t %%.10M %%.6D %%R"\n');

    % ---------------------------------------------------------------------
    otherwise
        error('async_transducer_pipeline:unknownPlatform', ...
            'Unknown platform ''%s''. Expected ''matlab'', ''slurm'', or ''qsub''.', platform);
end

end % async_transducer_pipeline

% =========================================================================
%% Stage runner
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

function p = make_preproc_params(base, affix_preproc, transducer_idx, ow)
    p = base;
    if nargin >= 3 && ~isempty(transducer_idx)
        % Per-transducer grid: align/crop to this transducer alone, reusing
        % the same single-transducer preproc_head path as solo/sequential mode.
        p.transducer = base.transducer(transducer_idx);
    end
    % Otherwise (shared mode) keep ALL transducers so grid_tissue_setup sizes the grid to contain all of them.
    % Preproc ignores drive level; clear the targets so a sweep vector does not
    % route this job into multi_isppa_pipeline (cache_signature excludes them).
    for ti = 1:numel(p.transducer)
        p.transducer(ti).target_isppa_wcm2 = NaN;
    end
    p.modules.run_source_setup       = 0;
    p.modules.run_acoustic_sims      = 0;
    p.modules.run_water_baseline     = 0;
    p.modules.run_heating_sims       = 0;
    p.modules.run_thermal_analysis   = 0;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 0;
    p.io.output_affix                = affix_preproc;
    p.io.preproc_affix               = affix_preproc;
    p.io.reuse_upstream_caches       = {};  % this stage builds the caches
    p                                = own_preproc_caches(p);
    p.io.overwrite_files             = ow.intermediates;
    p.io.overwrite_simnibs           = 0;
    p.io.save_grid_cache             = 1;  % override save_matrices: 0 from config
    p.io.save_medium_cache           = 1;  % cache medium too for acoustic sub-runs
    p.simulation.interactive         = 0;
    p.simulation.transducer_coupling = 'coherent';  % prevent re-entry
end

function p = make_acoustic_params(base, transducer_idx, affix, affix_preproc, ow)
    p = base;
    p.transducer                     = base.transducer(transducer_idx);
    % Clear target so the acoustic stage does not trigger multi-ISPPA mode.
    % Scaling is done per-sweep-point inside combine_async_intensity.  As a
    % consequence, ISPPA values from per-transducer acoustic runs are at the
    % simulated drive level — compare against combine_async_intensity log output
    % (combination_isppa_per_transducer_wcm2) for target-level per-transducer ISPPAs.
    p.transducer.target_isppa_wcm2   = NaN;
    p.modules.run_source_setup       = 1;
    p.modules.run_acoustic_sims      = 1;
    p.modules.run_water_baseline     = 1;  % needed for per-transducer baseline
    p.modules.run_heating_sims       = 0;
    p.modules.run_thermal_analysis   = 0;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 0;
    p.io.output_affix                = affix;
    p.io.preproc_affix               = affix_preproc;  % load shared grid sized for all transducers
    p                                = own_preproc_caches(p);
    % Tell load_grid_cache_file which transducer this sub-job is, so that on a
    % cache hit it copies this transducer's grid-space trans_pos/focus_pos
    % (set by grid_transducer_location during Stage 0) instead of leaving the
    % caller's input-space (e.g. T1-voxel) positions in place.
    p.io.async_transducer_idx        = transducer_idx;
    p.io.overwrite_files             = ow.intermediates;
    p.io.reuse_upstream_caches       = {'head', 'grid', 'medium'};  % built by Stage 0
    p.io.overwrite_simnibs           = 0;
    p.io.save_acoustic_matrices      = 1;
    p.io.save_kwave_medium           = 1;  % combine_async_intensity needs per-voxel rho/c
    p.io.save_thermal_matrices       = 0;
    p.simulation.interactive         = 0;
    p.simulation.transducer_coupling = 'coherent';  % prevent re-entry
end

function p = make_combine_params(base, files_in, file_out, targets_wcm2)
    p = base;
    % Targets travel in async_combine.targets_wcm2; clear the sweep vectors so
    % this job does not route into multi_isppa_pipeline.
    for ti = 1:numel(p.transducer)
        p.transducer(ti).target_isppa_wcm2 = NaN;
    end
    p.modules.run_source_setup       = 0;
    p.modules.run_acoustic_sims      = 0;
    p.modules.run_water_baseline     = 0;
    p.modules.run_heating_sims       = 0;
    p.modules.run_thermal_analysis   = 0;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 0;
    p.modules.combine_async          = 1;
    p.async_combine.files_in         = files_in;
    p.async_combine.file_out         = file_out;
    p.async_combine.targets_wcm2     = targets_wcm2;
    p.io.output_affix                = '';
    p.io.overwrite_files             = 'never';
    p.simulation.interactive         = 0;
    p.simulation.transducer_coupling = 'coherent';
end

function p = make_thermal_params(base, affix_combined, affix_thermal, affix_preproc, ow)
    p = base;
    % Clear all per-transducer targets — scaling was already done in combine.
    for ti = 1:numel(p.transducer)
        p.transducer(ti).target_isppa_wcm2 = NaN;
    end
    p.modules.run_source_setup       = 0;
    p.modules.run_acoustic_sims      = 0;
    p.modules.run_water_baseline     = 0;
    % An acoustic-only uncertainty variant (Trans-H/L) reports the combined
    % intensity without heating.
    heating = ~(isfield(base, 'uncertainty') && isfield(base.uncertainty, 'acoustic_only') && ...
                base.uncertainty.acoustic_only);
    p.modules.run_heating_sims       = heating;
    p.modules.run_thermal_analysis   = heating;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 1;
    p.io.acoustic_cache_affix        = affix_combined;
    p.io.output_affix                = affix_thermal;
    % The combined field is already on the full grid (combine_async_intensity
    % back-projects partial-FOV fields), so this job must not crop again.
    p.grid.acoustic_fov_diameter_mm  = [];
    p.grid.acoustic_fov_length_mm    = [];
    p.io.preproc_affix               = affix_preproc;  % reuse trans-01 preproc cache
    p                                = own_preproc_caches(p);
    p.io.overwrite_files             = ow.outputs;
    p.io.reuse_upstream_caches       = {'head', 'grid', 'medium'};  % the combined field was computed on these
    p.io.overwrite_simnibs           = 0;
    p.io.save_acoustic_matrices      = 0;
    p.io.save_thermal_matrices       = heating;
    p.simulation.interactive         = 0;
    p.simulation.transducer_coupling = 'coherent';
end

function p = make_water_params(base, transducer_idx, affix, ow)
    p = base;
    p.transducer                      = base.transducer(transducer_idx);
    p.transducer.target_isppa_wcm2    = NaN;   % no scaling; raw drive level
    p.simulation.medium               = 'water';
    p.simulation.transducer_coupling  = 'coherent';  % prevent re-entry
    p.modules.run_source_setup        = 1;
    p.modules.run_acoustic_sims       = 1;
    p.modules.run_water_baseline      = 0;
    p.modules.run_heating_sims        = 0;
    p.modules.run_thermal_analysis    = 0;
    p.modules.run_posthoc_water_sims  = 0;   % prevent recursive call
    p.modules.generate_report         = 1;
    p.io.output_affix                 = affix;
    p.io.overwrite_files              = ow.outputs;
    p.io.overwrite_simnibs            = 0;
    p.io.save_acoustic_matrices       = 0;
    p.io.save_kwave_medium            = 0;
    p.io.save_thermal_matrices        = 0;
    p.simulation.interactive          = 0;
end

function p = own_preproc_caches(p)
% The medium cache follows io.preproc_affix here: Stage 0 builds it per grid.
% Drop an io.medium_affix inherited from an enclosing pipeline (an uncertainty
% variant), whose medium Stage 0 already maps under the variant's preproc affix.
    if isfield(p.io, 'medium_affix')
        p.io = rmfield(p.io, 'medium_affix');
    end
end

function affix = preproc_affix_for(affix_preproc, ti)
% Resolve the preproc affix for transducer ti regardless of grid mode:
% a cell array (per_transducer mode) indexes by ti; a plain string (shared
% mode) is the same for every transducer.
    if iscell(affix_preproc)
        affix = affix_preproc{ti};
    else
        affix = affix_preproc;
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
    elseif isfield(p, 'hpc') && isfield(p.hpc, 'partition')
        p.hpc = rmfield(p.hpc, 'partition');
    end
end
