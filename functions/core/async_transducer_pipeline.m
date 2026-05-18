function async_transducer_pipeline(parameters, options)
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
%                  Affixes: _desc-t01, _desc-t02, …, _desc-tNN
%
%   Stage N+1…N+M  One COMBINE_ASYNC_INTENSITY call per sweep point j.
%                  Passes the j-th target for each transducer; writes a
%                  combined cache with affix _desc-asyncCombinedSjj (or
%                  _desc-asyncCombined when M=1).
%
%   Stage N+M+1…   One thermal job per sweep point, reading the matching
%   N+2M           combined cache.  Affixes: _desc-asyncThermalSjj (or
%                  _desc-asyncThermal when M=1).
%
% Use as:
%   async_transducer_pipeline(parameters)
%   async_transducer_pipeline(parameters, options)
%
% Input:
%   parameters - PRESTUS config with parameters.transducer as a struct array
%                (N ≥ 2) and simulation.transducer_coupling = 'async'.
%                Each transducer may carry a scalar or vector
%                target_isppa_wcm2.  Do NOT set io.output_affix.
%   options    - (optional) struct with HPC knobs:
%                .acoustic_timelimit   wall time per acoustic job ('06:00:00')
%                .combine_timelimit    wall time per combine job  ('00:30:00')
%                .thermal_timelimit    wall time per thermal job  ('04:00:00')
%                .acoustic_memorylimit / combine_memorylimit / thermal_memorylimit
%                .acoustic_partition  / combine_partition  / thermal_partition
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

% =========================================================================
%% Resolve options
% =========================================================================

if ~isfield(options, 'acoustic_timelimit');   options.acoustic_timelimit   = '06:00:00'; end
if ~isfield(options, 'combine_timelimit');    options.combine_timelimit    = '00:30:00'; end
if ~isfield(options, 'thermal_timelimit');    options.thermal_timelimit    = '04:00:00'; end
if ~isfield(options, 'acoustic_memorylimit'); options.acoustic_memorylimit = [];         end
if ~isfield(options, 'combine_memorylimit');  options.combine_memorylimit  = [];         end
if ~isfield(options, 'thermal_memorylimit');  options.thermal_memorylimit  = [];         end
if ~isfield(options, 'acoustic_partition');   options.acoustic_partition   = [];         end
if ~isfield(options, 'combine_partition');    options.combine_partition    = [];         end
if ~isfield(options, 'thermal_partition');    options.thermal_partition    = [];         end

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

% Build targets matrix [N × n_points]; scalars are broadcast.
targets_matrix = NaN(N, n_points);
for ti = 1:N
    if isscalar(target_vecs{ti})
        targets_matrix(ti, :) = target_vecs{ti};
    else
        targets_matrix(ti, :) = target_vecs{ti};
    end
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
    affix_per_t{ti} = sprintf('%s_desc-t%02d', base_affix, ti);
end

% Combined and thermal affixes — include sweep index when M > 1.
affix_combined = cell(1, n_points);
affix_thermal  = cell(1, n_points);
for j = 1:n_points
    if n_points == 1
        affix_combined{j} = [base_affix, '_desc-asyncCombined'];
        affix_thermal{j}  = [base_affix, '_desc-asyncThermal'];
    else
        affix_combined{j} = sprintf('%s_desc-asyncCombinedS%02d', base_affix, j);
        affix_thermal{j}  = sprintf('%s_desc-asyncThermalS%02d',  base_affix, j);
    end
end

% =========================================================================
%% Build parameter structs
% =========================================================================

p_acoustic = cell(1, N);
for ti = 1:N
    p_acoustic{ti} = make_acoustic_params(parameters, ti, affix_per_t{ti});
end

p_thermal = cell(1, n_points);
for j = 1:n_points
    p_thermal{j} = make_thermal_params(parameters, affix_combined{j}, affix_thermal{j});
end

% =========================================================================
%% Resolve file paths
% =========================================================================

subject_id = parameters.subject_id;
medium     = parameters.simulation.medium;
output_dir = get_output_dir(parameters);
cache_dir  = fullfile(output_dir, 'cache');

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
fprintf('╚══════════════════════════════════════════╝\n\n');

% =========================================================================
%% Execute
% =========================================================================

switch platform

    % ---------------------------------------------------------------------
    case 'matlab'
    % ---------------------------------------------------------------------

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

    % ---------------------------------------------------------------------
    case {'slurm', 'qsub'}
    % ---------------------------------------------------------------------

        subj = sprintf('sub-%03d', subject_id);

        % ── Stages 1…N: Acoustic (parallel) ──────────────────────────────
        job_ids_acoustic   = zeros(1, N);
        acoustic_submitted = false(1, N);

        for ti = 1:N
            fprintf('[Stage %d/%d] Acoustic — transducer %d\n', ti, N, ti);
            if isfile(files_acoustic{ti})
                fprintf('             Cache exists — skipping.\n');
            else
                pa = p_acoustic{ti};
                pa.hpc.timelimit = options.acoustic_timelimit;
                pa.hpc.job_name  = sprintf('PRESTUS-at%02d-ac_%s', ti, subj);
                pa = apply_memorylimit(pa, options.acoustic_memorylimit);
                pa = apply_partition(pa, options.acoustic_partition);
                job_ids_acoustic(ti)   = prestus_pipeline_start(pa);
                acoustic_submitted(ti) = true;
                fprintf('             → job %s\n', num2str(job_ids_acoustic(ti)));
            end
        end
        fprintf('\n');

        pending_acoustic = job_ids_acoustic(acoustic_submitted);

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
            if isfile(files_combined{j})
                fprintf('           Cache exists — skipping.\n');
            else
                p_comb = make_combine_params(parameters, files_acoustic, ...
                                             files_combined{j}, targets_matrix(:,j)');
                p_comb.hpc.timelimit = options.combine_timelimit;
                p_comb.hpc.job_name  = sprintf('PRESTUS-at-cb%02d_%s', j, subj);
                p_comb = apply_memorylimit(p_comb, options.combine_memorylimit);
                p_comb = apply_partition(p_comb, options.combine_partition);
                if ~isempty(pending_acoustic)
                    p_comb.hpc.depend_job_ids = pending_acoustic;
                end
                job_ids_combine(j)   = prestus_pipeline_start(p_comb);
                combine_submitted(j) = true;
                fprintf('           → job %s\n', num2str(job_ids_combine(j)));
            end
        end
        fprintf('\n');

        % ── Stages N+M+1…N+2M: Thermal per sweep point ───────────────────
        for j = 1:n_points
            stage_label = sprintf('Stage %d', N + n_points + j);
            if n_points > 1
                fprintf('[%s] Thermal — sweep %d/%d\n', stage_label, j, n_points);
            else
                fprintf('[%s] Thermal\n', stage_label);
            end
            if isfile(files_thermal{j})
                fprintf('           Output exists — skipping.\n');
            else
                pt = p_thermal{j};
                pt.hpc.timelimit = options.thermal_timelimit;
                pt.hpc.job_name  = sprintf('PRESTUS-at-th%02d_%s', j, subj);
                pt = apply_memorylimit(pt, options.thermal_memorylimit);
                pt = apply_partition(pt, options.thermal_partition);
                if combine_submitted(j)
                    pt.hpc.depend_job_id = job_ids_combine(j);
                end
                job_id_th = prestus_pipeline_start(pt);
                fprintf('           → job %s\n', num2str(job_id_th));
            end
        end
        fprintf('\n');

        fprintf('All jobs submitted. Monitor: squeue -u $USER\n');

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

function p = make_acoustic_params(base, transducer_idx, affix)
    p = base;
    p.transducer                     = base.transducer(transducer_idx);
    % Clear target vector so the acoustic stage doesn't trigger multi-ISPPA.
    p.transducer.target_isppa_wcm2   = NaN;
    p.modules.run_source_setup       = 1;
    p.modules.run_acoustic_sims      = 1;
    p.modules.run_water_baseline     = 1;  % needed for per-transducer baseline
    p.modules.run_heating_sims       = 0;
    p.modules.run_thermal_analysis   = 0;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 0;
    p.io.output_affix                = affix;
    p.io.overwrite_files             = 'never';
    p.io.overwrite_simnibs           = 0;
    p.io.save_acoustic_matrices      = 1;
    p.io.save_thermal_matrices       = 0;
    p.simulation.interactive         = 0;
    p.simulation.transducer_coupling = 'coherent';  % prevent re-entry
end

function p = make_combine_params(base, files_in, file_out, targets_wcm2)
    p = base;
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

function p = make_thermal_params(base, affix_combined, affix_thermal)
    p = base;
    % Clear all per-transducer targets — scaling was already done in combine.
    for ti = 1:numel(p.transducer)
        p.transducer(ti).target_isppa_wcm2 = NaN;
    end
    p.modules.run_source_setup       = 0;
    p.modules.run_acoustic_sims      = 0;
    p.modules.run_water_baseline     = 0;
    p.modules.run_heating_sims       = 1;
    p.modules.run_thermal_analysis   = 1;
    p.modules.run_posthoc_water_sims = 0;
    p.modules.generate_report        = 1;
    p.io.acoustic_cache_affix        = affix_combined;
    p.io.output_affix                = affix_thermal;
    p.io.preproc_affix               = '';
    p.io.overwrite_files             = 'never';
    p.io.overwrite_simnibs           = 0;
    p.io.save_acoustic_matrices      = 0;
    p.io.save_thermal_matrices       = 1;
    p.simulation.interactive         = 0;
    p.simulation.transducer_coupling = 'coherent';
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
