function plan = predict_outputs(parameters)
% PREDICT_OUTPUTS  Predict which files a PRESTUS run will write, before submission.
%
% Reads a resolved PRESTUS config (as returned by LOAD_PARAMETERS) and reports
% which of the four pipeline modes PRESTUS_PIPELINE_START dispatches to — the
% plain single-subject pipeline, ASYNC_TRANSDUCER_PIPELINE, MULTI_ISPPA_PIPELINE
% or UNCERTAINTY_PIPELINE — and the primary reportable artifacts (CSV/HTML) each
% one writes. Pure and read-only: no simulation runs, nothing is written to disk.
%
% Mode is detected with the SAME predicates PRESTUS_PIPELINE_START itself uses
% (IS_ASYNC_MODE, IS_MULTI_ISPPA_MODE, parameters.simulation.uncertainty), so this
% function cannot disagree with what actually runs. Dispatch PRECEDENCE differs
% by platform (see PRESTUS_PIPELINE_START): on 'matlab', async > multi-ISPPA >
% uncertainty; on 'slurm'/'qsub', uncertainty > async > multi-ISPPA. This only
% matters when uncertainty is competing with async or multi-ISPPA (an async
% sweep or a vector multi-ISPPA target on their own are unambiguous on every
% platform). parameters.platform = 'auto' (config_default.yaml's own default) is
% resolved to a concrete platform by hpc_detect_system() *inside*
% PRESTUS_PIPELINE_START, so a prediction made beforehand cannot know which way
% it will go; in that ambiguous+auto case, plan.mode is only a representative
% pick and plan.warnings names both possible outcomes explicitly.
%
% This only predicts the DEFAULT options of each sub-pipeline (e.g. the
% multi-ISPPA summary report, water validation): a non-default
% prestus_pipeline_start(parameters, options) call can change which of these
% files actually get written — see plan.warnings.
%
% cache/, debug/, log/, nii/ and img/ are written the same way in every mode
% (see RESOLVE_IO_DIRS) and are reported once in plan.common_dirs rather than
% repeated per stage/variant.
%
% Use as:
%   plan = predict_outputs(parameters)
%
% Input:
%   parameters - PRESTUS config as returned by load_parameters: subject_id,
%                path.sim, simulation.medium, transducer(:).target_isppa_wcm2
%                and (optionally) platform/modules.*/io.* should already be
%                resolved, as they are after load_parameters.
%
% Output:
%   plan - struct with fields:
%     .mode          'single' | 'async' | 'multi_isppa' | 'uncertainty'
%     .platform_used precedence this prediction assumed: 'matlab' | 'hpc' | 'auto'
%     .ambiguous     true only when uncertainty competes with async/multi-ISPPA
%                    (precedence between those two never actually varies by
%                    platform, so that combination alone is never ambiguous)
%     .output_dir    get_output_dir(parameters)
%     .artifacts     struct array (.name, .note) of primary output files,
%                    relative to output_dir
%     .common_dirs   cellstr of subfolders shared by every mode
%     .warnings      cellstr of caveats worth reading before submission
%     .jobs          cellstr describing the deployment structure: how many
%                    jobs/stages this run submits and what each one does
%
% See also: PRESTUS_PIPELINE_START, IS_ASYNC_MODE, IS_MULTI_ISPPA_MODE,
%           ASYNC_TRANSDUCER_PIPELINE, MULTI_ISPPA_PIPELINE, UNCERTAINTY_PIPELINE,
%           GET_OUTPUT_DIR, RESOLVE_IO_DIRS, PRINT_PREDICTED_OUTPUTS

arguments
    parameters (1,1) struct
end

subject_id = parameters.subject_id;
medium     = parameters.simulation.medium;
base_affix = '';
if isfield(parameters, 'io') && isfield(parameters.io, 'output_affix') && ...
        ~isempty(parameters.io.output_affix)
    base_affix = parameters.io.output_affix;
end

async_on       = is_async_mode(parameters);
multi_isppa_on = is_multi_isppa_mode(parameters);
uncertainty_on = isfield(parameters, 'simulation') && isfield(parameters.simulation, 'uncertainty') && ...
                 logical(parameters.simulation.uncertainty);

% platform='auto' (config_default.yaml's own default!) is resolved to a concrete
% platform by hpc_detect_system() *inside* PRESTUS_PIPELINE_START — by the time
% its dispatch switch runs, parameters.platform is never still 'auto'. Predicting
% from a config alone, before that call, we genuinely don't know which way 'auto'
% will go: keep it a distinct third state rather than guessing 'hpc'.
if isfield(parameters, 'platform') && strcmp(parameters.platform, 'matlab')
    platform_used = 'matlab';
elseif isfield(parameters, 'platform') && ismember(parameters.platform, {'slurm', 'qsub'})
    platform_used = 'hpc';
else
    platform_used = 'auto';
end

active = struct('async', async_on, 'multi_isppa', multi_isppa_on, 'uncertainty', uncertainty_on);
precedence_matlab = {'async', 'multi_isppa', 'uncertainty'};
precedence_hpc    = {'uncertainty', 'async', 'multi_isppa'};
mode_matlab = first_active(active, precedence_matlab);
mode_hpc    = first_active(active, precedence_hpc);

% A vector target_isppa_wcm2 under async coupling is NOT ambiguous: async mode
% absorbs it as its own internal sweep (see ASYNC_TRANSDUCER_PIPELINE's header),
% identically on every platform. The one platform-DEPENDENT case is uncertainty
% mode competing with async/multi-ISPPA — HPC checks uncertainty first,
% MATLAB checks it last (see PRESTUS_PIPELINE_START) — so only that
% combination is flagged, and only then does platform_used=='auto' matter.
ambiguous = uncertainty_on && (async_on || multi_isppa_on);
switch platform_used
    case 'matlab'
        mode = mode_matlab;
    case 'hpc'
        mode = mode_hpc;
    otherwise    % 'auto': representative pick; see the 'auto' warning below when ambiguous
        mode = mode_matlab;
end

warnings_ = {};
if ambiguous
    active_names = {};
    if async_on;       active_names{end+1} = 'async';       end
    if multi_isppa_on; active_names{end+1} = 'multi_isppa'; end
    if uncertainty_on; active_names{end+1} = 'uncertainty'; end
    base_msg = sprintf(['More than one pipeline trigger is active at once (%s). ' ...
        'PRESTUS_PIPELINE_START resolves this differently by platform: async > multi-ISPPA > ' ...
        'uncertainty on ''matlab'', but uncertainty > async > multi-ISPPA on ''slurm''/''qsub''.'], ...
        strjoin(active_names, ' + '));
    if strcmp(platform_used, 'auto')
        warnings_{end+1} = sprintf(['%s platform is ''auto'': PRESTUS decides at submission time ' ...
            '(hpc_detect_system), so this could run as ''%s'' mode (no scheduler found) or as ' ...
            '''%s'' mode (SLURM/qsub found). This preview shows ''%s'' mode; set platform ' ...
            'explicitly to pin down which one.'], base_msg, mode_matlab, mode_hpc, mode);
    else
        warnings_{end+1} = sprintf(['%s For platform_used=''%s'' this prediction assumes ''%s'' ' ...
            'mode runs; the others are ignored here, not combined.'], base_msg, platform_used, mode);
    end
end
if ~strcmp(mode, 'single')
    warnings_{end+1} = sprintf(['Assumes %s runs with its default options (see its own help ' ...
        'text); a non-default call to prestus_pipeline_start(parameters, options) can add to ' ...
        'or remove from this list.'], mode_function(mode));
end

plan = struct();
plan.mode          = mode;
plan.platform_used  = platform_used;
plan.ambiguous      = ambiguous;
plan.output_dir     = get_output_dir(parameters);
plan.common_dirs    = {'nii/', 'img/', 'log/', 'cache/'};
if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'debug') && ...
        parameters.simulation.debug
    plan.common_dirs{end+1} = 'debug/';
end
plan.warnings = warnings_;

artifacts = struct('name', {}, 'note', {});
jobs = {};
switch mode
    case 'single'
        artifacts(end+1) = art(sprintf('sub-%03d_%s%s.csv', subject_id, medium, base_affix), ...
            'metrics');
        artifacts(end+1) = art(sprintf('sub-%03d_%s_report%s.html', subject_id, medium, base_affix), ...
            'HTML report');
        if strcmp(medium, 'layered') && get_flag(parameters, 'modules', 'run_targeting_analysis', true)
            artifacts(end+1) = art(sprintf('sub-%03d_%s%s_targeting.csv', subject_id, medium, base_affix), ...
                'targeting accuracy');
        end
        jobs{end+1} = '1 job — the full pipeline, start to finish (prestus_pipeline.m)';

    case 'async'
        N = numel(parameters.transducer);
        n_points = async_sweep_points(parameters);
        run_water = get_flag(parameters, 'modules', 'run_posthoc_water_sims', true);
        if run_water
            for ti = 1:N
                artifacts(end+1) = art(sprintf('sub-%03d_water%s.csv', subject_id, ...
                    [base_affix, trans_suffix(ti, N)]), sprintf('free-water validation — transducer %d', ti));
            end
        end
        for j = 1:n_points
            affix = async_affix_thermal(base_affix, j, n_points);
            label = 'joint thermal + report';
            if n_points > 1
                label = sprintf('%s (sweep point %d/%d)', label, j, n_points);
            end
            artifacts(end+1) = art(sprintf('sub-%03d_%s%s.csv', subject_id, medium, affix), label);
            artifacts(end+1) = art(sprintf('sub-%03d_%s_report%s.html', subject_id, medium, affix), label);
        end
        async_grid_mode = 'shared';
        if isfield(parameters, 'simulation') && isfield(parameters.simulation, 'async_grid_mode') && ...
                ~isempty(parameters.simulation.async_grid_mode)
            async_grid_mode = parameters.simulation.async_grid_mode;
        end
        if strcmp(async_grid_mode, 'per_transducer')
            jobs{end+1} = sprintf(['%d preprocessing jobs — async_grid_mode=''per_transducer'', ' ...
                'one grid each'], N);
        else
            jobs{end+1} = ['1 shared preprocessing stage — async_grid_mode=''shared'', one grid ' ...
                'for all transducers'];
        end
        jobs{end+1} = sprintf('%d acoustic jobs — one per transducer', N);
        if run_water
            jobs{end+1} = sprintf('%d free-water validation jobs — one per transducer', N);
        end
        if n_points > 1
            jobs{end+1} = sprintf('%d combine job(s) — one per sweep point', n_points);
            jobs{end+1} = sprintf('%d joint thermal + report job(s) — one per sweep point', n_points);
        else
            jobs{end+1} = '1 combine job';
            jobs{end+1} = '1 joint thermal + report job';
        end

    case 'multi_isppa'
        targets = parameters.transducer(1).target_isppa_wcm2(:)';
        targets = targets(isfinite(targets));
        for ti = 1:numel(targets)
            affix = sprintf('%s_desc-isppa%03d', base_affix, round(targets(ti)));
            label = sprintf('%.0f W/cm^2 — thermal + report', targets(ti));
            artifacts(end+1) = art(sprintf('sub-%03d_%s%s.csv', subject_id, medium, affix), label);
            artifacts(end+1) = art(sprintf('sub-%03d_%s_report%s.html', subject_id, medium, affix), label);
        end
        artifacts(end+1) = art(sprintf('sub-%03d_%s_desc-multiIsppa_report%s.html', subject_id, medium, base_affix), ...
            'summary comparing all targets');
        jobs{end+1} = '1 acoustic + free-water baseline job (stage 1)';
        jobs{end+1} = sprintf('%d thermal jobs (stages 2–%d, one per target)', numel(targets), numel(targets) + 1);
        jobs{end+1} = '1 summary report job (final stage)';

    case 'uncertainty'
        variants = {'default', 'liberal', 'conservative'};
        for k = 1:numel(variants)
            v = variants{k};
            if strcmp(v, 'default')
                affix = base_affix;
            else
                affix = [base_affix, '_desc-', v];
            end
            artifacts(end+1) = art(sprintf('sub-%03d_%s%s.csv', subject_id, medium, affix), ...
                [v, ' variant metrics']);
            artifacts(end+1) = art(sprintf('sub-%03d_%s_report%s.html', subject_id, medium, affix), ...
                [v, ' variant report']);
        end
        artifacts(end+1) = art(sprintf('sub-%03d_%s_uncertainty_report.html', subject_id, medium), ...
            'combined uncertainty report (default/liberal/conservative)');
        jobs{end+1} = '1 preprocessing job (stage 1)';
        jobs{end+1} = '3 variant simulation jobs: default, liberal, conservative (stages 2-4)';
        jobs{end+1} = '1 uncertainty report job (stage 5)';
end
if ~strcmp(mode, 'single')
    switch platform_used
        case 'hpc'
            jobs{end+1} = 'On SLURM/qsub these submit as separate jobs with scheduler (afterok) dependencies.';
        case 'matlab'
            jobs{end+1} = 'In this MATLAB session these run sequentially, one after another.';
        otherwise
            jobs{end+1} = ['platform=''auto'': separate SLURM/qsub jobs with dependencies if a ' ...
                'scheduler is found at submission time, otherwise sequential stages in one MATLAB session.'];
    end
end
plan.artifacts = artifacts;
plan.jobs = jobs;
end

% =========================================================================
%% Helpers
% =========================================================================

function a = art(name, note)
a = struct('name', name, 'note', note);
end

function mode = first_active(active, precedence)
mode = 'single';
for k = 1:numel(precedence)
    if active.(precedence{k})
        mode = precedence{k};
        return;
    end
end
end

function tf = get_flag(parameters, section, field, default)
tf = default;
if isfield(parameters, section) && isfield(parameters.(section), field)
    tf = logical(parameters.(section).(field));
end
end

function n_points = async_sweep_points(parameters)
% Mirrors async_transducer_pipeline.m's own n_points derivation: the longest
% target_isppa_wcm2 vector across transducers (1 when every transducer carries
% a scalar or unset target).
n_points = 1;
for ti = 1:numel(parameters.transducer)
    t = parameters.transducer(ti);
    if isfield(t, 'target_isppa_wcm2') && any(isfinite(t.target_isppa_wcm2))
        n_points = max(n_points, numel(t.target_isppa_wcm2));
    end
end
end

function affix = async_affix_thermal(base_affix, j, n_points)
if n_points == 1
    affix = [base_affix, '_desc-asyncJoint'];
else
    affix = sprintf('%s_desc-asyncJointS%02d', base_affix, j);
end
end

function name = mode_function(mode)
switch mode
    case 'async';       name = 'async_transducer_pipeline';
    case 'multi_isppa';  name = 'multi_isppa_pipeline';
    case 'uncertainty';  name = 'uncertainty_pipeline';
    otherwise;           name = 'prestus_pipeline';
end
end
