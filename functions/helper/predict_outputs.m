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
% function cannot disagree with what actually runs. Precedence is the same on
% every platform: uncertainty > async > multi-ISPPA. Uncertainty combined with
% async or a multi-ISPPA sweep nests that pipeline inside each variant (see
% UNCERTAINTY_PIPELINE); plan.nested names it. A vector target under async
% coupling is async's own sweep.
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
%     .nested        uncertainty mode only: '' | 'async' | 'multi_isppa', the
%                    pipeline each variant runs; '' in every other mode
%     .platform_used 'matlab' | 'hpc' | 'auto' (only changes plan.jobs wording)
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

% parameters.platform = 'auto' (config_default.yaml's own default) is resolved
% inside PRESTUS_PIPELINE_START; it only changes how plan.jobs describes the run.
if isfield(parameters, 'platform') && strcmp(parameters.platform, 'matlab')
    platform_used = 'matlab';
elseif isfield(parameters, 'platform') && ismember(parameters.platform, {'slurm', 'qsub'})
    platform_used = 'hpc';
else
    platform_used = 'auto';
end

% Same order as PRESTUS_PIPELINE_START on every platform.
if uncertainty_on
    mode = 'uncertainty';
elseif async_on
    mode = 'async';
elseif multi_isppa_on
    mode = 'multi_isppa';
else
    mode = 'single';
end
nested = '';
if uncertainty_on && async_on
    nested = 'async';
elseif uncertainty_on && multi_isppa_on
    nested = 'multi_isppa';
end

warnings_ = {};
if ~strcmp(mode, 'single')
    warnings_{end+1} = sprintf(['Assumes %s runs with its default options (see its own help ' ...
        'text); a non-default call to prestus_pipeline_start(parameters, options) can add to ' ...
        'or remove from this list.'], mode_function(mode));
end

plan = struct();
plan.mode          = mode;
plan.nested         = nested;
plan.platform_used  = platform_used;
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
        variants = uncertainty_variants();
        % Output affix suffix of each sweep point an inner pipeline writes ('' when not nested).
        switch nested
            case 'multi_isppa'
                targets = parameters.transducer(1).target_isppa_wcm2(:)';
                targets = targets(isfinite(targets));
                point_suffixes = arrayfun(@(t) sprintf('_desc-isppa%03d', round(t)), targets, ...
                    'UniformOutput', false);
            case 'async'
                n_points = async_sweep_points(parameters);
                point_suffixes = arrayfun(@(j) async_affix_thermal('', j, n_points), 1:n_points, ...
                    'UniformOutput', false);
            otherwise
                point_suffixes = {''};
        end
        for k = 1:numel(variants)
            variant_affix = [base_affix, variants(k).desc];
            vlabel = [variants(k).label, ' variant'];
            if strcmp(variants(k).pair, 'trans')
                vlabel = [vlabel, ' (acoustic only)'];
            end
            for j = 1:numel(point_suffixes)
                affix = [variant_affix, point_suffixes{j}];
                label = vlabel;
                if numel(point_suffixes) > 1
                    label = sprintf('%s, sweep point %d/%d', label, j, numel(point_suffixes));
                end
                artifacts(end+1) = art(sprintf('sub-%03d_%s%s.csv', subject_id, medium, affix), ...
                    [label, ' metrics']); %#ok<AGROW>
                artifacts(end+1) = art(sprintf('sub-%03d_%s_report%s.html', subject_id, medium, affix), ...
                    [label, ' report']); %#ok<AGROW>
            end
            if strcmp(nested, 'multi_isppa')
                artifacts(end+1) = art(sprintf('sub-%03d_%s_desc-multiIsppa_report%s.html', ...
                    subject_id, medium, variant_affix), [vlabel, ': summary comparing all targets']); %#ok<AGROW>
            end
        end
        for j = 1:numel(point_suffixes)
            label = 'combined uncertainty report (default, Temp-H/L, Trans-H/L)';
            if numel(point_suffixes) > 1
                label = sprintf('%s, sweep point %d/%d', label, j, numel(point_suffixes));
            end
            artifacts(end+1) = art(sprintf('sub-%03d_%s_desc-uncertainty_report%s%s.html', ...
                subject_id, medium, base_affix, point_suffixes{j}), label); %#ok<AGROW>
        end
        switch nested
            case 'multi_isppa'
                jobs{end+1} = '1 preprocessing job (stage 1)';
                jobs{end+1} = sprintf(['5 multi-ISPPA runs, one per variant (default, Temp-H, Temp-L, ' ...
                    'Trans-H, Trans-L): default, Temp-H and Temp-L each 1 acoustic job, %d thermal jobs ' ...
                    'and 1 summary job; Trans-H and Trans-L run acoustic only, without thermal jobs'], ...
                    numel(point_suffixes));
                jobs{end+1} = sprintf('%d uncertainty report job(s), one per target (stage 7)', ...
                    numel(point_suffixes));
            case 'async'
                jobs{end+1} = ['5 async multi-transducer runs, one per variant (default, Temp-H, ' ...
                    'Temp-L, Trans-H, Trans-L), each with its own preprocessing, acoustic, combine and ' ...
                    'joint jobs; Trans-H and Trans-L run acoustic only'];
                jobs{end+1} = sprintf('%d uncertainty report job(s), one per sweep point (stage 7)', ...
                    numel(point_suffixes));
            otherwise
                jobs{end+1} = '1 preprocessing job (stage 1)';
                jobs{end+1} = ['5 variant simulation jobs: default, Temp-H, Temp-L, Trans-H, Trans-L ' ...
                    '(stages 2-6); Trans-H and Trans-L run acoustic only'];
                jobs{end+1} = '1 uncertainty report job (stage 7)';
        end
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
