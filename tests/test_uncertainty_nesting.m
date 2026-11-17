classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_uncertainty_nesting < matlab.unittest.TestCase
% TEST_UNCERTAINTY_NESTING  Job graph of uncertainty mode nested with a
% multi-ISPPA sweep or an async multi-transducer run (no simulation).
%
% prestus_pipeline_start and prestus_pipeline are replaced by recorders
% (fixtures written to a temp folder at the front of the path), so the
% chained pipelines build their full job graph without simulating or
% submitting anything. Each recorded job is what the real launcher would
% receive; the checks are on that graph:
%   - no recorded job re-enters a chained mode (e.g. async
%     preprocessing/combine jobs carrying sweep vectors),
%   - affixes, cache affixes and dependencies of every stage,
%   - one uncertainty report per sweep point, after every job it reads.
%
%   Run with:   results = runtests('tests/test_uncertainty_nesting.m');
%
% See also: UNCERTAINTY_PIPELINE, MULTI_ISPPA_PIPELINE, ASYNC_TRANSDUCER_PIPELINE

    properties
        TmpDir
        StubDir
        StartStubDir
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir  = tempname;
            tc.StubDir      = fullfile(tc.TmpDir, 'stubs');
            tc.StartStubDir = fullfile(tc.TmpDir, 'start_stub');
            mkdir(tc.StubDir);
            mkdir(tc.StartStubDir);
            write_recorders(tc.StubDir, tc.StartStubDir);
            addpath(tc.StubDir);
            addpath(tc.StartStubDir);
            setappdata(0, 'prestus_test_jobs', {});
            tc.addTeardown(@() rmdir(tc.TmpDir, 's'));
            tc.addTeardown(@() rmappdata_safe());
            tc.addTeardown(@() rmpath_safe(tc.StubDir));
            tc.addTeardown(@() rmpath_safe(tc.StartStubDir));
        end
    end

    methods (Test)

        function sweepHpcGraph(tc)
            p = base_parameters(tc.TmpDir, 'slurm');
            p.transducer(1).target_isppa_wcm2 = [10 30];
            out = uncertainty_pipeline(p);
            jobs = recorded();
            base = p.io.output_affix;
            [name, variants, heating] = variant_list(base);

            % 1 stage-1 + 5 x (acoustic + 2 per-point + summary) + 2 reports
            tc.verifyNumElements(jobs, 23);
            tc.verifyEqual(sort(out.all), sort(cellfun(@(j) j.id, jobs)), ...
                'every submitted job is returned to the caller');
            verify_no_reentry(tc, jobs);
            verify_unique_names(tc, jobs);

            s1 = find_jobs(jobs, 'u1-preproc');
            tc.verifyNumElements(s1, 1);
            tc.verifyTrue(all(isnan(s1{1}.targets{1})), 'stage 1 must not carry the sweep');

            all_variant_ids = s1{1}.id;
            for v = 1:numel(name)
                ac = find_jobs(jobs, sprintf('u%d-%s-mi1-acoustic', v + 1, name{v}));
                tc.verifyNumElements(ac, 1);
                tc.verifyEqual(ac{1}.depend, s1{1}.id, 'acoustic waits for stage 1');
                tc.verifyEqual(ac{1}.output_affix, variants{v});
                tc.verifyEqual(ac{1}.preproc_affix, '', 'acoustic reuses the stage-1 head/grid');
                tc.verifyEqual(ac{1}.medium_affix, variants{v}, 'medium cache is per variant');

                th = find_jobs(jobs, sprintf('u%d-%s-mi', v + 1, name{v}), 'thermal');
                tc.verifyNumElements(th, 2);
                for k = 1:2
                    tc.verifyEqual(th{k}.depend, ac{1}.id);
                    tc.verifyEqual(th{k}.preproc_affix, '');
                    tc.verifyEqual(th{k}.medium_affix, variants{v});
                    tc.verifyTrue(isscalar(th{k}.targets{1}));
                    tc.verifyEqual(logical(th{k}.modules.run_heating_sims), heating(v), ...
                        sprintf('%s: heating only for the default and Temp variants', name{v}));
                end
                tc.verifyEqual(sort(cellfun(@(j) j.output_affix, th, 'UniformOutput', false)), ...
                    sort({[variants{v} '_desc-isppa010'], [variants{v} '_desc-isppa030']}));
                sm = find_jobs(jobs, sprintf('u%d-%s-mi-summary', v + 1, name{v}));
                tc.verifyNumElements(sm, 1);
                all_variant_ids = [all_variant_ids, ac{1}.id, cellfun(@(j) j.id, th), sm{1}.id]; %#ok<AGROW>
            end

            reports = find_jobs(jobs, 'u7-report');
            tc.verifyNumElements(reports, 2);
            suffixes = sort(cellfun(@(j) j.report_suffix, reports, 'UniformOutput', false));
            tc.verifyEqual(suffixes, {[base '_desc-isppa010'], [base '_desc-isppa030']}, ...
                'report named after the default variant''s affix at that point');
            for k = 1:2
                r = reports{k};
                point = r.report_suffix(numel(base)+1:end);
                tc.verifyEqual(sort(r.depend), sort(all_variant_ids), ...
                    'each report waits for every variant job');
                tc.verifyEqual(r.report_affixes.default, r.report_suffix);
                for v = 2:numel(name)
                    tc.verifyEqual(r.report_affixes.(name{v}), [variants{v} point]);
                end
                tc.verifyEqual(r.modules.run_grid_setup, 0);
                tc.verifyTrue(all(isnan(r.targets{1})));
            end
        end

        function asyncHpcGraph(tc)
            p = base_parameters(tc.TmpDir, 'slurm');
            p.transducer(2) = p.transducer(1);
            p.simulation.transducer_coupling = 'async';
            p.transducer(1).target_isppa_wcm2 = [10 20];
            p.transducer(2).target_isppa_wcm2 = [10 20];
            uncertainty_pipeline(p);
            jobs = recorded();
            base = p.io.output_affix;
            [name, variants, heating] = variant_list(base);

            % No stage 1; 5 x (preproc + 2 acoustic + 2 combine + 2 joint) + 2 reports
            tc.verifyNumElements(jobs, 37);
            verify_no_reentry(tc, jobs);
            verify_unique_names(tc, jobs);
            tc.verifyEmpty(find_jobs(jobs, 'u1-preproc'));

            variant_ids = [];
            for v = 1:numel(name)
                vj = find_jobs(jobs, sprintf('u%d-%s-at', v + 1, name{v}));
                tc.verifyNumElements(vj, 7);
                pp = find_jobs(vj, '-pp_');
                tc.verifyNumElements(pp, 1);
                tc.verifyEqual(pp{1}.preproc_affix, [variants{v} '_desc-asyncPreproc'], ...
                    'async stage 0 builds its own caches');
                tc.verifyEmpty(pp{1}.reuse, 'stage 0 must not reuse upstream caches');
                % Combine jobs read the acoustic caches only; the others map or load a medium.
                for k = find(~contains(cellfun(@(j) j.name, vj, 'UniformOutput', false), '-cb'))
                    tc.verifyEqual(vj{k}.medium_affix, '<none>', ...
                        'async sub-jobs key the medium cache on their preproc affix');
                end
                th = find_jobs(vj, '-th');
                tc.verifyEqual(sort(cellfun(@(j) j.output_affix, th, 'UniformOutput', false)), ...
                    sort({[variants{v} '_desc-asyncJointS01'], [variants{v} '_desc-asyncJointS02']}));
                tc.verifyEqual(cellfun(@(j) logical(j.modules.run_heating_sims), th), ...
                    repmat(heating(v), 1, 2), sprintf('%s: joint jobs heat only for default and Temp', name{v}));
                variant_ids = [variant_ids, cellfun(@(j) j.id, vj)]; %#ok<AGROW>
            end

            reports = find_jobs(jobs, 'u7-report');
            tc.verifyNumElements(reports, 2);
            for k = 1:2
                tc.verifyEqual(sort(reports{k}.depend), sort(variant_ids));
                tc.verifyEqual(reports{k}.coupling, 'coherent');
            end
        end

        function sweepMatlabOrder(tc)
            % MATLAB platform: stages run in order, inner pipelines in place.
            p = base_parameters(tc.TmpDir, 'matlab');
            p.transducer(1).target_isppa_wcm2 = [10 30];
            uncertainty_pipeline(p);
            jobs = recorded();
            tc.verifyNumElements(jobs, 23);
            verify_no_reentry(tc, jobs);
            kinds = cellfun(@stage_kind, jobs, 'UniformOutput', false);
            tc.verifyEqual(kinds, [{'stage1'}, ...
                repmat({'acoustic', 'thermal', 'thermal', 'summary'}, 1, 3), ...
                repmat({'acoustic', 'point', 'point', 'summary'}, 1, 2), ...
                {'report', 'report'}]);
        end

        function matlabDispatchIsUncertaintyFirst(tc)
            % The real launcher on 'matlab' must route uncertainty + sweep to
            % uncertainty_pipeline (uncertainty is the outer mode).
            % telemetry_setup prompts once per user in an interactive session.
            tc.assumeTrue(batchStartupOptionUsed || ...
                isfile(fullfile(prefdir_prestus(), 'telemetry.json')), ...
                'telemetry consent not recorded: the launcher would prompt');
            rmpath(tc.StartStubDir);
            p = base_parameters(tc.TmpDir, 'matlab');
            p.transducer(1).target_isppa_wcm2 = [10 30];
            prestus_pipeline_start(p);
            jobs = recorded();
            tc.verifyNumElements(jobs, 23);
            tc.verifyEqual(stage_kind(jobs{1}), 'stage1');
            tc.verifyEqual(stage_kind(jobs{end}), 'report');
        end

        function sequentialNestedErrors(tc)
            p = base_parameters(tc.TmpDir, 'slurm');
            p.transducer(1).target_isppa_wcm2 = [10 30];
            opts.sequential_configs = struct('config_2', struct());
            tc.verifyError(@() uncertainty_pipeline(p, opts), ...
                'uncertainty_pipeline:nestedSequential');
            tc.verifyEmpty(recorded(), 'nothing may be submitted before the error');
        end

        function sequentialRunsChainOnReports(tc)
            % A follow-up run submits the five variants after run 1's
            % report, then its own report after them; only the heating
            % variants adopt run 1's end temperature.
            p = base_parameters(tc.TmpDir, 'slurm');
            seq = p;
            seq.simulation.uncertainty = false;
            uncertainty_pipeline(p, struct('sequential_configs', struct('config_2', seq)));
            jobs = recorded();
            [name, affixes, heating] = variant_list(p.io.output_affix);
            % stage 1 + 5 sims + report, then 5 sims + report
            tc.verifyNumElements(jobs, 13);
            verify_unique_names(tc, jobs);
            s1 = find_jobs(jobs, 'r1-u1-preproc');
            r1 = find_jobs(jobs, 'r1-u7-report');
            r2 = find_jobs(jobs, 'r2-u7-report');
            tc.assertNumElements(s1, 1); tc.assertNumElements(r1, 1); tc.assertNumElements(r2, 1);
            run1 = []; run2 = [];
            for v = 1:numel(name)
                a = find_jobs(jobs, sprintf('r1-u%d-sim-%s_', v + 1, name{v}));
                b = find_jobs(jobs, sprintf('r2-u%d-sim-%s_', v + 1, name{v}));
                tc.assertNumElements(a, 1, name{v}); tc.assertNumElements(b, 1, name{v});
                tc.verifyEqual(a{1}.depend, s1{1}.id, [name{v} ': run 1 waits for stage 1']);
                tc.verifyEqual(b{1}.depend, r1{1}.id, [name{v} ': run 2 waits for run 1''s report']);
                tc.verifyEqual(b{1}.output_affix, [affixes{v} '_seq2']);
                handoff = get_field(b{1}.p, {'io', 'adopted_heatmap'}, '');
                tc.verifyEqual(~isempty(handoff), heating(v), [name{v} ': heat handoff']);
                if heating(v)
                    tc.verifyTrue(endsWith(handoff, ['_T1w' affixes{v} '_heating_end.nii.gz']), name{v});
                end
                run1 = [run1, a{1}.id]; run2 = [run2, b{1}.id]; %#ok<AGROW>
            end
            tc.verifyEqual(r1{1}.depend, sort(run1), 'run 1 report waits for run 1 sims');
            tc.verifyEqual(r2{1}.depend, sort(run2), 'run 2 report waits for run 2 sims');
        end

        function coherentSweepErrorsBeforeSubmission(tc)
            p = base_parameters(tc.TmpDir, 'slurm');
            p.transducer(2) = p.transducer(1);
            p.transducer(1).target_isppa_wcm2 = [10 30];
            tc.verifyError(@() uncertainty_pipeline(p), 'prestus:isppa:coherentSweep');
            tc.verifyEmpty(recorded());
        end

        function plainUncertaintyUnchanged(tc)
            p = base_parameters(tc.TmpDir, 'slurm');
            uncertainty_pipeline(p);
            jobs = recorded();
            tc.verifyNumElements(jobs, 7);
            verify_no_reentry(tc, jobs);
            r = find_jobs(jobs, 'u7-report');
            tc.verifyNumElements(r, 1);
            tc.verifyEqual(r{1}.report_suffix, p.io.output_affix, ...
                'report named after the base affix');
            tc.verifyEqual(r{1}.modules.run_grid_setup, p.modules.run_grid_setup);
        end

        function plainVariantsMediumAndHeating(tc)
            % One job per variant: its medium file is merged in, and only the
            % default and Temp variants run heating.
            p = base_parameters(tc.TmpDir, 'slurm');
            uncertainty_pipeline(p);
            jobs = recorded();
            [name, variants, heating] = variant_list(p.io.output_affix);
            folder = fullfile(get_prestus_path(), 'config', 'uncertainty');
            for v = 1:numel(name)
                sim = find_jobs(jobs, sprintf('u%d-sim-%s_', v + 1, name{v}));
                tc.assertNumElements(sim, 1, name{v});
                tc.verifyEqual(sim{1}.output_affix, variants{v});
                tc.verifyEqual(logical(sim{1}.modules.run_heating_sims), heating(v), name{v});
                tc.verifyEqual(logical(sim{1}.modules.run_thermal_analysis), heating(v), name{v});
                if v > 1
                    file = yaml.loadFile(fullfile(folder, ['config_medium_' name{v} '.yaml']), ...
                        'ConvertToArray', true);
                    tc.verifyEqual(sim{1}.p.medium_properties.skull.sound_speed, ...
                        file.medium_properties.skull.sound_speed, name{v});
                    tc.verifyEqual(sim{1}.p.medium_properties.brain.alpha_power, 1, name{v});
                end
            end
        end

        function pctOnlyInDefaultVariant(tc)
            % With a pseudo-CT, the bound variants keep the shared head
            % geometry (pct stays enabled) but model a homogeneous skull.
            p = base_parameters(tc.TmpDir, 'slurm');
            p.pct.enabled = 1;
            uncertainty_pipeline(p);
            jobs = recorded();
            [name, ~, ~] = variant_list(p.io.output_affix);
            for v = 1:numel(name)
                sim = find_jobs(jobs, sprintf('u%d-sim-%s_', v + 1, name{v}));
                tc.assertNumElements(sim, 1, name{v});
                tc.verifyEqual(sim{1}.p.pct.enabled, 1, name{v});
                tc.verifyEqual(pct_skull_mapped(sim{1}.p), v == 1, name{v});
            end
            r = find_jobs(jobs, 'u7-report');
            tc.verifyEqual(sort(fieldnames(r{1}.p.uncertainty.bound_skull)), sort(name(2:end))');
        end

        function legacyOptionsError(tc)
            p = base_parameters(tc.TmpDir, 'slurm');
            tc.verifyError(@() uncertainty_pipeline(p, struct('liberal_config', 'x.yaml')), ...
                'uncertainty_pipeline:legacyVariant');
            tc.verifyError(@() uncertainty_pipeline(p, struct('affixes', ...
                struct('default', '', 'liberal', '_l', 'conservative', '_c'))), ...
                'uncertainty_pipeline:legacyVariant');
            tc.verifyError(@() uncertainty_pipeline(p, struct('temp_high_config', 'x.yaml')), ...
                'uncertainty_pipeline:variantFile');
            tc.verifyEmpty(recorded(), 'nothing may be submitted before the error');
        end

        function predictOutputsNested(tc)
            p = base_parameters(tc.TmpDir, 'matlab');
            p.transducer(1).target_isppa_wcm2 = [10 30];
            plan = predict_outputs(p);
            tc.verifyEqual(plan.mode, 'uncertainty');
            tc.verifyEqual(plan.nested, 'multi_isppa');
            names = {plan.artifacts.name};
            tc.verifyNumElements(names, 27);
            tc.verifyTrue(any(strcmp(names, sprintf( ...
                'sub-001_%s_desc-uncertainty_report%s_desc-isppa030.html', ...
                p.simulation.medium, p.io.output_affix))));

            p.platform = 'slurm';
            tc.verifyEqual(predict_outputs(p).mode, 'uncertainty', ...
                'same precedence on every platform');
        end

    end
end

% ======================================================================
%% Fixtures
% ======================================================================

function [names, affixes, heating] = variant_list(base)
% Variant names, output affixes and heating flags in run order.
    v = uncertainty_variants();
    names   = {v.name};
    affixes = cellfun(@(d) [base d], {v.desc}, 'UniformOutput', false);
    heating = [v.heating];
end

function p = base_parameters(tmp, platform)
    p = make_minimal_parameters(fullfile(tmp, 'sims'));
    p.platform = platform;
    p.simulation.medium      = 'layered';
    p.simulation.uncertainty = true;
    p.simulation.transducer_coupling = 'coherent';
    p.modules.run_heating_sims     = 1;
    p.modules.run_thermal_analysis = 1;
    p.transducer(1).target_isppa_wcm2 = 30;
end

function write_recorders(stub_dir, start_stub_dir)
% Recorders for the two entry points the chained pipelines call per stage.
    body = { ...
        'function [job_id, jobs] = prestus_pipeline_start(parameters, options) %#ok<INUSD>'
        '    job_id = record_test_job(parameters);'
        '    jobs = struct();'
        'end'};
    write_lines(fullfile(start_stub_dir, 'prestus_pipeline_start.m'), body);
    body = { ...
        'function prestus_pipeline(parameters, options) %#ok<INUSD>'
        '    record_test_job(parameters);'
        'end'};
    write_lines(fullfile(stub_dir, 'prestus_pipeline.m'), body);
    body = { ...
        'function id = record_test_job(p)'
        '    jobs = getappdata(0, ''prestus_test_jobs'');'
        '    id = 1000 + numel(jobs);'
        '    j.id = id;'
        '    j.p = p;'
        '    jobs{end+1} = j;'
        '    setappdata(0, ''prestus_test_jobs'', jobs);'
        'end'};
    write_lines(fullfile(stub_dir, 'record_test_job.m'), body);
end

function write_lines(file, lines)
    fid = fopen(file, 'w');
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end

function rmpath_safe(d)
    if contains(path, d)
        rmpath(d);
    end
end

function rmappdata_safe()
    if isappdata(0, 'prestus_test_jobs')
        rmappdata(0, 'prestus_test_jobs');
    end
end

function jobs = recorded()
% Recorded jobs as flat structs with the fields the checks use.
    raw = getappdata(0, 'prestus_test_jobs');
    jobs = cell(1, numel(raw));
    for k = 1:numel(raw)
        p = raw{k}.p;
        j = struct();
        j.id            = raw{k}.id;
        j.p             = p;
        j.name          = get_field(p, {'hpc', 'job_name'}, '');
        j.output_affix  = get_field(p, {'io', 'output_affix'}, '');
        j.preproc_affix = get_field(p, {'io', 'preproc_affix'}, '<none>');
        j.medium_affix  = get_field(p, {'io', 'medium_affix'}, '<none>');
        j.reuse         = get_field(p, {'io', 'reuse_upstream_caches'}, {});
        j.coupling      = get_field(p, {'simulation', 'transducer_coupling'}, '');
        j.modules       = p.modules;
        j.report_suffix = get_field(p, {'uncertainty', 'report_suffix'}, '<none>');
        j.report_affixes = get_field(p, {'uncertainty', 'affixes'}, struct());
        j.targets       = arrayfun(@(t) t.target_isppa_wcm2, p.transducer, 'UniformOutput', false);
        dep = [get_field(p, {'hpc', 'depend_job_ids'}, []), get_field(p, {'hpc', 'depend_job_id'}, [])];
        j.depend        = unique(dep(:)');
        jobs{k} = j;
    end
end


function out = find_jobs(jobs, varargin)
% Jobs whose name contains every pattern.
    keep = cellfun(@(j) all(cellfun(@(pat) contains(j.name, pat), varargin)), jobs);
    out = jobs(keep);
end

function kind = stage_kind(j)
    m = j.modules;
    if isfield(m, 'uncertainty_report') && m.uncertainty_report
        kind = 'report';
    elseif isfield(m, 'multi_isppa_report') && m.multi_isppa_report
        kind = 'summary';
    elseif m.run_acoustic_sims
        kind = 'acoustic';
    elseif m.run_heating_sims
        kind = 'thermal';
    elseif isfield(m, 'generate_report') && m.generate_report
        kind = 'point';     % acoustic-only per-point job (Trans variants)
    else
        kind = 'stage1';
    end
end

function verify_no_reentry(tc, jobs)
% The launcher re-dispatches each job: none may trigger a chained mode again.
    for k = 1:numel(jobs)
        p = jobs{k}.p;
        tc.verifyFalse(is_multi_isppa_mode(p), sprintf('job %s re-enters multi-ISPPA', jobs{k}.name));
        tc.verifyFalse(is_async_mode(p), sprintf('job %s re-enters async', jobs{k}.name));
        tc.verifyFalse(isfield(p.simulation, 'uncertainty') && logical(p.simulation.uncertainty), ...
            sprintf('job %s re-enters uncertainty', jobs{k}.name));
    end
end

function verify_unique_names(tc, jobs)
    names = cellfun(@(j) j.name, jobs, 'UniformOutput', false);
    tc.verifyEqual(numel(unique(names)), numel(names), 'job names must be unique');
end
