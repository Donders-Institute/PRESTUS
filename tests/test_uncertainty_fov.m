classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_uncertainty_fov < matlab.unittest.TestCase
% TEST_UNCERTAINTY_FOV  Uncertainty mode with the bound variants in a partial
% acoustic FOV (simulation.uncertainty_fov = 'partial'); no simulation.
%
%   - Mismatch factors and rescaling (UNCERTAINTY_FOV_MISMATCH,
%     UNCERTAINTY_FOV_RESCALE) for every metric kind, with guards and controls.
%   - Job graph: prestus_pipeline_start and prestus_pipeline are replaced by
%     recorders, so UNCERTAINTY_PIPELINE builds its stages without simulating.
%   - Report: GENERATE_UNCERTAINTY_REPORT rescales synthetic variant tables
%     and writes the factor table.
%
%   Run with:   results = runtests('tests/test_uncertainty_fov.m');
%
% See also: UNCERTAINTY_PIPELINE, UNCERTAINTY_FOV_MISMATCH, UNCERTAINTY_FOV_RESCALE

    properties
        TmpDir
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir = tempname;
            stub_dir  = fullfile(tc.TmpDir, 'stubs');
            mkdir(stub_dir);
            write_recorders(stub_dir);
            addpath(stub_dir);
            setappdata(0, 'prestus_test_jobs', {});
            tc.addTeardown(@() rmdir(tc.TmpDir, 's'));
            tc.addTeardown(@() rmappdata_safe());
            tc.addTeardown(@() rmpath_safe(stub_dir));
        end
    end

    methods (Test)

        % ---- factors and rescaling --------------------------------------

        function everyMetricKindIsClassified(tc)
            ratio = {'Isppa', 'Isppa_brain', 'Isppa_skull', 'Isppa_skin', 'Isppa_after_exitplane', ...
                     'Ispta', 'Ipa_target', 'MI', 'MI_brain', 'MI_tc', 'Psptp', 'Psptp_brain', ...
                     'Ptp_target', 'riseT_brain', 'riseT_skull', 'riseT_skin', ...
                     'CEM43_brain', 'CEM43_skull', 'CEM43_skin'};
            rise  = {'maxT', 'maxT_brain', 'maxT_skull', 'maxT_skin', 'endT_brain'};
            none  = {'IsppaScale', 'real_focal_distance_bowl_mm', 'real_focal_distance_ep_mm', ...
                     'focal_vol_6dB_mm3', 'focal_fwhm_axial_mm', 'halfmax_ISPPA_volume_brain_mm3', ...
                     'aircav_flag', 'subject_id'};
            for n = ratio, tc.verifyEqual(uncertainty_fov_metric_kind(n{1}), 'ratio', n{1}); end
            for n = rise,  tc.verifyEqual(uncertainty_fov_metric_kind(n{1}), 'rise',  n{1}); end
            for n = none,  tc.verifyEqual(uncertainty_fov_metric_kind(n{1}), 'none',  n{1}); end
        end

        function factorsPerKind(tc)
            [full, fov] = default_tables();
            m = uncertainty_fov_mismatch(full, fov, struct('brain', 37, 'skull', 36));
            f = @(name) m.factor(strcmp(m.metric, name));
            tc.verifyEqual(f('Isppa_brain'), 10 / 8, 'AbsTol', 1e-12);
            tc.verifyEqual(f('MI_brain'), 0.5 / 0.4, 'AbsTol', 1e-12);
            tc.verifyEqual(f('riseT_brain'), 1.0 / 0.8, 'AbsTol', 1e-12);
            % absolute temperature: the factor of the rise above the tissue baseline
            tc.verifyEqual(f('maxT_brain'), (38.0 - 37) / (37.8 - 37), 'AbsTol', 1e-9);
            tc.verifyEqual(f('maxT_skull'), (38.0 - 36) / (37.0 - 36), 'AbsTol', 1e-9);
            tc.verifyEqual(m.baseline(strcmp(m.metric, 'maxT')), 37, 'no tissue suffix: 37');
            tc.verifyTrue(isnan(f('real_focal_distance_bowl_mm')), 'geometry has no factor');
            tc.verifyEqual(m.diff_pct(strcmp(m.metric, 'Isppa_brain')), -20, 'AbsTol', 1e-9);
        end

        function rescaleMapsTheDefaultOntoItself(tc)
            % The partial-FOV default, rescaled, must give the whole-grid default
            % for every rescaled metric: the defining property of the factors.
            [full, fov] = default_tables();
            m = uncertainty_fov_mismatch(full, fov, struct('brain', 37, 'skull', 36));
            out = uncertainty_fov_rescale(fov, m);
            for r = 1:height(m)
                name = m.metric{r};
                if strcmp(m.kind{r}, 'none')
                    tc.verifyEqual(out.(name), fov.(name), [name ' must stay unscaled']);
                else
                    tc.verifyEqual(out.(name), full.(name), 'AbsTol', 1e-9, name);
                end
            end
        end

        function rescaleAppliesFactorsToAVariant(tc)
            [full, fov] = default_tables();
            m = uncertainty_fov_mismatch(full, fov, struct('brain', 37, 'skull', 36));
            v = fov;
            v.Isppa_brain = 16; v.maxT_brain = 38.6; v.real_focal_distance_bowl_mm = 61;
            v.other_column = 5;
            out = uncertainty_fov_rescale(v, m);
            tc.verifyEqual(out.Isppa_brain, 16 * 1.25, 'AbsTol', 1e-12);
            tc.verifyEqual(out.maxT_brain, 37 + 1.6 * 1.25, 'AbsTol', 1e-9);
            tc.verifyEqual(out.real_focal_distance_bowl_mm, 61);
            tc.verifyEqual(out.other_column, 5, 'columns outside the factor table are untouched');
        end

        function identicalDefaultsLeaveVariantsUnchanged(tc)
            % Negative control: no mismatch, no change.
            full = default_tables();
            m = uncertainty_fov_mismatch(full, full, 37);
            scaled = m.factor(~strcmp(m.kind, 'none'));
            tc.verifyEqual(scaled, ones(size(scaled)), 'AbsTol', 1e-12);
            v = full; v.Isppa_brain = 3; v.maxT_brain = 39;
            tc.verifyEqual(uncertainty_fov_rescale(v, m), v);
        end

        function invalidValuesGiveNoFactor(tc)
            [full, fov] = default_tables();
            fov.Isppa_brain = 0;        % zero partial-FOV value
            fov.MI_brain    = NaN;      % missing
            fov.maxT_brain  = 37;       % no rise in the partial FOV
            full.CEM43_brain = Inf;
            m = uncertainty_fov_mismatch(full, fov, 37);
            for name = {'Isppa_brain', 'MI_brain', 'maxT_brain', 'CEM43_brain'}
                tc.verifyTrue(isnan(m.factor(strcmp(m.metric, name{1}))), name{1});
            end
            v = fov; v.Isppa_brain = 12;
            out = uncertainty_fov_rescale(v, m);
            tc.verifyEqual(out.Isppa_brain, 12, 'a metric without a factor stays as simulated');
        end

        % ---- job graph ---------------------------------------------------

        function matlabStagesAndGrids(tc)
            for configured = [false true]
                setappdata(0, 'prestus_test_jobs', {});
                p = base_parameters(tc.TmpDir, 'matlab');
                p.simulation.uncertainty_fov = 'partial';
                if configured
                    p.grid.acoustic_fov_length_mm   = 110;
                    p.grid.acoustic_fov_diameter_mm = 70;
                    want = {110, 70};
                else
                    want = {'auto', 'auto'};
                end
                uncertainty_pipeline(p);
                jobs = recorded();
                % stage 1, five variants, the partial-FOV default, the report
                tc.verifyNumElements(jobs, 8);
                names = {uncertainty_variants().name};
                for v = 1:numel(names)
                    g = jobs{v + 1}.p.grid;
                    if strcmp(names{v}, 'default')
                        tc.verifyEmpty(g.acoustic_fov_length_mm, 'default runs on the whole grid');
                        tc.verifyEmpty(g.acoustic_fov_diameter_mm, 'default runs on the whole grid');
                    else
                        tc.verifyEqual({g.acoustic_fov_length_mm, g.acoustic_fov_diameter_mm}, want, names{v});
                    end
                end
                ref = jobs{7};
                tc.verifyEqual(ref.p.io.output_affix, [p.io.output_affix '_desc-defaultfov']);
                tc.verifyEqual({ref.p.grid.acoustic_fov_length_mm, ref.p.grid.acoustic_fov_diameter_mm}, want);
                tc.verifyEqual(get_field(ref.p, {'medium_properties'}, []), ...
                    get_field(jobs{2}.p, {'medium_properties'}, []), 'the reference uses the default medium');
                tc.verifyEqual(logical(ref.p.modules.run_heating_sims), true);
                tc.verifyEqual(ref.p.io.preproc_affix, '', 'reference reuses the stage-1 head/grid');
                tc.verifyEqual(ref.p.io.medium_affix, ref.p.io.output_affix);
                rep = jobs{8};
                tc.verifyEqual(rep.p.modules.uncertainty_report, 1);
                tc.verifyEqual(rep.p.uncertainty.fov_reference_affix, ref.p.io.output_affix);
                for k = 1:numel(jobs)
                    tc.verifyFalse(logical(jobs{k}.p.simulation.uncertainty), 'no job re-enters uncertainty mode');
                end
            end
        end

        function slurmGraph(tc)
            p = base_parameters(tc.TmpDir, 'slurm');
            p.simulation.uncertainty_fov = 'partial';
            out = uncertainty_pipeline(p);
            jobs = recorded();
            tc.verifyNumElements(jobs, 8);
            tc.verifyEqual(sort(out.all), sort(cellfun(@(j) j.id, jobs)));
            name = cellfun(@(j) j.p.hpc.job_name, jobs, 'UniformOutput', false);
            tc.verifyNumElements(unique(name), 8, 'job names are unique');
            s1  = jobs{contains(name, 'u1-preproc')};
            ref = jobs{contains(name, 'u7-sim-default_fov')};
            rep = jobs{contains(name, 'u8-report')};
            sims = jobs(contains(name, '-sim-'));
            tc.verifyNumElements(sims, 6);
            for k = 1:numel(sims)
                tc.verifyEqual(sims{k}.p.hpc.depend_job_id, s1.id, 'simulations wait for stage 1');
            end
            tc.verifyEqual(sort(rep.p.hpc.depend_job_ids), sort(cellfun(@(j) j.id, sims)), ...
                'the report waits for all six simulations');
            tc.verifyEqual(rep.p.uncertainty.fov_reference_affix, ref.p.io.output_affix);
        end

        function configModeIsUnchanged(tc)
            % Negative control: without 'partial', no extra run and the
            % configured FOV reaches every variant untouched.
            for mode = {'', 'config'}
                setappdata(0, 'prestus_test_jobs', {});
                p = base_parameters(tc.TmpDir, 'matlab');
                if ~isempty(mode{1}), p.simulation.uncertainty_fov = mode{1}; end
                p.grid.acoustic_fov_length_mm   = 110;
                p.grid.acoustic_fov_diameter_mm = 70;
                uncertainty_pipeline(p);
                jobs = recorded();
                tc.verifyNumElements(jobs, 7);
                for k = 2:6
                    tc.verifyEqual(jobs{k}.p.grid.acoustic_fov_length_mm, 110);
                    tc.verifyEqual(jobs{k}.p.grid.acoustic_fov_diameter_mm, 70);
                end
                tc.verifyFalse(isfield(jobs{7}.p.uncertainty, 'fov_reference_affix'));
            end
        end

        function customReferenceAffix(tc)
            p = base_parameters(tc.TmpDir, 'matlab');
            p.simulation.uncertainty_fov = 'partial';
            options.affixes = struct('default', '_a', 'temp_high', '_b', 'temp_low', '_c', ...
                'trans_high', '_d', 'trans_low', '_e', 'default_fov', '_ref');
            uncertainty_pipeline(p, options);
            jobs = recorded();
            tc.verifyEqual(jobs{7}.p.io.output_affix, '_ref');
            tc.verifyEqual(jobs{8}.p.uncertainty.fov_reference_affix, '_ref');
        end

        function unsupportedCombinationsError(tc)
            p = base_parameters(tc.TmpDir, 'matlab');
            p.simulation.uncertainty_fov = 'partial';
            q = p; q.transducer(1).target_isppa_wcm2 = [10 30];
            tc.verifyError(@() uncertainty_pipeline(q), 'uncertainty_pipeline:fovNested');
            q = p; q.simulation.transducer_coupling = 'async';       % async needs two stances
            q.transducer(2) = q.transducer(1);
            q.transducer(2).trans_pos = q.transducer(1).trans_pos + [4 0 0];
            tc.verifyError(@() uncertainty_pipeline(q), 'uncertainty_pipeline:fovNested');
            options.sequential_configs.config_2 = p;
            tc.verifyError(@() uncertainty_pipeline(p, options), 'uncertainty_pipeline:fovSequential');
            q = p; q.simulation.uncertainty_fov = 'full';
            tc.verifyError(@() uncertainty_pipeline(q), 'uncertainty_pipeline:fovMode');
            tc.verifyEmpty(recorded(), 'nothing is run or submitted before the error');
        end

        % ---- report ------------------------------------------------------

        function reportRescalesBoundVariants(tc)
            p = base_parameters(tc.TmpDir, 'matlab');
            p = resolve_io_dirs(p);
            out_dir = p.io.dir_output;
            if ~isfolder(out_dir), mkdir(out_dir); end
            [full, fov] = default_tables();
            variants = uncertainty_variants();
            affixes  = struct();
            csv = @(affix) fullfile(out_dir, sprintf('sub-%03d_%s%s.csv', p.subject_id, p.simulation.medium, affix));
            for v = 1:numel(variants)
                affixes.(variants(v).name) = variants(v).desc;
                if strcmp(variants(v).name, 'default')
                    writetable(full, csv(variants(v).desc));
                else
                    t = fov; t.Isppa_brain = 8 + v;      % partial-FOV values of this variant
                    writetable(t, csv(variants(v).desc));
                end
            end
            writetable(fov, csv('_desc-defaultfov'));
            p.uncertainty.fov_reference_affix = '_desc-defaultfov';
            p.uncertainty.report_suffix = '';
            report = generate_uncertainty_report(p, affixes);
            tc.assertTrue(isfile(report), 'report written');
            html = fileread(report);
            tc.verifyTrue(contains(html, 'Partial-FOV Rescaling'));

            mm_file = fullfile(out_dir, sprintf('sub-%03d_%s_desc-uncertainty_fovmismatch.csv', ...
                p.subject_id, p.simulation.medium));
            tc.assertTrue(isfile(mm_file), 'factor table written');
            mm = readtable(mm_file);
            row = strcmp(mm.metric, 'Isppa_brain');
            tc.verifyEqual(mm.factor(row), 1.25, 'AbsTol', 1e-9);
            for v = 2:numel(variants)
                tc.verifyEqual(mm.([variants(v).name '_fov'])(row), 8 + v, 'AbsTol', 1e-9);
                tc.verifyEqual(mm.([variants(v).name '_scaled'])(row), (8 + v) * 1.25, 'AbsTol', 1e-9);
            end
            % the variant CSVs keep their partial-FOV values
            t = readtable(csv(variants(2).desc));
            tc.verifyEqual(t.Isppa_brain, 10, 'AbsTol', 1e-9);

            % Negative control: without the reference affix nothing is rescaled
            p.uncertainty = rmfield(p.uncertainty, 'fov_reference_affix');
            delete(mm_file);
            html = fileread(generate_uncertainty_report(p, affixes));
            tc.verifyFalse(contains(html, 'Partial-FOV Rescaling'));
            tc.verifyFalse(isfile(mm_file));
        end

        function reportWithoutReferenceStaysUnscaled(tc)
            p = base_parameters(tc.TmpDir, 'matlab');
            p = resolve_io_dirs(p);
            if ~isfolder(p.io.dir_output), mkdir(p.io.dir_output); end
            full = default_tables();
            writetable(full, fullfile(p.io.dir_output, sprintf('sub-%03d_%s.csv', p.subject_id, p.simulation.medium)));
            p.uncertainty.fov_reference_affix = '_desc-defaultfov';   % its CSV does not exist
            html = fileread(generate_uncertainty_report(p));
            tc.verifyTrue(contains(html, 'shown unscaled'));
        end

    end
end

% =========================================================================

function [full, fov] = default_tables()
% Metrics of the default medium on the whole grid and in the partial FOV.
    full = table(1, 10, 0.5, 1.0, 38.0, 38.0, 38.2, 0.02, 60, ...
        'VariableNames', {'subject_id', 'Isppa_brain', 'MI_brain', 'riseT_brain', 'maxT_brain', ...
                          'maxT_skull', 'maxT', 'CEM43_brain', 'real_focal_distance_bowl_mm'});
    fov = full;
    fov.Isppa_brain = 8;   fov.MI_brain = 0.4;  fov.riseT_brain = 0.8;
    fov.maxT_brain  = 37.8; fov.maxT_skull = 37.0; fov.maxT = 38.0;
    fov.CEM43_brain = 0.01; fov.real_focal_distance_bowl_mm = 59;
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

function write_recorders(stub_dir)
% Recorders for the two entry points UNCERTAINTY_PIPELINE calls per stage.
    write_lines(fullfile(stub_dir, 'prestus_pipeline_start.m'), { ...
        'function [job_id, jobs] = prestus_pipeline_start(parameters, options) %#ok<INUSD>'
        '    job_id = record_test_job(parameters);'
        '    jobs = struct();'
        'end'});
    write_lines(fullfile(stub_dir, 'prestus_pipeline.m'), { ...
        'function prestus_pipeline(parameters, options) %#ok<INUSD>'
        '    record_test_job(parameters);'
        'end'});
    write_lines(fullfile(stub_dir, 'record_test_job.m'), { ...
        'function id = record_test_job(p)'
        '    jobs = getappdata(0, ''prestus_test_jobs'');'
        '    id = 1000 + numel(jobs);'
        '    j.id = id;'
        '    j.p = p;'
        '    jobs{end+1} = j;'
        '    setappdata(0, ''prestus_test_jobs'', jobs);'
        'end'});
end

function write_lines(file, lines)
    fid = fopen(file, 'w');
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end

function jobs = recorded()
    jobs = getappdata(0, 'prestus_test_jobs');
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
