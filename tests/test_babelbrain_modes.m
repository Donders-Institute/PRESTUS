classdef test_babelbrain_modes < matlab.unittest.TestCase
% TEST_BABELBRAIN_MODES  PRESTUS run modes combined with acoustic_backend /
% thermal_backend 'babelbrain', checked without BabelBrain, a GPU or data.
%
% Covers the modes listed as untested with BabelBrain in doc_babelbrain.md
% (Limitations and unsupported PRESTUS modes):
%   - multi_isppa: every derived job keeps the BabelBrain backends; with
%     BabelBrain thermal the per-target thermal jobs run one after another
%     (they share one BabelBrain thermal output file).
%   - uncertainty: each variant's medium_properties reach the material JSON
%     babelbrain_wrapper passes to run_babelbrain_case.py; with ct_mapping
%     the JSON is still passed next to --use_ct.
%   - post-hoc water: the re-entry run keeps the BabelBrain backends and sets
%     simulation.medium 'water' (babelbrain_wrapper then passes --water_only).
%   - water / phantom medium: the label map write_babelbrain_input writes.
%   - Single preset: --aperture and --focal_length reach run_babelbrain_case.py.
%
% prestus_pipeline_start, prestus_pipeline and run_command are replaced by
% recorders (fixtures written to a temp folder at the front of the path), so
% the pipelines build their job graphs and babelbrain_wrapper runs its
% dry-run path (babelbrain.dry_run) up to the Python call, which is recorded
% instead of executed. What BabelBrain then does with the files is checked
% on the cluster by v05/code/smoke/smoke_babelbrain_modes.m.
%
%   Run with:   results = runtests('tests/test_babelbrain_modes.m');
%
% See also: MULTI_ISPPA_PIPELINE, UNCERTAINTY_PIPELINE,
%           POSTHOC_WATER_PARAMETERS, BABELBRAIN_WRAPPER, WRITE_BABELBRAIN_INPUT

    properties
        TmpDir
        StubDir
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir  = tempname;
            tc.StubDir = fullfile(tc.TmpDir, 'stubs');
            mkdir(tc.StubDir);
            write_recorders(tc.StubDir);
            addpath(tc.StubDir);
            setappdata(0, 'prestus_test_jobs', {});
            setappdata(0, 'prestus_test_cmds', {});
            tc.addTeardown(@() rmdir(tc.TmpDir, 's'));
            tc.addTeardown(@() rmappdata_safe());
            tc.addTeardown(@() rmpath_safe(tc.StubDir));
        end
    end

    methods (Test, TestTags = {'babelbrain'})

        %% multi_isppa

        function multiIsppaJobsKeepBabelBrain(tc)
            p = bb_parameters(tc.TmpDir, 'slurm');
            p.transducer(1).target_isppa_wcm2 = [10 30];
            out = multi_isppa_pipeline(p);
            jobs = recorded();
            tc.verifyNumElements(jobs, 4, 'acoustic + 2 thermal + summary');
            for k = 1:numel(jobs)
                tc.verifyEqual(jobs{k}.p.acoustic_backend, 'babelbrain', jobs{k}.name);
                tc.verifyEqual(jobs{k}.p.thermal_backend,  'babelbrain', jobs{k}.name);
                tc.verifyFalse(is_multi_isppa_mode(jobs{k}.p), 'no job may re-enter the sweep');
            end
            th = find_jobs(jobs, 'thermal');
            tc.verifyNumElements(th, 2);
            tc.verifyEqual(cellfun(@(j) j.p.transducer(1).target_isppa_wcm2, th), [10 30]);
            tc.verifyTrue(all(cellfun(@(j) ~j.p.modules.run_acoustic_sims, th)), ...
                'thermal jobs read the stage-1 acoustic cache (and its DataForSim.h5)');
            tc.verifyEqual(sort(out.all), sort(cellfun(@(j) j.id, jobs)));
        end

        function multiIsppaBabelBrainThermalSerialised(tc)
            % BabelBrain names its thermal output after the acoustic file only
            % (<stem>-ThermalField_AllCombinations.h5), so two targets running at
            % once would overwrite each other's result.
            p = bb_parameters(tc.TmpDir, 'slurm');
            p.transducer(1).target_isppa_wcm2 = [10 20 30];
            multi_isppa_pipeline(p);
            jobs = recorded();
            ac = find_jobs(jobs, 'acoustic');
            th = find_jobs(jobs, 'thermal');
            tc.verifyNumElements(th, 3);
            tc.verifyEqual(th{1}.depend, ac{1}.id);
            tc.verifyEqual(sort(th{2}.depend), sort([ac{1}.id, th{1}.id]));
            tc.verifyEqual(sort(th{3}.depend), sort([ac{1}.id, th{2}.id]));
        end

        function multiIsppaKwaveThermalStaysParallel(tc)
            p = bb_parameters(tc.TmpDir, 'slurm');
            p.thermal_backend = 'kwave';
            p.transducer(1).target_isppa_wcm2 = [10 30];
            multi_isppa_pipeline(p);
            jobs = recorded();
            ac = find_jobs(jobs, 'acoustic');
            th = find_jobs(jobs, 'thermal');
            for k = 1:numel(th)
                tc.verifyEqual(th{k}.depend, ac{1}.id, 'k-Wave thermal jobs wait for stage 1 only');
            end
        end

        function multiIsppaMatlabKeepsBabelBrain(tc)
            p = bb_parameters(tc.TmpDir, 'matlab');
            p.transducer(1).target_isppa_wcm2 = [10 30];
            multi_isppa_pipeline(p, struct('generate_summary_report', false));
            jobs = recorded();
            tc.verifyNumElements(jobs, 3);
            tc.verifyTrue(jobs{1}.p.modules.run_acoustic_sims == 1);
            for k = 1:numel(jobs)
                tc.verifyEqual(jobs{k}.p.acoustic_backend, 'babelbrain');
                tc.verifyEqual(jobs{k}.p.thermal_backend,  'babelbrain');
            end
        end

        %% uncertainty

        function uncertaintyVariantsReachMaterialJson(tc)
            p = bb_parameters(tc.TmpDir, 'slurm');
            p.simulation.uncertainty = true;
            uncertainty_pipeline(p);
            jobs = recorded();
            base = p.io.output_affix;
            names = {'default', 'liberal', 'conservative'};
            affixes = {base, [base '_desc-liberal'], [base '_desc-conservative']};
            repo = fileparts(fileparts(mfilename('fullpath')));
            expected = {p.medium_properties, ...
                yaml_medium(fullfile(repo, 'config', 'uncertainty', 'config_medium_liberal.yaml')), ...
                yaml_medium(fullfile(repo, 'config', 'uncertainty', 'config_medium_conservative.yaml'))};

            [masks, fp, tp] = phantom_masks(p, {'skin', 'skull_cortical', 'skull_trabecular', 'brain'});
            c_cort = zeros(1, 3);
            for v = 1:3
                sims = jobs(cellfun(@(j) strcmp(j.output_affix, affixes{v}) && ...
                    j.p.modules.run_acoustic_sims == 1, jobs));
                tc.assertNumElements(sims, 1, names{v});
                pv = sims{1}.p;
                tc.verifyEqual(pv.acoustic_backend, 'babelbrain', names{v});
                tc.verifyEqual(pv.thermal_backend,  'babelbrain', names{v});

                [js, argv] = run_wrapper(tc, pv, masks, fp, tp, []);
                for t = {'brain', 'skin', 'skull_cortical', 'skull_trabecular'}
                    for f = {'sound_speed', 'density', 'alpha_coeff', 'alpha_power'}
                        if isfield(expected{v}, t{1}) && isfield(expected{v}.(t{1}), f{1})
                            tc.verifyEqual(js.(t{1}).(f{1}), expected{v}.(t{1}).(f{1}), 'RelTol', 1e-12, ...
                                sprintf('%s: %s.%s in the material JSON', names{v}, t{1}, f{1}));
                        end
                    end
                end
                tc.verifyTrue(has_flag(argv, '--medium_properties_json'), ...
                    'the JSON is passed to run_babelbrain_case.py');
                tc.verifyFalse(has_flag(argv, '--use_ct'));
                c_cort(v) = js.skull_cortical.sound_speed;
            end
            tc.verifyEqual(numel(unique(c_cort)), 3, 'variants differ in cortical sound speed');
        end

        function uncertaintyWithCtMapping(tc)
            % With ct_mapping, BabelBrain maps bone from HU; the JSON is still
            % passed (its Skin/Brain entries remain BabelBrain's soft-tissue
            % materials on the CT path, BabelIntegrationBASE.py).
            p = bb_parameters(tc.TmpDir, 'matlab');
            p.babelbrain.ct_mapping = 'Webb-Marsac';
            [masks, fp, tp] = phantom_masks(p, {'skin', 'skull', 'brain'});
            pct = zeros(size(masks));
            pct(masks == find(strcmp(fieldnames(p.medium_properties), 'skull'))) = 1200;
            [js, argv] = run_wrapper(tc, p, masks, fp, tp, pct);
            tc.verifyTrue(has_flag(argv, '--use_ct'));
            k = find(strcmp(argv, '--ct_mapping'));
            tc.verifyEqual(argv{k + 1}, 'Webb-Marsac');
            tc.verifyTrue(has_flag(argv, '--medium_properties_json'));
            tc.verifyEqual(js.brain.sound_speed, p.medium_properties.brain.sound_speed);
            tc.verifyTrue(isfile(fullfile(p.io.dir_cache, 'sub-001', ...
                sprintf('%s_250kHz_%dPPW_CT.nii.gz', target_name(p), ppw_of(argv)))), 'CT.nii.gz written');
        end

        %% post-hoc water

        function posthocWaterKeepsBabelBrain(tc)
            p = bb_parameters(tc.TmpDir, 'slurm');
            p.modules.run_posthoc_water_sims = 1;
            p.modules.run_heating_sims = 1;
            p.grid.dims = [40 40 60];
            w = posthoc_water_parameters(p, 1);
            tc.verifyEqual(w.simulation.medium, 'water');
            tc.verifyEqual(w.acoustic_backend, 'babelbrain');
            tc.verifyEqual(w.thermal_backend,  'babelbrain');
            tc.verifyEqual(w.modules.run_posthoc_water_sims, 0, 'no further re-entry');
            tc.verifyEqual(w.modules.run_heating_sims, 0);
            tc.verifyEqual(w.grid.default_dims, p.grid.dims, 'same grid as the run');
            tc.verifyEqual(w.io.output_affix, p.io.output_affix, ...
                'same BabelBrain target name as the run (and its free-water baseline)');
            tc.verifyEqual(w.babelbrain, p.babelbrain);
        end

        function posthocWaterWrapperRunsWaterOnly(tc)
            p = bb_parameters(tc.TmpDir, 'matlab');
            p.grid.dims = [40 40 60];
            w = posthoc_water_parameters(p, 1);
            masks = ones(w.grid.dims);   % grid_tissue_setup's water mask
            [~, argv, nii] = run_wrapper(tc, w, masks, [20 20 45], [20 20 1], []);
            tc.verifyTrue(has_flag(argv, '--water_only'));
            lab = niftiread(nii);
            tc.verifyEqual(unique(lab(:))', uint8([0 5]), 'water-only label map: background + focal marker');
        end

        %% single-element preset geometry

        function singlePresetPassesGeometry(tc)
            % BabelIntegrationSingle.RUN_SIM.RunCases reads Aperture and
            % FocalLength from its kwargs, and Babel_SingleTx/default.yaml has
            % neither: the wrapper passes both from transducer.annular.
            p = bb_parameters(tc.TmpDir, 'matlab');
            p.babelbrain.tx_system = 'Single';
            p.grid.dims = [40 40 60];
            [masks, fp, tp] = phantom_masks(p, {'skin', 'skull_cortical', 'brain'});
            [~, argv] = run_wrapper(tc, p, masks, fp, tp, []);
            a = p.transducer(1).annular;
            tc.verifyEqual(str2double(argv{find(strcmp(argv, '--aperture')) + 1}), ...
                max(a.elem_od_mm) / 1000, 'AbsTol', 1e-6);
            tc.verifyEqual(str2double(argv{find(strcmp(argv, '--focal_length')) + 1}), ...
                a.curv_radius_mm / 1000, 'AbsTol', 1e-6);

            p.babelbrain.tx_system = 'R15473';   % ring presets keep their yaml geometry
            [~, argv] = run_wrapper(tc, p, masks, fp, tp, []);
            tc.verifyFalse(has_flag(argv, '--aperture'));
            tc.verifyFalse(has_flag(argv, '--focal_length'));
        end

        %% water / phantom label maps

        function waterMaskLabels(tc)
            p = bb_parameters(tc.TmpDir, 'matlab');
            for masks = {ones(40, 40, 60), zeros(40, 40, 60)}
                lab = write_labels(tc, p, masks{1}, [20 20 45], [20 20 1]);
                tc.verifyEqual(unique(lab(:))', uint8([0 5]));
                tc.verifyEqual(nnz(lab == 5), 1);
            end
        end

        function phantomMaskLabels(tc)
            % Each tissue keeps its own BabelBrain label (1 skin, 2 cortical,
            % 3 trabecular, 4 brain), the materials run_babelbrain_case.py
            % injects; non-background labels are not all treated as skull.
            p = bb_parameters(tc.TmpDir, 'matlab');
            [masks, fp, tp] = phantom_masks(p, {'skin', 'skull_cortical', 'brain'});
            lab = write_labels(tc, p, masks, fp, tp);
            tc.verifyEqual(unique(lab(:))', uint8([0 1 2 4 5]));
            tc.verifyEqual(nnz(lab == 5), 1);
            % depth order along the beam after the pre-flip: brain, bone, skin, water
            prof = squeeze(lab(round(end/2), round(end/2), :))';
            prof = prof(prof ~= 5);
            order = prof([true, diff(double(prof)) ~= 0]);
            tc.verifyEqual(order, uint8([4 2 1 0]));

            [masks, fp, tp] = phantom_masks(p, {'skin', 'skull_cortical', 'skull_trabecular', 'skull_cortical', 'brain'});
            lab = write_labels(tc, p, masks, fp, tp);
            tc.verifyEqual(unique(lab(:))', uint8([0 1 2 3 4 5]));

            [masks, fp, tp] = phantom_masks(p, {'skull'});
            lab = write_labels(tc, p, masks, fp, tp);
            tc.verifyEqual(unique(lab(:))', uint8([0 2 5]), 'undifferentiated skull -> cortical');
        end

        function phantomUnknownTissueBecomesWater(tc)
            % A medium_properties tissue other than water/skin/skull*/brain
            % has no BabelBrain label: it is written as background (water).
            p = bb_parameters(tc.TmpDir, 'matlab');
            p.medium_properties.gel = p.medium_properties.brain;
            [masks, fp, tp] = phantom_masks(p, {'gel'});
            lab = write_labels(tc, p, masks, fp, tp);
            tc.verifyEqual(unique(lab(:))', uint8([0 5]));
        end

    end
end

% ======================================================================
%% Fixtures
% ======================================================================

function p = bb_parameters(tmp, platform)
    p = make_minimal_parameters(fullfile(tmp, 'sims'));
    p.platform = platform;
    p.simulation.medium = 'layered';
    p.simulation.uncertainty = false;
    p.modules.run_heating_sims     = 1;
    p.modules.run_thermal_analysis = 1;
    p.transducer(1).target_isppa_wcm2 = 30;
    p.acoustic_backend = 'babelbrain';
    p.thermal_backend  = 'babelbrain';
    % babelbrain_wrapper only checks that <root>/TranscranialModeling exists
    bb_root = fullfile(tmp, 'BabelBrain');
    mkdir(fullfile(bb_root, 'TranscranialModeling'));
    p.babelbrain.root              = bb_root;
    p.babelbrain.python_env        = '';
    p.babelbrain.tx_system         = 'R15473';
    p.babelbrain.computing_backend = 4;
    p.babelbrain.device_name       = 'CPU';
    p.babelbrain.ppw               = 3;      % coarse native grid: small test volumes
    p.babelbrain.focal_calibration = false;
    p.babelbrain.dry_run           = true;
    p.babelbrain.ct_mapping        = '';
    p.io.dir_cache = fullfile(tmp, 'cache');
    % no explicit phases: geometric focusing, no weights file
    p.transducer(1).annular = rmfield(p.transducer(1).annular, 'elem_phase_rad');
end

function [masks, fp, tp] = phantom_masks(p, layers)
% 40 x 40 x 60 grid, beam +z from the transducer at z = 1; water up to z = 10,
% then 4-voxel slabs of the given tissues, the last one filling to the end.
    names = fieldnames(p.medium_properties);
    masks = repmat(find(strcmp(names, 'water')), 40, 40, 60);
    z = 11;
    for k = 1:numel(layers)
        z_end = 60;
        if k < numel(layers); z_end = z + 3; end
        masks(:, :, z:z_end) = find(strcmp(names, layers{k}));
        z = z_end + 1;
    end
    fp = [20 20 45];
    tp = [20 20 1];
end

function lab = write_labels(tc, p, masks, fp, tp)
    planimg = struct('t1_header', [], 'transf', [], 'inv_transf', [], 'origin_ras_mm', [0 0 0]);
    out = fullfile(tc.TmpDir, 'labels');
    nii = write_babelbrain_input(masks, fp, p, planimg, out, sprintf('t%d', randi(1e6)), 0, [], tp);
    lab = niftiread(nii);
end

function [js, argv, nii] = run_wrapper(tc, p, masks, fp, tp, pct)
% babelbrain_wrapper dry run with run_command recorded; returns the material
% JSON it wrote and the argv of the BabelBrain call.
    setappdata(0, 'prestus_test_cmds', {});
    planimg = struct('t1_header', [], 'transf', [], 'inv_transf', [], 'origin_ras_mm', [0 0 0]);
    sd = babelbrain_wrapper(p, masks, struct(), fp, tp, planimg, struct(), pct);
    tc.verifyEqual(sd.babelbrain_source, 'dry_run');
    cmds = getappdata(0, 'prestus_test_cmds');
    is_case = cellfun(@(c) any(contains(c, 'run_babelbrain_case.py')) && ~any(strcmp(c, '--get_version')), cmds);
    tc.assertEqual(nnz(is_case), 1, 'one BabelBrain case call');
    argv = cmds{is_case};
    k = find(strcmp(argv, '--medium_properties_json'));
    tc.assertNotEmpty(k, 'material JSON flag');
    js = jsondecode(fileread(argv{k + 1}));
    tgt = target_name(p);
    nii = fullfile(p.io.dir_cache, 'sub-001', ...
        sprintf('%s_250kHz_%dPPW_BabelViscoInput.nii.gz', tgt, ppw_of(argv)));
    tc.verifyTrue(isfile(nii));
end

function t = target_name(p)
    t = strrep(p.io.output_affix, '_', '');
end

function n = ppw_of(argv)
    n = str2double(argv{find(strcmp(argv, '--ppw')) + 1});
end

function tf = has_flag(argv, flag)
    tf = any(strcmp(argv, flag));
end

function mp = yaml_medium(file)
    y = yaml.loadFile(file, 'ConvertToArray', true);
    mp = y.medium_properties;
end

function write_recorders(stub_dir)
% Recorders for the pipeline entry points and for external commands.
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
    write_lines(fullfile(stub_dir, 'run_command.m'), { ...
        'function [status, cmdout] = run_command(argv, varargin) %#ok<INUSD>'
        '    cmds = getappdata(0, ''prestus_test_cmds'');'
        '    cmds{end+1} = argv;'
        '    setappdata(0, ''prestus_test_cmds'', cmds);'
        '    status = 0;'
        '    cmdout = '''';'
        'end'});
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
    for k = {'prestus_test_jobs', 'prestus_test_cmds'}
        if isappdata(0, k{1})
            rmappdata(0, k{1});
        end
    end
end

function jobs = recorded()
    raw = getappdata(0, 'prestus_test_jobs');
    jobs = cell(1, numel(raw));
    for k = 1:numel(raw)
        p = raw{k}.p;
        j = struct();
        j.id = raw{k}.id;
        j.p  = p;
        j.name = get_field(p, {'hpc', 'job_name'}, '');
        j.output_affix = get_field(p, {'io', 'output_affix'}, '');
        dep = [get_field(p, {'hpc', 'depend_job_ids'}, []), get_field(p, {'hpc', 'depend_job_id'}, [])];
        j.depend = unique(dep(:)');
        jobs{k} = j;
    end
end

function v = get_field(s, path, default)
    v = default;
    for k = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{k}); return; end
        s = s.(path{k});
    end
    v = s;
end

function out = find_jobs(jobs, varargin)
    keep = cellfun(@(j) all(cellfun(@(pat) contains(j.name, pat), varargin)), jobs);
    out = jobs(keep);
end
