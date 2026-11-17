classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_placement_report < matlab.unittest.TestCase
% TEST_PLACEMENT_REPORT  Tests for generate_placement_report and the PlanTUS
%                        criteria evaluation (plantus_evaluation.py).
%
%   'report'      MATLAB only: report written / skipped / disabled, PlanTUS
%                 section without an evaluation.
%   'evaluation'  Runs plantus_evaluation.py on a synthetic spherical head
%                 and renders its JSON in the report. Needs Python 3 with
%                 numpy, scipy and matplotlib: PRESTUS_TEST_PYTHON, or the
%                 first python3/python on the PATH that has them. Skipped
%                 otherwise.
%
% Run:
%   results = runtests('tests/test_placement_report.m');

    properties
        TmpDir (1,:) char
    end

    methods (TestMethodSetup)
        function create_tmp_dir(tc)
            tc.TmpDir = tempname();
            mkdir(tc.TmpDir);
        end
    end

    methods (TestMethodTeardown)
        function remove_tmp_dir(tc)
            if isfolder(tc.TmpDir); rmdir(tc.TmpDir, 's'); end
        end
    end

    methods (Test, TestTags = {'report'})

        function test_manual_report_written(tc)
            p = make_params(tc.TmpDir, 'manual');
            r = generate_placement_report(p);
            tc.verifyTrue(isfile(r), 'Report file not written.');
            html = fileread(r);
            tc.verifySubstring(html, 'PRESTUS Placement Report');
            tc.verifySubstring(html, '[30, 30, 50]');          % trans_pos
            tc.verifyTrue(isfile(fullfile(p.io.dir_img, 'sub-001_placement_manual.png')), ...
                'T1 overlay not drawn for manual mode.');
            tc.verifyEmpty(strfind(html, 'class="error-notice"'), 'A report section failed.');
        end

        function test_flags_are_plain_tinted_boxes(tc)
            % every .pl-flag rule: tint only, no coloured left accent bar
            html = fileread(generate_placement_report(make_params(tc.TmpDir, 'manual')));
            rules = regexp(html, '\.pl-flag[^{]*\{[^}]*\}', 'match');
            tc.assertNotEmpty(rules, 'No .pl-flag rules found.');
            for k = 1:numel(rules)
                tc.verifyEmpty(regexp(rules{k}, 'border-(left|color)', 'once'), rules{k});
            end
        end

        function test_no_output_dir_skips_quietly(tc)
            p = make_params(tc.TmpDir, 'manual');
            p = rmfield(p, 'io');
            r = tc.verifyWarningFree(@() generate_placement_report(p));
            tc.verifyEmpty(r);
        end

        function test_report_disabled(tc)
            p = make_params(tc.TmpDir, 'manual');
            p.placement.report = false;
            preproc_transducer_placement(p);
            tc.verifyEmpty(dir(fullfile(p.io.dir_reports, '*.html')));
        end

        function test_report_written_by_placement_stage(tc)
            p = make_params(tc.TmpDir, 'manual');
            preproc_transducer_placement(p);
            tc.verifyTrue(isfile(fullfile(p.io.dir_reports, 'sub-001_placement_manual.html')));
        end

        function test_plantus_without_evaluation_shows_placeholder(tc)
            p = make_params(tc.TmpDir, 'plantus');
            html = fileread(generate_placement_report(p, 'plantus'));
            tc.verifySubstring(html, 'No PlanTUS evaluation found');
        end
    end

    methods (Test, TestTags = {'evaluation'})

        function test_plantus_evaluation_end_to_end(tc)
            py = find_python_with_matplotlib();
            tc.assumeNotEmpty(py, 'No Python 3 with numpy/scipy/matplotlib — skipping.');

            p = make_params(tc.TmpDir, 'plantus');
            p.path.sim = tc.TmpDir;                   % placement runs: <sim>/placement/...
            eval_dir = placement_run_dir(p, 'plantus');
            mkdir(eval_dir);
            metrics = fullfile(eval_dir, 'plantus_metrics.mat');
            write_sphere_metrics(metrics);

            evalpy = fullfile(prestus_root(), 'functions', 'neuronav', 'plantus_evaluation.py');
            [status, out] = system(sprintf('"%s" "%s" "%s"', py, evalpy, metrics));
            tc.assertEqual(status, 0, out);

            ev = jsondecode(fileread(fullfile(eval_dir, 'plantus_evaluation.json')));
            tc.verifyEqual(ev.n_vertices, numel(load(metrics).skin_target_distances));
            tc.verifyNumElements(ev.sweep, 12);                % configured, equal, 5 heavy, 5 only
            tc.verifyTrue(isfile(fullfile(eval_dir, 'plantus_criteria_maps.png')));
            codes = arrayfun(@(f) f.code, ev.flags, 'UniformOutput', false);
            tc.verifyTrue(ismember('intersection_sparse', codes), ...
                'Small target should flag the intersection criterion.');

            html = fileread(generate_placement_report(p, 'plantus'));
            tc.verifySubstring(html, 'Criteria at the placed spot');
            tc.verifySubstring(html, 'Weight sweep');
            tc.verifySubstring(html, 'intersection is non-zero');
            tc.verifyEmpty(strfind(html, 'class="error-notice"'), 'A report section failed.');
            tc.verifyGreaterThanOrEqual(numel(strfind(html, 'data:image/png')), 3, ...
                'Expected T1 overlay and both evaluation figures.');
        end
    end
end

% ── Fixtures ──────────────────────────────────────────────────────────────

function p = make_params(root, mode)
    anat = fullfile(root, 'anat');
    mkdir(anat);
    vol = zeros(64, 64, 64, 'single');
    vol(12:52, 12:52, 12:58) = 100;
    niftiwrite(vol, fullfile(anat, 'sub-001_T1w.nii'));

    p = struct();
    p.subject_id = 1;
    p.path.anat = anat;
    p.path.t1_pattern = 'sub-%03d_T1w.nii';
    p.path.seg = fullfile(root, 'seg');
    p.io.dir_output = fullfile(root, 'out');
    p.io.dir_img = fullfile(root, 'out', 'img');
    p.io.dir_reports = fullfile(root, 'out');
    p.io.output_affix = '';
    p.io.deface_plots = 0;
    p.placement.mode = mode;
    if strcmp(mode, 'plantus')
        p.placement.plantus.target_name = 'test_target';
        p.placement.plantus.target_ras = [5 0 0];
    end
    p.transducer(1).trans_pos = [30 30 50];
    p.transducer(1).focus_pos = [32 32 32];
end

function write_sphere_metrics(path)
% Spherical "head" (r = 70 mm) with a small spherical target (r = 3 mm).
    [X, Y, Z] = meshgrid(-80:2:80);
    fv = isosurface(X, Y, Z, sqrt(X.^2 + Y.^2 + Z.^2), 70);
    v = fv.vertices;
    n = v ./ vecnorm(v, 2, 2);                           % outward normals
    c = [20 0 0];
    to_c = c - v;
    d = vecnorm(to_c, 2, 2);
    ang = acosd(max(-1, min(1, sum(-n .* to_c, 2) ./ d)));
    t = sum(-n .* to_c, 2);
    perp2 = d.^2 - t.^2;
    inter = 2 * sqrt(max(0, 9 - perp2)) .* (t > 0 & perp2 < 9);
    rng(0);
    S = struct();
    S.skin_coordinates = v;
    S.skin_normals = n;
    S.skin_faces = int32(fv.faces - 1);                  % 0-based, as the launchers write
    S.target_center = c;
    S.skin_target_distances = d;
    S.skin_target_angles = ang;
    S.skin_target_intersections = inter;
    S.skin_skull_angles = 20 * rand(size(d));
    th = 4 + 2 * rand(size(d));
    th(v(:, 3) < -40) = NaN;
    S.skull_thickness = th;
    S.avoidance = double(v(:, 3) > -50);
    S.score = zeros(size(d));
    S.best_vertex = 0;
    [~, sel] = min(vecnorm(v - [70 0 0], 2, 2));
    S.selected_vertex = sel - 1;
    S.selection_mode = 'auto';
    S.min_distance = 30;
    S.max_distance = 110;
    S.optimal_distance = 60;
    S.weights = [0.2 0.2 0.2 0.2 0.2];
    save(path, '-struct', 'S', '-v7');
end

function r = prestus_root()
    r = fileparts(fileparts(mfilename('fullpath')));
end

function py = find_python_with_matplotlib()
    py = '';
    cands = {getenv('PRESTUS_TEST_PYTHON'), 'python3', 'python'};
    for i = 1:numel(cands)
        c = cands{i};
        if isempty(c); continue; end
        if ~isfile(c)
            [s, out] = system(sprintf('which %s 2>/dev/null', c));
            if s ~= 0 || isempty(strtrim(out)); continue; end
            c = strtrim(out);
        end
        [s, ~] = system(sprintf('"%s" -c "import numpy, scipy, matplotlib" 2>/dev/null', c));
        if s == 0; py = c; return; end
    end
end
