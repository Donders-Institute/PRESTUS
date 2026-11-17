classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_debug_reports < matlab.unittest.TestCase
% TEST_DEBUG_REPORTS  Debug settings, debug manifest and report, head-model QC records and report
%
%   Synthetic data only. Covers DEBUG_ON / DEBUG_STAGES (all stages and
%   kinds), DEBUG_REGISTER / DEBUG_SAVE_FIGURE, GENERATE_DEBUG_REPORT,
%   SEG_QC through the instrumented head steps, GENERATE_SEGMENTATION_REPORT
%   and REPORT_LINK.
%
%   Run with:   results = runtests('tests/test_debug_reports.m');

    properties
        tmp
    end

    methods (TestMethodSetup)
        function make_tmp(tc)
            tc.tmp = tempname; mkdir(tc.tmp);
            tc.addTeardown(@() rmdir(tc.tmp, 's'));
            tc.addTeardown(@() seg_qc('end'));
        end
    end

    methods (Access = private)
        function p = run_params(tc)
            p = struct();
            p.subject_id = 1;
            p.simulation.medium = 'layered';
            p.simulation.debug = 1;
            p.io.dir_output = tc.tmp;
            p.io.output_affix = '';
            p = resolve_io_dirs(p);
            p.modules.generate_report = 0;
        end

        function [p, seg, c] = head(tc)
            % Concentric synthetic head: brain, CSF, skull, skin, an eye
            p = tc.run_params();
            p.io.dir_cache = fullfile(tc.tmp, 'cache'); mkdir(p.io.dir_cache);
            p.grid.resolution_mm = 1;
            p.grid.mode = 'focal_axis';
            p.headmodel.smooth_method = 'gaussian';
            p.headmodel.smooth_fwhm_mm = 1;
            p.headmodel.smooth_threshold_skull = 0.5;
            p.headmodel.smooth_threshold_other = 0.5;
            p.headmodel.skull_fill_method = 'imclose';
            p.headmodel.csf_expansion = 3;
            p.pct.enabled = 0;
            p.layers.water = [0, 3, 6, 9, 10];
            p.layers.brain = [1, 2];
            p.layers.skin  = 5;
            p.layers.skull = 4;
            p.medium_properties = struct('water', struct(), 'brain', struct(), ...
                'skin', struct(), 'skull', struct());
            p.simulation.debug = 0;
            dims = [33, 33, 33];
            [X, Y, Z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
            c = round((dims + 1) / 2);
            R = sqrt((X-c(1)).^2 + (Y-c(2)).^2 + (Z-c(3)).^2);
            seg = zeros(dims, 'uint8');
            seg(R < 7) = 1;
            seg(R >= 7 & R < 8) = 3;
            seg(R >= 8 & R <= 10) = 4;
            seg(R > 10 & R <= 12) = 5;
            seg(R >= 8 & R <= 10 & abs(X - c(1)) <= 1 & abs(Y - c(2)) <= 1 & Z > c(3)) = 3;  % skull gap
            seg(abs(X - c(1) - 9) <= 1 & abs(Y - c(2)) <= 1 & abs(Z - c(3)) <= 1) = 6;      % eye
        end

        function steps = record_head(tc, p, seg, c)
            seg_qc('begin', p);
            mm = preproc_medium_mask(seg, p);
            seg_qc('add_mapping', p, seg, mm, c, 'test grid');
            mm = skull_fill_holes(p, mm, c, seg);
            preproc_crop_eCSF(p, mm, seg, c, c);
            steps = seg_qc('steps');
            tc.assertNotEmpty(steps);
        end
    end

    % ------------------------------------------------------------------ %
    %% Debug settings
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'debug'})

        function test_master_switch_gates_every_stage_and_kind(tc)
            [~, known] = debug_stages(struct());
            kinds = {'', 'figures', 'volumes', 'other'};
            p = tc.run_params();
            for m = [0 1]
                p.simulation.debug = m;
                for s = known
                    for k = kinds
                        tc.verifyEqual(debug_on(p, s{1}, k{1}), m == 1, ...
                            sprintf('debug=%d stage=%s kind=%s', m, s{1}, k{1}));
                    end
                end
            end
        end

        function test_stage_selection_and_kinds(tc)
            [~, known] = debug_stages(struct());
            p = tc.run_params();
            for i = 1:numel(known)
                p.simulation.debug_stages = known(i);
                for j = 1:numel(known)
                    tc.verifyEqual(debug_on(p, known{j}), i == j, sprintf('%s selected, %s asked', known{i}, known{j}));
                end
            end
            p.simulation.debug_stages = 'all';
            p.simulation.debug_figures = 0;
            tc.verifyFalse(debug_on(p, 'preproc', 'figures'));
            tc.verifyTrue(debug_on(p, 'preproc', 'volumes'));
            p.simulation.debug_volumes = 0;
            tc.verifyFalse(debug_on(p, 'medium', 'volumes'));
            p.simulation.debug_stages = {'all'};
            tc.verifyEqual(debug_stages(p), known, 'the Nexus writes all as a one-item list');
        end

        function test_unknown_stage_or_kind_errors(tc)
            p = tc.run_params();
            p.simulation.debug_stages = {'preproc', 'medum'};
            tc.verifyError(@() debug_stages(p), 'PRESTUS:debug:unknownStage');
            tc.verifyError(@() debug_on(p, 'preproc'), 'PRESTUS:debug:unknownStage');
            p.simulation.debug_stages = 'all';
            tc.verifyError(@() debug_on(p, 'preproc', 'pngs'), 'PRESTUS:debug:unknownKind');
            tc.verifyError(@() debug_on(p, 'preproc', 'fig_files'), 'PRESTUS:debug:unknownKind', ...
                '.fig copies are no longer a debug output');
        end

        function test_save_figure_registers_one_png(tc)
            p = tc.run_params();
            h = figure('Visible', 'off'); plot(1:3);
            f = debug_save_figure(p, 'preproc', h, 'sub-001_test', 'A caption');
            tc.verifyEqual(numel(f), 1);
            tc.verifyTrue(isfile(f{1}));
            tc.verifyFalse(isgraphics(h), 'the figure is closed');
            tc.verifyEmpty(dir(fullfile(p.io.dir_debug_preproc, '*.fig')), 'no .fig files');
            r = jsondecode(strtrim(fileread(fullfile(p.io.dir_debug, 'manifest.jsonl'))));
            tc.verifyEqual(r.file, 'preproc/sub-001_test.png', 'stored relative to debug/');
            tc.verifyEqual(r.kind, 'figure');
            tc.verifyEqual(r.group, '', 'a single view has no group');
            tc.verifySubstring(r.source, 'test_save_figure_registers_one_png', ...
                'the caller of debug_save_figure, not debug_save_figure itself');
        end

        function test_save_figure_views_of_a_3d_figure(tc)
            % every view is one PNG of one group, every axes and light follows
            % the camera, and the figure is closed once at the end
            p = tc.run_params();
            views = {'beam', 'side', 'oblique'};
            h = figure('Visible', 'off');
            ax1 = subplot(1, 2, 1); [X, Y, Z] = sphere(8); surf(ax1, X, Y, Z); camlight(ax1);
            ax2 = subplot(1, 2, 2); scatter3(ax2, X(:), Y(:), Z(:));
            f = debug_save_figure(p, 'source', h, 'sub-001_mask', 'Mask', 'Views', views);
            tc.verifyEqual(numel(f), numel(views));
            for i = 1:numel(views)
                tc.verifyTrue(isfile(fullfile(p.io.dir_debug_source, sprintf('sub-001_mask_view-%s.png', views{i}))), views{i});
            end
            tc.verifyFalse(isgraphics(h), 'the figure is closed');
            lines = splitlines(strtrim(fileread(fullfile(p.io.dir_debug, 'manifest.jsonl'))));
            recs = cellfun(@jsondecode, lines, 'UniformOutput', false);
            recs = [recs{:}];
            tc.verifyEqual(numel(recs), numel(views));
            tc.verifyTrue(all(strcmp({recs.group}, 'sub-001_mask')), 'one group');
            tc.verifyEqual(numel(unique({recs.caption})), numel(views), 'each view has its own caption');

            % last view applied to both axes
            h = figure('Visible', 'off');
            ax1 = subplot(1, 2, 1); surf(ax1, X, Y, Z);
            ax2 = subplot(1, 2, 2); surf(ax2, X, Y, Z);
            debug_save_figure(p, 'source', h, 'sub-001_two', 'Two', 'Views', {'side'}, 'Close', false);
            [az1, el1] = view(ax1); [az2, el2] = view(ax2);
            tc.verifyEqual([az1 el1; az2 el2], [0 0; 0 0]);
            close(h);
        end

        function test_save_figure_unknown_view_errors(tc)
            p = tc.run_params();
            h = figure('Visible', 'off');
            tc.addTeardown(@() close(h));
            tc.verifyError(@() debug_save_figure(p, 'source', h, 'x', 'X', 'Views', {'top'}), ...
                'PRESTUS:debug:unknownView');
        end

        function test_register_resolves_niftiwrite_names_and_needs_debug_dir(tc)
            p = tc.run_params();
            mkdir(p.io.dir_debug_medium);
            f = fullfile(p.io.dir_debug_medium, 'matrix_x.nii');
            niftiwrite(zeros(3, 3, 3, 'single'), f, 'Compressed', true);   % writes matrix_x.nii.gz
            debug_register(p, 'medium', f, 'X');
            r = jsondecode(strtrim(fileread(fullfile(p.io.dir_debug, 'manifest.jsonl'))));
            tc.verifyEqual(r.file, 'medium/matrix_x.nii.gz');
            tc.verifyEqual(r.kind, 'volume');
            q = rmfield(p, 'io');
            tc.verifyWarningFree(@() debug_register(q, 'medium', f, 'X'), 'no io.dir_debug: nothing to do');
        end
    end

    % ------------------------------------------------------------------ %
    %% Debug report
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'debug', 'report'})

        function test_debug_report_lists_registered_and_unregistered_files(tc)
            p = tc.run_params();
            h = figure('Visible', 'off'); plot(1:3);
            debug_save_figure(p, 'preproc', h, 'sub-001_fig', 'Preproc caption');
            mkdir(p.io.dir_debug_medium);
            f = fullfile(p.io.dir_debug_medium, 'matrix_density');
            niftiwrite(zeros(3, 3, 3, 'single'), f, 'Compressed', true);
            debug_register(p, 'medium', f, 'Density caption', 'volume');
            mkdir(p.io.dir_debug_source);
            fid = fopen(fullfile(p.io.dir_debug_source, 'stray.txt'), 'w'); fprintf(fid, 'x'); fclose(fid);

            rp = generate_debug_report(p);
            tc.assertTrue(isfile(rp));
            html = fileread(rp);
            tc.verifySubstring(html, 'Preproc caption');
            tc.verifySubstring(html, 'src="debug/preproc/sub-001_fig.png"');
            tc.verifySubstring(html, 'matrix_density.nii.gz');
            tc.verifySubstring(html, 'Density caption');
            tc.verifySubstring(html, 'Unregistered files (1)');
            tc.verifySubstring(html, 'source/stray.txt');
            tc.verifyEmpty(strfind(html, 'Run error'));
        end

        function test_debug_report_shows_views_side_by_side(tc)
            p = tc.run_params();
            h = figure('Visible', 'off'); surf(peaks(8));
            debug_save_figure(p, 'source', h, 'sub-001_mask', 'Mask', 'Views', {'beam', 'side', 'oblique'});
            h = figure('Visible', 'off'); plot(1:3);
            debug_save_figure(p, 'source', h, 'sub-001_flat', 'Flat');
            html = fileread(generate_debug_report(p));
            row = regexp(html, '<div class="debug-views">(.*?)</div>', 'tokens', 'once');
            tc.assertNotEmpty(row, 'the views share one row');
            tc.verifyEqual(numel(strfind(row{1}, '<img')), 3);
            tc.verifyEmpty(strfind(row{1}, 'sub-001_flat'), 'a single-view figure is not in the row');
            tc.verifySubstring(html, 'sub-001_flat.png');
        end

        function test_debug_report_after_error_and_switches(tc)
            p = tc.run_params();
            ME = MException('PRESTUS:test:boom', 'Something failed');
            html = fileread(generate_debug_report(p, ME));
            tc.verifySubstring(html, 'Run error');
            tc.verifySubstring(html, 'Something failed');
            p.report.debug = 0;
            tc.verifyEqual(generate_debug_report(p), '');
            p.report.debug = 1; p.simulation.debug = 0;
            tc.verifyEqual(generate_debug_report(p), '', 'no debug report without debug mode');
        end

        function test_report_link(tc)
            tc.verifyEqual(report_link('/a/b/c/x.png', '/a/b'), 'c/x.png');
            tc.verifyEqual(report_link('/a/x.html', '/a/b'), '../x.html');
            tc.verifyEqual(report_link('/z/x.html', '/a/b'), 'file:///z/x.html');
        end
    end

    % ------------------------------------------------------------------ %
    %% Head-model QC records and segmentation report
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'segqc'})

        function test_steps_recorded_in_order(tc)
            [p, seg, c] = tc.head();
            % eyes mapped into skin, so the eyes-to-water step has voxels to change
            % (with the default layers the eye is water already)
            p.layers.water = [0, 3, 9, 10];
            p.layers.skin  = [5, 6];
            steps = tc.record_head(p, seg, c);
            tc.verifyEqual({steps.id}, {'segmentation', 'smoothing', 'skull_fill', 'eyes', 'csf_crop'});
            tc.verifyEqual(steps(1).label_space, 'charm');
            tc.verifyTrue(all(strcmp({steps(2:end).label_space}, 'medium')));
            tc.verifyGreaterThan(steps(3).n_changed, 0, 'the skull gap is filled');
            tc.verifyGreaterThan(steps(4).n_changed, 0, 'the eye becomes water');
        end

        function test_every_step_conserves_voxels(tc)
            % for every step with a before state: counts sum to the volume,
            % transitions sum to the changed voxels, and net change per
            % tissue equals what moved in minus what moved out
            [p, seg, c] = tc.head();
            steps = tc.record_head(p, seg, c);
            for s = steps(~arrayfun(@(x) isempty(x.counts_before), steps))
                n = prod(s.size);
                tc.verifyEqual(sum(s.counts_before), n, s.id);
                tc.verifyEqual(sum(s.counts_after), n, s.id);
                tc.verifyEqual(sum(s.transitions(:)), s.n_changed, s.id);
                tc.verifyEqual(trace(s.transitions), 0, [s.id ': only changes are counted']);
                net = s.counts_after - s.counts_before;
                tc.verifyEqual(net, sum(s.transitions, 1) - sum(s.transitions, 2)', s.id);
            end
        end

        function test_no_records_without_session_or_when_off(tc)
            [p, seg, c] = tc.head();
            seg_qc('end');
            mm = preproc_medium_mask(seg, p);
            seg_qc('add_mapping', p, seg, mm, c, 'test grid');
            skull_fill_holes(p, mm, c, seg);
            tc.verifyEmpty(seg_qc('steps'), 'no session: nothing recorded');
            p.report.segmentation = 0;
            seg_qc('begin', p);
            seg_qc('add_mapping', p, seg, mm, c, 'test grid');
            tc.verifyEmpty(seg_qc('steps'), 'report.segmentation 0: nothing recorded');
            p.report.segmentation = 1; p.simulation.medium = 'water';
            seg_qc('begin', p);
            seg_qc('add_mapping', p, seg, mm, c, 'test grid');
            tc.verifyEmpty(seg_qc('steps'), 'non-layered medium: nothing recorded');
        end

        function test_segmentation_report_round_trips_slices(tc)
            [p, seg, c] = tc.head();
            steps = tc.record_head(p, seg, c);
            seg_qc('save', p);
            seg_qc('end');
            tc.assertTrue(isfile(seg_qc('file', p)));
            rp = generate_segmentation_report(p);
            tc.assertTrue(isfile(rp));
            html = fileread(rp);
            for s = steps
                tc.verifySubstring(html, s.title);
            end
            tok = regexp(html, '<script type="application/json" id="segqc-data">(.*?)</script>', 'tokens', 'once');
            tc.assertNotEmpty(tok);
            D = jsondecode(strrep(tok{1}, '<\/', '</'));
            for i = 1:numel(steps)
                S = pick(D.steps, i);
                for pl = 1:3
                    P = pick(S.planes, pl);
                    tc.verifyEqual(unrle(P.after, P.w * P.h), ...
                        reshape(fliplr(steps(i).slices_after{pl}), 1, []), ...
                        sprintf('step %d plane %d', i, pl));
                end
            end
        end

        function test_segmentation_report_without_records_and_switch(tc)
            [p, ~, ~] = tc.head();
            html = fileread(generate_segmentation_report(p));
            tc.verifySubstring(html, 'No postprocessing records');
            p.report.segmentation = 0;
            tc.verifyEqual(generate_segmentation_report(p), '');
            p.report.segmentation = 1; p.simulation.medium = 'water';
            tc.verifyEqual(generate_segmentation_report(p), '', 'water runs have no segmentation report');
        end
    end
end

function x = pick(a, i)
% jsondecode gives a struct array for uniform objects, a cell otherwise
    if iscell(a); x = a{i}; else; x = a(i); end
end

function out = unrle(e, n)
% MATLAB twin of the report's JavaScript decoder
    v = matlab.net.base64decode(e.v);
    l = double(typecast(matlab.net.base64decode(e.n), 'uint16'));
    out = repelem(v, l);
    assert(numel(out) == n);
end
