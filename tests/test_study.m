classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_study < matlab.unittest.TestCase
% TEST_STUDY  Unit tests for the study (multi-subject) setup and group analysis.
%
%   Covers the study manifest reader, key/value config overrides (also via
%   load_parameters), the run table (affixes, duplicate guard), a dry run of
%   prestus_study_start, the per-run targeting analysis on synthetic NIfTIs,
%   the group MNI maps, and the study-mode group report. No MRI / SimNIBS /
%   GPU needed.
%
%   Run with:  results = runtests('tests/test_study.m');

    properties
        tmp   % per-test temporary folder
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.tmp = fullfile(tempdir, 'prestus_study_unittest');
            if exist(tc.tmp, 'dir'); rmdir(tc.tmp, 's'); end
            mkdir(tc.tmp);
            ws = warning('off', 'prestus:toolboxDefault');
            ws(2) = warning('off', 'prestus:noSimNIBS');
            tc.addTeardown(@() warning(ws));
        end
    end

    methods (TestMethodTeardown)
        function teardown(tc)
            close all force
            if exist(tc.tmp, 'dir'); rmdir(tc.tmp, 's'); end
        end
    end

    methods (Test, TestTags = {'study'})

        % ── manifest ────────────────────────────────────────────────────
        function test_manifest_parses_reserved_and_key_columns(tc)
            f = write_text(fullfile(tc.tmp, 'm.tsv'), { ...
                '# a comment line', ...
                sprintf('subject_id\tcondition\tinclude\tplacement.localite.target\ttransducer.target_isppa_wcm2\tanalysis.target_mni_mm\tnotes'), ...
                sprintf('1\tM1\t1\tL_M1\t30\t[-38, -20, 60]\tfirst'), ...
                sprintf('sub-003\tM1\tyes\tL_M1\t\t\t'), ...
                sprintf('3\tPUL\t0\tL_PUL\t20\t\t')});
            [m, info] = study_read_manifest(f);
            tc.verifyEqual(height(m), 3);
            tc.verifyEqual(m.subject_id', [1 3 3]);
            tc.verifyEqual(m.condition', ["M1" "M1" "PUL"]);
            tc.verifyEqual(m.include', [true true false]);
            tc.verifyEqual(info.override_keys, {'placement.localite.target', ...
                'transducer.target_isppa_wcm2', 'analysis.target_mni_mm'});
            ov = m.overrides{1};
            tc.verifyEqual(ov(:, 1)', info.override_keys);
            tc.verifyEqual(ov{1, 2}, 'L_M1');
            tc.verifyEqual(ov{2, 2}, 30);
            tc.verifyEqual(ov{3, 2}, [-38 -20 60]);
            tc.verifyEqual(size(m.overrides{2}, 1), 1, 'empty cells must keep the study value');
        end

        function test_manifest_csv_quotes_and_config_path(tc)
            write_text(fullfile(tc.tmp, 'row.yaml'), {'grid:', '  resolution_mm: 1'});
            f = write_text(fullfile(tc.tmp, 'm.csv'), { ...
                'subject_id,config,placement.mni.focus_pos_mm', ...
                '2,row.yaml,"[10, 20, 30]"'});
            m = study_read_manifest(f);
            tc.verifyEqual(m.overrides{1}{1, 2}, [10 20 30]);
            tc.verifyEqual(char(m.config(1)), fullfile(tc.tmp, 'row.yaml'));
        end

        function test_manifest_rejects_bad_input(tc)
            f = write_text(fullfile(tc.tmp, 'a.tsv'), {sprintf('subject\tcondition'), sprintf('1\tx')});
            tc.verifyError(@() study_read_manifest(f), 'study_read_manifest:noSubjectColumn');
            f = write_text(fullfile(tc.tmp, 'b.tsv'), {sprintf('subject_id\tnot a key'), sprintf('1\tx')});
            tc.verifyError(@() study_read_manifest(f), 'study_read_manifest:badColumn');
            f = write_text(fullfile(tc.tmp, 'c.tsv'), {sprintf('subject_id'), sprintf('abc')});
            tc.verifyError(@() study_read_manifest(f), 'study_read_manifest:badSubject');
        end

        % ── overrides ───────────────────────────────────────────────────
        function test_overrides_nested_and_indexed(tc)
            p = struct('a', 1, 'b', struct('c', 2));
            p.transducer = struct('serial', {'X', 'Y'}, 'focus_pos', {[], []});
            p = apply_config_overrides(p, {'b.c', 5; 'transducer(2).serial', 'Z'});
            tc.verifyEqual(p.b.c, 5);
            tc.verifyEqual(p.transducer(2).serial, 'Z');
            tc.verifyEqual(p.transducer(1).serial, 'X');
            tc.verifyError(@() apply_config_overrides(p, {'transducer.serial', 'Q'}), ...
                'apply_config_overrides:ambiguous');
            tc.verifyError(@() apply_config_overrides(p, {'transducer(3).serial', 'Q'}), ...
                'apply_config_overrides:indexRange');
            tc.verifyWarning(@() apply_config_overrides(p, {'b.typo', 1}), ...
                'apply_config_overrides:unknownKey');
        end

        function test_load_parameters_applies_overrides_before_derivation(tc)
            p = load_parameters({'subject_id', 7; 'placement.localite.target', 'L_MD'; ...
                                 'io.output_affix', 'cond A'});
            tc.verifyEqual(p.subject_id, 7);
            tc.verifyEqual(p.placement.localite.target, 'L_MD');
            tc.verifyEqual(p.io.output_affix, '_desc-cond_A', 'override must be sanitised like the YAML value');
        end

        % ── run table / study setup ─────────────────────────────────────
        function test_run_table_affixes_and_duplicates(tc)
            f = write_text(fullfile(tc.tmp, 'm.tsv'), { ...
                sprintf('subject_id\tcondition\tio.output_affix'), ...
                sprintf('1\tM1\t'), sprintf('1\tPUL\t'), sprintf('2\t\tcustom'), sprintf('3\t\t')});
            m = study_read_manifest(f);
            p = struct('io', struct('output_affix', '_desc-liberal'), ...
                       'simulation', struct('medium', 'layered'), 'path', struct('sim', tc.tmp));
            runs = study_run_table(p, m);
            tc.verifyEqual(runs.affix', ["_desc-liberal-M1" "_desc-liberal-PUL" "_desc-custom" "_desc-liberal"]);
            tc.verifyEqual(runs.label(1), "sub-001 · M1");

            f = write_text(fullfile(tc.tmp, 'dup.tsv'), {sprintf('subject_id\tcondition'), ...
                sprintf('1\tM1'), sprintf('1\tM1')});
            tc.verifyError(@() study_run_table(p, study_read_manifest(f)), 'study_run_table:duplicateRun');
        end

        function test_study_dry_run_resolves_every_row(tc)
            cfg = make_study(tc.tmp, { ...
                sprintf('subject_id\tcondition\ttransducer.target_isppa_wcm2\tplacement.mode\tplacement.heuristic.target_name'), ...
                sprintf('1\tA\t30\theuristic\tM1'), ...
                sprintf('2\tA\t20\theuristic\tM1'), ...
                sprintf('2\tB\t\tmanual\t')});
            runs = prestus_study_start(cfg, 'DryRun', true, 'Platform', 'matlab');
            tc.verifyEqual(height(runs), 3);
            tc.verifyEqual(runs.status', repmat("planned", 1, 3));
            tc.verifyEqual(runs.target_isppa_wcm2', ["30" "20" "15"]);
            tc.verifyEqual(runs.placement_mode', ["heuristic" "heuristic" "manual"]);
            tc.verifyEqual(runs.target(1), "M1");
            plans = dir(fullfile(tc.tmp, 'sim', 'group', 'study', 'log', 'plan_*.tsv'));
            tc.verifyNotEmpty(plans, 'dry run must write a plan');
            tc.verifyFalse(isfolder(fullfile(tc.tmp, 'sim', 'sub-001')), 'dry run must not run anything');

            % Selection options
            runs = prestus_study_start(cfg, 'DryRun', true, 'Platform', 'matlab', 'Conditions', "B");
            tc.verifyEqual(height(runs), 1);
            tc.verifyEqual(runs.subject_id, 2);
        end

        function test_row_parameters_merge_order(tc)
            write_text(fullfile(tc.tmp, 'row.yaml'), {'transducer:', '  target_isppa_wcm2: 10', ...
                'grid:', '  resolution_mm: 1'});
            cfg = make_study(tc.tmp, {sprintf('subject_id\tconfig\ttransducer.target_isppa_wcm2'), ...
                sprintf('4\trow.yaml\t25'), sprintf('5\trow.yaml\t')});
            s = study_setup(cfg);
            p4 = study_row_parameters(s.config, s.runs(1, :));
            p5 = study_row_parameters(s.config, s.runs(2, :));
            tc.verifyEqual(p4.subject_id, 4);
            tc.verifyEqual(p4.transducer(1).target_isppa_wcm2, 25, 'manifest column beats row config');
            tc.verifyEqual(p5.transducer(1).target_isppa_wcm2, 10, 'row config beats study config');
            tc.verifyEqual(p5.grid.resolution_mm, 1);
        end

        % ── targeting ───────────────────────────────────────────────────
        function test_targeting_analysis_synthetic(tc)
            sub_dir = fullfile(tc.tmp, 'sim', 'sub-001');
            nii_dir = fullfile(sub_dir, 'nii'); mkdir(nii_dir);
            % T1w grid: 1 mm, origin at voxel 0 = (0,0,0) mm. Focal blob centred at voxel (20,20,30).
            T = eye(4);
            [mask, I] = blob([40 40 50], [20 20 30], 3);
            write_nii(fullfile(nii_dir, 'sub-001_layered_T1w_desc-A_focal_mask_3dB'), uint8(mask), T);
            write_nii(fullfile(nii_dir, 'sub-001_layered_T1w_desc-A_Ipa'), single(I), T);
            % MNI grid: 1 mm, x flipped (LAS storage) with origin at +20 mm
            Tm = eye(4); Tm(1, 1) = -1; Tm(4, 1:3) = [20 -20 -25];
            write_nii(fullfile(nii_dir, 'sub-001_layered_MNI_desc-A_focal_mask_3dB'), uint8(mask), Tm);
            write_nii(fullfile(nii_dir, 'sub-001_layered_MNI_desc-A_Ipa'), single(I), Tm);
            % ROI on a 2 mm grid covering x in [-2, 2] (MNI mm around the blob centre)
            Tr = eye(4) * 2; Tr(4, :) = [-20 -20 -20 1];
            roi = false(20, 20, 20); roi(9:11, 9:12, 11:14) = true;
            write_nii(fullfile(tc.tmp, 'roi'), uint8(roi), Tr);

            p = struct('subject_id', 1, 'simulation', struct('medium', 'layered'), ...
                'io', struct('output_affix', '_desc-A', 'dir_output', sub_dir), ...
                'placement', struct('mode', 'manual'), ...
                'analysis', struct('target_mni_mm', [0 0 5], 'target_roi_mni', fullfile(tc.tmp, 'roi.nii.gz')), ...
                'study', struct('condition', 'A'));
            % planned focus 3 mm shallower along +z (beam from z=0 towards +z)
            t = targeting_analysis(p, 'PlannedFocus', [20 20 27], 'PlannedTrans', [20 20 1]);
            tc.verifyTrue(isfile(fullfile(sub_dir, 'sub-001_layered_desc-A_targeting.csv')));
            tc.verifyEqual([t.focal_centroid_ras_x t.focal_centroid_ras_y t.focal_centroid_ras_z], ...
                [19 19 29], 'AbsTol', 1e-9);
            tc.verifyEqual(t.focus_error_mm, 3, 'AbsTol', 1e-9);
            tc.verifyEqual(t.focus_error_axial_mm, 3, 'AbsTol', 1e-9, '+ = deeper than planned');
            tc.verifyEqual(t.focus_error_lateral_mm, 0, 'AbsTol', 1e-9);
            % MNI: voxel 19 (0-based) on x -> 20 - 19 = 1 mm; y -> -1; z -> 4
            tc.verifyEqual([t.focal_centroid_mni_x t.focal_centroid_mni_y t.focal_centroid_mni_z], ...
                [1 -1 4], 'AbsTol', 1e-9);
            tc.verifyEqual(t.target_mni_source, "analysis.target_mni_mm");
            tc.verifyEqual(t.target_error_mni_mm, sqrt(1 + 1 + 1), 'AbsTol', 1e-9);
            tc.verifyGreaterThan(t.roi_focal_fraction, 0);
            tc.verifyLessThanOrEqual(t.roi_focal_fraction, 1);
            tc.verifyGreaterThan(t.roi_mean_intensity, 0);
        end

        % ── group maps + study report ───────────────────────────────────
        function test_group_maps_and_study_report(tc)
            cfg = make_study(tc.tmp, [{sprintf('subject_id\tcondition\ttransducer.target_isppa_wcm2')}, ...
                arrayfun(@(k) sprintf('%d\tA\t30', k), 1:6, 'UniformOutput', false), ...
                arrayfun(@(k) sprintf('%d\tB\t20', k), 1:3, 'UniformOutput', false), ...
                {sprintf('9\tB\t20')}]);   % sub-009 has no outputs
            s = study_setup(cfg);
            Tm = eye(4); Tm(1, 1) = -1; Tm(4, 1:3) = [20 -20 -25];
            for r = 1:height(s.runs) - 1
                f = study_run_files(s.runs.subject_id(r), 'layered', s.runs.affix(r), s.runs.sim_dir(r));
                if ~isfolder(f.nii_dir); mkdir(f.nii_dir); end
                c = [20 20 30] + [r 0 0] * (s.runs.condition(r) == "B");
                [mask, I] = blob([40 40 50], c, 3);
                write_nii(regexprep(f.mni('focal_mask_3dB'), '\.nii\.gz$', ''), uint8(mask), Tm);
                write_nii(regexprep(f.mni('Ipa'), '\.nii\.gz$', ''), single(I), Tm);
                isppa = 20 + r;
                if r == 6; isppa = 80; end              % one outlier in condition A
                row = synthetic_csv_row(s.runs.subject_id(r), isppa);
                writetable(row, f.csv);
                tg = table(s.runs.subject_id(r), string(s.runs.condition(r)), "manual", "L_M1", 1.5 + r / 10, ...
                    0.5, 1.2, 1 - c(1) + 19, -1, 4, 0, 0, 5, ...
                    'VariableNames', {'subject_id', 'condition', 'placement_mode', 'target_label', ...
                    'focus_error_mm', 'focus_error_axial_mm', 'focus_error_lateral_mm', ...
                    'focal_centroid_mni_x', 'focal_centroid_mni_y', 'focal_centroid_mni_z', ...
                    'target_mni_x', 'target_mni_y', 'target_mni_z'});
                tg.target_error_mni_mm = 3 + r / 10;
                writetable(tg, f.targeting_csv);
            end

            maps = create_group_maps(s.runs(s.runs.include, :), s.out_dir);
            tc.verifyEqual(numel(maps), 2);
            tc.verifyEqual([maps.n], [6 3]);
            tc.verifyEqual(maps(2).n_missing, 1);
            ov = niftiread(maps(1).files.overlap);
            tc.verifyEqual(max(ov(:)), single(1), 'identical masks in A overlap fully');
            ovB = niftiread(maps(2).files.overlap);
            tc.verifyLessThan(max(ovB(:)), single(1.0001));
            tc.verifyTrue(isfile(maps(1).figure), 'overlap figure missing');
            tc.verifyTrue(isfield(maps(1).files, 'intensity_mean'));

            rp = prestus_group_report_start(cfg);
            tc.verifyTrue(isfile(rp), 'no study report produced');
            tc.verifyEqual(rp, fullfile(s.out_dir, 'group_study_report.html'));
            txt = fileread(rp);
            tc.verifyEqual(numel(strfind(txt, 'Section failed')), 0, 'a section builder errored');
            for id = {'id="design"', 'id="compare"', 'id="targeting"', 'id="maps"', 'id="outliers"'}
                tc.verifyNotEmpty(strfind(txt, id{1}), sprintf('%s section missing', id{1}));
            end
            tc.verifyNotEmpty(strfind(txt, 'sub-009 · B'), 'missing run should be listed');
            tc.verifyNotEmpty(strfind(txt, 'value="condition" selected'), 'group-by condition missing');
            tc.verifyTrue(isfile(fullfile(s.out_dir, 'group_study_runs.csv')), 'tidy export missing');
            exported = readtable(fullfile(s.out_dir, 'group_study_runs.csv'));
            tc.verifyEqual(height(exported), 9);
            tc.verifyTrue(ismember('focus_error_mm', exported.Properties.VariableNames));
            % exactly one outlier row (Isppa of sub-006 in A); roster flags it
            tc.verifyNotEmpty(regexp(txt, '<td><a href="#row-6">sub-006 · A</a></td><td>A</td><td>Isppa_brain', 'once'));
        end

        function test_legacy_subject_list_still_works(tc)
            sim = fullfile(tc.tmp, 'sim');
            for id = 1:2
                d = fullfile(sim, sprintf('sub-%03d', id)); mkdir(d);
                writetable(synthetic_csv_row(id, 20 + id), fullfile(d, sprintf('sub-%03d_layered.csv', id)));
            end
            p = struct('path', struct('sim', sim), 'simulation', struct('medium', 'layered'), ...
                       'io', struct('output_affix', ''), 'modules', struct('run_heating_sims', 0));
            rp = generate_group_report(p);
            tc.verifyEqual(rp, fullfile(sim, 'group_layered_report.html'));
            txt = fileread(rp);
            tc.verifyEqual(numel(strfind(txt, 'Section failed')), 0);
            tc.verifyEmpty(strfind(txt, 'id="compare"'), 'no comparison with a single condition');
        end
    end

    % study_read_manifest
    methods (Test, TestTags = {'study'})

        function test_manifest_keeps_labels_and_serials_as_text(tc)
            f = fullfile(tc.tmp, 'manifest.tsv');
            fid = fopen(f, 'w');
            fprintf(fid, '# comment line\n');
            fprintf(fid, 'subject_id\tplacement.localite.target\ttransducer.serial\ttransducer.target_isppa_wcm2\n');
            fprintf(fid, '1\tno\t001\t30\n');
            fclose(fid);
            m = study_read_manifest(f);
            ov = m.overrides{1};
            tc.verifyEqual(ov{1, 2}, 'no');
            tc.verifyEqual(ov{2, 2}, '001');
            tc.verifyEqual(ov{3, 2}, 30);
        end

        function test_manifest_error_reports_file_line(tc)
            f = fullfile(tc.tmp, 'manifest.tsv');
            fid = fopen(f, 'w');
            fprintf(fid, '# comment\n\nsubject_id\tcondition\n1\tA\textra\n');
            fclose(fid);
            try
                study_read_manifest(f);
                tc.verifyFail('Expected tooManyCells error');
            catch ME
                tc.verifyEqual(ME.identifier, 'study_read_manifest:tooManyCells');
                tc.verifySubstring(ME.message, 'line 4');
            end
        end

    end
end

% =========================================================================
%  Local helpers
% =========================================================================
function f = write_text(f, lines)
    fid = fopen(f, 'w');
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end

function cfg = make_study(tmp, manifest_lines)
    write_text(fullfile(tmp, 'study.tsv'), manifest_lines);
    cfg = fullfile(tmp, 'config_study.yaml');
    write_text(cfg, { ...
        'path:', sprintf('  sim: ''%s''', fullfile(tmp, 'sim')), ...
        'simulation:', '  medium: layered', ...
        'transducer:', '  target_isppa_wcm2: 15', ...
        'study:', '  manifest: study.tsv'});
end

function [mask, I] = blob(sz, c, r)
    [x, y, z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
    d2 = (x - c(1)).^2 + (y - c(2)).^2 + (z - c(3)).^2;
    mask = d2 <= r^2;
    I = 30 * exp(-d2 / (2 * r^2));
end

function write_nii(stem, vol, T)
    niftiwrite(vol, stem, 'Compressed', true);
    info = niftiinfo([stem '.nii.gz']);
    info.Transform = affine3d(T);
    info.TransformName = 'Sform';
    niftiwrite(vol, stem, info, 'Compressed', true);
end

function row = synthetic_csv_row(id, isppa_brain)
    row = table();
    row.subject_id = id; row.freq_Hz = 5e5;
    row.Isppa = isppa_brain + 1; row.Ipa_target = isppa_brain - 1;
    row.real_focal_distance_mm = 64; row.real_focal_distance_ep_mm = 40;
    row.MI_tc = 0.7; row.MI_brain = 0.7; row.MI_skull = 0.9; row.MI_skin = 0.4; row.TIC = 0.8;
    row.Isppa_brain = isppa_brain; row.Isppa_skull = 9; row.Isppa_skin = 8;
    row.focal_vol_3dB_mm3 = 110;
end
