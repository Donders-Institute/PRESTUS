classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_group < matlab.unittest.TestCase
% TEST_GROUP  Group-level helpers (functions/group) on synthetic inputs.
%
% Covers CHANGEM_VECTORIZED (label relabelling), DISCOVER_GROUP_SUBJECTS
% (finding finished subjects by their per-subject CSV in a temporary
% path.sim) and STUDY_RUN_DESCRIPTORS (per-run settings for the study plan
% and group report). The study table / manifest functions are covered by
% test_study. CREATE_GROUP_MNI_PLOTS and COMBINE_PLOTS_BY_SUFFIX (figures,
% image montage from real subject outputs) are not covered.
%
%   results = runtests('tests/test_group.m');

    properties
        Dir
    end


    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tempname;
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
        end
    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'discover_group_subjects'})

        function test_discover_finds_subjects_with_csv(tc)
            p = group_params(tc.Dir);
            touch_csv(tc.Dir, 3, 'layered', '');
            touch_csv(tc.Dir, 1, 'layered', '');
            mkdir(fullfile(tc.Dir, 'sub-002'));                  % no CSV: not finished
            touch_csv(tc.Dir, 4, 'water', '');                   % other medium
            tc.verifyEqual(discover_group_subjects(p), [1 3]);
        end

        function test_discover_respects_output_affix(tc)
            p = group_params(tc.Dir);
            p.io.output_affix = '_desc-A';
            touch_csv(tc.Dir, 1, 'layered', '');
            touch_csv(tc.Dir, 2, 'layered', '_desc-A');
            tc.verifyEqual(discover_group_subjects(p), 2);
        end

        function test_discover_ignores_non_subject_folders(tc)
            p = group_params(tc.Dir);
            touch_csv(tc.Dir, 5, 'layered', '');
            mkdir(fullfile(tc.Dir, 'sub-abc'));
            mkdir(fullfile(tc.Dir, 'sub-007_old'));
            fclose(fopen(fullfile(tc.Dir, 'sub-008'), 'w'));    % a file, not a folder
            tc.verifyEqual(discover_group_subjects(p), 5);
        end

        function test_discover_empty_and_missing_paths(tc)
            p = group_params(tc.Dir);
            tc.verifyEmpty(discover_group_subjects(p));
            p.path.sim = fullfile(tc.Dir, 'does_not_exist');
            tc.verifyWarning(@() discover_group_subjects(p), 'discover_group_subjects:simPathMissing');
            p.path.sim = '';
            tc.verifyWarning(@() discover_group_subjects(p), 'discover_group_subjects:noSimPath');
        end

    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'study_run_descriptors'})

        function test_descriptors_localite(tc)
            p = struct();
            p.placement = struct('mode', 'localite', 'localite', struct('target', 'M1'));
            p.transducer = struct('name', {'CTX-500', ''}, 'focal_distance_ep', {45, []}, ...
                'target_isppa_wcm2', {[10 20], 5});
            p.study = struct('condition', "active");
            d = study_run_descriptors(p);
            tc.verifyEqual(d.placement_mode, "localite");
            tc.verifyEqual(d.target_label, "M1");
            tc.verifyEqual(d.transducer_name, "CTX-500");
            tc.verifyEqual(d.n_transducers, 2);
            tc.verifyEqual(d.focal_distance_ep_mm, [45 NaN]);
            tc.verifyEqual(d.target_isppa_wcm2, [10 20]);
            tc.verifyEqual(d.condition, "active");
        end

        function test_descriptors_localite_series_fallback(tc)
            p.placement = struct('mode', 'localite', 'localite', struct('target', '', 'position', [3 4]));
            d = study_run_descriptors(p);
            tc.verifyEqual(d.target_label, "series 3");
        end

        function test_descriptors_mni_and_serial_fallback(tc)
            p.placement = struct('mode', 'mni', 'mni', struct('focus_pos_mm', [-38 -20 60]));
            p.transducer = struct('name', '', 'serial', 'SN123');
            d = study_run_descriptors(p);
            tc.verifyEqual(d.target_label, "MNI [-38 -20 60]");
            tc.verifyEqual(d.transducer_name, "SN123");
        end

        function test_descriptors_empty_parameters(tc)
            d = study_run_descriptors(struct());
            tc.verifyEqual(d.placement_mode, "");
            tc.verifyEqual(d.target_label, "");
            tc.verifyEqual(d.n_transducers, 0);
            tc.verifyEqual(d.condition, "");
        end

    end

end

% =========================================================================

function p = group_params(root)
    p = struct('path', struct('sim', root), ...
               'simulation', struct('medium', 'layered'), ...
               'io', struct('output_affix', ''));
end

function touch_csv(root, id, medium, affix)
    d = fullfile(root, sprintf('sub-%03d', id));
    if ~isfolder(d); mkdir(d); end
    fclose(fopen(fullfile(d, sprintf('sub-%03d_%s%s.csv', id, medium, affix)), 'w'));
end
