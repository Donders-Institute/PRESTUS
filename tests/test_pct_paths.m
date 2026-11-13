classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_pct_paths < matlab.unittest.TestCase
% TEST_PCT_PATHS  Pseudo-CT location: GET_PCT_DIR and FIND_PSEUDOCT_DIR.
%
% The pseudo-CT goes to <path.pct>/sub-NNN, or <parent of path.seg>/pct/sub-NNN
% when path.pct is empty, never into the SimNIBS m2m folder. A pseudo-CT that
% older versions wrote into m2m_sub-NNN is still found, with a warning.
%
% Run with:   results = runtests('tests/test_pct_paths.m');

    properties
        TmpDir
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
            tc.addTeardown(@() rmdir(tc.TmpDir, 's'));
        end
    end

    methods (Test)
        function test_default_is_next_to_seg(tc)
            p = base_params(tc.TmpDir);
            tc.verifyEqual(get_pct_dir(p), fullfile(tc.TmpDir, 'pct', 'sub-009'));
            p.path.seg = [p.path.seg filesep];        % trailing separator
            tc.verifyEqual(get_pct_dir(p), fullfile(tc.TmpDir, 'pct', 'sub-009'));
        end

        function test_path_pct_wins(tc)
            p = base_params(tc.TmpDir);
            p.path.pct = fullfile(tc.TmpDir, 'pct_kosciessa');
            tc.verifyEqual(get_pct_dir(p), fullfile(tc.TmpDir, 'pct_kosciessa', 'sub-009'));
        end

        function test_missing_paths_error(tc)
            p = struct('subject_id', 9, 'path', struct());
            tc.verifyError(@() get_pct_dir(p), 'get_pct_dir:missingPath');
        end

        function test_find_new_location(tc)
            p = base_params(tc.TmpDir);
            touch(fullfile(get_pct_dir(p), 'pseudoCT.nii.gz'));
            tc.verifyEqual(tc.verifyWarningFree(@() find_pseudoct_dir(p)), get_pct_dir(p));
        end

        function test_find_legacy_m2m_location_warns(tc)
            p = base_params(tc.TmpDir);
            legacy = fullfile(p.path.seg, 'm2m_sub-009');
            touch(fullfile(legacy, 'pseudoCT.nii.gz'));
            tc.verifyWarning(@() find_pseudoct_dir(p), 'PRESTUS:pct:legacyLocation');
            w = warning('off', 'PRESTUS:pct:legacyLocation');
            tc.addTeardown(@() warning(w));
            tc.verifyEqual(find_pseudoct_dir(p), legacy);
            % The new location wins once the file exists there
            touch(fullfile(get_pct_dir(p), 'pseudoCT.nii.gz'));
            tc.verifyEqual(find_pseudoct_dir(p), get_pct_dir(p));
        end

        function test_find_none(tc)
            p = base_params(tc.TmpDir);
            tc.verifyEqual(find_pseudoct_dir(p), '');
        end
    end
end

function p = base_params(root)
    p = struct();
    p.subject_id = 9;
    p.path.seg = fullfile(root, 'simnibs');
    p.path.pct = '';
end

function touch(f)
    if ~isfolder(fileparts(f)); mkdir(fileparts(f)); end
    fclose(fopen(f, 'w'));
end
