classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_get_t1_mni_file < matlab.unittest.TestCase
% TEST_GET_T1_MNI_FILE  The T1-in-MNI image lives next to the segmentation, not in m2m.
%
% GET_T1_MNI_FILE resolves <parent of path.seg>/mni/sub-NNN/; SIMNIBS_T1_TO_MNI
% checks and writes that file and never writes into the m2m folder.
%
%   results = runtests('tests/test_get_t1_mni_file.m');

    properties
        Root (1,:) char
    end

    methods (TestMethodSetup)
        function make_root(tc)
            tc.Root = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
        end
    end

    methods (Test, TestTags = {'paths'})

        function test_path_is_outside_seg(tc)
            seg = fullfile(tc.Root, 'seg');
            for s = {seg, [seg filesep]}                    % with and without trailing separator
                [f, f_m2m] = get_t1_mni_file(s{1}, 7);
                tc.verifyEqual(f, fullfile(tc.Root, 'mni', 'sub-007', 'T1_to_MNI_post-hoc.nii.gz'));
                tc.verifyFalse(startsWith(f, seg), 'T1-in-MNI path is inside path.seg.');
                tc.verifyEqual(f_m2m, fullfile(seg, 'm2m_sub-007', 'toMNI', 'T1_to_MNI_post-hoc.nii.gz'));
            end
        end

        function test_existing_file_skips_without_m2m_writes(tc)
            % before: the step checked m2m/toMNI, found nothing, and registered
            % into the m2m folder
            [p, m2m] = tc.subject();
            f = get_t1_mni_file(p.path.seg, p.subject_id);
            mkdir(fileparts(f)); fclose(fopen(f, 'w'));
            p.io.overwrite_files = 'never';
            tc.verifyWarningFree(@() simnibs_t1_to_mni(p, true, m2m));
            tc.verifyFalse(isfolder(fullfile(m2m, 'toMNI')), 'Wrote into the m2m folder.');
        end

        function test_registration_targets_derived_folder(tc)
            % no warp field: the native warp fails after the output folder is made
            [p, m2m] = tc.subject();
            p.io.overwrite_files = 'always';
            p.analysis.mni_warp_method = 'native';
            tc.verifyError(@() simnibs_t1_to_mni(p, true, m2m), ?MException);
            tc.verifyTrue(isfolder(fileparts(get_t1_mni_file(p.path.seg, p.subject_id))));
            d = dir(m2m);
            tc.verifyEqual(setdiff({d.name}, {'.', '..'}), {'T1.nii.gz'}, 'm2m folder changed.');
        end

        function test_non_layered_is_skipped(tc)
            % negative control: no registration, no folder
            [p, m2m] = tc.subject();
            simnibs_t1_to_mni(p, false, m2m);
            tc.verifyFalse(isfolder(fullfile(tc.Root, 'mni')));
        end
    end

    methods (Access = private)
        function [p, m2m] = subject(tc)
            p = struct('subject_id', 3);
            p.path.seg = fullfile(tc.Root, 'seg');
            p.io = struct('save_MNI', true, 'overwrite_files', 'never');
            p.simulation.interactive = 0;
            m2m = fullfile(p.path.seg, 'm2m_sub-003');
            mkdir(m2m);
            fclose(fopen(fullfile(m2m, 'T1.nii.gz'), 'w'));
        end
    end
end
