classdef test_deface_m2m < matlab.unittest.TestCase
% TEST_DEFACE_M2M  Unit tests for preproc_deface_m2m.
%
%   Tests use a temporary fake m2m folder and a stub pydeface binary —
%   no real MRI data or SimNIBS installation required.
%
%   Run with:   results = runtests('tests/test_deface_m2m.m');

    properties
        TmpDir       % root temp directory for this test run
        M2mDir       % fake m2m_sub-001/ folder
        BinDir       % directory containing the stub pydeface binary
        Params       % minimal parameters struct
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);

            % Minimal fake m2m folder with one anatomical NIfTI
            tc.M2mDir = fullfile(tc.TmpDir, 'm2m_sub-001');
            mkdir(tc.M2mDir);
            fclose(fopen(fullfile(tc.M2mDir, 'T1.nii.gz'), 'w'));
            fclose(fopen(fullfile(tc.M2mDir, 'sub-001.msh'), 'w'));

            % Stub pydeface binary that succeeds (exit 0, no-op)
            tc.BinDir = fullfile(tc.TmpDir, 'bin');
            mkdir(tc.BinDir);
            stub = fullfile(tc.BinDir, 'pydeface');
            fid = fopen(stub, 'w');
            fprintf(fid, '#!/bin/sh\nexit 0\n');
            fclose(fid);
            system(sprintf('chmod +x "%s"', stub));

            tc.Params = struct('startup', struct('simnibs_bin_path', tc.BinDir));
        end
    end

    methods (TestMethodTeardown)
        function cleanup(tc)
            if isfolder(tc.TmpDir)
                rmdir(tc.TmpDir, 's');
            end
        end
    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'deface'})

        function test_creates_defaced_folder(tc)
            out = preproc_deface_m2m(tc.M2mDir, tc.Params);
            expected = [tc.M2mDir, '_defaced'];
            tc.verifyEqual(out, expected);
            tc.verifyTrue(isfolder(expected), '_defaced folder should be created');
        end

        function test_writes_sentinel(tc)
            preproc_deface_m2m(tc.M2mDir, tc.Params);
            sentinel = fullfile([tc.M2mDir, '_defaced'], '.defaced');
            tc.verifyTrue(isfile(sentinel), '.defaced sentinel should exist');
        end

        function test_copies_non_anat_files(tc)
            preproc_deface_m2m(tc.M2mDir, tc.Params);
            mesh_copy = fullfile([tc.M2mDir, '_defaced'], 'sub-001.msh');
            tc.verifyTrue(isfile(mesh_copy), 'Non-anatomical files should be copied unchanged');
        end

        function test_idempotent_on_second_call(tc)
            % First call creates the folder and sentinel
            preproc_deface_m2m(tc.M2mDir, tc.Params);
            sentinel = fullfile([tc.M2mDir, '_defaced'], '.defaced');
            t1 = dir(sentinel);

            % Second call should return immediately without re-running
            out2 = preproc_deface_m2m(tc.M2mDir, tc.Params);
            t2   = dir(sentinel);
            tc.verifyEqual(out2, [tc.M2mDir, '_defaced']);
            tc.verifyEqual(t1.datenum, t2.datenum, 'Sentinel mtime should not change on re-run');
        end

        function test_skips_when_no_simnibs_bin_path(tc)
            p = struct('startup', struct('simnibs_bin_path', ''));
            w = warning('off', 'all');
            out = preproc_deface_m2m(tc.M2mDir, p);
            warning(w);
            tc.verifyEmpty(out, 'Return value should be empty when simnibs_bin_path is unset');
            defaced = [tc.M2mDir, '_defaced'];
            tc.verifyFalse(isfolder(defaced), '_defaced folder should not be created when bin_path is empty');
        end

        function test_skips_when_pydeface_missing(tc)
            p = struct('startup', struct('simnibs_bin_path', fullfile(tc.TmpDir, 'nonexistent')));
            w = warning('off', 'all');
            out = preproc_deface_m2m(tc.M2mDir, p);
            warning(w);
            tc.verifyEmpty(out, 'Return value should be empty when pydeface is missing');
            defaced = [tc.M2mDir, '_defaced'];
            tc.verifyFalse(isfolder(defaced), '_defaced folder should not be created when pydeface is missing');
        end

    end

end
