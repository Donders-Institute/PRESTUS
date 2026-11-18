classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_hpc_scratch < matlab.unittest.TestCase
% TEST_HPC_SCRATCH  Node-local scratch for the k-Wave .h5 files of a SLURM job.
%
% HPC_SCRATCH_SETUP creates <root>/prestus_<jobid>/sub-NNN<affix>/ only with
% hpc.tmp_gb > 0 inside a SLURM job, finds the root from hpc.scratch_dir or
% the environment, and its cleanup removes the folder again; a nested call
% that reuses the folder does not remove it. The SLURM environment is faked
% with environment variables and restored after each test.
%
%   results = runtests('tests/test_hpc_scratch.m');

    properties
        tmp
    end

    methods (TestMethodSetup)
        function fake_slurm(tc)
            tc.tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            names = {'SLURM_JOB_ID', 'SLURM_JOB_NAME', 'SLURM_TMPDIR', 'PRESTUS_SCRATCH', 'PRESTUS_SCRATCH_BASE', 'TMPDIR', 'USER'};
            for k = 1:numel(names)
                old = getenv(names{k});
                tc.addTeardown(@() setenv(names{k}, old));
                setenv(names{k}, '');
            end
            setenv('SLURM_JOB_ID', '4242');
            setenv('SLURM_JOB_NAME', 'PS_sub-001');
            setenv('USER', 'tester');
        end
    end

    methods (Test)

        function test_folder_created_and_removed(tc)
            setenv('PRESTUS_SCRATCH', tc.tmp);
            [p, cleanup] = hpc_scratch_setup(params(50, 'auto'));
            expected = fullfile(tc.tmp, 'prestus_4242', 'sub-001_x');
            tc.verifyEqual(p.io.dir_scratch, expected);
            tc.verifyTrue(isfolder(expected));
            fid = fopen(fullfile(expected, 'left_behind.h5'), 'w'); fclose(fid);
            clear cleanup
            tc.verifyFalse(isfolder(expected));
            tc.verifyFalse(isfolder(fullfile(tc.tmp, 'prestus_4242')));
        end

        function test_off_without_tmp_gb_or_outside_slurm(tc)
            setenv('PRESTUS_SCRATCH', tc.tmp);
            for gb = {0, [], 'absent'}
                p_in = params(gb{1}, 'auto');
                if strcmp(gb{1}, 'absent'), p_in.hpc = rmfield(p_in.hpc, 'tmp_gb'); end
                [p, cleanup] = hpc_scratch_setup(p_in);
                tc.verifyFalse(isfield(p.io, 'dir_scratch'));
                tc.verifyEmpty(cleanup);
            end
            setenv('SLURM_JOB_ID', '');
            [p, cleanup] = hpc_scratch_setup(params(50, 'auto'));
            tc.verifyFalse(isfield(p.io, 'dir_scratch'));
            tc.verifyEmpty(cleanup);
            tc.verifyEmpty(dir(fullfile(tc.tmp, 'prestus_*')));
        end

        function test_root_order_and_created_job_folder(tc)
            % no candidate exists, but <base>/$USER does: <base>/$USER/$SLURM_JOB_NAME is
            % created and removed again; an explicit scratch_dir wins over 'auto'
            setenv('PRESTUS_SCRATCH_BASE', tc.tmp);
            mkdir(fullfile(tc.tmp, 'tester'));
            [p, cleanup] = hpc_scratch_setup(params(50, 'auto'));
            job_root = fullfile(tc.tmp, 'tester', 'PS_sub-001');
            tc.verifyTrue(startsWith(p.io.dir_scratch, job_root));
            clear cleanup
            tc.verifyFalse(isfolder(job_root));

            explicit = fullfile(tc.tmp, 'own_${SLURM_JOB_ID}');
            [p, cleanup] = hpc_scratch_setup(params(50, explicit)); %#ok<ASGLU>
            tc.verifyTrue(startsWith(p.io.dir_scratch, fullfile(tc.tmp, 'own_4242')));
        end

        function test_unset_variable_never_collapses(tc)
            % a path naming an unset variable is skipped; with nothing else found,
            % the run warns and keeps the default .h5 location
            [p, cleanup] = tc.verifyWarning(@() hpc_scratch_setup(params(50, '/x/${NOT_SET_ANYWHERE}')), ...
                'prestus:scratch');
            tc.verifyFalse(isfield(p.io, 'dir_scratch'));
            tc.verifyEmpty(cleanup);
        end

        function test_nested_call_does_not_remove(tc)
            setenv('PRESTUS_SCRATCH', tc.tmp);
            [p, outer] = hpc_scratch_setup(params(50, 'auto'));
            [~, inner] = hpc_scratch_setup(params(50, 'auto'));
            tc.verifyEmpty(inner);
            tc.verifyTrue(isfolder(p.io.dir_scratch));
            clear outer
            tc.verifyFalse(isfolder(p.io.dir_scratch));
        end

    end
end

function p = params(tmp_gb, scratch_dir)
    p = struct('subject_id', 1);
    p.hpc = struct('tmp_gb', tmp_gb, 'scratch_dir', scratch_dir);
    p.io = struct('output_affix', '_x');
end
