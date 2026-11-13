classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_hpc_submit < matlab.unittest.TestCase
% TEST_HPC_SUBMIT  Job IDs from sbatch output (PARSE_SBATCH_JOB_ID).
%
% Regression: the first number anywhere in sbatch's output was taken as the
% job ID, so a warning line such as 'sbatch: warning: 2 GPUs ...' gave job 2.
%
%   results = runtests('tests/test_hpc_submit.m');


    methods (Test, TestTags = {'hpc_submit'})

        function test_parsable_plain(tc)
            tc.verifyEqual(parse_sbatch_job_id(sprintf('4815162\n')), 4815162);
        end

        function test_parsable_with_cluster(tc)
            tc.verifyEqual(parse_sbatch_job_id('4815162;snellius'), 4815162);
        end

        function test_warnings_before_job_id(tc)
            out = sprintf(['sbatch: warning: requested 2 GPUs on partition gpu40g\n' ...
                           'sbatch: note: 128 GB memory\n4815162\n']);
            tc.verifyEqual(parse_sbatch_job_id(out), 4815162);
        end

        function test_classic_output(tc)
            tc.verifyEqual(parse_sbatch_job_id('Submitted batch job 271828'), 271828);
        end

        function test_classic_output_after_warning(tc)
            out = sprintf('sbatch: warning: 2 GPUs\nSubmitted batch job 271828\n');
            tc.verifyEqual(parse_sbatch_job_id(out), 271828);
        end

        function test_no_id_is_nan(tc)
            tc.verifyTrue(isnan(parse_sbatch_job_id('sbatch: error: invalid partition')));
            tc.verifyTrue(isnan(parse_sbatch_job_id('')));
        end

        function test_string_input(tc)
            tc.verifyEqual(parse_sbatch_job_id("123"), 123);
        end

    end


    % Job scripts
    methods (Test, TestTags = {'hpc_submit'})

        function test_qsub_script_walltime_and_cpu_request(tc)
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            p = hpc_params();
            p.simulation.code_type = 'matlab_cpu';
            try %#ok<TRYNC>
                hpc_submit_job('qsub', 'job_m', p, tmp);
            end
            txt = fileread(only_file(tc, tmp, 'temp_qsub_*.sh'));
            tc.verifySubstring(txt, 'walltime=04:00:00');
            tc.verifyEmpty(strfind(txt, 'gpus='));
        end

        function test_snellius_timelimit_check(tc)
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            p = hpc_params();
            p.hpc.name = 'snellius';
            p.hpc.partition = 'gpu_a100';
            p.snellius.gpu_a100 = struct('memorylimit', 120, 'cores', 18, ...
                'timelimit', '120:00:00', 'n_gpu', 1);
            try %#ok<TRYNC>
                hpc_submit_job('slurm', 'job_m', p, tmp);
            end
            tc.verifySubstring(fileread(only_file(tc, tmp, 'temp_slurm_*.sh')), '--time=04:00:00');

            p.hpc.timelimit = '121:00:00';
            tc.verifyError(@() hpc_submit_job('slurm', 'job_m', p, tmp), ?MException);
        end

    end
end

function p = hpc_params()
    p.subject_id = 1;
    p.hpc = struct('name', 'default', 'partition', '', 'gpu', '', 'reservation', '', ...
        'timelimit', '04:00:00', 'memorylimit', 20, 'job_prefix', 'PS');
    p.simulation.code_type = 'matlab_gpu';
    p.io.output_affix = '';
end

function f = only_file(tc, folder, pattern)
    d = dir(fullfile(folder, pattern));
    tc.assertNumElements(d, 1, sprintf('Expected one %s', pattern));
    f = fullfile(d.folder, d.name);
end
