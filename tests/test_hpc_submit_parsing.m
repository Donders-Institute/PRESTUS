classdef test_hpc_submit_parsing < matlab.unittest.TestCase
% TEST_HPC_SUBMIT_PARSING  Job IDs from sbatch output (PARSE_SBATCH_JOB_ID).
%
% Regression: the first number anywhere in sbatch's output was taken as the
% job ID, so a warning line such as 'sbatch: warning: 2 GPUs ...' gave job 2.
%
%   results = runtests('tests/test_hpc_submit_parsing.m');

    methods (TestClassSetup)
        function add_paths(~)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(repo, 'functions')));
        end
    end

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

end
