classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_hpc_submit < matlab.unittest.TestCase
% TEST_HPC_SUBMIT  Job IDs from sbatch output, job scripts, parameter validation, and job polling.
%
% PARSE_SBATCH_JOB_ID takes the job ID, not the first number in sbatch's
% output (a warning line such as 'sbatch: warning: 2 GPUs ...' precedes it).
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

        function test_job_logs_are_per_job(tc)
            % Each job writes its own log. SLURM expands %j in #SBATCH
            % directives; PBS expands nothing in #PBS directives, so the
            % qsub script must not name its logs there with %j.
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            p = hpc_params();
            try %#ok<TRYNC>
                hpc_submit_job('qsub', 'job_m', p, tmp);
            end
            txt = fileread(only_file(tc, tmp, 'temp_qsub_*.sh'));
            tc.verifyEmpty(regexp(txt, '#PBS -[oe] ', 'once'), 'PBS -o/-e with a fixed name is shared by all jobs');
            tc.verifyEmpty(strfind(txt, '%j'));
            tc.verifySubstring(txt, 'sub-001_qsub_output_${PBS_JOBID}.log');
            tc.verifySubstring(txt, 'sub-001_qsub_error_${PBS_JOBID}.log');
            % The redirect precedes the MATLAB call
            tc.verifyLessThan(strfind(txt, 'exec >'), strfind(txt, 'matlab -batch'));

            try %#ok<TRYNC>
                hpc_submit_job('slurm', 'job_m', p, tmp);
            end
            txt = fileread(only_file(tc, tmp, 'temp_slurm_*.sh'));
            tc.verifySubstring(txt, '#SBATCH --output=sub-001_slurm_output_%j.log');
            tc.verifySubstring(txt, '#SBATCH --error=sub-001_slurm_error_%j.log');
        end

        function test_tmp_request_follows_tmp_gb(tc)
            % --tmp only with hpc.tmp_gb > 0; 0, empty, or no key requests nothing
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            cases = {50, '#SBATCH --tmp=50G'; 12.5, '#SBATCH --tmp=13G'; 0, ''; [], ''; 'absent', ''};
            for k = 1:size(cases, 1)
                p = hpc_params();
                if ~strcmp(cases{k, 1}, 'absent'), p.hpc.tmp_gb = cases{k, 1}; end
                delete(fullfile(tmp, 'temp_*.sh'));
                try %#ok<TRYNC>
                    hpc_submit_job('slurm', 'job_m', p, tmp);
                end
                txt = fileread(only_file(tc, tmp, 'temp_slurm_*.sh'));
                if isempty(cases{k, 2})
                    tc.verifyEmpty(strfind(txt, '--tmp'), sprintf('case %d', k));
                else
                    tc.verifySubstring(txt, cases{k, 2}, sprintf('case %d', k));
                end
            end
        end

        function test_job_inherits_run_type(tc)
            % A run type set in this session reaches both job scripts
            % (sbatch/qsub start from MATLAB's launch environment), before
            % MATLAB starts; deployment and invalid values write nothing.
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            orig = getenv('PRESTUS_RUN_TYPE');
            restore = onCleanup(@() setenv('PRESTUS_RUN_TYPE', orig)); %#ok<NASGU>
            p = hpc_params();
            cases = {'development', 'export PRESTUS_RUN_TYPE=development'; ...
                     'technical',   'export PRESTUS_RUN_TYPE=technical'; ...
                     '',            ''; ...
                     'rm -rf x',    ''};
            for k = 1:size(cases, 1)
                setenv('PRESTUS_RUN_TYPE', cases{k, 1});
                for sched = {'slurm', 'qsub'}
                    delete(fullfile(tmp, 'temp_*.sh'));
                    try %#ok<TRYNC>
                        hpc_submit_job(sched{1}, 'job_m', p, tmp);
                    end
                    txt = fileread(only_file(tc, tmp, ['temp_' sched{1} '_*.sh']));
                    label = sprintf('%s, PRESTUS_RUN_TYPE=''%s''', sched{1}, cases{k, 1});
                    if isempty(cases{k, 2})
                        tc.verifyEmpty(strfind(txt, 'PRESTUS_RUN_TYPE'), label);
                    else
                        tc.verifySubstring(txt, cases{k, 2}, label);
                        tc.verifyLessThan(strfind(txt, cases{k, 2}), strfind(txt, 'matlab -batch'), label);
                    end
                end
            end
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

        function test_gpu_request_follows_the_solvers(tc)
            % The GPU request comes from the solvers that run and their compute
            % backends, not from simulation.code_type alone.
            cases = { ...
                % acoustic, thermal, heating, code_type, jwave.use_gpu, bb backend, GPU?
                'kwave',      'kwave',      0, 'matlab_gpu', [],    'auto', true; ...
                'kwave',      'kwave',      0, 'matlab_cpu', [],    'auto', false; ...
                'babelbrain', 'babelbrain', 1, 'matlab_cpu', [],    'auto', true; ...
                'babelbrain', 'babelbrain', 1, 'matlab_gpu', [],    'cpu',  false; ...
                'babelbrain', 'kwave',      1, 'matlab_gpu', [],    'cpu',  true; ...   % k-Wave heating
                'jwave',      'kwave',      0, 'matlab_cpu', [],    'auto', true; ...   % auto-detect
                'jwave',      'kwave',      0, 'matlab_gpu', false, 'auto', false};
            for k = 1:size(cases, 1)
                tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
                p = hpc_params();
                [p.acoustic_backend, p.thermal_backend] = cases{k, 1:2};
                p.modules.run_heating_sims = cases{k, 3};
                p.simulation.code_type = cases{k, 4};
                p.jwave.use_gpu = cases{k, 5};
                p.babelbrain.computing_backend = cases{k, 6};
                try %#ok<TRYNC>
                    hpc_submit_job('slurm', 'job_m', p, tmp);
                end
                txt = fileread(only_file(tc, tmp, 'temp_slurm_*.sh'));
                tc.verifyEqual(contains(txt, '--partition=gpu'), cases{k, 7}, ...
                    sprintf('case %d: %s/%s, code_type %s', k, cases{k, 1}, cases{k, 2}, cases{k, 4}));
            end
        end

        function test_gpu_request_with_yaml_nulls(tc)
            % A YAML '~' loads as a yaml.Null object, which isempty does not
            % catch; every field the GPU decision reads must treat it as unset
            % (jwave.use_gpu: ~ is the toolbox default and stopped submission).
            cases = { ...
                % acoustic, thermal, null field, GPU?
                'jwave',      'kwave',      {'jwave', 'use_gpu'},                 true; ...   % auto-detect
                'babelbrain', 'babelbrain', {'babelbrain', 'computing_backend'},  true; ...   % default 'auto'
                'kwave',      'kwave',      {'modules', 'run_heating_sims'},      true; ...   % k-Wave acoustic on GPU
                'kwave',      'kwave',      {'time_reversal', 'enabled'},         true; ...
                'jwave',      'kwave',      {'simulation', 'code_type'},          true};      % jwave still auto-detects
            for k = 1:size(cases, 1)
                tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
                p = hpc_params();
                [p.acoustic_backend, p.thermal_backend] = cases{k, 1:2};
                p.simulation.code_type = 'matlab_gpu';
                p.modules.run_heating_sims = 1;
                p.jwave.use_gpu = [];
                p.babelbrain.computing_backend = 'auto';
                f = cases{k, 3};
                p.(f{1}).(f{2}) = yaml.Null;
                try %#ok<TRYNC>
                    hpc_submit_job('slurm', 'job_m', p, tmp);
                end
                txt = fileread(only_file(tc, tmp, 'temp_slurm_*.sh'));
                tc.verifyEqual(contains(txt, '--partition=gpu'), cases{k, 4}, ...
                    sprintf('case %d: %s.%s = yaml.Null', k, f{1}, f{2}));
            end
        end

    end


    % Parameter validation
    methods (Test, TestTags = {'hpc_validate_parameters'})

        function test_validate_returns_non_interactive(tc)
            for platform = {'slurm', 'qsub'}
                for interactive = [0 1]
                    p = hpc_params();
                    p.simulation.interactive = interactive;
                    p.io.overwrite_files = 'always';
                    [~, out] = evalc('hpc_validate_parameters(p, platform{1})');
                    tc.verifyFalse(logical(out.simulation.interactive));
                    tc.verifyEqual(rmfield(out, 'simulation'), rmfield(p, 'simulation'));
                end
            end
        end

        function test_validate_without_output(tc)
            % A caller that ignores the output still validates.
            p = hpc_params();
            p.simulation.interactive = 0;
            p.io.overwrite_files = 'never';
            tc.verifyWarningFree(@() hpc_validate_parameters(p, 'slurm'));
        end

        function test_validate_rejects_ask(tc)
            % Negative control: 'ask' would block an unattended job.
            p = hpc_params();
            p.simulation.interactive = 0;
            p.io.overwrite_files = 'ask';
            tc.verifyError(@() hpc_validate_parameters(p, 'slurm'), ?MException);
        end

    end


    % Job polling
    methods (Test, TestTags = {'hpc_wait_for_completion'})

        function test_qsub_wait_state_is_per_call(tc)
            % qstat stub: the job is never listed. Each call must first wait
            % for the job to appear, so with one check every call times out;
            % no state carries over from an earlier call.
            stubs = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            write_function(stubs, 'run_command', ...
                'function [status, out] = run_command(varargin)', 'status = 1; out = '''';');
            write_function(stubs, 'pause', 'function pause(varargin)', '');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(stubs));
            rehash;
            tc.assertEqual(which('run_command'), fullfile(stubs, 'run_command.m'));
            for k = 1:3
                out = evalc('hpc_wait_for_completion(''123.server'', ''qsub'', 1)');
                tc.verifySubstring(out, 'to appear', sprintf('call %d', k));
                tc.verifyEmpty(strfind(out, 'COMPLETE'), sprintf('call %d', k));
            end
            % Negative control: two checks see the job gone after it was waited for
            out = evalc('hpc_wait_for_completion(''123.server'', ''qsub'', 2)');
            tc.verifySubstring(out, 'COMPLETE');
        end

    end
end

function write_function(folder, name, signature, body)
    fid = fopen(fullfile(folder, [name '.m']), 'w');
    fprintf(fid, '%s\n', signature, body, 'end');
    fclose(fid);
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
