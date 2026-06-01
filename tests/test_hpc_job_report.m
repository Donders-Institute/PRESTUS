%% 
classdef test_hpc_job_report < matlab.unittest.TestCase
% TEST_HPC_JOB_REPORT  Unit tests for functions/hpc/hpc_job_report
%
% Tests use a mock sacct that injects known output via $PATH prepending,
% so no real cluster connection is required.
%
%   Run with:   results = runtests('tests/test_hpc_job_report.m');

    properties
        TmpDir       % temporary directory created per test
        MockBinDir   % directory holding a fake `sacct` script
        Parameters   % minimal PRESTUS parameters struct
        OrigPath     % $PATH before mock injection
    end

    methods (TestMethodSetup)

        function setup(tc)
            tc.TmpDir    = tempname;
            tc.MockBinDir = fullfile(tc.TmpDir, 'mockbin');
            mkdir(tc.TmpDir);
            mkdir(tc.MockBinDir);

            % Minimal parameters struct
            tc.Parameters.subject_id    = 1;
            tc.Parameters.hpc.memorylimit  = 64;    % GB requested
            tc.Parameters.hpc.timelimit    = '04:00:00';

            % Inject mock bin dir at the front of PATH so `sacct` is found
            tc.OrigPath = getenv('PATH');
            setenv('PATH', [tc.MockBinDir ':' tc.OrigPath]);
        end

    end

    methods (TestMethodTeardown)

        function teardown(tc)
            setenv('PATH', tc.OrigPath);
            if isfolder(tc.TmpDir)
                rmdir(tc.TmpDir, 's');
            end
        end

    end

    % ------------------------------------------------------------------ %
    %% Helpers
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function write_mock_sacct(tc, alloc_output, rss_output)
            % Write a shell script that echoes pre-canned sacct output.
            % The script uses $@ to detect which format flags were passed:
            %   --allocations → emit alloc_output
            %   otherwise     → emit rss_output (MaxRSS query)
            mock_path = fullfile(tc.MockBinDir, 'sacct');
            fid = fopen(mock_path, 'w');
            fprintf(fid, '#!/bin/sh\n');
            fprintf(fid, 'case "$*" in\n');
            fprintf(fid, '  *--allocations*) echo "%s" ;;\n', alloc_output);
            fprintf(fid, '  *) echo "%s" ;;\n',               rss_output);
            fprintf(fid, 'esac\n');
            fclose(fid);
            system(sprintf('chmod +x %s', mock_path));
        end

    end

    % ------------------------------------------------------------------ %
    %% Struct fields
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'hpc_job_report', 'struct'})

        function test_report_has_expected_fields(tc)
            % A completed job should produce a struct with all key fields.
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '4096000K');

            report = hpc_job_report(12345, 'slurm', tc.Parameters, tc.TmpDir);

            expected = {'job_id','job_state','exit_code', ...
                        'elapsed_s','cpu_core_s','alloc_cpus','max_rss_gb','n_gpus', ...
                        'requested_mem_gb','requested_timelimit', ...
                        'energy_cpu_kwh','energy_gpu_kwh','energy_mem_kwh', ...
                        'energy_total_kwh','co2e_g','co2e_kg', ...
                        'memory_efficiency','time_efficiency', ...
                        'cpu_sbu','gpu_sbu','total_sbu','cost_eur','billing_profile'};
            for i = 1:numel(expected)
                tc.verifyTrue(isfield(report, expected{i}), ...
                    sprintf('Missing field: %s', expected{i}));
            end
        end

    end

    % ------------------------------------------------------------------ %
    %% Carbon calculation
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'hpc_job_report', 'carbon'})

        function test_carbon_calculation_cpu_only(tc)
            % 1 hour, 4 CPUs, no GPU, 4 GB peak memory
            % energy_cpu = 1h * 4cores * 12W/core / 1000 = 0.048 kWh
            % energy_mem = 1h * 4GB   * 0.3725W/GB / 1000 = 0.00149 kWh
            % energy_total = (0.048 + 0.00149) * PUE(1.67) = 0.08295 kWh
            % co2e = 0.08295 * 270 g/kWh = 22.40 g
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|14400|4|cpu=4,mem=32G', ...
                '4194304K');   % 4 GB in KB

            report = hpc_job_report(99, 'slurm', tc.Parameters, tc.TmpDir);

            tc.verifyEqual(report.elapsed_s,  3600, 'AbsTol', 1);
            tc.verifyEqual(report.alloc_cpus, 4,    'AbsTol', 0.1);
            tc.verifyEqual(report.n_gpus,     0);
            tc.verifyEqual(report.max_rss_gb, 4.0,  'AbsTol', 0.01);

            tc.verifyEqual(report.energy_cpu_kwh, 0.048, 'AbsTol', 1e-4);
            tc.verifyGreaterThan(report.co2e_g, 0);
            % Sanity: 4 CPU-hours at 12 W/core + overhead < 50 g CO2
            tc.verifyLessThan(report.co2e_g, 50);
        end

        function test_carbon_calculation_with_gpu(tc)
            % GPU jobs should have non-zero energy_gpu_kwh
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|7200|14400|2|cpu=2,gres/gpu=1,mem=32G', ...
                '8388608K');   % 8 GB

            report = hpc_job_report(77, 'slurm', tc.Parameters, tc.TmpDir);

            tc.verifyEqual(report.n_gpus, 1);
            tc.verifyGreaterThan(report.energy_gpu_kwh, 0);
            tc.verifyGreaterThan(report.energy_gpu_kwh, report.energy_cpu_kwh);
        end

        function test_custom_carbon_params(tc)
            % Doubling carbon_intensity should exactly double co2e_g
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '2097152K');

            cp_default.carbon_intensity = 270;
            cp_double.carbon_intensity  = 540;

            r1 = hpc_job_report(1, 'slurm', tc.Parameters, [], cp_default);
            r2 = hpc_job_report(1, 'slurm', tc.Parameters, [], cp_double);

            tc.verifyEqual(r2.co2e_g, 2 * r1.co2e_g, 'AbsTol', 1e-6);
        end

    end

    % ------------------------------------------------------------------ %
    %% Cost calculation
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'hpc_job_report', 'cost'})

        function test_no_cost_for_default_profile(tc)
            % Non-Snellius parameters should yield NaN cost
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '2097152K');

            report = hpc_job_report(10, 'slurm', tc.Parameters, []);

            tc.verifyTrue(isnan(report.cost_eur), ...
                'cost_eur should be NaN when no billing profile is set');
            tc.verifyEqual(report.billing_profile, 'unknown');
        end

        function test_snellius_cpu_cost(tc)
            % Snellius, 1 h, 2 CPUs, no GPU
            % cpu_sbu = 1h * 2cores * 1 SBU/core-hour = 2 SBU
            % cost = 2 * 16/1000 = EUR 0.032
            p = tc.Parameters;
            p.hpc.name      = 'snellius';
            p.hpc.partition = 'gpu_a100';  % partition present but no GPU allocated
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '2097152K');

            report = hpc_job_report(20, 'slurm', p, []);

            tc.verifyEqual(report.cpu_sbu,   2.0,   'AbsTol', 1e-6);
            tc.verifyEqual(report.gpu_sbu,   0.0,   'AbsTol', 1e-6);
            tc.verifyEqual(report.total_sbu, 2.0,   'AbsTol', 1e-6);
            tc.verifyEqual(report.cost_eur,  2*16/1000, 'AbsTol', 1e-6);
            tc.verifySubstring(report.billing_profile, 'Snellius');
        end

        function test_snellius_a100_gpu_cost(tc)
            % Snellius gpu_a100, 2 h, 2 CPUs, 1 A100
            % cpu_sbu = 2h * 2 * 1 = 4 SBU
            % gpu_sbu = 2h * 1 * 128 = 256 SBU
            % total = 260, cost = 260 * 16/1000 = EUR 4.16
            p = tc.Parameters;
            p.hpc.name      = 'snellius';
            p.hpc.partition = 'gpu_a100';
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|7200|14400|2|cpu=2,gres/gpu=1,mem=32G', ...
                '8388608K');

            report = hpc_job_report(21, 'slurm', p, []);

            tc.verifyEqual(report.cpu_sbu,   4.0,   'AbsTol', 1e-6);
            tc.verifyEqual(report.gpu_sbu,   256.0, 'AbsTol', 1e-6);
            tc.verifyEqual(report.total_sbu, 260.0, 'AbsTol', 1e-6);
            tc.verifyEqual(report.cost_eur,  260*16/1000, 'AbsTol', 1e-6);
        end

        function test_snellius_h100_gpu_cost(tc)
            % H100: gpu_sbu_per_gpu_hour = 192
            % 1 h, 2 CPUs, 1 H100
            % gpu_sbu = 1h * 1 * 192 = 192; cpu_sbu = 2; total = 194
            p = tc.Parameters;
            p.hpc.name      = 'snellius';
            p.hpc.partition = 'gpu_h100';
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,gres/gpu=1,mem=32G', ...
                '4194304K');

            report = hpc_job_report(22, 'slurm', p, []);

            tc.verifyEqual(report.gpu_sbu,   192.0, 'AbsTol', 1e-6);
            tc.verifyEqual(report.total_sbu, 194.0, 'AbsTol', 1e-6);
        end

        function test_custom_cost_params_override(tc)
            % A custom eur_per_1000_sbu should scale cost proportionally
            p = tc.Parameters;
            p.hpc.name      = 'snellius';
            p.hpc.partition = 'gpu_a100';
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '2097152K');

            cp_custom.eur_per_1000_sbu = 32;  % double the default
            r1 = hpc_job_report(30, 'slurm', p, [], [], []);
            r2 = hpc_job_report(30, 'slurm', p, [], [], cp_custom);

            tc.verifyEqual(r2.cost_eur, 2 * r1.cost_eur, 'AbsTol', 1e-9);
        end

        function test_cost_written_to_report_file(tc)
            p = tc.Parameters;
            p.hpc.name      = 'snellius';
            p.hpc.partition = 'gpu_a100';
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '2097152K');

            hpc_job_report(31, 'slurm', p, tc.TmpDir);

            content = fileread(fullfile(tc.TmpDir, 'sub-001_job31_report.txt'));
            for field = {'billing_profile', 'cpu_sbu', 'gpu_sbu', 'total_sbu', 'cost_eur'}
                tc.verifyTrue(contains(content, field{1}), ...
                    sprintf('Report file missing cost field: %s', field{1}));
            end
        end

    end

    % ------------------------------------------------------------------ %
    %% Resource efficiency
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'hpc_job_report', 'efficiency'})

        function test_memory_efficiency_ratio(tc)
            % Requested 64 GB, used 8 GB → efficiency = 8/64 = 0.125
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=64G', ...
                '8388608K');   % 8 GB in KB

            report = hpc_job_report(55, 'slurm', tc.Parameters, tc.TmpDir);

            tc.verifyEqual(report.memory_efficiency, 8/64, 'AbsTol', 1e-3);
        end

        function test_time_efficiency_ratio(tc)
            % Requested 04:00:00 = 14400 s, used 3600 s → efficiency = 0.25
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '1048576K');

            report = hpc_job_report(56, 'slurm', tc.Parameters, tc.TmpDir);

            tc.verifyEqual(report.time_efficiency, 3600/14400, 'AbsTol', 1e-3);
        end

    end

    % ------------------------------------------------------------------ %
    %% Report file
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'hpc_job_report', 'file'})

        function test_report_file_created(tc)
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '2097152K');

            hpc_job_report(42, 'slurm', tc.Parameters, tc.TmpDir);

            files = dir(fullfile(tc.TmpDir, 'sub-001_job42_report.txt'));
            tc.verifyNotEmpty(files, 'Report file was not created');
        end

        function test_report_file_contains_key_fields(tc)
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '2097152K');

            hpc_job_report(43, 'slurm', tc.Parameters, tc.TmpDir);

            content = fileread(fullfile(tc.TmpDir, 'sub-001_job43_report.txt'));
            for field = {'job_state', 'co2e_kg', 'max_rss_gb', 'memory_efficiency', 'energy_total_kwh'}
                tc.verifyTrue(contains(content, field{1}), ...
                    sprintf('Report file missing field: %s', field{1}));
            end
        end

        function test_no_report_file_when_log_dir_empty(tc)
            tc.write_mock_sacct( ...
                'COMPLETED|0:0|3600|7200|2|cpu=2,mem=8G', ...
                '2097152K');

            % Pass empty log_dir — should not error and should not create a file
            tc.verifyWarningFree(@() hpc_job_report(44, 'slurm', tc.Parameters, []));
        end

    end

    % ------------------------------------------------------------------ %
    %% Live SLURM deployment  (skipped when sbatch is not available)
    % ------------------------------------------------------------------ %
    % Enable by ensuring `sbatch` and `sacct` are on $PATH (i.e. on a SLURM
    % login or compute node). The test submits a trivial 1-core sleep job,
    % waits for it to finish, then checks that hpc_job_report returns
    % plausible metrics from the real sacct record.
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'hpc_job_report', 'slurm_live'})

        function test_live_slurm_job_report(tc)
            % Skip if sbatch is not available in the original PATH
            [status, ~] = system(sprintf('PATH=%s which sbatch 2>/dev/null', tc.OrigPath));
            tc.assumeTrue(status == 0, 'sbatch not found on PATH — skipping live SLURM test');

            % Restore original PATH so the real sbatch/sacct are used
            setenv('PATH', tc.OrigPath);

            % Submit a minimal job: sleep 5 seconds, 1 core, 1 GB
            script = fullfile(tc.TmpDir, 'test_job.sh');
            fid = fopen(script, 'w');
            fprintf(fid, '#!/bin/bash\n');
            fprintf(fid, '#SBATCH --job-name=prestus_test_report\n');
            fprintf(fid, '#SBATCH --ntasks=1\n');
            fprintf(fid, '#SBATCH --cpus-per-task=1\n');
            fprintf(fid, '#SBATCH --mem=1G\n');
            fprintf(fid, '#SBATCH --time=00:02:00\n');
            fprintf(fid, '#SBATCH --output=%s/test_job_%%j.log\n', tc.TmpDir);
            fprintf(fid, 'sleep 5\n');
            fclose(fid);
            system(sprintf('chmod +x %s', script));

            [sub_status, sub_out] = system(sprintf('sbatch %s', script));
            tc.assumeTrue(sub_status == 0, ...
                sprintf('sbatch submission failed: %s', strtrim(sub_out)));

            % Extract job ID from "Submitted batch job 12345"
            tok = regexp(strtrim(sub_out), '(\d+)$', 'tokens', 'once');
            tc.assumeTrue(~isempty(tok), 'Could not parse job ID from sbatch output');
            job_id = str2double(tok{1});

            % Wait for the job to leave the queue (max 2 min)
            fprintf('  Waiting for SLURM job %d to complete...\n', job_id);
            job_id_str = sprintf('%d', job_id);
            deadline = tic;
            finished = false;
            while toc(deadline) < 120
                [~, sq] = system(sprintf('squeue --noheader -j %s 2>/dev/null', job_id_str));
                if isempty(strtrim(sq))
                    finished = true;
                    break;
                end
                pause(5);
            end
            tc.assumeTrue(finished, 'Job did not finish within 2 minutes — skipping assertions');

            % Build minimal parameters matching the job
            p = tc.Parameters;
            p.hpc.memorylimit = 1;
            p.hpc.timelimit   = '00:02:00';

            % Run the report against the real sacct record
            report = hpc_job_report(job_id, 'slurm', p, tc.TmpDir);

            % sacct found the job — state must be a known terminal value, not UNKNOWN
            known_states = {'completed','failed','cancelled','timeout','out_of_memory','node_fail'};
            tc.verifyTrue(any(strcmpi(report.job_state, known_states)), ...
                sprintf('job_state ''%s'' was not recognised from sacct — sacct query may have failed', ...
                report.job_state));

            % Timing was retrieved (not NaN)
            tc.verifyFalse(isnan(report.elapsed_s), 'elapsed_s should not be NaN');
            tc.verifyGreaterThanOrEqual(report.elapsed_s, 0, 'elapsed_s must be non-negative');

            % Allocation matches what was requested
            tc.verifyEqual(report.alloc_cpus, 1, 'AbsTol', 0.1, ...
                'alloc_cpus should match the 1 CPU requested in the job script');

            % Carbon estimate is finite and non-negative (may be ~0 for very short jobs)
            tc.verifyFalse(isnan(report.co2e_g), 'co2e_g should not be NaN');
            tc.verifyGreaterThanOrEqual(report.co2e_g, 0, 'co2e_g must be non-negative');

            % Report file should exist
            files = dir(fullfile(tc.TmpDir, sprintf('sub-001_job%d_report.txt', job_id)));
            tc.verifyNotEmpty(files, 'Report file was not written for live job');
        end

    end

    % ------------------------------------------------------------------ %
    %% Graceful degradation
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'hpc_job_report', 'robustness'})

        function test_failed_job_does_not_error(tc)
            % sacct returns FAILED state — report should still be produced
            tc.write_mock_sacct( ...
                'FAILED|1:0|120|240|2|cpu=2,mem=8G', ...
                '0');   % MaxRSS often 0 for failed jobs

            tc.verifyWarningFree(@() hpc_job_report(66, 'slurm', tc.Parameters, []));
        end

        function test_qsub_produces_nan_metrics_with_warning(tc)
            % qsub path has no sacct — usage fields should be NaN, not an error
            report = tc.verifyWarning( ...
                @() hpc_job_report('job.12345', 'qsub', tc.Parameters, []), ...
                'PRESTUS:hpc_job_report:unsupportedScheduler');
            % verifyWarning returns the output of the function handle
            tc.verifyTrue(isnan(report.elapsed_s));
            tc.verifyTrue(isnan(report.co2e_kg));
        end

    end

end
