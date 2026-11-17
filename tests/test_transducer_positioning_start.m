classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_transducer_positioning_start < matlab.unittest.TestCase
% TEST_TRANSDUCER_POSITIONING_START  Platform dispatch of TRANSDUCER_POSITIONING_START.
%
% TRANSDUCER_POSITIONING, HPC_DETECT_SYSTEM, HPC_VALIDATE_PARAMETERS,
% HPC_SETUP_TEMP_FILES, HPC_SUBMIT_JOB and HPC_WAIT_FOR_COMPLETION are
% shadowed by recording stubs, so nothing is computed or submitted.
%
%   results = runtests('tests/test_transducer_positioning_start.m');

    properties
        StubDir
    end

    methods (TestMethodSetup)
        function install_stubs(tc)
            tc.StubDir = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            d = tc.StubDir;
            write_stub(d, 'transducer_positioning', 'function transducer_positioning(varargin)', '');
            write_stub(d, 'hpc_detect_system', 'function p = hpc_detect_system(varargin)', 'p = ''slurm'';');
            write_stub(d, 'hpc_validate_parameters', 'function parameters = hpc_validate_parameters(varargin)', ...
                'parameters = varargin{1}; parameters.simulation.interactive = false;');
            write_stub(d, 'hpc_setup_temp_files', ...
                'function [log_dir, prestus_path, temp_data_path, temp_m_path, temp_m_file] = hpc_setup_temp_files(varargin)', ...
                sprintf(['log_dir = ''%s''; prestus_path = log_dir; temp_m_file = ''job_m''; ' ...
                         'temp_data_path = fullfile(log_dir, ''job.mat''); temp_m_path = fullfile(log_dir, ''job_m.m'');'], d));
            write_stub(d, 'hpc_submit_job', 'function [job_id, parameters] = hpc_submit_job(varargin)', ...
                'job_id = 4242; parameters = varargin{3};');
            write_stub(d, 'hpc_wait_for_completion', 'function hpc_wait_for_completion(varargin)', '');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(d));
            rehash;
        end
    end

    methods (Test, TestTags = {'transducer_positioning_start'})

        function test_matlab_platform_without_hpc_runs_locally(tc)
            % platform 'matlab' runs in MATLAB even when parameters has no
            % hpc struct (no scheduler detection).
            p = struct('platform', 'matlab');
            evalc('transducer_positioning_start(p, struct(), ''M1'', struct())');
            tc.verifyTrue(tc.called('transducer_positioning'));
            tc.verifyFalse(tc.called('hpc_detect_system'));
            tc.verifyFalse(tc.called('hpc_submit_job'));
        end

        function test_auto_platform_detects(tc)
            % Negative control: 'auto' does detect the scheduler.
            p = hpc_params('auto', 0);
            evalc('transducer_positioning_start(p, struct(), ''M1'', struct())');
            tc.verifyTrue(tc.called('hpc_detect_system'));
            tc.verifyTrue(tc.called('hpc_submit_job'));
            tc.verifyFalse(tc.called('transducer_positioning'));
        end

        function test_hpc_platforms_wait_and_report(tc)
            % Both schedulers submit; waiting prints the job ID.
            for platform = {'slurm', 'qsub'}
                p = hpc_params(platform{1}, 1);
                out = evalc('transducer_positioning_start(p, struct(), ''M1'', struct())');
                tc.verifySubstring(out, 'Job 4242 completed', platform{1});
                tc.verifyTrue(tc.called('hpc_wait_for_completion'), platform{1});
                sub = load(fullfile(tc.StubDir, 'hpc_submit_job.mat'));
                tc.verifyEqual(sub.args{1}, platform{1});
                tc.verifyEqual(sub.args{3}.hpc.job_prefix, 'TP');
                % The validated parameters (interactive off) reach the job
                saved = load(fullfile(tc.StubDir, 'job.mat'), 'parameters');
                tc.verifyFalse(logical(saved.parameters.simulation.interactive), platform{1});
                delete(fullfile(tc.StubDir, '*.mat'));
            end
        end

    end

    methods (Access = private)
        function tf = called(tc, name)
            tf = isfile(fullfile(tc.StubDir, [name '.mat']));
        end
    end

end

% =========================================================================

function p = hpc_params(platform, wait)
    p.platform = platform;
    p.subject_id = 1;
    p.simulation.interactive = 1;
    p.io.overwrite_files = 'always';
    p.hpc = struct('wait_for_job', wait, 'max_wait_checks', 1, ...
        'memorylimit', 20, 'timelimit', '01:00:00', 'job_prefix', 'PS');
end

function write_stub(dir_, name, signature, body)
    fid = fopen(fullfile(dir_, [name '.m']), 'w');
    fprintf(fid, '%s\n', signature, ...
        'args = varargin; %#ok<NASGU>', ...
        sprintf('save(fullfile(fileparts(mfilename(''fullpath'')), ''%s.mat''), ''args'');', name), ...
        body, 'end');
    fclose(fid);
end
