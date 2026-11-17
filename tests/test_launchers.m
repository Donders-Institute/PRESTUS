classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_launchers < matlab.unittest.TestCase
% TEST_LAUNCHERS  Dispatch logic of CALIBRATION_PIPELINE_START and a headless-safe PRESTUS_GUI smoke test.
%
% CALIBRATION_PIPELINE_START: unknown platform errors; 'matlab' calls
% CALIBRATION_TRANSDUCER directly; 'slurm' serialises the arguments (forcing
% platform 'matlab' and interactive 0 inside the job, job prefix 'CAL') and
% writes a bootstrap script. CALIBRATION_TRANSDUCER, HPC_VALIDATE_PARAMETERS
% and HPC_SUBMIT_JOB are shadowed by recording stubs, so nothing is simulated
% or submitted.
%
% PRESTUS_GUI: builds the window and closes it. Skipped when MATLAB cannot
% show figures (no display, e.g. CI on Linux). Telemetry consent is
% pre-recorded in a temporary user.home so no consent dialog can block.
%
%   results = runtests('tests/test_launchers.m');

    properties
        Dir
        StubDir
    end


    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tempname;
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
        end
    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'calibration_pipeline_start'})

        function test_unknown_platform_errors(tc)
            p = struct('platform', 'lsf');
            try
                evalc('calibration_pipeline_start(''p.csv'', ''CTX_500'', 30, 50, p, 1)');
                tc.verifyFail('no error for an unknown platform');
            catch ME
                tc.verifySubstring(ME.message, 'unknown platform ''lsf''');
            end
        end

        function test_matlab_platform_calls_calibration_directly(tc)
            tc.install_stubs();
            p = struct('platform', 'matlab');
            [~, job_id] = evalc('calibration_pipeline_start(''prof.csv'', ''CTX_500'', [20 30], 55, p, 7)');
            tc.verifyEmpty(job_id);
            rec = load(fullfile(tc.StubDir, 'calibration_transducer.mat'));
            tc.verifyEqual(rec.args{1}, 'prof.csv');
            tc.verifyEqual(rec.args{2}, 'CTX_500');
            tc.verifyEqual(rec.args{3}, [20 30]);
            tc.verifyEqual(rec.args{4}, 55);
            tc.verifyEqual(rec.args{6}, 7);
            tc.verifyFalse(isfile(fullfile(tc.StubDir, 'hpc_submit_job.mat')));
        end

        function test_slurm_platform_serialises_job(tc)
            tc.install_stubs();
            p = struct('platform', 'slurm', ...
                'simulation', struct('interactive', 1), ...
                'hpc', struct('wait_for_job', 0), ...
                'calibration', struct('path_output', tc.Dir));
            [~, job_id] = evalc('calibration_pipeline_start(''prof.csv'', ''CTX_500'', 30, 55, p, 7)');
            tc.verifyEqual(job_id, 4242);
            sub = load(fullfile(tc.StubDir, 'hpc_submit_job.mat'));
            tc.verifyEqual(sub.args{1}, 'slurm');
            log_dir = fullfile(tc.Dir, 'sub-007', 'log_hpc');
            tc.verifyEqual(sub.args{4}, log_dir);
            % Bootstrap script and payload
            script = fullfile(log_dir, [sub.args{2} '.m']);
            tc.assertTrue(isfile(script));
            txt = fileread(script);
            tc.verifySubstring(txt, 'calibration_transducer(cal_args.profile_empirical');
            mat = regexp(txt, 'load\(''([^'']+)''\)', 'tokens', 'once');
            tc.assertNotEmpty(mat);
            payload = load(mat{1});
            tc.verifyEqual(payload.parameters.platform, 'matlab');          % no nested submission
            tc.verifyEqual(payload.parameters.simulation.interactive, 0);
            tc.verifyEqual(payload.parameters.subject_id, 7);
            tc.verifyEqual(payload.parameters.hpc.job_prefix, 'CAL');
            tc.verifyEqual(payload.cal_args.equipment_name, 'CTX_500');
            tc.verifyEqual(payload.cal_args.desired_intensity, 30);
            % The job's parameters as submitted carry the HPC prefix too
            tc.verifyEqual(sub.args{3}.hpc.job_prefix, 'CAL');
        end

    end

    % ------------------------------------------------------------------ %
    methods (Access = private)

        function install_stubs(tc)
            % Recording stubs: each saves its inputs to <stubdir>/<name>.mat
            tc.StubDir = tempname;
            mkdir(tc.StubDir);
            tc.addTeardown(@() rmdir(tc.StubDir, 's'));
            write_stub(tc.StubDir, 'calibration_transducer', ...
                'function calibration_transducer(varargin)', '');
            write_stub(tc.StubDir, 'hpc_validate_parameters', ...
                'function parameters = hpc_validate_parameters(varargin)', ...
                'parameters = varargin{1};');
            write_stub(tc.StubDir, 'hpc_submit_job', ...
                'function [job_id, parameters] = hpc_submit_job(varargin)', ...
                'job_id = 4242; parameters = varargin{3};');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(tc.StubDir));
            rehash;
            for n = {'calibration_transducer', 'hpc_validate_parameters', 'hpc_submit_job'}
                tc.assertEqual(which(n{1}), fullfile(tc.StubDir, [n{1} '.m']), ...
                    sprintf('stub %s is not first on the path', n{1}));
            end
        end

    end

end

% =========================================================================

function write_stub(dir_, name, signature, body)
    fid = fopen(fullfile(dir_, [name '.m']), 'w');
    fprintf(fid, '%s\n', signature, ...
        'args = varargin; %#ok<NASGU>', ...
        sprintf('save(fullfile(fileparts(mfilename(''fullpath'')), ''%s.mat''), ''args'');', name), ...
        body, 'end');
    fclose(fid);
end
