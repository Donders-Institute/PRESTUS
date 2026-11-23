classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'water'}) test_telemetry_pipeline < matlab.unittest.TestCase
% TEST_TELEMETRY_PIPELINE  PRESTUS_PIPELINE sends its telemetry events.
%
% A small 2D free-water run, opted in, with webwrite shadowed by a recording
% stub (no request leaves the machine): the run sends run_start, then run_end
% with the same run_id and status success; a run whose acoustic stage fails
% sends run_error with that run_id, the error identifier and failed_stage
% 'acoustic'. The JVM's user.home is redirected, so the real ~/.prestus is
% never read or written. What the server accepts is TEST_TELEMETRY_LIVE.
%
%   results = runtests('tests/test_telemetry_pipeline.m');

    properties
        SimPath
        Parameters
        WebLog
    end

    methods (TestMethodSetup)
        function setup(tc)
            home = tempname;
            mkdir(home);
            tc.addTeardown(@() rmdir(home, 's'));
            orig_home = char(java.lang.System.getProperty('user.home'));
            java.lang.System.setProperty('user.home', home);
            tc.addTeardown(@() java.lang.System.setProperty('user.home', orig_home));
            tc.assertEqual(prefdir_prestus(), fullfile(home, '.prestus'), ...
                'user.home redirection failed; refusing to touch the real ~/.prestus');
            for v = {'SLURM_JOB_ID', 'PBS_JOBID', 'PRESTUS_RUN_TYPE', 'PRESTUS_TELEMETRY', ...
                     'PRESTUS_TEST_WEBWRITE_LOG'}
                old = getenv(v{1});
                tc.addTeardown(@() setenv(v{1}, old));
            end
            setenv('SLURM_JOB_ID', '');
            setenv('PBS_JOBID', '');
            setenv('PRESTUS_TELEMETRY', '');
            setenv('PRESTUS_RUN_TYPE', 'technical');

            stub_dir = tc.install_stubs({'webwrite', { ...
                'function out = webwrite(url, data, varargin)', ...
                '% Test stub: record the request instead of sending it.', ...
                'fid = fopen(getenv(''PRESTUS_TEST_WEBWRITE_LOG''), ''a'');', ...
                'fprintf(fid, ''%s\n'', data);', ...
                'fclose(fid);', ...
                'out = [];', ...
                'end'}});
            tc.WebLog = fullfile(stub_dir, 'requests.log');
            setenv('PRESTUS_TEST_WEBWRITE_LOG', tc.WebLog);
            tc.assertEqual(which('webwrite'), fullfile(stub_dir, 'webwrite.m'), ...
                'webwrite stub is not first on the path; refusing to risk a real request');
            evalc('telemetry_set_consent(true)');

            tc.SimPath = tempname;
            mkdir(tc.SimPath);
            tc.addTeardown(@() rmdir(tc.SimPath, 's'));
            p = make_minimal_parameters(tc.SimPath);
            p.grid.axisymmetric = 1;                 % 2D: seconds instead of minutes
            p.grid.default_dims = [73, 128];         % [2 x radial, axial], as test_physics_water
            p.transducer.trans_pos = [37, 5];
            p.transducer.focus_pos = [37, 64];
            tc.Parameters = p;
        end
    end

    methods (Test, TestTags = {'water', 'pipeline', 'telemetry'})

        function test_run_sends_start_then_end(tc)
            prestus_pipeline(tc.Parameters);
            ev = tc.events();
            tc.assertNotEmpty(ev, 'an opted-in run sent nothing');
            tc.verifyEqual(ev{1}.event, 'run_start', 'first event');
            tc.verifyEqual(ev{end}.event, 'run_end', 'last event');
            tc.verifyEqual(ev{end}.status, 'success');
            tc.verifyEqual(ev{end}.run_id, ev{1}.run_id, 'run_end does not pair with run_start');
            tc.verifyTrue(all(cellfun(@(e) strcmp(e.run_type, 'technical'), ev)), ...
                'an event lost the run type');
            tc.verifyFalse(any(cellfun(@(e) strcmp(e.event, 'run_error'), ev)));
        end

        function test_failed_stage_is_reported(tc)
            % The acoustic stage throws; the run fails and says where.
            tc.install_stubs({'acoustic_stage', { ...
                'function varargout = acoustic_stage(varargin)', ...
                'error(''test:telemetry:injected'', ''injected acoustic failure'');', ...
                'end'}});
            tc.verifyError(@() prestus_pipeline(tc.Parameters), 'test:telemetry:injected');
            ev = tc.events();
            tc.assertNotEmpty(ev, 'a failed run sent nothing');
            names = cellfun(@(e) e.event, ev, 'UniformOutput', false);
            tc.verifyEqual(names{1}, 'run_start');
            tc.verifyFalse(any(strcmp(names, 'run_end')), 'a failed run sent run_end');
            k = find(strcmp(names, 'run_error'), 1, 'last');
            tc.assertNotEmpty(k, 'no run_error event');
            err = ev{k};
            tc.verifyEqual(err.status, 'error');
            tc.verifyEqual(err.error_id, 'test:telemetry:injected');
            tc.verifyEqual(err.failed_stage, 'acoustic');
            tc.verifyEqual(err.run_id, ev{1}.run_id, 'run_error does not pair with run_start');
        end

    end

    methods (Access = private)

        function ev = events(tc)
            % The payloads the stub recorded, in order.
            ev = {};
            if ~isfile(tc.WebLog)
                return
            end
            lines = strsplit(strtrim(fileread(tc.WebLog)), newline);
            ev = cellfun(@jsondecode, lines, 'UniformOutput', false);
        end

        function d = install_stubs(tc, stubs)
            % stubs: n-by-2 cell {name, {lines}} in a fresh folder first on the path.
            d = tempname;
            mkdir(d);
            tc.addTeardown(@() rmdir(d, 's'));
            for k = 1:size(stubs, 1)
                fid = fopen(fullfile(d, [stubs{k, 1} '.m']), 'w');
                fprintf(fid, '%s\n', stubs{k, 2}{:});
                fclose(fid);
            end
            tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'MATLAB:dispatcher:nameConflict'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(d));
            rehash;
        end

    end

end
