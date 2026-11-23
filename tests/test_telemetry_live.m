classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'telemetry_live'}) test_telemetry_live < matlab.unittest.TestCase
% TEST_TELEMETRY_LIVE  The telemetry server accepts what TRACK_USAGE sends.
%
% TRACK_USAGE builds the payloads of run_start, run_end (with simulation
% statistics and skull medians) and run_error (with failed_stage) from
% artificial parameters; webwrite is shadowed by a stub that records URL and
% body instead of sending. The test then posts each recorded body to the
% recorded URL itself (TRACK_USAGE ignores replies, so a rejected event would
% otherwise go unnoticed) and requires a 2xx reply without an error. No
% simulation runs.
%
% Needs network access; each run adds 4 rows to the real events table: three
% with run type 'technical' and one 'development' (the run type the server
% once rejected), both of which usage statistics leave out. No 'deployment'
% event is sent, as that would count as real usage. Run it weekly and before
% a release, not on every push:
%
%   run_all_tests('telemetry_live')
%
% The pipeline's own telemetry calls are TEST_TELEMETRY_PIPELINE (no network).

    properties
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
            for v = {'PRESTUS_RUN_TYPE', 'PRESTUS_TELEMETRY', 'PRESTUS_TEST_WEBWRITE_LOG'}
                old = getenv(v{1});
                tc.addTeardown(@() setenv(v{1}, old));
            end
            setenv('PRESTUS_TELEMETRY', '');

            d = tempname;
            mkdir(d);
            tc.addTeardown(@() rmdir(d, 's'));
            fid = fopen(fullfile(d, 'webwrite.m'), 'w');
            fprintf(fid, '%s\n', ...
                'function out = webwrite(url, data, varargin)', ...
                '% Test stub: record the request instead of sending it.', ...
                'fid = fopen(getenv(''PRESTUS_TEST_WEBWRITE_LOG''), ''a'');', ...
                'fprintf(fid, ''%s\n%s\n'', url, data);', ...
                'fclose(fid);', ...
                'out = [];', ...
                'end');
            fclose(fid);
            tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'MATLAB:dispatcher:nameConflict'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(d));
            rehash;
            tc.WebLog = fullfile(d, 'requests.log');
            setenv('PRESTUS_TEST_WEBWRITE_LOG', tc.WebLog);
            tc.assertEqual(which('webwrite'), fullfile(d, 'webwrite.m'), ...
                'webwrite stub is not first on the path');
            evalc('telemetry_set_consent(true)');    % current notice: skull medians included
        end
    end

    methods (Test, TestTags = {'telemetry'})

        function test_server_accepts_every_event(tc)
            setenv('PRESTUS_RUN_TYPE', 'technical');
            p = artificial_parameters();
            rid = generate_run_id();
            stats = struct('acoustic_dims', [240 71], 'acoustic_nt', 2871, 'acoustic_s', 12.5, ...
                'thermal_steps', 880, 'stage_acoustic_s', 12.5, 'peak_ram_gb', 3.2, ...
                'skull_c_median', 2550, 'skull_rho_median', 1850, 'skull_atten_median', 4.2, ...
                'skull_path_c_median', 2600, 'skull_path_rho_median', 1900, 'skull_path_atten_median', 4.5);
            track_usage('run_start', p, struct('run_id', rid));
            track_usage('run_end', p, struct('run_id', rid, 'status', 'success', ...
                'duration_s', 61.5, 'sim_stats', stats));
            track_usage('run_error', p, struct('run_id', generate_run_id(), 'status', 'error', ...
                'duration_s', 7.25, 'error_id', 'test:telemetry:live', 'failed_stage', 'acoustic', ...
                'sim_stats', struct('stage_preproc_s', 3.5)));
            setenv('PRESTUS_RUN_TYPE', 'development');
            track_usage('run_end', p, struct('run_id', generate_run_id(), 'status', 'success', ...
                'duration_s', 1));

            req = recorded_requests(tc.WebLog);
            tc.assertNumElements(req, 4, 'track_usage did not build four requests');
            sent = cellfun(@(r) jsondecode(r.body), req, 'UniformOutput', false);
            events = cellfun(@(s) s.event, sent, 'UniformOutput', false);
            run_types = cellfun(@(s) s.run_type, sent, 'UniformOutput', false);
            tc.verifyEqual(events, {'run_start', 'run_end', 'run_error', 'run_end'});
            tc.verifyEqual(run_types, {'technical', 'technical', 'technical', 'development'});
            tc.verifyTrue(isfield(sent{2}.sim_stats, 'skull_c_median'), 'skull medians not in the payload');
            for k = 1:numel(req)
                [ok, detail] = post_json(req{k}.url, req{k}.body);
                tc.verifyTrue(ok, sprintf('server rejected %s (%s): %s', ...
                    events{k}, run_types{k}, detail));
            end
        end

    end

end

% =========================================================================

function req = recorded_requests(weblog)
% {struct(url, body)} from the stub's log (alternating URL and body lines).
    req = {};
    if ~isfile(weblog)
        return
    end
    lines = strsplit(strtrim(fileread(weblog)), newline);
    for k = 1:2:numel(lines) - 1
        req{end + 1} = struct('url', strtrim(lines{k}), 'body', lines{k + 1}); %#ok<AGROW>
    end
end

function [ok, detail] = post_json(url, body)
% POST body (JSON text) as sent by TRACK_USAGE; ok for a 2xx reply whose body
% reports no error. matlab.net.http is not shadowed by the webwrite stub.
    import matlab.net.http.*
    msg = MessageBody();
    msg.Payload = unicode2native(body, 'UTF-8');
    request = RequestMessage('POST', field.ContentTypeField('application/json'), msg);
    try
        reply = request.send(matlab.net.URI(url), HTTPOptions('ConnectTimeout', 20));
    catch me
        ok = false;
        detail = me.message;
        return
    end
    code = double(reply.StatusCode);
    data = reply.Body.Data;
    if isstruct(data)
        detail = sprintf('HTTP %d %s', code, jsonencode(data));
        has_error = isfield(data, 'error') && ~isempty(data.error);
    else
        detail = sprintf('HTTP %d %s', code, strtrim(char(string(data))));
        has_error = false;
    end
    ok = code >= 200 && code < 300 && ~has_error;
end

function p = artificial_parameters()
    p = struct();
    p.subject_id = 1;
    p.platform = 'slurm';
    p.simulation = struct('medium', 'layered', 'code_type', 'cpp_gpu', 'precision', 'single');
    p.path = struct('sim', tempdir);
    p.hpc = struct('name', 'test', 'timelimit', '04:00:00', 'memorylimit', 32, ...
        'gpu', 'nvidia_a100', 'partition', 'gpu');
    p.grid = struct('resolution_mm', 0.5);
    p.transducer = struct('type', 'annular', 'freq_hz', 500e3, 'annular', struct('elem_n', 4));
    p.modules = struct('run_acoustic_sims', 1, 'run_heating_sims', 1);
    p.layers = struct('brain', 1, 'skull', 4);
end
