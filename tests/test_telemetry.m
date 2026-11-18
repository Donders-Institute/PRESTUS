classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_telemetry < matlab.unittest.TestCase
% TEST_TELEMETRY  Opt-in telemetry: consent handling, non-interactive path, no network without consent.
%
% Covers TELEMETRY_SETUP (recorded decision -> silent; SLURM/PBS job or
% `matlab -batch` session -> notice only, nothing recorded, no prompt;
% audit item B1), TELEMETRY_IS_NONINTERACTIVE, TELEMETRY_IS_HEADLESS_CLUSTER
% and TELEMETRY_RUN_TYPE (environment variables), TELEMETRY_ENV_OVERRIDE
% (PRESTUS_TELEMETRY off/on/invalid, precedence over a recorded decision),
% TELEMETRY_SET_CONSENT / TELEMETRY_SETUP_RESET, PREFDIR_PRESTUS,
% GENERATE_RUN_ID, and TRACK_USAGE (sends nothing and writes no ID unless
% opted in; when opted in, one POST whose payload carries no paths or
% subject IDs). The JVM's user.home is redirected to a temporary folder, so
% the real ~/.prestus is never read or written, and webwrite is shadowed by
% a recording stub, so no request leaves the machine. input and uiconfirm
% are shadowed by stubs that throw, so a prompt fails the test instead of
% blocking it. batchStartupOptionUsed is shadowed where a test depends on
% it, so results do not depend on whether the suite itself runs under
% `matlab -batch`.
%
% Not covered: answering the interactive input() prompt and the GUI
% uiconfirm dialog; both need a user.
%
%   results = runtests('tests/test_telemetry.m');

    properties
        Home        % redirected user.home
        StubDir     % folder holding the webwrite stub
        WebLog      % file the stub appends requests to
    end


    methods (TestMethodSetup)
        function isolate(tc)
            tc.Home = tempname;
            mkdir(tc.Home);
            tc.addTeardown(@() rmdir(tc.Home, 's'));
            orig_home = char(java.lang.System.getProperty('user.home'));
            java.lang.System.setProperty('user.home', tc.Home);
            tc.addTeardown(@() java.lang.System.setProperty('user.home', orig_home));
            tc.assertEqual(prefdir_prestus(), fullfile(tc.Home, '.prestus'), ...
                'user.home redirection failed; refusing to touch the real ~/.prestus');
            for v = {'SLURM_JOB_ID', 'PBS_JOBID', 'PRESTUS_RUN_TYPE', ...
                     'PRESTUS_TEST_WEBWRITE_LOG', 'PRESTUS_TELEMETRY', 'PRESTUS_TEST_BATCH'}
                tc.preserve_env(v{1});
            end
            setenv('SLURM_JOB_ID', '');
            setenv('PBS_JOBID', '');
            setenv('PRESTUS_TELEMETRY', '');     % the caller's setting must not leak in
            clear('telemetry_env_override');     % reset its warn-once memory
            tc.addTeardown(@() clear('telemetry_env_override'));
            tc.install_prompt_traps();
        end
    end

    methods (Test, TestTags = {'telemetry'})

        function test_prefdir_under_redirected_home(tc)
            d = prefdir_prestus();
            tc.verifyEqual(d, fullfile(tc.Home, '.prestus'));
            tc.verifyTrue(isfolder(d));
        end

        function test_recorded_decision_is_silent(tc)
            telemetry_set_consent_quiet(false);
            f = fullfile(prefdir_prestus(), 'telemetry.json');
            before = fileread(f);
            out = evalc('telemetry_setup()');
            tc.verifyEmpty(strtrim(out));
            tc.verifyEqual(fileread(f), before);
        end

        function test_slurm_job_prints_notice_records_nothing(tc)
            % B1: a SLURM job must neither prompt nor record a decision
            setenv('SLURM_JOB_ID', '12345');
            out = evalc('telemetry_setup()');
            tc.verifySubstring(out, 'Non-interactive session');
            tc.verifySubstring(out, 'no data sent');
            tc.verifyFalse(isfile(fullfile(prefdir_prestus(), 'telemetry.json')));
        end

        function test_pbs_job_prints_notice_records_nothing(tc)
            setenv('PBS_JOBID', '678.pbs');
            out = evalc('telemetry_setup()');
            tc.verifySubstring(out, 'Non-interactive session');
            tc.verifyFalse(isfile(fullfile(prefdir_prestus(), 'telemetry.json')));
        end

        function test_batch_session_prints_notice_never_prompts(tc)
            % matlab -batch (cron, CI, local script), with or without a
            % display: same as a SLURM job, never a prompt or dialog.
            tc.install_batch_stub();
            setenv('PRESTUS_TEST_BATCH', '1');
            out = evalc('telemetry_setup()');
            tc.verifySubstring(out, 'Non-interactive session');
            tc.verifySubstring(out, 'no data sent');
            tc.verifyFalse(isfile(fullfile(prefdir_prestus(), 'telemetry.json')));
        end

        function test_headless_cluster_detection(tc)
            tc.verifyFalse(telemetry_is_headless_cluster());
            setenv('SLURM_JOB_ID', '1');
            tc.verifyTrue(telemetry_is_headless_cluster());
            setenv('SLURM_JOB_ID', '');
            setenv('PBS_JOBID', '1');
            tc.verifyTrue(telemetry_is_headless_cluster());
        end

        function test_noninteractive_detection(tc)
            tc.install_batch_stub();
            setenv('PRESTUS_TEST_BATCH', '1');
            tc.verifyTrue(telemetry_is_noninteractive(), 'batch session');
            setenv('PRESTUS_TEST_BATCH', '0');
            setenv('SLURM_JOB_ID', '1');
            tc.verifyTrue(telemetry_is_noninteractive(), 'SLURM job, not batch');
            % Not batch, no scheduler: depends on desktop/display, not asserted.
        end

        function test_run_type_from_environment(tc)
            setenv('PRESTUS_RUN_TYPE', '');
            tc.verifyEqual(telemetry_run_type(), 'deployment');
            setenv('PRESTUS_RUN_TYPE', ' Technical ');
            tc.verifyEqual(telemetry_run_type(), 'technical');
            setenv('PRESTUS_RUN_TYPE', 'production');
            tc.verifyEqual(telemetry_run_type(), 'deployment');
        end

        function test_env_override_values(tc)
            vals = {'off', 'off'; '0', 'off'; 'FALSE', 'off'; ' Off ', 'off'; ...
                    'on', 'on'; '1', 'on'; 'True', 'on'; '', ''; '  ', ''};
            for k = 1:size(vals, 1)
                setenv('PRESTUS_TELEMETRY', vals{k, 1});
                tc.verifyEqual(telemetry_env_override(), vals{k, 2}, ...
                    sprintf('PRESTUS_TELEMETRY=''%s''', vals{k, 1}));
            end
        end

        function test_env_off_no_prompt_no_decision_file(tc)
            % Even in a session that would otherwise reach the interactive
            % path (stubbed non-batch; the input trap would throw).
            tc.install_batch_stub();
            setenv('PRESTUS_TEST_BATCH', '0');
            setenv('PRESTUS_TELEMETRY', 'off');
            out = evalc('telemetry_setup()');
            tc.verifySubstring(out, 'Telemetry off for this process');
            tc.verifyFalse(contains(out, 'What IS collected'), 'notice printed despite override');
            tc.verifyFalse(isfile(fullfile(prefdir_prestus(), 'telemetry.json')));
        end

        function test_env_off_overrides_recorded_opt_in(tc)
            tc.install_webwrite_stub();
            telemetry_set_consent_quiet(true);
            f = fullfile(prefdir_prestus(), 'telemetry.json');
            before = fileread(f);
            setenv('PRESTUS_TELEMETRY', '0');
            evalc('telemetry_setup()');
            track_usage('run_end', track_params(), struct('status', 'success'));
            tc.verifyFalse(isfile(tc.WebLog), 'webwrite was called with PRESTUS_TELEMETRY=0');
            tc.verifyFalse(isfile(fullfile(prefdir_prestus(), 'uuid.txt')));
            tc.verifyEqual(fileread(f), before, 'recorded decision was changed');
        end

        function test_env_on_sends_without_recording(tc)
            tc.install_webwrite_stub();
            tc.install_batch_stub();
            setenv('PRESTUS_TEST_BATCH', '0');
            setenv('PRESTUS_TELEMETRY', 'TRUE');
            out = evalc('telemetry_setup()');
            tc.verifySubstring(out, 'Telemetry on for this process');
            tc.verifyFalse(isfile(fullfile(prefdir_prestus(), 'telemetry.json')), ...
                'PRESTUS_TELEMETRY=on wrote a decision file');
            track_usage('run_end', track_params(), struct('status', 'success'));
            tc.assertTrue(isfile(tc.WebLog), 'PRESTUS_TELEMETRY=on sent nothing');
            tc.verifyNumElements(strsplit(strtrim(fileread(tc.WebLog)), newline), 2);
        end

        function test_env_on_never_overrides_recorded_opt_out(tc)
            % A recorded opt-out is the user's explicit choice; a process-wide
            % variable (e.g. set by a site module) must not turn sending back on.
            tc.install_webwrite_stub();
            telemetry_set_consent_quiet(false);
            f = fullfile(prefdir_prestus(), 'telemetry.json');
            before = fileread(f);
            setenv('PRESTUS_TELEMETRY', '1');
            track_usage('run_end', track_params(), struct('status', 'success'));
            tc.verifyFalse(isfile(tc.WebLog), 'PRESTUS_TELEMETRY=1 overrode a recorded opt-out');
            out = evalc('telemetry_setup()');
            tc.verifySubstring(out, 'PRESTUS_TELEMETRY=on ignored');
            tc.verifyEqual(fileread(f), before, 'recorded decision was changed');
        end

        function test_env_invalid_warns_and_is_ignored(tc)
            tc.install_webwrite_stub();
            telemetry_set_consent_quiet(false);
            setenv('PRESTUS_TELEMETRY', 'maybe');
            tc.verifyWarning(@() track_usage('run_end', track_params(), struct()), ...
                'prestus:telemetry:invalidEnv');
            tc.verifyFalse(isfile(tc.WebLog), 'invalid value overrode the recorded opt-out');
            tc.verifyWarningFree(@() track_usage('run_end', track_params(), struct()), ...
                'warning repeated for the same value');
            % With no recorded decision it behaves as unset: batch -> notice only.
            evalc('telemetry_setup_reset()');
            tc.install_batch_stub();
            setenv('PRESTUS_TEST_BATCH', '1');
            out = evalc('telemetry_setup()');
            tc.verifySubstring(out, 'Non-interactive session');
            tc.verifyFalse(isfile(fullfile(prefdir_prestus(), 'telemetry.json')));
        end

        function test_set_consent_and_reset(tc)
            f = fullfile(prefdir_prestus(), 'telemetry.json');
            telemetry_set_consent_quiet(true);
            cfg = jsondecode(fileread(f));
            tc.verifyTrue(cfg.opt_in);
            tc.verifyMatches(cfg.decided_on, '^\d{4}-\d{2}-\d{2}$');
            telemetry_set_consent_quiet(false);
            cfg = jsondecode(fileread(f));
            tc.verifyFalse(cfg.opt_in);
            evalc('telemetry_setup_reset()');
            tc.verifyFalse(isfile(f));
        end

        function test_run_id_is_uuid_v4(tc)
            a = generate_run_id();
            tc.verifyMatches(a, '^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$');
            tc.verifyNotEqual(a, generate_run_id());
        end

        function test_track_usage_sends_nothing_without_decision(tc)
            tc.install_webwrite_stub();
            track_usage('run_end', track_params(), struct('status', 'success'));
            tc.verifyFalse(isfile(tc.WebLog), 'webwrite was called without consent');
            tc.verifyFalse(isfile(fullfile(prefdir_prestus(), 'uuid.txt')));
        end

        function test_track_usage_sends_nothing_after_opt_out(tc)
            tc.install_webwrite_stub();
            telemetry_set_consent_quiet(false);
            track_usage('run_end', track_params(), struct('status', 'success'));
            tc.verifyFalse(isfile(tc.WebLog), 'webwrite was called after opt-out');
            tc.verifyFalse(isfile(fullfile(prefdir_prestus(), 'uuid.txt')));
        end

        function test_track_usage_payload_after_opt_in(tc)
            tc.install_webwrite_stub();
            telemetry_set_consent_quiet(true);
            setenv('PRESTUS_RUN_TYPE', 'technical');
            p = track_params();
            track_usage('run_end', p, struct('run_id', 'abc', 'status', 'success', 'duration_s', 1));
            tc.assertTrue(isfile(tc.WebLog), 'opted-in run sent nothing');
            lines = strsplit(strtrim(fileread(tc.WebLog)), newline);
            tc.verifyNumElements(lines, 2);                       % one request: URL + body
            payload = jsondecode(lines{2});
            tc.verifyEqual(payload.event, 'run_end');
            tc.verifyEqual(payload.run_type, 'technical');
            tc.verifyEqual(payload.status, 'success');
            tc.verifyFalse(contains(lines{2}, p.path.sim), 'payload contains a path');
            tc.verifyFalse(isfield(payload, 'subject_id'));
            tc.verifyTrue(isfile(fullfile(prefdir_prestus(), 'uuid.txt')));
        end

        function test_track_usage_sim_stats_allow_list(tc)
            tc.install_webwrite_stub();
            telemetry_set_consent_quiet(true);
            rec = struct('acoustic_dims', [240 71], 'acoustic_nt', 2871, 'acoustic_s', 12.5, ...
                'thermal_steps', 880, 'host', 'dccn-c091', 'timestamp', '2026-11-16T10:00:00Z', ...
                'hpc_partition', 'gpu40g', 'stage_preproc_s', 214.5, 'stage_note_s', 'text');
            track_usage('run_end', track_params(), struct('status', 'success', 'sim_stats', rec));
            lines = strsplit(strtrim(fileread(tc.WebLog)), newline);
            payload = jsondecode(lines{2});
            tc.verifyEqual(payload.sim_stats.acoustic_dims(:)', [240 71]);
            tc.verifyEqual(payload.sim_stats.acoustic_nt, 2871);
            tc.verifyEqual(payload.sim_stats.thermal_steps, 880);
            tc.verifyEqual(payload.sim_stats.stage_preproc_s, 214.5);
            tc.verifyFalse(contains(lines{2}, 'dccn-c091'), 'payload contains the host name');
            tc.verifyEqual(sort(fieldnames(payload.sim_stats))', ...
                sort({'acoustic_dims', 'acoustic_nt', 'acoustic_s', 'stage_preproc_s', 'thermal_steps'}));
        end

        function test_track_usage_skull_medians_need_current_notice(tc)
            % Skull properties: medians only, and only with consent to notice 2 or later.
            tc.install_webwrite_stub();
            rec = struct('acoustic_nt', 2871, 'skull_c_median', 2550, 'skull_c_p10', 2200, ...
                'skull_rho_median', 1850, 'skull_atten_median', 4.2, 'skull_path_mm', [6.5 7], ...
                'skull_path_c_median', 2600, 'skull_path_rho_median', 1900, 'skull_path_atten_median', 4.5);
            medians = {'skull_atten_median', 'skull_c_median', 'skull_path_atten_median', ...
                'skull_path_c_median', 'skull_path_rho_median', 'skull_rho_median'};

            telemetry_set_consent_quiet(true);               % records the current notice
            tc.verifyEqual(read_json(fullfile(prefdir_prestus(), 'telemetry.json')).notice, ...
                telemetry_notice_version());
            track_usage('run_end', track_params(), struct('status', 'success', 'sim_stats', rec));
            sent = last_payload(tc.WebLog).sim_stats;
            tc.verifyEqual(sort(fieldnames(sent))', sort([{'acoustic_nt'}, medians]));
            tc.verifyEqual(sent.skull_c_median, 2550);

            % consent recorded before notice 2 (no 'notice' field): no skull fields
            fid = fopen(fullfile(prefdir_prestus(), 'telemetry.json'), 'w');
            fprintf(fid, '{"opt_in":true,"decided_on":"2026-01-01","prestus_ver":"x"}\n');
            fclose(fid);
            track_usage('run_end', track_params(), struct('status', 'success', 'sim_stats', rec));
            sent = last_payload(tc.WebLog).sim_stats;
            tc.verifyEqual(fieldnames(sent)', {'acoustic_nt'});
        end

    end

    methods (Access = private)

        function preserve_env(tc, name)
            old = getenv(name);
            tc.addTeardown(@() setenv(name, old));
        end

        function install_webwrite_stub(tc)
            tc.StubDir = tc.install_stubs({'webwrite', { ...
                'function out = webwrite(url, data, varargin)', ...
                '% Test stub: record the request instead of sending it.', ...
                'fid = fopen(getenv(''PRESTUS_TEST_WEBWRITE_LOG''), ''a'');', ...
                'fprintf(fid, ''%s\n%s\n'', url, data);', ...
                'fclose(fid);', ...
                'out = [];', ...
                'end'}});
            tc.WebLog = fullfile(tc.StubDir, 'requests.log');
            setenv('PRESTUS_TEST_WEBWRITE_LOG', tc.WebLog);
            tc.assertEqual(which('webwrite'), fullfile(tc.StubDir, 'webwrite.m'), ...
                'webwrite stub is not first on the path; refusing to risk a real request');
        end

        function install_prompt_traps(tc)
            % input and uiconfirm throw, so a prompt fails fast instead of
            % blocking the test run.
            d = tc.install_stubs({ ...
                'input', {'function varargout = input(varargin)', ...
                    'error(''test:telemetry:prompted'', ''input() prompt reached'');', 'end'}; ...
                'uiconfirm', {'function varargout = uiconfirm(varargin)', ...
                    'error(''test:telemetry:prompted'', ''uiconfirm dialog reached'');', 'end'}});
            tc.assertEqual(which('input'), fullfile(d, 'input.m'), 'input trap not active');
        end

        function install_batch_stub(tc)
            % batchStartupOptionUsed returns PRESTUS_TEST_BATCH == '1'.
            d = tc.install_stubs({'batchStartupOptionUsed', { ...
                'function tf = batchStartupOptionUsed()', ...
                'tf = strcmp(getenv(''PRESTUS_TEST_BATCH''), ''1'');', 'end'}});
            tc.assertEqual(which('batchStartupOptionUsed'), ...
                fullfile(d, 'batchStartupOptionUsed.m'), 'batchStartupOptionUsed stub not active');
        end

        function d = install_stubs(tc, stubs)
            % stubs: n-by-2 cell {name, {lines}}; written to a fresh folder
            % that is first on the path for this test only.
            d = tempname;
            mkdir(d);
            tc.addTeardown(@() rmdir(d, 's'));
            for k = 1:size(stubs, 1)
                fid = fopen(fullfile(d, [stubs{k, 1} '.m']), 'w');
                fprintf(fid, '%s\n', stubs{k, 2}{:});
                fclose(fid);
            end
            % Shadowing a built-in (input, batchStartupOptionUsed) warns.
            tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'MATLAB:dispatcher:nameConflict'));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(d));
            rehash;
        end

    end

end

% =========================================================================

function telemetry_set_consent_quiet(tf)
    evalc('telemetry_set_consent(tf)');
end

function s = read_json(f)
    s = jsondecode(fileread(f));
end

function p = last_payload(weblog)
    lines = strsplit(strtrim(fileread(weblog)), newline);
    p = jsondecode(lines{end});
end

function p = track_params()
    p = struct();
    p.subject_id = 42;
    p.platform = 'matlab';
    p.simulation = struct('medium', 'layered', 'code_type', 'matlab_cpu', 'precision', 'single');
    p.path = struct('sim', '/private/study/sub-042/secret_path');
    p.transducer = struct('type', 'annular', 'freq_hz', 500e3, 'annular', struct('elem_n', 4));
    p.modules = struct('run_acoustic_sims', 1, 'run_heating_sims', 0);
    p.layers = struct('brain', 1, 'skull', 4);
end
