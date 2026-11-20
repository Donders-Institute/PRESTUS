classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_run_command < matlab.unittest.TestCase
% TEST_RUN_COMMAND  Unit tests for run_command.m / run_shell_script.m.
%
% Context: run_command() replaces system(cmd_string) with
% java.lang.ProcessBuilder(argv) for the "invoke a program with these
% arguments" call sites in PRESTUS, removing the shell entirely rather
% than defending against it (see shell_quote_arg.m for the injection
% class this eliminates and functions/helper/run_command.m's header for
% full background).
%
% These tests exercise the real JVM/OS process machinery (no mocking) to
% prove the [status, cmdout] contract matches system()'s, and that shell
% metacharacters in an argument are structurally inert — there is no
% shell in the loop for them to reach.
%
% Run with:   results = runtests('tests/test_run_command.m');

    properties
        TmpDir = ''   % stays '' when setup skips the test (Windows)
    end

    methods (TestMethodSetup)
        function setup(tc)
            if ispc
                tc.assumeFail('run_command POSIX-shell-free tests assume a Unix toolset (echo/sh/cat/false).');
            end
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
        end
    end

    methods (TestMethodTeardown)
        function cleanup(tc)
            if isfolder(tc.TmpDir)
                rmdir(tc.TmpDir, 's');
            end
        end
    end

    methods (Test, TestTags = {'security'})

        function test_basic_argv_and_status(tc)
            [status, out] = run_command({'echo', 'hello world'});
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), 'hello world');
        end

        function test_argument_with_spaces_is_one_argument(tc)
            % A single argv element containing spaces must reach the
            % program as ONE argument, not be word-split the way a shell
            % would split an unquoted "%s" substitution.
            [status, out] = run_command({'echo', 'one arg with spaces'});
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), 'one arg with spaces');
        end

        function test_command_substitution_payload_is_inert(tc)
            % The canonical BabelBrain-report payload shape: with a shell
            % in the loop this would run `touch <marker>`; with no shell,
            % it is just argument text.
            marker = fullfile(tc.TmpDir, 'PWNED_cmdsub');
            payload = sprintf('$(touch %s)', marker);
            [status, out] = run_command({'echo', payload});
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), payload);
            tc.verifyFalse(isfile(marker), ...
                '$(...) in an argv element must not execute without a shell');
        end

        function test_backtick_payload_is_inert(tc)
            marker = fullfile(tc.TmpDir, 'PWNED_backtick');
            payload = sprintf('`touch %s`', marker);
            [status, out] = run_command({'echo', payload});
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), payload);
            tc.verifyFalse(isfile(marker));
        end

        function test_semicolon_pipe_ampersand_payload_is_inert(tc)
            marker = fullfile(tc.TmpDir, 'PWNED_chain');
            payload = sprintf('x; touch %s; y | touch %s && touch %s', marker, marker, marker);
            [status, out] = run_command({'echo', payload});
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), payload);
            tc.verifyFalse(isfile(marker));
        end

        function test_nonzero_exit_status(tc)
            [status, ~] = run_command({'false'});
            tc.verifyNotEqual(status, 0);
        end

        function test_combined_stdout_stderr(tc)
            [status, out] = run_command({'sh', '-c', 'echo out1; echo err1 1>&2'});
            tc.verifyEqual(status, 0);
            tc.verifyTrue(contains(out, 'out1'));
            tc.verifyTrue(contains(out, 'err1'));
        end

        function test_missing_binary_returns_status_127(tc)
            [status, out] = run_command({'this_binary_does_not_exist_prestus_test'});
            tc.verifyEqual(status, 127);
            tc.verifyTrue(contains(out, 'this_binary_does_not_exist_prestus_test'));
        end

        function test_working_dir_option(tc)
            [status, out] = run_command({'pwd'}, 'WorkingDir', tc.TmpDir);
            tc.verifyEqual(status, 0);
            % Resolve symlinks (e.g. macOS /tmp -> /private/tmp) on both
            % sides before comparing.
            tc.verifyEqual(realpath(strtrim(out)), realpath(tc.TmpDir));
        end

        function test_env_option(tc)
            [status, out] = run_command({'sh', '-c', 'echo $PRESTUS_TEST_VAR'}, ...
                'Env', struct('PRESTUS_TEST_VAR', 'abc123'));
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), 'abc123');
        end

        function test_stdin_option(tc)
            [status, out] = run_command({'cat'}, 'Stdin', 'piped text');
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), 'piped text');
        end

        function test_sees_live_setenv_path_update(tc)
            % Regression test: java.lang.ProcessBuilder's default inherited
            % environment is a snapshot taken at JVM startup and does NOT
            % see later MATLAB setenv() changes (a real JVM caching quirk,
            % confirmed against this MATLAB build). PRESTUS relies on
            % setenv('PATH', ...) both in tests (mocking a binary) and in
            % real HPC/tool setup, so run_command must resolve PATH
            % lookups against the CURRENT MATLAB environment, not a stale
            % one from when the JVM first started.
            orig_path = getenv('PATH');
            cleanup = onCleanup(@() setenv('PATH', orig_path)); %#ok<NASGU>

            bin_dir = fullfile(tc.TmpDir, 'bin');
            mkdir(bin_dir);
            stub = fullfile(bin_dir, 'prestus_test_stub_tool');
            fid = fopen(stub, 'w');
            fprintf(fid, '#!/bin/sh\necho STUB_RAN\n');
            fclose(fid);
            run_command({'chmod', '+x', stub});

            setenv('PATH', [bin_dir pathsep orig_path]);

            [status, out] = run_command({'prestus_test_stub_tool'});
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), 'STUB_RAN');
        end

        function test_background_option_returns_immediately(tc)
            marker = fullfile(tc.TmpDir, 'bg_done.txt');
            t0 = tic;
            [status, out] = run_command({'sh', '-c', sprintf('sleep 1 && touch %s', marker)}, ...
                'Background', true);
            elapsed = toc(t0);
            tc.verifyLessThan(elapsed, 0.9, 'Background=true must not block for the child''s sleep 1');
            tc.verifyTrue(isnan(status));
            tc.verifyEqual(out, '');
            tc.verifyFalse(isfile(marker), 'child should not have finished yet');

            t1 = tic;
            while toc(t1) < 5 && ~isfile(marker)
                pause(0.1);
            end
            tc.verifyTrue(isfile(marker), 'backgrounded process should eventually complete on its own');
        end

        function test_run_shell_script_supports_shell_semantics(tc)
            % Module-load-style chaining genuinely needs a real shell.
            [status, out] = run_shell_script('export FOO=bar && echo $FOO');
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), 'bar');
        end

        function test_run_shell_script_applies_shell_exactly_once(tc)
            % A value already shell_quote_arg()-quoted inside the script
            % text must survive intact -- run_shell_script must not
            % re-tokenize/re-quote the script through an extra outer
            % shell layer (which is exactly what system(script) would do).
            marker = fullfile(tc.TmpDir, 'PWNED_rss');
            payload = sprintf('$(touch %s)', marker);
            q = shell_quote_arg(payload);
            [status, out] = run_shell_script(sprintf('echo %s', q));
            tc.verifyEqual(status, 0);
            tc.verifyEqual(strtrim(out), payload);
            tc.verifyFalse(isfile(marker));
        end

    end
end

function p = realpath(p)
    p = char(java.io.File(p).getCanonicalPath());
end
