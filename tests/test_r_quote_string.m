classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_r_quote_string < matlab.unittest.TestCase
% TEST_R_QUOTE_STRING  Unit tests for r_quote_string (R string-literal injection defense).
%
% Context: position_transducer_tus_entry.m generates a small R script by
% sprintf-substituting mask/target paths into readNIfTI("%s")-style R
% string literals, then runs it with Rscript. A path containing a literal
% " closes that R string early, so anything after it is parsed as R code
% — e.g. a path ending in
%   foo.nii.gz"); system("touch PWNED"); x <- readNIfTI("
% would have R itself execute an arbitrary shell command. r_quote_string
% escapes \ and " the way R does, so the substituted value cannot break
% out of its string literal.
%
% Where Rscript is available on PATH, these tests actually run the
% generated R source and check for the sentinel side effect, proving the
% fix end-to-end rather than just checking the escaping algorithm.
%
% Run with:   results = runtests('tests/test_r_quote_string.m');

    properties
        TmpDir
        HaveR
        RscriptBin
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
            [status, out] = system('which Rscript 2>/dev/null');
            tc.HaveR = (status == 0) && ~isempty(strtrim(out));
            if tc.HaveR
                tc.RscriptBin = strtrim(out);
            end
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

        function test_escapes_double_quote(tc)
            tc.verifyEqual(r_quote_string('a"b'), 'a\"b');
        end

        function test_escapes_backslash(tc)
            tc.verifyEqual(r_quote_string('a\b'), 'a\\b');
        end

        function test_escapes_backslash_before_quote(tc)
            % Backslash must be escaped first, or a payload ending in \"
            % would produce \\" (an escaped backslash followed by an
            % unescaped, string-terminating quote) instead of \\\".
            tc.verifyEqual(r_quote_string('a\"b'), 'a\\\"b');
        end

        function test_plain_string_unchanged(tc)
            tc.verifyEqual(r_quote_string('plain/path.nii.gz'), 'plain/path.nii.gz');
        end

        function test_readNIfTI_breakout_payload_stays_inert_via_Rscript(tc)
            if ~tc.HaveR
                tc.assumeFail('Rscript not found on PATH; skipping end-to-end R-injection test.');
            end
            sentinel = fullfile(tc.TmpDir, 'PWNED_R');
            payload = sprintf('x.nii.gz"); system("touch %s"); y<-readNIfTI("', sentinel);

            src = sprintf('z <- readNIfTI("%s")\ncat("ran\\n")\n', r_quote_string(payload));
            script = fullfile(tc.TmpDir, 'poc.R');
            fid = fopen(script, 'w');
            fprintf(fid, '%s', src);
            fclose(fid);

            % readNIfTI will fail (the "file" does not exist) but that is
            % expected and irrelevant here — R must fail while parsing the
            % quoted string as ONE argument, never having executed system().
            system(sprintf('%s --vanilla %s', shell_quote_arg(tc.RscriptBin), shell_quote_arg(script)));

            tc.verifyFalse(isfile(sentinel), ...
                'a " inside a mask path must not break out of the R string literal and run system()');
        end

    end
end
