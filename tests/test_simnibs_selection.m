classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_simnibs_selection < matlab.unittest.TestCase
% TEST_SIMNIBS_SELECTION  SimNIBS install version and existing-segmentation checks.
%
%   simnibs_install_info reads the version and charm.ini of a fake install
%   (conda layout and installer wrapper); segmentation_check_existing
%   compares a fake m2m folder with it. No SimNIBS installation is run.
%   Mirrors nexus/tests/test_simnibs.py.
%
%   Run with:   results = runtests('tests/test_simnibs_selection.m');

    properties
        TmpDir
    end

    properties (Constant)
        Ini = sprintf('# standard settings\n[preprocess]\ndenoise = false\n\n[samseg]\nthreads = 8\n')
    end

    methods (TestMethodSetup)
        function setup(tc)
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

    methods (Access = private)
        function bin = conda_env(tc, name, version)
            env = fullfile(tc.TmpDir, name);
            site = fullfile(env, 'lib', 'python3.11', 'site-packages');
            mkdir(fullfile(site, 'simnibs'));
            mkdir(fullfile(site, ['simnibs-' version '.dist-info']));
            write_text(fullfile(site, 'simnibs', 'charm.ini'), tc.Ini);
            bin = fullfile(env, 'bin');
            mkdir(bin);
            write_text(fullfile(bin, 'charm'), sprintf('#!/usr/bin/env python\nfrom simnibs.cli.charm import main\n'));
        end

        function m2m = m2m_folder(tc, version, settings)
            m2m = fullfile(tc.TmpDir, 'seg', 'm2m_sub-001');
            mkdir(m2m);
            write_text(fullfile(m2m, 'final_tissues.nii.gz'), '');
            write_text(fullfile(m2m, 'charm_log.html'), sprintf('INFO: simnibs version %s\n', version));
            write_text(fullfile(m2m, 'settings.ini'), settings);
        end
    end

    methods (Test)
        function conda_install_version(tc)
            info = simnibs_install_info(tc.conda_env('simnibs_v4.6.0', '4.6.0'));
            tc.verifyEqual(info.version, '4.6.0');
            tc.verifyEqual(info.major, 4);
            tc.verifyTrue(endsWith(info.default_settings, 'charm.ini'));
        end

        function installer_wrapper_version(tc)
            root = fullfile(tc.TmpDir, 'SimNIBS-4.0');
            site = fullfile(root, 'simnibs_env', 'lib', 'python3.9', 'site-packages');
            mkdir(fullfile(site, 'simnibs', 'cli'));
            mkdir(fullfile(site, 'simnibs-4.0.0.dist-info'));
            mkdir(fullfile(root, 'bin'));
            write_text(fullfile(root, 'bin', 'charm'), sprintf( ...
                '#! /bin/bash -e\n"%s/simnibs_env/bin/python" -E -u "%s/simnibs/cli/charm.py" "$@"\n', root, site));
            info = simnibs_install_info(fullfile(root, 'bin'));
            tc.verifyEqual(info.version, '4.0.0');
        end

        function linked_paths_follow_the_install(tc)
            % hpc.ld_library_path and the PlanTUS Python follow startup.simnibs_bin_path while empty
            bin = tc.conda_env('simnibs_v4.6.0', '4.6.0');
            env = fileparts(bin);
            py = fullfile(env, 'bin', 'python3.11');
            write_text(py, '');
            write_text(fullfile(bin, 'charm'), sprintf('#!%s\nfrom simnibs.cli.charm import main\n', py));
            lib = fullfile(env, 'lib', 'python3.11', 'site-packages', 'simnibs', 'external', 'lib', 'linux');
            mkdir(lib);
            info = simnibs_install_info(bin);
            tc.verifyEqual(info.python, py);
            tc.verifyEqual(info.env_root, env);
            tc.verifyEqual(info.ld_library_path, lib);

            p = struct('startup', struct('simnibs_bin_path', bin), 'hpc', struct('ld_library_path', ''));
            if ispc
                tc.verifyEmpty(simnibs_ld_library_path(p));
                return
            end
            tc.verifyEqual(simnibs_ld_library_path(p), lib);
            old = getenv('LD_LIBRARY_PATH');
            tc.addTeardown(@() setenv('LD_LIBRARY_PATH', old));
            setenv('LD_LIBRARY_PATH', '/usr/lib/x');
            env_over = simnibs_env_overrides(p);
            tc.verifyEqual(env_over.LD_LIBRARY_PATH, [lib ':/usr/lib/x']);   % prefixed, not replaced
            p.hpc.ld_library_path = '/opt/gcc/7.2.0/lib64';                 % a set value wins
            tc.verifyEqual(simnibs_ld_library_path(p), '/opt/gcc/7.2.0/lib64');
        end

        function no_charm_is_unknown(tc)
            info = simnibs_install_info(tc.TmpDir);
            tc.verifyEmpty(info.version);
            tc.verifyTrue(isnan(info.major));
        end

        function matching_segmentation_has_no_differences(tc)
            bin = tc.conda_env('simnibs_v4.6.0', '4.6.0');
            % comments and spacing do not count
            m2m = tc.m2m_folder('4.6.0', sprintf('[preprocess]\ndenoise=false  # off\n[samseg]\n threads =  8\n'));
            p = struct('startup', struct('simnibs_bin_path', bin));
            tc.verifyEmpty(segmentation_check_existing(m2m, p));
        end

        function differing_version_and_settings(tc)
            bin = tc.conda_env('simnibs_v4.6.0', '4.6.0');
            m2m = tc.m2m_folder('4.0.0', strrep(tc.Ini, 'false', 'true'));
            p = struct('startup', struct('simnibs_bin_path', bin));
            diffs = segmentation_check_existing(m2m, p);
            tc.verifyNumElements(diffs, 2);
            tc.verifyTrue(contains(diffs{1}, '4.0.0'));
            tc.verifyTrue(contains(diffs{2}, 'settings'));

            % the same settings in segmentation.settings_file: only the version differs
            custom = fullfile(tc.TmpDir, 'custom.ini');
            write_text(custom, strrep(tc.Ini, 'false', 'true'));
            p.segmentation = struct('settings_file', custom);
            diffs = segmentation_check_existing(m2m, p);
            tc.verifyNumElements(diffs, 1);
        end
    end
end


function write_text(file, text)
    fid = fopen(file, 'w');
    fprintf(fid, '%s', text);
    fclose(fid);
end
