classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_install_check < matlab.unittest.TestCase
% TEST_INSTALL_CHECK  Installation report and Nexus Python discovery.
%
% Covers PRESTUS_CHECK_INSTALL (report fields, required parts of this
% checkout), CONDA_ENV_DIRS and FIND_NEXUS_PYTHON (an unusable
% $PRESTUS_NEXUS_PYTHON is tried and skipped). No data, no Nexus needed.
%
%   results = runtests('tests/test_install_check.m');

    methods (Test)

        function test_report_fields(tc)
            c = prestus_check_install(Quiet=true, Nexus=false);
            tc.verifyNotEmpty(c);
            tc.verifyEqual(sort(fieldnames(c))', ...
                sort({'name', 'required', 'ok', 'detail', 'fix', 'submodule'}));
            tc.verifyTrue(all(cellfun(@isempty, {c([c.ok]).fix})), ...
                'a check that passed has no fix');
        end

        function test_required_parts_of_this_checkout(tc)
            % the test run itself needs PRESTUS on the path and k-Wave
            c = prestus_check_install(Quiet=true, Nexus=false);
            names = {c.name};
            tc.verifyTrue(c(strcmp(names, 'PRESTUS on path')).ok);
            tc.verifyTrue(c(strcmp(names, 'k-Wave')).ok);
            tc.verifyTrue(all(cellfun(@isempty, {c.submodule})));
        end

        function test_conda_env_dirs_exist(tc)
            d = conda_env_dirs();
            tc.verifyClass(d, 'cell');
            tc.verifyTrue(all(cellfun(@isfolder, d)));
        end

        function test_unusable_python_is_skipped(tc)
            old = getenv('PRESTUS_NEXUS_PYTHON');
            tc.addTeardown(@() setenv('PRESTUS_NEXUS_PYTHON', old));
            setenv('PRESTUS_NEXUS_PYTHON', 'prestus-no-such-python');
            [python, tried] = find_nexus_python();
            tc.verifyNotEqual(python, 'prestus-no-such-python');
            tc.verifyEqual(tried(1).python, 'prestus-no-such-python');
            tc.verifyFalse(tried(1).ok);
        end
    end
end
