classdef test_suite_wiring < matlab.unittest.TestCase
% TEST_SUITE_WIRING  Every test file is run by RUN_ALL_TESTS or listed as standalone.
%
% Test files that no level of run_all_tests picks up are never run and go
% stale unnoticed (in 2026-10, five such files were found, one with failing
% tests). A new tests/test_*.m must be added to a level in run_all_tests.m,
% or to its list of standalone suites.
%
%   results = runtests('tests/test_suite_wiring.m');

    methods (Test, TestTags = {'wiring'})

        function test_every_test_file_is_referenced(tc)
            here = fileparts(mfilename('fullpath'));
            runner = fileread(fullfile(here, 'run_all_tests.m'));
            files = dir(fullfile(here, 'test_*.m'));
            missing = {};
            for i = 1:numel(files)
                if ~contains(runner, files(i).name)
                    missing{end+1} = files(i).name; %#ok<AGROW>
                end
            end
            tc.verifyEmpty(missing, sprintf(['Not run by run_all_tests and not listed as ' ...
                'standalone: %s'], strjoin(missing, ', ')));
        end

    end

end
