classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_suite_wiring < matlab.unittest.TestCase
% TEST_SUITE_WIRING  Every test carries exactly one level tag of RUN_ALL_TESTS.
%
% RUN_ALL_TESTS selects tests by level tag (PRESTUS_TEST_LEVELS). A test
% without a level tag is never run and goes stale unnoticed; a test with two
% runs at the wrong level. Tag a new test class with its level, e.g.
%
%   classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_x < matlab.unittest.TestCase
%
% or tag each methods block when one file spans several levels.
%
%   results = runtests('tests/test_suite_wiring.m');

    methods (Test, TestTags = {'wiring'})

        function test_every_test_has_one_level_tag(tc)
            here = fileparts(mfilename('fullpath'));
            levels = prestus_test_levels();
            tag_lists = values(levels);
            level_tags = unique([tag_lists{:}, {'standalone'}]);
            suite = matlab.unittest.TestSuite.fromFolder(here);
            n = arrayfun(@(t) numel(intersect(t.Tags, level_tags)), suite);
            bad = unique(regexprep({suite(n ~= 1).Name}, '\(.*\)$', ''));
            tc.verifyEmpty(bad, sprintf('Tests without exactly one level tag (%s): %s', ...
                strjoin(level_tags, ', '), strjoin(bad, ', ')));
        end

        function test_every_test_file_is_a_test_class(tc)
            % fromFolder silently skips a test_*.m that is not a TestCase class
            here = fileparts(mfilename('fullpath'));
            files = dir(fullfile(here, 'test_*.m'));
            suite = matlab.unittest.TestSuite.fromFolder(here);
            classes = unique({suite.TestParentName});
            tc.verifyEmpty(setdiff(erase({files.name}, '.m'), classes), ...
                'test_*.m files that hold no tests');
        end

    end

end
