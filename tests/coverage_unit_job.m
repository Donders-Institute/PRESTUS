function results = coverage_unit_job()
% COVERAGE_UNIT_JOB  Run the unit level with line coverage of functions/ (Cobertura XML).
%
% Runs the tests RUN_ALL_TESTS('unit') runs, with a CodeCoveragePlugin over
% functions/ (subfolders included), and writes tests/coverage/unit_coverage.xml
% (Cobertura format; tests/coverage/ is gitignored). tests/docs/gen_test_map.py
% then adds a line-coverage column to its local copy of the test map
% (tests/coverage/tests_function_map_coverage.md).
%
% It is a function, not a script, so TestSuite.fromFolder does not take it for
% a script-based test. It errors after writing the XML when a test fails.
%
% Use as (cluster: in an sbatch job, never on the login node):
%   cd <PRESTUS>/tests; matlab -batch coverage_unit_job
%   python ../tests/docs/gen_test_map.py
%
% Output:
%   results - TestResult array
%
% See also: RUN_ALL_TESTS, PRESTUS_TEST_LEVELS

    import matlab.unittest.TestRunner
    import matlab.unittest.TestSuite
    import matlab.unittest.selectors.HasTag
    import matlab.unittest.plugins.CodeCoveragePlugin
    import matlab.unittest.plugins.codecoverage.CoberturaFormat

    here = fileparts(mfilename('fullpath'));
    root = fileparts(here);
    out_dir = fullfile(here, 'coverage');
    if ~isfolder(out_dir); mkdir(out_dir); end
    xml_file = fullfile(out_dir, 'unit_coverage.xml');

    % same run type as run_all_tests: no usage telemetry from test runs
    orig_run_type = getenv('PRESTUS_RUN_TYPE');
    setenv('PRESTUS_RUN_TYPE', 'technical');
    restore_run_type = onCleanup(@() setenv('PRESTUS_RUN_TYPE', orig_run_type)); %#ok<NASGU>

    levels = prestus_test_levels();
    tags = levels('unit');
    selector = HasTag(tags{1});
    for k = 2:numel(tags)
        selector = selector | HasTag(tags{k});
    end
    suite = TestSuite.fromFolder(here, selector);
    fprintf('coverage_unit_job: %d tests, coverage over %s\n', numel(suite), fullfile(root, 'functions'));

    runner = TestRunner.withTextOutput;
    runner.addPlugin(CodeCoveragePlugin.forFolder(fullfile(root, 'functions'), ...
        'IncludingSubfolders', true, 'Producing', CoberturaFormat(xml_file)));
    results = runner.run(suite);

    fprintf('coverage_unit_job: wrote %s\n', xml_file);
    fprintf('coverage_unit_job: %d passed, %d failed, %d incomplete\n', ...
        nnz([results.Passed]), nnz([results.Failed]), nnz([results.Incomplete]));
    if any([results.Failed])
        error('coverage_unit_job:failed', '%d unit test(s) failed; coverage was still written.', ...
            nnz([results.Failed]));
    end
end
