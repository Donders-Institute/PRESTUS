%% RUN_ALL_TESTS  Entry point for the PRESTUS test suite.
%
%   Usage (from the tests/ folder, or with it on the path):
%
%     run_all_tests              % unit tests only (default)
%     run_all_tests('unit')      % same as above
%     run_all_tests('water')     % unit + water pipeline (no MRI/SimNIBS needed)
%     run_all_tests('head')      % water + full placement + smoke_config + smoke_head
%     run_all_tests('acoustic')  % head + smoke_acoustic
%     run_all_tests('all')       % full suite including thermal
%     run_all_tests('physics')   % unit + free-water physics benchmarks (2D axisymmetric,
%                                % O'Neil / linearity / ISPPA targeting; minutes on CPU)
%     run_all_tests('physics_3d') % 3D O'Neil benchmark (GPU/HPC; PRESTUS_PHYSICS_CODE_TYPE)
%     run_all_tests('itrusst')   % ITRUSST BM1 benchmark vs. reference fields
%                                % (needs PRESTUS_ITRUSST_DATA; 3D, GPU/HPC)
%     run_all_tests('slurm_live') % submit a real SLURM job and check hpc_job_report output
%                                 % (skipped automatically when sbatch is not on PATH)
%
%   Demo-data levels (require PRESTUS_DEMO_DATA to point to prestus_testdata/):
%     run_all_tests('demo_inputs')   % verify all demo input files exist
%     run_all_tests('demo_config')   % load all three placement-mode configs
%     run_all_tests('demo_localite') % parse Localite XML -> voxel positions
%     run_all_tests('demo_pct')      % verify pCT data is present/readable
%     run_all_tests('demo_head')     % head preprocessing on real sub-009 data
%     run_all_tests('demo_acoustic') % full acoustic pipeline on sub-009
%     run_all_tests('demo_thermal')  % thermal pipeline on cached sub-009 acoustics
%     run_all_tests('demo_golden')   % compare demo key metrics with the stored reference
%
%   Environment variables (for integration levels >= 'head'):
%     PRESTUS_TEST_DATA   — path to folder with demo subject data & m2m_* output
%     PRESTUS_DEMO_CONFIG — path to study config YAML  (default: config_tutorial.yaml)
%
%   Environment variables (for demo levels):
%     PRESTUS_DEMO_DATA   — path to prestus_testdata root folder
%
%   Tests are selected by level tag (PRESTUS_TEST_LEVELS): every test class in
%   tests/ is found automatically. Failures print the full test diagnostics.
%   The function errors on any failure, making it suitable for CI:
%     matlab -batch "run_all_tests"
%
%   results = run_all_tests(...) also returns the TestResult array.
%
%   Resources: every test and test class is timed and its memory recorded
%   (PRESTUS_TEST_RESOURCES). The table is written to
%   <dir>/prestus_test_resources_<level>_<SLURM job id or time stamp>.csv, where
%   <dir> is PRESTUS_TEST_RESOURCES_DIR (default: tempdir); the run prints
%   the path and the slowest and most memory-hungry tests.
%
% See also: PRESTUS_TEST_LEVELS, PRESTUS_PATHS, TEST_SUITE_WIRING

function results = run_all_tests(level)

    arguments
        level (1,:) char = 'unit'
    end

    import matlab.unittest.TestRunner
    import matlab.unittest.TestSuite
    import matlab.unittest.selectors.HasTag
    import matlab.unittest.plugins.DiagnosticsRecordingPlugin

    clear functions %#ok<CLFUNC>

    orig_run_type = getenv('PRESTUS_RUN_TYPE');
    setenv('PRESTUS_RUN_TYPE', 'technical');
    restore_run_type = onCleanup(@() setenv('PRESTUS_RUN_TYPE', orig_run_type));

    levels = prestus_test_levels();
    level = lower(level);
    if ~isKey(levels, level)
        error('run_all_tests:unknownLevel', 'Unknown level ''%s''. Use one of: %s', ...
              level, strjoin(levels.keys, ' | '));
    end

    tags = levels(level);
    selector = HasTag(tags{1});
    for k = 2:numel(tags)
        selector = selector | HasTag(tags{k});
    end
    suite = TestSuite.fromFolder(fileparts(mfilename('fullpath')), selector);

    prev_warn = warning('off', 'prestus:toolboxDefault');
    warning('off', 'prestus:noSimNIBS');
    restore_warn = onCleanup(@() warning(prev_warn));

    runner = TestRunner.withNoPlugins();
    runner.addPlugin(DiagnosticsRecordingPlugin());
    resources = prestus_test_resources();
    runner.addPlugin(resources);
    results = runner.run(suite);

    print_summary(results, level);
    write_resources(resources.Rows, results, level);

    n_fail = sum([results.Failed]);
    if n_fail > 0
        error('run_all_tests:failures', '%d test(s) failed.', n_fail);
    end
    if nargout == 0
        clear results
    end
end

function print_summary(results, level)
    fprintf('\n========================================\n');
    fprintf('  PRESTUS test suite — level: %s\n', upper(level));
    fprintf('========================================\n');
    for r = reshape(results, 1, [])
        if r.Passed
            fprintf('  [PASS]  %s\n', r.Name);
        elseif r.Failed
            fprintf('  [FAIL]  %s\n', r.Name);
            for d = r.Details.DiagnosticRecord
                if ~strcmp(d.Event, 'DiagnosticLogged')
                    fprintf('%s\n', indent(d.Report));
                end
            end
        else
            fprintf('  [SKIP]  %s\n', r.Name);
            reason = r.Details.DiagnosticRecord;
            reason = reason(strcmp({reason.Event}, 'AssumptionFailed'));
            if ~isempty(reason)
                fprintf('          -> %s\n', first_line(reason(1)));
            end
        end
    end
    fprintf('----------------------------------------\n');
    fprintf('  Passed: %d   Failed: %d   Skipped: %d\n', sum([results.Passed]), ...
            sum([results.Failed]), sum([results.Incomplete] & ~[results.Failed]));
    fprintf('========================================\n\n');
end

function write_resources(rows, results, level)
% Join the outcome to each test row, write the CSV, print the heaviest tests
    if isempty(rows); return; end
    T = struct2table(rows, 'AsArray', true);
    T.kind = string(T.kind); T.name = string(T.name);
    T.outcome = repmat("", height(T), 1);
    [is_test, at] = ismember(T.name, string({results.Name}));
    for k = find(is_test)'
        r = results(at(k));
        if r.Failed
            T.outcome(k) = "fail";
        elseif r.Incomplete
            T.outcome(k) = "skip";
        else
            T.outcome(k) = "pass";
        end
    end
    T.cores_busy = T.cpu_s ./ max(T.wall_s, eps);

    out_dir = getenv('PRESTUS_TEST_RESOURCES_DIR');
    if isempty(out_dir); out_dir = tempdir; end
    run_id = getenv('SLURM_JOB_ID');
    if isempty(run_id); run_id = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')); end
    f = fullfile(out_dir, sprintf('prestus_test_resources_%s_%s.csv', level, run_id));
    try
        if ~isfolder(out_dir); mkdir(out_dir); end
        writetable(T, f);
    catch err
        warning('run_all_tests:resources', 'Resource table not written (%s): %s', f, err.message);
        f = '';
    end

    tests = T(T.kind == "test", :);
    fprintf('Resources (%d tests): %.1f min wall, %.1f min CPU, MATLAB peak %s MB\n', ...
        height(tests), sum(tests.wall_s) / 60, sum(tests.cpu_s) / 60, num2str(max(T.peak_mb), '%.0f'));
    print_top(tests, 'wall_s', 'Slowest tests', '%7.1f s');
    if any(isfinite(tests.peak_delta_mb))
        print_top(tests, 'peak_delta_mb', 'Most memory above the test''s start', '%7.0f MB');
    end
    if ~isempty(f); fprintf('Resource table: %s\n\n', f); end
end

function print_top(T, col, title, fmt)
    [~, order] = sort(T.(col), 'descend', 'MissingPlacement', 'last');
    fprintf('%s:\n', title);
    for k = order(1:min(10, numel(order)))'
        fprintf(['  ' fmt '  %s\n'], T.(col)(k), T.name(k));
    end
end

function s = indent(txt)
    s = regexprep(strtrim(txt), '(^|\n)', '$1          ');
end

function s = first_line(rec)
    % The user-supplied assumption message, else the first report line
    txt = {rec.TestDiagnosticResults.DiagnosticText};
    lines = strtrim(splitlines(string(strjoin(txt, newline))));
    if isempty(lines) || all(lines == "")
        lines = strtrim(splitlines(string(rec.Report)));
    end
    s = char(lines(find(lines ~= "", 1)));
end
