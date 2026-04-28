%% RUN_ALL_TESTS  Entry point for the PRESTUS test suite.
%
%   Usage (from repo root or tests/ folder):
%
%     run_all_tests            % unit tests only (default)
%     run_all_tests('unit')    % same as above
%     run_all_tests('water')   % unit + water pipeline (no MRI/SimNIBS needed)
%     run_all_tests('head')    % water + smoke_config + smoke_head
%     run_all_tests('acoustic')% head + smoke_acoustic
%     run_all_tests('all')     % full suite including thermal
%
%   Environment variables (for integration levels >= 'head'):
%     PRESTUS_TEST_DATA   — path to folder with demo subject data & m2m_* output
%     PRESTUS_DEMO_CONFIG — path to study config YAML  (default: tutorial_config.yaml)
%
%   The script exits with a non-zero status on any failure, making it
%   suitable for CI:
%     matlab -batch "run_all_tests"

function run_all_tests(level)

    if nargin < 1
        level = 'unit';
    end

    % ---- Resolve paths ------------------------------------------------
    here      = fileparts(mfilename('fullpath'));
    repo_root = fileparts(here);
    addpath(genpath(fullfile(repo_root, 'functions')));
    addpath(fullfile(repo_root, 'configs'));
    addpath(fullfile(here, 'fixtures'));

    % ---- Select test files and tags by level --------------------------
    unit_files = {
        fullfile(here, 'test_helper.m')
        fullfile(here, 'test_thermal_parameters.m')
        fullfile(here, 'test_transform.m')
        fullfile(here, 'test_load_parameters.m')
        fullfile(here, 'test_head_preprocessing.m')
    };

    integration_file  = fullfile(here, 'test_integration_pipeline.m');
    water_test_file   = fullfile(here, 'test_integration_water.m');
    isppa_test_file   = fullfile(here, 'test_isppa_rescaling.m');
    placement_test_file = fullfile(here, 'test_transducer_placement.m');

    switch lower(level)
        case 'unit'
            suites = [testsuite(unit_files), ...
                      testsuite(placement_test_file, 'Tag', 'manual')];
        case 'water'
            suites = [testsuite(unit_files), ...
                      testsuite(placement_test_file, 'Tag', 'manual'), ...
                      testsuite(water_test_file), ...
                      testsuite(isppa_test_file)];
        case 'head'
            suites = [testsuite(unit_files), ...
                      testsuite(placement_test_file), ...
                      testsuite(water_test_file), ...
                      testsuite(isppa_test_file), ...
                      testsuite(integration_file, 'Tag', 'smoke_config'), ...
                      testsuite(integration_file, 'Tag', 'smoke_head')];
        case 'acoustic'
            suites = [testsuite(unit_files), ...
                      testsuite(placement_test_file), ...
                      testsuite(water_test_file), ...
                      testsuite(isppa_test_file), ...
                      testsuite(integration_file, 'Tag', 'smoke_config'), ...
                      testsuite(integration_file, 'Tag', 'smoke_head'), ...
                      testsuite(integration_file, 'Tag', 'smoke_acoustic')];
        case 'all'
            suites = [testsuite(unit_files), ...
                      testsuite(placement_test_file), ...
                      testsuite(water_test_file), ...
                      testsuite(isppa_test_file), ...
                      testsuite(integration_file)];
        otherwise
            error('Unknown level ''%s''. Use: unit | head | acoustic | all', level);
    end

    % ---- Run ----------------------------------------------------------
    runner  = matlab.unittest.TestRunner.withTextOutput();
    results = runner.run(suites);

    % ---- Summary ------------------------------------------------------
    n_pass = sum([results.Passed]);
    n_fail = sum([results.Failed]);
    n_skip = sum([results.Incomplete]);

    fprintf('\n========================================\n');
    fprintf('  PRESTUS test suite — level: %s\n', upper(level));
    fprintf('  Passed: %d   Failed: %d   Skipped: %d\n', n_pass, n_fail, n_skip);
    fprintf('========================================\n\n');

    if n_fail > 0
        disp(table(results));
        error('run_all_tests:failures', '%d test(s) failed.', n_fail);
    end

end
