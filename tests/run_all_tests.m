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
%     PRESTUS_DEMO_CONFIG — path to study config YAML  (default: config_tutorial.yaml)
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
    addpath(fullfile(repo_root, 'config'));
    addpath(genpath(fullfile(repo_root, 'external')));
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
    % Suppress warnings that are expected/irrelevant during testing
    prev_warn = warning('off', 'prestus:toolboxDefault');
    warning('off', 'prestus:noSimNIBS');

    runner  = matlab.unittest.TestRunner.withNoPlugins();
    results = runner.run(suites);

    warning(prev_warn);

    % ---- Per-test output ----------------------------------------------
    descriptions = build_descriptions();

    fprintf('\n========================================\n');
    fprintf('  PRESTUS test suite — level: %s\n', upper(level));
    fprintf('========================================\n');

    for i = 1:numel(results)
        r    = results(i);
        name = r.Name;

        % Strip "classdef/method" → "method" for display
        parts = strsplit(name, '/');
        method = parts{end};

        % Look up human-readable description; fall back to method name
        if isKey(descriptions, method)
            label = descriptions(method);
        else
            label = method;
        end

        if r.Passed
            status = 'PASS';
        elseif r.Incomplete
            status = 'SKIP';
        else
            status = 'FAIL';
        end

        fprintf('  [%s]  %s\n', status, label);

        if r.Failed && isfield(r.Details, 'DiagnosticRecord') && ~isempty(r.Details.DiagnosticRecord)
            for d = r.Details.DiagnosticRecord
                if ~isempty(d.Exception)
                    fprintf('         -> %s\n', d.Exception.message);
                end
            end
        end
    end

    % ---- Summary ------------------------------------------------------
    n_pass = sum([results.Passed]);
    n_fail = sum([results.Failed]);
    n_skip = sum([results.Incomplete]);

    fprintf('----------------------------------------\n');
    fprintf('  Passed: %d   Failed: %d   Skipped: %d\n', n_pass, n_fail, n_skip);
    fprintf('========================================\n\n');

    if n_fail > 0
        error('run_all_tests:failures', '%d test(s) failed.', n_fail);
    end

end

% =========================================================================
function m = build_descriptions()
% Returns a containers.Map of test method name → human-readable description.
    keys = {
        % --- test_load_parameters ---
        'test_default_config_loads',                'Default config loads without error'
        'test_top_level_keys_present',              'Default config contains all required top-level sections'
        'test_simulation_defaults',                 'simulation.medium / code_type / precision / interactive defaults are correct'
        'test_modules_defaults',                    'Module enable flags match documented defaults'
        'test_grid_defaults',                       'Grid resolution, PML size, and kWaveArray defaults are correct'
        'test_hpc_defaults',                        'HPC wait_checks and wait_for_job defaults are correct'
        'test_study_config_overrides_default',      'Study-config override replaces default value (medium -> water)'
        'test_override_does_not_clobber_unrelated_fields', 'Override does not disturb unrelated fields (grid resolution)'
        'test_output_affix_sanitized',              'output_affix with illegal chars is sanitized to alphanumeric/underscore'
        'test_output_affix_valid_passthrough',      'Valid output_affix passes through unchanged'
        % --- test_helper ---
        'test_round_if_integer_exact',              'round_if_integer: exact integer value rounds correctly'
        'test_round_if_integer_within_tolerance',   'round_if_integer: near-integer within tolerance is accepted'
        'test_round_if_integer_array',              'round_if_integer: element-wise rounding works for [1.0, 2.0]'
        'test_round_if_integer_non_integer_errors', 'round_if_integer: non-integer input raises an error'
        'test_find_min_factor_power_of_two',        'find_min_factor: prefers 128 (power of 2) over primes in [120,130]'
        'test_find_min_factor_single_element_range','find_min_factor: single-value range returns that value'
        'test_find_min_factor_returns_value_in_range','find_min_factor: result falls within the requested range'
        'test_get_crop_dims_single_voxel',          'get_crop_dims: single-voxel mask with margin gives correct bounds'
        'test_get_crop_dims_cuboid_object',         'get_crop_dims: cuboid mask with margin gives correct min/max corners'
        'test_get_crop_dims_grid_size_consistent',  'get_crop_dims: returned size equals max - min + 1'
        'test_masked_max_3d_known_location',        'masked_max_3d: finds correct value and voxel index in full mask'
        'test_masked_max_3d_mask_excludes_true_max','masked_max_3d: mask excludes true global max, returns in-mask max'
        'test_charm_seg_labels_returns_struct',     'charm_seg_labels: returns a struct'
        'test_charm_seg_labels_expected_fields',    'charm_seg_labels: has all expected tissue fields'
        'test_charm_seg_labels_values',             'charm_seg_labels: numeric label values are correct'
        'test_flhm_symmetric_gaussian',             'get_flhm_center_position: Gaussian at 0 gives centre = 0'
        'test_flhm_shifted_gaussian',               'get_flhm_center_position: Gaussian at 10 gives centre = 10'
        'test_cast_struct_single_level',            'cast_struct: top-level numeric fields cast to target type'
        'test_cast_struct_nested',                  'cast_struct: nested numeric fields cast to target type'
        'test_cast_struct_preserves_non_numeric',   'cast_struct: non-numeric fields (char) are left unchanged'
        'test_get_xyz_mesh_size',                   'get_xyz_mesh: output has N_vox rows and 3 columns'
        'test_get_xyz_mesh_range',                  'get_xyz_mesh: coordinate ranges span [1, dim] for each axis'
        'test_zip_fields_order',                    'zip_fields: interleaves fieldnames and values in order'
        'test_subset_fields_copies_requested',      'subset_fields: copies only the requested fields'
        % --- test_thermal_parameters ---
        'test_duty_cycle',                          'thermal_parameters: duty cycle = pd / pri'
        'test_prf',                                 'thermal_parameters: pulse repetition frequency = 1 / pri'
        'test_pulses_per_train',                    'thermal_parameters: pulse count per train = ptd / pri'
        'test_ptri_repetitions',                    'thermal_parameters: PTRI repetition count = ptrd / ptri'
        'test_on_off_steps_sum_to_pri',             'thermal_parameters: ON + OFF step durations sum to PRI'
        'test_non_integer_pulses_per_train_errors', 'thermal_parameters: non-integer pulse count raises an error'
        % --- test_transform ---
        'test_ras_to_grid_identity_affine',         'ras_to_grid: identity affine leaves RAS coords unchanged'
        'test_ras_to_grid_translation',             'ras_to_grid: translation affine applies correct grid offset'
        'test_ras_to_grid_column_input',            'ras_to_grid: column-vector RAS input returns row output'
        'test_radial_expand_output_size',           'radialExpand2DTo3D: output volume is [2Nr x 2Nr x C]'
        'test_radial_expand_symmetry',              'radialExpand2DTo3D: output is symmetric (transposition invariant)'
        % --- test_transducer_placement ---
        'test_manual_mode_is_passthrough',          'placement (manual): trans_pos and focus_pos are unchanged'
        'test_default_mode_is_manual',              'placement (no field): defaults to manual passthrough'
        'test_unknown_mode_errors',                 'placement: unknown mode string raises an error'
        'test_localite_missing_file_errors',        'placement (localite): missing file/path raises an error'
        'test_localite_matrix_to_positions_identity','localite_matrix_to_positions: identity translation round-trips correctly'
        'test_localite_xml_parse_returns_valid_positions','position_transducer_localite: real XML returns finite non-zero positions'
        % --- test_isppa_rescaling ---
        'test_multi_isppa_runs_without_error',      'multi_isppa_pipeline: completes for two targets without error'
        'test_output_csv_created_per_target',       'multi_isppa_pipeline: one output CSV written per ISPPA target'
        'test_intensity_scales_with_target',        'multi_isppa_pipeline: focal Isppa ratio matches target ratio (+-20%)'
        'test_each_output_isppa_positive',          'multi_isppa_pipeline: focal Isppa is positive for each target'
        'test_acoustic_cache_reused',               'multi_isppa_pipeline: acoustic cache .mat is written by Stage 1'
        'test_single_target_passthrough',           'multi_isppa_pipeline: scalar target_isppa raises an error'
    };

    m = containers.Map(keys(:,1), keys(:,2));
end
