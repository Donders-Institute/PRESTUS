classdef test_integration_pipeline < matlab.unittest.TestCase
% TEST_INTEGRATION_PIPELINE  Smoke/integration tests for full and partial pipeline runs.
%
%   These tests require demo input data (Ernie template or equivalent).
%   They are intentionally coarse: they verify that the pipeline reaches
%   expected output checkpoints without crashing, rather than checking
%   numerical accuracy.
%
%   Test levels (controlled by TestTags):
%     'smoke_config'    — config loading only (seconds; no I/O)
%     'smoke_head'      — head preprocessing up to medium masks (~1–5 min)
%     'smoke_acoustic'  — full acoustic simulation (minutes–hours; requires GPU/CPU)
%     'smoke_thermal'   — thermal simulation (minutes–hours; requires acoustic outputs)
%
%   Run a specific level:
%     runtests('tests/test_integration_pipeline.m', 'Tag', 'smoke_config')
%
%   Run all integration tests:
%     runtests('tests/test_integration_pipeline.m')
%
%   Prerequisites:
%     - Set environment variable PRESTUS_TEST_DATA to the folder containing
%       the demo subject's T1/T2 and SimNIBS m2m_* output.
%     - Set PRESTUS_DEMO_CONFIG to the path of the demo study config YAML.
%       (Defaults to configs/tutorial_config.yaml if not set.)

    properties (Constant)
        DEMO_SUBJECT_ID = 1
    end

    properties
        params       % loaded parameter struct for this subject
        repo_root    % absolute path to repo root
        data_root    % path to test data (from env var or skipped)
    end

    methods (TestClassSetup)
        function resolve_paths(tc)
            test_dir      = fileparts(mfilename('fullpath'));
            tc.repo_root  = fileparts(test_dir);
            addpath(genpath(fullfile(tc.repo_root, 'functions')));
            addpath(fullfile(tc.repo_root, 'configs'));
        end
    end

    methods (TestMethodSetup)
        function load_demo_params(tc)
            % Resolve demo data root — skip gracefully if not available
            data_env = getenv('PRESTUS_TEST_DATA');
            if isempty(data_env) || ~isfolder(data_env)
                tc.assumeFail('PRESTUS_TEST_DATA not set or not a valid folder. Skipping integration tests.');
            end
            tc.data_root = data_env;

            cfg_env = getenv('PRESTUS_DEMO_CONFIG');
            if isempty(cfg_env)
                cfg_env = fullfile(tc.repo_root, 'configs', 'tutorial_config.yaml');
            end
            tc.assumeTrue(isfile(cfg_env), ...
                sprintf('Demo config not found: %s', cfg_env));

            cfg_struct = yaml.loadFile(cfg_env, 'ConvertToArray', true);
            cfg_struct.simulation.debug = 0;
            tc.params = load_parameters(cfg_struct);
            tc.params.subject_id             = tc.DEMO_SUBJECT_ID;
            tc.params.simulation.interactive = 0;
            tc.params.io.overwrite_files    = 'always';
            tc.params.platform              = 'matlab';
        end
    end

    % ------------------------------------------------------------------ %
    %% Level 1: config loading (no data needed — prerequisite always runs)
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'smoke_config'})

        function test_tutorial_config_loads(tc)
            cfg = fullfile(tc.repo_root, 'configs', 'tutorial_config.yaml');
            tc.assumeTrue(isfile(cfg));
            % Merge the YAML config with a debug-suppressing struct so that
            % load_parameters does not print the parameter summary.
            cfg_struct = yaml.loadFile(cfg, 'ConvertToArray', true);
            cfg_struct.simulation.debug = 0;
            p = load_parameters(cfg_struct);
            tc.verifyClass(p, 'struct');
            tc.verifyTrue(isfield(p, 'transducer'));
        end

    end

    % ------------------------------------------------------------------ %
    %% Level 2: head preprocessing (segmentation + medium masks)
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'smoke_head'})

        function test_segmentation_file_exists(tc)
            % Verify the expected charm output is present for this subject
            m2m = fullfile(tc.data_root, ...
                sprintf('m2m_sub-%03d', tc.DEMO_SUBJECT_ID), 'final_tissues.nii.gz');
            tc.verifyTrue(isfile(m2m), ...
                'SimNIBS charm output not found — run segmentation first.');
        end

        function test_preproc_head_runs(tc)
            % Run up to and including head preprocessing; check output struct
            p = tc.params;
            p.modules.run_medium_setup  = 0;
            p.modules.run_source_setup  = 0;
            p.modules.run_acoustic_sims = 0;
            p.modules.run_acoustic_analysis = 0;
            p.modules.run_heating_sims  = 0;
            p.modules.run_nifti_creation = 0;
            p.modules.generate_report   = 0;
            p.modules.run_posthoc_water_sims = 0;

            % Call preproc directly to avoid full pipeline path setup
            [p_out, medium_masks, segmentation, ~, ~] = preproc_head(p);
            tc.verifyNotEmpty(medium_masks);
            tc.verifyNotEmpty(segmentation);
            tc.verifyTrue(isfield(p_out, 'grid'));
            tc.verifyTrue(isfield(p_out.grid, 'dims'));
        end

        function test_medium_masks_label_range(tc)
            p = tc.params;
            [~, medium_masks, ~, ~, ~] = preproc_head(p);
            n_layers = numel(fieldnames(p.medium_properties));
            tc.verifyGreaterThanOrEqual(min(medium_masks(:)), 0);
            tc.verifyLessThanOrEqual(max(medium_masks(:)), n_layers);
        end

    end

    % ------------------------------------------------------------------ %
    %% Level 3: full acoustic pipeline
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'smoke_acoustic'})

        function test_acoustic_output_file_created(tc)
            p = tc.params;
            p.modules.run_heating_sims       = 0;
            p.modules.run_thermal_analysis   = 0;
            p.modules.run_posthoc_water_sims = 0;
            p.modules.generate_report        = 0;

            prestus_pipeline(p);

            expected_file = fullfile(p.path.sim, ...
                sprintf('sub-%03d', tc.DEMO_SUBJECT_ID), ...
                sprintf('sub-%03d_%s%s_intensity_orig_coord.nii.gz', ...
                    tc.DEMO_SUBJECT_ID, p.simulation.medium, p.io.output_affix));
            tc.verifyTrue(isfile(expected_file), ...
                'Acoustic output NIfTI not created.');
        end

        function test_output_csv_created(tc)
            p = tc.params;
            p.modules.run_heating_sims       = 0;
            p.modules.run_posthoc_water_sims = 0;
            p.modules.generate_report        = 0;

            prestus_pipeline(p);

            csv_pattern = fullfile(p.path.sim, ...
                sprintf('sub-%03d', tc.DEMO_SUBJECT_ID), ...
                sprintf('sub-%03d_%s_output_table*.csv', ...
                    tc.DEMO_SUBJECT_ID, p.simulation.medium));
            hits = dir(csv_pattern);
            tc.verifyGreaterThan(numel(hits), 0, 'Output CSV not found.');
        end

    end

    % ------------------------------------------------------------------ %
    %% Level 4: thermal pipeline (requires completed acoustic run)
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'smoke_thermal'})

        function test_thermal_output_file_created(tc)
            p = tc.params;
            p.modules.run_acoustic_sims      = 0;  % reuse cached acoustic output
            p.modules.run_heating_sims       = 1;
            p.modules.run_posthoc_water_sims = 0;
            p.modules.generate_report        = 0;
            p.io.overwrite_files             = 'never'; % preserve acoustics

            prestus_pipeline(p);

            expected_file = fullfile(p.path.sim, ...
                sprintf('sub-%03d', tc.DEMO_SUBJECT_ID), ...
                sprintf('sub-%03d_%s%s_cem43_orig_coord.nii.gz', ...
                    tc.DEMO_SUBJECT_ID, p.simulation.medium, p.io.output_affix));
            tc.verifyTrue(isfile(expected_file), ...
                'Thermal (CEM43) output NIfTI not created.');
        end

    end

end
