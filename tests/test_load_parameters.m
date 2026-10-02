classdef test_load_parameters < matlab.unittest.TestCase
% TEST_LOAD_PARAMETERS  Unit tests for load_parameters() config merging.
%
%   These tests operate on config_default.yaml only (no SimNIBS/k-Wave needed).
%   Run with:   results = runtests('tests/test_load_parameters.m');

    methods (TestClassSetup)
        function add_prestus_path(tc)
            % Resolve repo root relative to this test file and add functions/
            test_dir   = fileparts(mfilename('fullpath'));
            repo_root  = fileparts(test_dir);
            addpath(genpath(fullfile(repo_root, 'functions')));
            addpath(genpath(fullfile(repo_root, 'external')));
            % load_parameters expects to find config_default.yaml on the path
            addpath(fullfile(repo_root, 'config'));
        end
    end

    methods (Access = private)
        function override = quiet(~, override)
            % Suppress parameter summary output in tests that don't test it.
            if nargin < 2; override = struct(); end
            override.simulation.debug = 0;
        end
    end

    methods (Test, TestTags = {'defaults'})

        function test_default_config_loads(tc)
            % Must not error; interactive guard bypassed by non-desktop mode
            p = load_parameters(tc.quiet());
            tc.verifyClass(p, 'struct');
        end

        function test_top_level_keys_present(tc)
            p = load_parameters(tc.quiet());
            expected = {'simulation','path','io','modules','hpc','grid', ...
                        'headmodel','layers','medium_properties','timing','thermal','analysis'};
            for i = 1:numel(expected)
                tc.verifyTrue(isfield(p, expected{i}), ...
                    sprintf('Missing top-level key: %s', expected{i}));
            end
        end

        function test_simulation_defaults(tc)
            p = load_parameters(tc.quiet());
            tc.verifyTrue(ischar(p.simulation.medium) || isstring(p.simulation.medium), ...
                'simulation.medium should be a string');
            tc.verifyTrue(contains(string(p.simulation.medium), 'layered'), ...
                sprintf('Expected medium to contain ''layered'', got: %s', p.simulation.medium));
            tc.verifyTrue(contains(string(p.simulation.code_type), 'matlab_gpu'), ...
                sprintf('Expected code_type to contain ''matlab_gpu'', got: %s', p.simulation.code_type));
            % precision may be 'single' or "'single'" depending on YAML parser
            tc.verifyTrue(contains(string(p.simulation.precision), 'single'), ...
                sprintf('Expected precision to contain ''single'', got: %s', p.simulation.precision));
            % interactive may be logical false or numeric 0
            tc.verifyTrue(~p.simulation.interactive);
        end

        function test_modules_defaults(tc)
            p = load_parameters(tc.quiet());
            tc.verifyEqual(p.modules.run_grid_setup,            1);
            tc.verifyEqual(p.modules.run_heating_sims,          0);
            tc.verifyEqual(p.modules.segmentation_only,         0);
            tc.verifyEqual(p.modules.generate_report,           1);
            tc.verifyEqual(p.modules.run_transducer_placement,  1);
            tc.verifyEqual(string(p.modules.run_water_baseline), "auto");
        end

        function test_grid_defaults(tc)
            p = load_parameters(tc.quiet());
            tc.verifyEqual(p.grid.resolution_mm,     0.5);
            tc.verifyEqual(p.grid.pml_size,          "auto");
            tc.verifyEqual(p.grid.use_kWaveArray,    1);
        end

        function test_hpc_defaults(tc)
            p = load_parameters(tc.quiet());
            tc.verifyEqual(p.hpc.max_wait_checks, 540);
            tc.verifyFalse(p.hpc.wait_for_job);
        end

    end

    methods (Test, TestTags = {'merging'})

        function test_study_config_overrides_default(tc)
            % Supply a minimal struct override
            override = tc.quiet();
            override.simulation.medium = 'water';
            override.simulation.interactive = 0;
            p = load_parameters(override);
            tc.verifyEqual(p.simulation.medium, 'water');
        end

        function test_override_does_not_clobber_unrelated_fields(tc)
            override = tc.quiet();
            override.simulation.medium = 'water';
            override.simulation.interactive = 0;
            p = load_parameters(override);
            % grid fields untouched
            tc.verifyEqual(p.grid.resolution_mm, 0.5);
        end

    end

    methods (Test, TestTags = {'sanitization'})

        function test_output_affix_sanitized(tc)
            override = tc.quiet();
            override.io.output_affix = 'bad/chars!here';
            override.simulation.interactive = 0;
            p = load_parameters(override);
            % Illegal chars are replaced and value is wrapped in _desc- entity
            tc.verifyEqual(p.io.output_affix, '_desc-bad_chars_here');
        end

        function test_output_affix_valid_passthrough(tc)
            override = tc.quiet();
            override.io.output_affix = 'run_01';
            override.simulation.interactive = 0;
            p = load_parameters(override);
            % Valid affix is wrapped in _desc- entity
            tc.verifyEqual(p.io.output_affix, '_desc-run_01');
        end

    end

    methods (Test, TestTags = {'backends'})

        function test_unknown_thermal_backend_errors(tc)
            override = tc.quiet();
            override.simulation.interactive = 0;
            override.thermal_backend = 'fdtd';
            tc.verifyError(@() load_parameters(override), 'prestus:unknownThermalBackend');
        end

        function test_babelbrain_thermal_requires_babelbrain_acoustics(tc)
            override = tc.quiet();
            override.simulation.interactive = 0;
            override.acoustic_backend = 'kwave';
            override.thermal_backend  = 'babelbrain';
            tc.verifyError(@() load_parameters(override), ...
                'prestus:thermalBackendRequiresBabelBrainAcoustics');
        end

        function test_babelbrain_thermal_with_babelbrain_acoustics_ok(tc)
            override = tc.quiet();
            override.simulation.interactive = 0;
            override.acoustic_backend = 'babelbrain';
            override.thermal_backend  = 'babelbrain';
            p = load_parameters(override);
            tc.verifyEqual(p.thermal_backend, 'babelbrain');
        end

        function test_computing_backend_name_normalised(tc)
            override = tc.quiet();
            override.simulation.interactive = 0;
            override.babelbrain.computing_backend = 'metal';
            p = load_parameters(override);
            tc.verifyEqual(p.babelbrain.computing_backend, 3);
        end

        function test_computing_backend_auto_kept(tc)
            p = load_parameters(tc.quiet());
            tc.verifyEqual(char(p.babelbrain.computing_backend), 'auto');
        end

    end

end
