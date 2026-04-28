classdef test_isppa_rescaling < matlab.unittest.TestCase
% TEST_ISPPA_RESCALING  Tests for multi-ISPPA pipeline and ISPPA scaling.
%
%   Tests the multi_isppa_pipeline dispatcher and the pressure-scaling logic
%   (apply_isppa_scaling) using a fast water simulation. No MRI or SimNIBS
%   required. Expected runtime: ~5–15 min on CPU (two thermal stages).
%
%   Run with:
%     results = runtests('tests/test_isppa_rescaling.m');

    properties
        SimPath
        Parameters
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.SimPath = fullfile(tempdir, sprintf('prestus_isppa_test_%s', ...
                                  datestr(now, 'yyyymmdd_HHMMSS')));
            mkdir(tc.SimPath);
            tc.Parameters = make_isppa_parameters(tc.SimPath);
        end
    end

    methods (TestMethodTeardown)
        function teardown(tc)
            if exist(tc.SimPath, 'dir')
                rmdir(tc.SimPath, 's');
            end
        end
    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'isppa', 'pipeline'})

        function test_multi_isppa_runs_without_error(tc)
            % multi_isppa_pipeline should complete for two targets on water
            tc.verifyWarningFree(@() multi_isppa_pipeline(tc.Parameters));
        end

        function test_output_csv_created_per_target(tc)
            % One output CSV should be written per ISPPA target
            multi_isppa_pipeline(tc.Parameters);
            targets = tc.Parameters.calibration.target_isppa_wcm2;
            for ti = 1:numel(targets)
                affix   = sprintf('_isppa%03d', round(targets(ti)));
                pattern = fullfile(tc.SimPath, '**', ...
                    sprintf('*output_table%s.csv', affix));
                files   = dir(pattern);
                tc.verifyNotEmpty(files, ...
                    sprintf('Expected output CSV for target %.0f W/cm² (affix %s)', ...
                            targets(ti), affix));
            end
        end

        function test_intensity_scales_with_target(tc)
            % The focal Isppa in each output table should be approximately
            % proportional to the target (linear regime: I ∝ p², p ∝ scaling factor).
            multi_isppa_pipeline(tc.Parameters);
            targets = tc.Parameters.calibration.target_isppa_wcm2;

            isppa_vals = zeros(1, numel(targets));
            for ti = 1:numel(targets)
                affix   = sprintf('_isppa%03d', round(targets(ti)));
                pattern = fullfile(tc.SimPath, '**', ...
                    sprintf('*output_table%s.csv', affix));
                files   = dir(pattern);
                tc.assertNotEmpty(files, ...
                    sprintf('Missing CSV for target %.0f W/cm²', targets(ti)));
                tbl = readtable(fullfile(files(1).folder, files(1).name));
                col = tbl{:, contains(tbl.Properties.VariableNames, 'Isppa', 'IgnoreCase', true)};
                isppa_vals(ti) = col(1);
            end

            % Ratio of outputs should match ratio of targets to within 20%
            expected_ratio = targets(2) / targets(1);
            actual_ratio   = isppa_vals(2) / isppa_vals(1);
            tc.verifyEqual(actual_ratio, expected_ratio, 'RelTol', 0.20, ...
                sprintf(['Isppa ratio (%.2f) should match target ratio (%.2f) ' ...
                         'to within 20%% (linear regime).'], ...
                         actual_ratio, expected_ratio));
        end

        function test_each_output_isppa_positive(tc)
            % Both thermal outputs should have positive focal intensity
            multi_isppa_pipeline(tc.Parameters);
            targets = tc.Parameters.calibration.target_isppa_wcm2;
            for ti = 1:numel(targets)
                affix   = sprintf('_isppa%03d', round(targets(ti)));
                pattern = fullfile(tc.SimPath, '**', ...
                    sprintf('*output_table%s.csv', affix));
                files   = dir(pattern);
                tc.assertNotEmpty(files);
                tbl = readtable(fullfile(files(1).folder, files(1).name));
                col = tbl{:, contains(tbl.Properties.VariableNames, 'Isppa', 'IgnoreCase', true)};
                tc.verifyGreaterThan(col(1), 0, ...
                    sprintf('Isppa should be positive for target %.0f W/cm²', targets(ti)));
            end
        end

        function test_acoustic_cache_reused(tc)
            % Stage 1 acoustic cache file should exist after the run
            multi_isppa_pipeline(tc.Parameters);
            p    = tc.Parameters;
            subj = p.subject_id;
            med  = p.simulation.medium;
            cache_pat = fullfile(tc.SimPath, '**', ...
                sprintf('sub-%03d_%s_results*.mat', subj, med));
            files = dir(cache_pat);
            tc.verifyNotEmpty(files, 'Acoustic cache .mat should be saved by Stage 1');
        end

        function test_single_target_passthrough(tc)
            % A scalar target_isppa_wcm2 should raise an error from multi_isppa_pipeline
            % (it requires at least two targets — scalar use goes via prestus_pipeline).
            p = tc.Parameters;
            p.calibration.target_isppa_wcm2 = 10;
            tc.verifyError(@() multi_isppa_pipeline(p), ...
                'multi_isppa_pipeline:*');
        end

    end

end

% =========================================================================
%% Fixture helper
% =========================================================================

function parameters = make_isppa_parameters(sim_path)
% Build minimal water parameters extended for multi-ISPPA testing.
    parameters = make_minimal_parameters(sim_path);

    % Enable heating and thermal analysis so thermal stages produce output tables
    parameters.modules.run_heating_sims     = 1;
    parameters.modules.run_thermal_analysis = 1;
    parameters.modules.run_water_baseline   = 1;
    parameters.modules.generate_report      = 0;   % skip HTML for speed

    % Two ISPPA targets — triggers multi_isppa_pipeline dispatch
    parameters.calibration.target_isppa_wcm2 = [5, 20];

    % Overwrite to keep the test hermetic
    parameters.io.overwrite_files  = 'always';
    parameters.io.output_affix     = '';   % managed by multi_isppa_pipeline
    parameters.io.save_acoustic_matrices = 1;
end
