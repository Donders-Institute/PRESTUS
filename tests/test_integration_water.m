classdef test_integration_water < matlab.unittest.TestCase
% TEST_INTEGRATION_WATER  End-to-end water pipeline integration test.
%
%   Runs the full prestus_pipeline with simulation.medium = 'water' using a
%   single-element bowl transducer on a small grid. No MRI, no SimNIBS, no GPU.
%   Expected runtime: ~2–5 min on CPU.
%
%   Run with:
%     results = runtests('tests/test_integration_water.m');

    properties
        SimPath   % temporary output directory
        Parameters
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.SimPath = fullfile(tempdir, sprintf('prestus_water_test_%s', ...
                                  datestr(now, 'yyyymmdd_HHMMSS')));
            mkdir(tc.SimPath);
            tc.Parameters = make_minimal_parameters(tc.SimPath);
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
    methods (Test, TestTags = {'water', 'pipeline'})

        function test_pipeline_runs_without_error(tc)
            % Pipeline should complete without throwing
            tc.verifyWarningFree(@() prestus_pipeline(tc.Parameters));
        end

        function test_output_csv_created(tc)
            prestus_pipeline(tc.Parameters);
            csv_files = dir(fullfile(tc.SimPath, '**', '*.csv'));
            tc.verifyNotEmpty(csv_files, 'Expected at least one CSV output table');
        end

        function test_nifti_output_created(tc)
            prestus_pipeline(tc.Parameters);
            nii_files = dir(fullfile(tc.SimPath, '**', '*.nii*'));
            tc.verifyNotEmpty(nii_files, 'Expected at least one NIfTI output file');
        end

        function test_intensity_positive(tc)
            % Read the output CSV and check Isppa > 0
            prestus_pipeline(tc.Parameters);
            csv_files = dir(fullfile(tc.SimPath, '**', '*output_table*.csv'));
            tc.verifyNotEmpty(csv_files);
            tbl = readtable(fullfile(csv_files(1).folder, csv_files(1).name));
            intensity_col = tbl{:, contains(tbl.Properties.VariableNames, 'Isppa', 'IgnoreCase', true)};
            tc.verifyGreaterThan(intensity_col(1), 0, 'Isppa (I_sppa) should be positive');
        end

        function test_axisymmetric_water_runs(tc)
            % 2D axisymmetric should also complete without error
            p = tc.Parameters;
            p.grid.axisymmetric  = 1;
            p.grid.default_dims  = [1, 72, 128];   % 2D: first dim = 1
            p.io.output_affix    = '_test_2d';
            tc.verifyWarningFree(@() prestus_pipeline(p));
        end

    end

end
