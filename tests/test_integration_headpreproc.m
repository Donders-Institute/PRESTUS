classdef test_integration_headpreproc < matlab.unittest.TestCase
% TEST_INTEGRATION_HEADPREPROC  Integration tests for head preprocessing functions.
%
%   Tests preproc_medium_mask and skull processing using a synthetic
%   concentric-shell segmentation. No MRI files needed.
%
%   Run with:
%     results = runtests('tests/test_integration_headpreproc.m');

    properties
        Parameters
        Segmentation
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.Parameters   = make_minimal_parameters();
            tc.Segmentation = make_synthetic_segmentation([64, 64, 80]);

            % Switch to layered for head preprocessing tests
            tc.Parameters.simulation.medium = 'layered';

            % Minimal headmodel settings
            tc.Parameters.headmodel.smooth_method          = 'gaussian';
            tc.Parameters.headmodel.smooth_fwhm_mm         = 1;
            tc.Parameters.headmodel.smooth_threshold_skull = 0.5;
            tc.Parameters.headmodel.smooth_threshold_other = 0.5;
            tc.Parameters.headmodel.smooth_properties      = false;
            tc.Parameters.headmodel.skull_fill_method      = 'imclose';  % faster than rubberwrap for tests

            % Single-skull layer setup (no cortical/trabecular split)
            tc.Parameters.layers.water  = [0, 3, 6, 9, 10];
            tc.Parameters.layers.brain  = [1, 2];
            tc.Parameters.layers.skin   = [5];
            tc.Parameters.layers.skull  = [4];
        end
    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'headpreproc', 'medium_mask'})

        function test_medium_mask_returns_correct_size(tc)
            masks = preproc_medium_mask(tc.Segmentation, tc.Parameters);
            tc.verifyEqual(size(masks), size(tc.Segmentation));
        end

        function test_medium_mask_has_expected_labels(tc)
            masks = preproc_medium_mask(tc.Segmentation, tc.Parameters);
            unique_labels = unique(masks);
            % Should contain at least water (0 or 1) and other tissues
            tc.verifyGreaterThanOrEqual(numel(unique_labels), 2);
        end

        function test_medium_mask_no_unlabelled_interior(tc)
            % All voxels inside the skin shell should be assigned a non-water label
            masks = preproc_medium_mask(tc.Segmentation, tc.Parameters);
            medium_labels = fieldnames(tc.Parameters.medium_properties);
            i_water = find(strcmp(medium_labels, 'water'));
            interior = tc.Segmentation > 0;   % inside skin sphere
            interior_masks = masks(interior);
            frac_water = mean(interior_masks == i_water);
            % Allow some water voxels at boundaries (smoothing), but not majority
            tc.verifyLessThan(frac_water, 0.5, ...
                'More than 50%% of interior voxels assigned as water — check layer mapping');
        end

    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'headpreproc', 'skull'})

        function test_skull_fill_holes_runs(tc)
            % skull_fill_holes should not error on synthetic skull
            skull_mask = ismember(tc.Segmentation, tc.Parameters.layers.skull);
            tc.verifyWarningFree(@() skull_fill_holes(skull_mask, tc.Parameters));
        end

        function test_skull_fill_holes_output_size(tc)
            skull_mask = ismember(tc.Segmentation, tc.Parameters.layers.skull);
            filled = skull_fill_holes(skull_mask, tc.Parameters);
            tc.verifyEqual(size(filled), size(skull_mask));
        end

        function test_skull_fill_holes_non_empty(tc)
            skull_mask = ismember(tc.Segmentation, tc.Parameters.layers.skull);
            filled = skull_fill_holes(skull_mask, tc.Parameters);
            tc.verifyGreaterThan(sum(filled(:)), 0, 'Filled skull should contain non-zero voxels');
        end

    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'headpreproc', 'crop'})

        function test_get_crop_dims_on_segmentation(tc)
            brain_mask = tc.Segmentation > 0;
            margin = 5;
            [mn, mx, sz] = get_crop_dims(brain_mask, margin);
            tc.verifyEqual(sz, mx - mn + 1);
            tc.verifyGreaterThan(min(sz), 0);
        end

    end

end
