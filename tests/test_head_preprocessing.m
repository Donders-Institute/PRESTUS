classdef test_head_preprocessing < matlab.unittest.TestCase
% TEST_HEAD_PREPROCESSING  Unit tests for head preprocessing helpers.
%
%   Tests functions that operate on synthetic segmentation volumes —
%   no SimNIBS or real MRI data required.
%
%   Run with:   results = runtests('tests/test_head_preprocessing.m');

    properties
        params   % minimal parameters struct
    end

    methods (TestMethodSetup)
        function build_params(tc)
            p = struct();
            p.grid.resolution_mm = 0.5;
            p.headmodel.smooth_method          = 'gaussian';
            p.headmodel.smooth_fwhm_mm         = 1;
            p.headmodel.smooth_threshold_skull = 0.5;
            p.headmodel.smooth_threshold_other = 0.5;
            p.headmodel.smooth_properties      = false;
            p.headmodel.head_pad_mm            = 0;
            p.headmodel.csf_expansion          = 5;
            p.headmodel.skull_fill_method      = 'rubberwrap';
            p.headmodel.skull_wrap_radius      = 3;
            p.headmodel.skull_wrap_visualize   = 0;
            p.pct.enabled                      = 0;
            % Fields required by skull_fill_holes
            p.subject_id                       = 1;
            p.simulation.medium                = 'layered';
            p.simulation.debug                 = 0;
            p.io.output_affix                  = '';
            p.io.output_dir                    = tempdir;
            p.io.debug_dir                     = tempdir;
            % Layers matching charm labels
            p.layers.water            = [0, 3, 6, 9, 10];
            p.layers.brain            = [1, 2];
            p.layers.skin             = [5];
            p.layers.skull            = [4];
            % Medium property fields — order matters: index used as label ID
            p.medium_properties.water = struct();
            p.medium_properties.brain = struct();
            p.medium_properties.skin  = struct();
            p.medium_properties.skull = struct();   % skull must be present; index = 4
            tc.params = p;
        end
    end

    % ------------------------------------------------------------------ %
    %% get_crop_dims with realistic head-shaped mask
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'crop'})

        function test_crop_dims_respects_margin(tc)
            vol = zeros(40, 40, 40);
            vol(10:30, 10:30, 10:30) = 1;
            margin = 5;
            [mn, mx, ~] = get_crop_dims(vol, margin);
            tc.verifyEqual(mn, [5, 5, 5]);
            tc.verifyEqual(mx, [35, 35, 35]);
        end

    end

    % ------------------------------------------------------------------ %
    %% preproc_medium_mask — synthetic segmentation
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'medium_mask'})

        function test_medium_mask_labels_in_range(tc)
            % Synthetic 10x10x10 segmentation: brain inside, skull shell, skin outer
            seg = ones(10, 10, 10);          % label 1 = wm everywhere
            seg([1,10],:,:) = 5;             % skin at borders
            seg([2,9],:,:)  = 4;             % skull just inside skin
            masks = preproc_medium_mask(seg, tc.params);
            n_layers = numel(fieldnames(tc.params.medium_properties));
            tc.verifyGreaterThanOrEqual(min(masks(:)), 0);
            tc.verifyLessThanOrEqual(max(masks(:)), n_layers);
        end

        function test_medium_mask_brain_label_assigned(tc)
            seg = zeros(10, 10, 10);
            seg(4:7, 4:7, 4:7) = 1;   % wm in centre
            masks = preproc_medium_mask(seg, tc.params);
            layer_names = fieldnames(tc.params.medium_properties);
            brain_id = find(strcmp(layer_names, 'brain'));
            tc.verifyTrue(any(masks(:) == brain_id), ...
                'Brain label should be assigned for wm voxels');
        end

        function test_medium_mask_size_preserved(tc)
            seg = uint8(ones(12, 14, 16));
            masks = preproc_medium_mask(seg, tc.params);
            tc.verifyEqual(size(masks), [12, 14, 16]);
        end

    end

    % ------------------------------------------------------------------ %
    %% skull_fill_holes — synthetic hollow skull
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'skull_fill'})

        function test_skull_fill_closes_small_hole(tc)
            tc.assumeTrue(license('test', 'Image_Toolbox'), ...
                'Image Processing Toolbox not available — skipping skull_fill_holes test');

            % Build a synthetic segmentation with a skull shell and a punched hole
            dims = [31, 31, 31];
            seg  = zeros(dims, 'uint8');
            [X, Y, Z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
            c = (dims+1)/2;
            R = sqrt((X-c(1)).^2 + (Y-c(2)).^2 + (Z-c(3)).^2);
            seg(R >= 8 & R <= 10) = 4;   % skull shell (charm label 4)
            seg(R < 8)            = 1;   % brain inside

            % Build medium_masks matching medium_properties field order
            skull_id = find(strcmp(fieldnames(tc.params.medium_properties), 'skull'));
            medium_masks = zeros(dims);
            medium_masks(seg == 4) = skull_id;

            % Punch a hole and count skull voxels before fill
            medium_masks(1:5, 1:5, 1:5) = 0;
            n_before = sum(medium_masks(:) > 0);

            p = tc.params;
            p.headmodel.skull_fill_method = 'imclose';
            [filled, ~] = skull_fill_holes(p, medium_masks, round(c), seg);
            tc.verifyGreaterThanOrEqual(sum(filled(:) > 0), n_before);
        end

    end

end
