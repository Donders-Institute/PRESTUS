classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_head_preprocessing < matlab.unittest.TestCase
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
            tc.assumeTrue(logical(license('test', 'Image_Toolbox')), ...
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

        function test_skull_fill_adds_cortical_bone(tc)
            % with a differentiated skull, filled voxels become cortical, not trabecular
            dims = [31, 31, 31];
            [X, Y, Z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
            c = (dims+1)/2;
            R = sqrt((X-c(1)).^2 + (Y-c(2)).^2 + (Z-c(3)).^2);
            seg = zeros(dims, 'uint8');
            seg(R >= 8 & R <= 10) = 7;   % cortical shell (charm label 7)
            seg(R < 8)            = 1;

            p = tc.params;
            p.headmodel.skull_fill_method = 'imclose';
            p.layers = rmfield(p.layers, 'skull');
            p.layers.skull_trabecular = 8;
            p.layers.skull_cortical   = 7;
            p.medium_properties.skull_trabecular = struct();
            p.medium_properties.skull_cortical   = struct();
            labels = fieldnames(p.medium_properties);
            cortical_i   = find(strcmp(labels, 'skull_cortical'));
            trabecular_i = find(strcmp(labels, 'skull_trabecular'));

            medium_masks = zeros(dims);
            medium_masks(seg == 7) = cortical_i;
            hole = R >= 8 & R <= 10 & abs(X - c(1)) <= 1 & abs(Y - c(2)) <= 1 & Z > c(3);
            medium_masks(hole) = 0;

            filled = skull_fill_holes(p, medium_masks, round(c), seg);
            added  = filled ~= medium_masks;
            tc.assertTrue(any(added(:)), 'the hole is filled');
            tc.verifyTrue(all(filled(added) == cortical_i), 'filled voxels are cortical');
            tc.verifyFalse(any(filled(:) == trabecular_i), 'no trabecular voxels are added');
        end

        function test_rubber_wrap_closes_gap_between_brain_and_skin(tc)
            % a skull gap where brain meets skin is part of the wrapped skull
            tc.assumeTrue(logical(license('test', 'Image_Toolbox')), ...
                'Image Processing Toolbox not available (niftiwrite)');
            % Flat layers: a curved voxel shell is too thin near its pole for
            % the wrap's closing to bridge a hole
            dims = [30, 30, 24];
            [X, Y, Z] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
            seg = zeros(dims, 'uint8');
            seg(:, :, 1:9)   = 1;   % brain
            seg(:, :, 10:14) = 4;   % skull
            seg(:, :, 15:18) = 5;   % skin
            % 3x3 gap through the skull: brain below, one empty voxel, skin above
            col = abs(X - 15) <= 1 & abs(Y - 15) <= 1;
            seg(col & Z >= 10 & Z <= 12) = 1;
            seg(col & Z == 13) = 0;
            seg(col & Z == 14) = 5;
            gap = col & Z == 13;

            labels = fieldnames(tc.params.medium_properties);
            medium_masks = zeros(dims, 'single');
            medium_masks(seg == 1) = find(strcmp(labels, 'brain'));
            medium_masks(seg == 5) = find(strcmp(labels, 'skin'));
            medium_masks(seg == 4) = find(strcmp(labels, 'skull'));

            % skull_rubber_wrap reads the T1 header from path.anat
            tmp = tempname; mkdir(tmp);
            tc.addTeardown(@() rmdir(tmp, 's'));
            niftiwrite(zeros(dims, 'single'), fullfile(tmp, 'T1.nii'));
            p = tc.params;
            p.path.anat = tmp; p.path.t1_pattern = 'T1.nii'; p.subject_id = '';

            wrapped = skull_rubber_wrap(p, seg == 4, medium_masks, seg);
            tc.verifyTrue(all(wrapped(gap)), 'gap voxels between brain and skin are skull');
            tc.verifyFalse(any(wrapped(medium_masks == find(strcmp(labels, 'brain')))), 'no brain voxel becomes skull');
        end

    end


    % smooth_img: replicate padding
    methods (Test, TestTags = {'smooth_img'})

        function test_gaussian_constant_volume_edges_preserved(tc)
            img = ones(16, 16, 16);
            out = smooth_img(img, 3, 1, 0, "gaussian");
            tc.verifyEqual(size(out), size(img));
            tc.verifyEqual(out, img, 'AbsTol', 1e-10);
        end

        function test_box_constant_volume_edges_preserved(tc)
            img = ones(16, 16, 16);
            out = smooth_img(img, 3, 1, 0, "box");
            tc.verifyEqual(size(out), size(img));
            tc.verifyEqual(out, img, 'AbsTol', 1e-10);
        end

        function test_binary_mask_touching_border_survives_threshold(tc)
            % A slab filling the lower half of the volume touches the border;
            % zero padding used to erode it at the faces after thresholding.
            img = zeros(16, 16, 16);
            img(:, :, 1:8) = 1;
            out = smooth_img(img, 3, 1, 0.5, "gaussian");
            tc.verifyTrue(all(out(:, :, 1:6) == 1, 'all'));
        end

    end
end
