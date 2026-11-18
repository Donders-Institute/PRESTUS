classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_medium_mask_unclaimed < matlab.unittest.TestCase
% TEST_MEDIUM_MASK_UNCLAIMED  Voxels that no smoothed layer keeps are water in
% the medium masks (as MEDIUM_SETUP simulates them), and unlabelled voxels do
% not widen the grid crop. Synthetic volumes; makeBowl (k-Wave) for the crop.
%
%   results = runtests('tests/test_medium_mask_unclaimed.m');

    methods (Test)
        function test_smoothing_gaps_become_water(tc)
            p = mask_params();
            seg = zeros(16, 16, 16);
            seg(5:12, 5:12, 5:12) = 1;                         % wm cube in "water"
            names   = fieldnames(p.medium_properties);
            i_water = find(strcmp(names, 'water'));
            i_brain = find(strcmp(names, 'brain'));
            brain_s = smooth_img(seg == 1, p.headmodel.smooth_fwhm_mm, p.grid.resolution_mm, ...
                p.headmodel.smooth_threshold_other, p.headmodel.smooth_method) ~= 0;
            water_s = smooth_img(seg ~= 1, p.headmodel.smooth_fwhm_mm, p.grid.resolution_mm, ...
                p.headmodel.smooth_threshold_other, p.headmodel.smooth_method) ~= 0;
            tc.assumeTrue(any(~brain_s(:) & ~water_s(:)), 'the threshold leaves no gap between the layers');

            masks = preproc_medium_mask(seg, p);
            tc.verifyFalse(any(masks(:) == 0), 'no voxel stays unlabelled');
            expected = repmat(single(i_water), size(seg));
            expected(brain_s) = i_brain;                       % negative control: kept tissue unchanged
            tc.verifyEqual(masks, expected);
        end

        function test_unlabelled_voxels_do_not_widen_the_crop(tc)
            [p, masks, seg, bone, tp, fp, i_water, i_brain] = crop_setup();
            ref = masks;
            ref(masks == 0) = i_water;                         % the same head, fully labelled
            [m_ref, ~, ~, ~, p_ref, tp_ref] = preproc_crop_grid(p, ref, seg, bone, tp, fp);
            [m_gap, ~, ~, ~, p_gap, tp_gap] = preproc_crop_grid(p, masks, seg, bone, tp, fp);
            tc.verifyEqual(p_gap.grid.dims, p_ref.grid.dims);
            tc.verifyEqual(tp_gap, tp_ref);
            tc.verifyEqual(m_gap, m_ref);
            tc.verifyFalse(any(m_gap(:) == 0));
            % Negative control: the crop still holds the whole brain and the transducer
            tc.verifyEqual(nnz(m_ref == i_brain), nnz(masks == i_brain));
            tc.verifyGreaterThanOrEqual(tp_ref, [1 1 1]);
        end
    end
end

function p = mask_params()
% Layered head settings with a high threshold, so the smoothed layers leave a gap
    p = struct();
    p.grid.resolution_mm = 0.5;
    p.headmodel.smooth_method          = 'gaussian';
    p.headmodel.smooth_fwhm_mm         = 1;
    p.headmodel.smooth_threshold_skull = 0.9;
    p.headmodel.smooth_threshold_other = 0.9;
    p.pct.enabled = 0;
    p.layers.water = [0, 3, 6, 9, 10];
    p.layers.brain = [1, 2];
    p.medium_properties.water = struct();
    p.medium_properties.brain = struct();
end

function [p, masks, seg, bone, tp, fp, i_water, i_brain] = crop_setup()
% 48^3 voxels of 1 mm: CSF block (kept by the CSF crop), brain cube, and an
% unlabelled slab beyond the brain on +x; transducer on -x
    p = load_parameters();
    p.grid.resolution_mm      = 1;
    p.grid.max_expand         = 0;          % no FFT enlargement: compare raw crop bounds
    p.headmodel.head_pad_mm   = 0;
    p.headmodel.csf_expansion = 10;
    p.transducer = struct('type', 'annular', 'annular', struct('elem_n', 1, ...
        'elem_od_mm', 9, 'elem_id_mm', 0, 'curv_radius_mm', 12));
    names   = fieldnames(p.medium_properties);
    i_water = find(strcmp(names, 'water'));
    i_brain = find(strcmp(names, 'brain'));
    n = 48;
    seg = zeros(n, n, n);
    seg(16:32, 16:32, 16:32) = getidx(charm_seg_labels(), 'csf');
    masks = repmat(single(i_water), n, n, n);
    masks(20:28, 20:28, 20:28) = i_brain;
    masks(29:36, 20:28, 20:28) = 0;           % unlabelled, inside the CSF crop
    bone = zeros(n, n, n);
    tp = [8 24 24];
    fp = [24 24 24];
end
