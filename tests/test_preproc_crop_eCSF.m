classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_preproc_crop_eCSF < matlab.unittest.TestCase
% TEST_PREPROC_CROP_ECSF  CSF dilation reaches headmodel.csf_expansion on both sides of each axis.
%
% A single CSF voxel at the centre of a cube; every voxel within
% csf_expansion (Chebyshev distance, in voxels) keeps its label and every
% voxel beyond it becomes water. Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_preproc_crop_eCSF.m');

    properties (TestParameter)
        % [csf_expansion mm, grid resolution mm, expected reach in voxels]
        setup = struct('half_mm_grid', [2 0.5 4], 'coarse', [3 1 3], ...
                       'rounds_up', [2.6 1 3], 'rounds_down', [2.4 1 2], 'zero', [0 1 0]);
    end

    methods (Test)
        function test_reach_is_symmetric_and_exact(tc, setup)
            [masks, water, r, c] = run_crop(setup);
            for ax = 1:3
                for side = [-1 1]
                    inside  = c; inside(ax)  = c(ax) + side * r;
                    outside = c; outside(ax) = c(ax) + side * (r + 1);
                    tc.verifyNotEqual(masks(inside(1), inside(2), inside(3)), water, ...
                        sprintf('axis %d, side %+d: voxel at the reach must be kept', ax, side));
                    tc.verifyEqual(masks(outside(1), outside(2), outside(3)), water, ...
                        sprintf('axis %d, side %+d: voxel beyond the reach must be water', ax, side));
                end
            end
        end

        function test_kept_region_is_a_cube(tc, setup)
            [masks, water, r] = run_crop(setup);
            tc.verifyEqual(nnz(masks ~= water), (2*r + 1)^3);
        end
    end

    methods (Test)
        function test_matches_convolution_with_cube(tc)
            % scattered CSF voxels, some on the array faces: the kept region equals a
            % convolution with the (2r+1)^3 cube, for reaches below and above the array size
            p = load_parameters();
            p.grid.resolution_mm = 1;
            names = fieldnames(p.medium_properties);
            water = find(strcmp(names, 'water'));
            other = find(~strcmp(names, 'water'), 1);
            rng(7);
            seg = zeros(23, 19, 17);
            seg(randperm(numel(seg), 12)) = getidx(charm_seg_labels(), 'csf');
            seg(1, 5, 9) = getidx(charm_seg_labels(), 'csf');
            seg(23, 19, 17) = getidx(charm_seg_labels(), 'csf');
            for r = [0, 1, 4, 30]
                p.headmodel.csf_expansion = r;
                masks = preproc_crop_eCSF(p, other * ones(size(seg)), seg, [12 10 9]);
                ref = convn(double(seg == getidx(charm_seg_labels(), 'csf')), ones(2*r + 1, 2*r + 1, 2*r + 1), 'same') > 0;
                tc.verifyEqual(masks ~= water, ref, sprintf('reach %d voxels', r));
            end
        end
    end
end

function [masks, water, r, c] = run_crop(setup)
    p = load_parameters();
    p.headmodel.csf_expansion = setup(1);
    p.grid.resolution_mm      = setup(2);
    r = setup(3);
    n = 2*r + 5;                         % room for one voxel beyond the reach on each side
    c = [1 1 1] * (r + 3);
    seg = zeros(n, n, n);
    seg(c(1), c(2), c(3)) = getidx(charm_seg_labels(), 'csf');
    names = fieldnames(p.medium_properties);
    water = find(strcmp(names, 'water'));
    other = find(~strcmp(names, 'water'), 1);
    masks = preproc_crop_eCSF(p, other * ones(n, n, n), seg, c);
end
