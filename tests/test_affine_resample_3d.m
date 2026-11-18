classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_affine_resample_3d < matlab.unittest.TestCase
% TEST_AFFINE_RESAMPLE_3D  Unit tests for the interpolation methods of affine_resample_3d
%
%   Run with:   results = runtests('tests/test_affine_resample_3d.m');
%
%   A quadratic field sampled half-way between voxels separates the methods:
%   cubic convolution reproduces quadratics exactly, linear interpolation
%   overestimates x^2 by r*(1-r) = 0.25 per quadratic axis, and nearest
%   returns a neighbouring voxel value.

    properties (TestParameter)
        nC = {1, 2};   % 3D and 4D (per-channel) code paths
    end

    methods (Static)
        function [vol, M, out_size, x, y, z] = quadratic_case(nC)
            % Input 10^3 voxels, f = x^2 + 2*y^2 + 3*z^2 (+ channel offset).
            % Output voxel i samples input index i + 2.5 (interior: all
            % cubic neighbours lie within 1..10, so boundary handling is not tested).
            [x, y, z] = ndgrid(1:10);
            vol = zeros([10 10 10 nC], 'single');
            for c = 1:nC
                vol(:,:,:,c) = x.^2 + 2*y.^2 + 3*z.^2 + 100*(c-1);
            end
            M = eye(4);
            M(4, 1:3) = -2.5;             % [out 1] = [in 1] * M  ->  in = out + 2.5
            out_size = [5 5 5];
            [x, y, z] = ndgrid((1:5) + 2.5);
        end
    end

    methods (Test)

        function test_cubic_reproduces_quadratic(tc, nC)
            [vol, M, out_size, x, y, z] = test_affine_resample_3d.quadratic_case(nC);
            out = affine_resample_3d(vol, M, out_size, 'cubic', 0);
            for c = 1:nC
                expected = single(x.^2 + 2*y.^2 + 3*z.^2 + 100*(c-1));
                tc.verifyEqual(out(:,:,:,c), expected, 'RelTol', single(1e-5), ...
                    'cubic must reproduce a quadratic field at interior points');
            end
        end

        function test_cubic_differs_from_linear(tc)
            [vol, M, out_size] = test_affine_resample_3d.quadratic_case(1);
            outC = affine_resample_3d(vol, M, out_size, 'cubic', 0);
            outL = affine_resample_3d(vol, M, out_size, 'linear', 0);
            tc.verifyEqual(double(outL - outC), 1.5 * ones(out_size), 'AbsTol', 1e-3, ...
                'linear minus cubic should be 0.25*(1+2+3) for this quadratic');
        end

        function test_linear_unchanged(tc, nC)
            % Negative control: linear keeps its known bias of r*(1-r) per axis
            [vol, M, out_size, x, y, z] = test_affine_resample_3d.quadratic_case(nC);
            out = affine_resample_3d(vol, M, out_size, 'linear', 0);
            for c = 1:nC
                expected = single(x.^2 + 2*y.^2 + 3*z.^2 + 1.5 + 100*(c-1));
                tc.verifyEqual(out(:,:,:,c), expected, 'RelTol', single(1e-5));
            end
        end

        function test_linear_exact_on_trilinear_field(tc)
            % Negative control: linear is exact for a trilinear field
            [x, y, z] = ndgrid(1:10);
            vol = single(x + 2*y - 3*z + x.*y.*z / 10);
            M = eye(4);  M(4, 1:3) = [-2.3 -1.7 -2.5];
            [xo, yo, zo] = ndgrid((1:5) + 2.3, (1:5) + 1.7, (1:5) + 2.5);
            out = affine_resample_3d(vol, M, [5 5 5], 'linear', 0);
            tc.verifyEqual(out, single(xo + 2*yo - 3*zo + xo.*yo.*zo / 10), 'AbsTol', single(1e-3));
        end

        function test_nearest_unchanged(tc)
            % Negative control: nearest returns input voxel values only
            [vol, ~, out_size] = test_affine_resample_3d.quadratic_case(1);
            M = eye(4);  M(4, 1:3) = -2.2;   % in = out + 2.2 -> nearest index out + 2
            out = affine_resample_3d(vol, M, out_size, 'nearest', 0);
            tc.verifyEqual(out, vol(3:7, 3:7, 3:7));
        end

        function test_fill_out_of_bounds_all_methods(tc)
            % Out-of-grid voxels get the fill value for every method
            vol = ones(6, 6, 6, 'single');
            M = eye(4);  M(4, 1:3) = 20;     % output lies entirely outside the input
            for method = {'nearest', 'linear', 'cubic'}
                out = affine_resample_3d(vol, M, [3 3 3], method{1}, -7);
                tc.verifyEqual(out, -7 * ones(3, 3, 3, 'single'), method{1});
            end
        end
    end
end
