classdef test_transform < matlab.unittest.TestCase
% TEST_TRANSFORM  Unit tests for functions/transform/
%
%   Run with:   results = runtests('tests/test_transform.m');

    methods (Test, TestTags = {'ras_to_grid'})

        function test_ras_to_grid_identity_affine(tc)
            % With identity affine, RAS == grid (homogeneous coords)
            hdr.Transform.T = eye(4);
            grid = ras_to_grid([5, 6, 7], hdr);
            tc.verifyEqual(grid, [5, 6, 7]);
        end

        function test_ras_to_grid_translation(tc)
            % Affine = shift by [10,20,30]
            T = eye(4);
            T(1:3, 4) = [10; 20; 30];
            hdr.Transform.T = T';  % niftiinfo stores T transposed
            grid = ras_to_grid([15, 26, 37], hdr);
            tc.verifyEqual(grid, [5, 6, 7]);
        end

        function test_ras_to_grid_column_input(tc)
            hdr.Transform.T = eye(4);
            grid = ras_to_grid([1; 2; 3], hdr);
            tc.verifyEqual(grid, [1, 2, 3]);
        end

    end

    methods (Test, TestTags = {'axisymmetric'})

        function test_radial_expand_output_size(tc)
            % radialExpand2DTo3D(data2D) with [C x Nr] input
            % returns [2*Nr x 2*Nr x C]
            C = 10; Nr = 20;
            data2d = rand(C, Nr);
            vol3d = radialExpand2DTo3D(data2d);
            tc.verifyEqual(size(vol3d), [2*Nr, 2*Nr, C]);
        end

        function test_radial_expand_symmetry(tc)
            % Output should be symmetric: vol3d(i,j,c) == vol3d(j,i,c)
            C = 2; Nr = 10;
            data2d = rand(C, Nr);
            vol3d = radialExpand2DTo3D(data2d);
            tc.verifyEqual(vol3d(:,:,1), vol3d(:,:,1)', 'AbsTol', 1e-10);
        end

    end

end
