classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_localite_export < matlab.unittest.TestCase
% TEST_LOCALITE_EXPORT  Heuristic-placement coordinates match the Localite planning image.
%
% The heuristic placement exports positions for Localite in the coordinates
% of the planning image CANONICAL_AFFINE_TRANSFORM writes (GRID_TO_LOCALITE_MM).
% These tests check, for a T1 with a rotated and shifted affine, that the
% exported coordinate of a voxel is the world coordinate the planning image
% assigns to that voxel, and that it maps back to the same voxel. Whether
% Localite displays the NIfTI world coordinate is outside these tests.
%
%   results = runtests('tests/test_localite_export.m');

    properties
        Dir
        T1          % niftiinfo of the synthetic T1
        Planning    % niftiinfo of the planning image
        Img
    end

    properties (TestParameter)
        voxel = struct('first', [1 1 1], 'inner', [11 21 31], 'last', [20 30 40], ...
                       'edge_mixed', [7 3 25]);
    end

    methods (TestMethodSetup)
        function make_images(tc)
            tc.Dir = tempname;
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));

            tc.Img = uint16(reshape(1:20*30*40, [20 30 40]));   % every voxel unique
            f_t1 = fullfile(tc.Dir, 't1.nii');
            niftiwrite(tc.Img, f_t1);
            info = niftiinfo(f_t1);

            % Oblique T1: anisotropic voxels, rotation about z, non-zero origin
            pixdim = [0.9 1.0 1.2];
            a = deg2rad(30);
            R = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
            T = eye(4);
            T(1:3, 1:3) = diag(pixdim) * R';            % row-vector convention
            T(4, 1:3) = [-90 -110 -70];
            info.PixelDimensions = pixdim;
            info.Transform = affine3d(T);
            info.TransformName = 'Sform';   % else niftiwrite drops Transform
            niftiwrite(tc.Img, f_t1, info);
            tc.T1 = niftiinfo(f_t1);
            tc.assertEqual(tc.T1.Transform.T, T, 'AbsTol', 1e-5, ...
                'Fixture T1 did not keep its oblique affine');

            f_plan = fullfile(tc.Dir, 'planning.nii');
            evalc('canonical_affine_transform(f_t1, f_plan);');
            tc.Planning = niftiinfo(f_plan);
        end
    end

    methods (Test, TestTags = {'placement'})

        function test_export_is_planning_world_coordinate(tc, voxel)
            exported = grid_to_localite_mm(voxel, tc.T1);
            world = grid_to_ras(voxel, tc.Planning);
            tc.verifyEqual(exported, world, 'AbsTol', 0.005 + 1e-9);
        end

        function test_export_maps_back_to_same_voxel(tc, voxel)
            exported = grid_to_localite_mm(voxel, tc.T1);
            tc.verifyEqual(ras_to_grid(exported, tc.Planning), voxel);
        end

        function test_planning_image_keeps_voxel_order(tc)
            % Same voxel index = same anatomy in T1 and planning image
            tc.verifyEqual(niftiread(tc.Planning), tc.Img);
        end

        function test_known_voxel(tc)
            % Voxel (11,21,31) at 0.9 x 1.0 x 1.2 mm lies at (9, 20, 36) mm
            tc.verifyEqual(grid_to_localite_mm([11 21 31], tc.T1), [9 20 36], 'AbsTol', 1e-9);
        end

        function test_planning_affine_is_zero_origin_diagonal(tc)
            tc.verifyEqual(tc.Planning.Transform.T, diag([0.9 1.0 1.2 1]), 'AbsTol', 1e-6);
        end
    end
end
