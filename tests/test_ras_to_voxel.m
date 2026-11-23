classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_ras_to_voxel < matlab.unittest.TestCase
% TEST_RAS_TO_VOXEL  RAS mm to 1-based voxel indices: one implementation, several points.
%
% RAS_TO_VOXEL is the single place that inverts a NIfTI affine; RAS_TO_GRID
% and the neuronav conversions call it. Checked on axis-aligned, oblique,
% and flipped affines, against GRID_TO_RAS (its inverse) and RAS_TO_GRID.
%
%   results = runtests('tests/test_ras_to_voxel.m');

    properties (TestParameter)
        affine = struct( ...
            'identity', eye(4), ...
            'scaled_shifted', [0.5 0 0 0; 0 0.75 0 0; 0 0 1.25 0; -90 -120 -70 1], ...
            'oblique', oblique_affine(), ...
            'flipped_x', [-1 0 0 0; 0 1 0 0; 0 0 1 0; 90 -126 -72 1]);
    end

    methods (Test)

        function test_inverts_grid_to_ras(tc, affine)
            hdr = struct('Transform', struct('T', affine));
            vox = [1 1 1; 11.25 21.5 31.75; 64 32 16];
            ras = grid_to_ras(vox, hdr);
            tc.verifyEqual(ras_to_voxel(ras, hdr), vox, 'AbsTol', 1e-9);
            tc.verifyEqual(ras_to_voxel(ras, affine), vox, 'AbsTol', 1e-9);
        end

        function test_round_matches_ras_to_grid(tc, affine)
            hdr = struct('Transform', struct('T', affine));
            ras = grid_to_ras([11.3 21.6 31.2; 5.49 6.51 7.5], hdr);
            for i = 1:size(ras, 1)
                tc.verifyEqual(ras_to_voxel(ras(i,:), hdr, Round=true), ras_to_grid(ras(i,:), hdr));
            end
        end

        function test_ras_to_grid_rejects_several_points(tc)
            hdr = struct('Transform', struct('T', eye(4)));
            tc.verifyError(@() ras_to_grid([1 2 3; 4 5 6; 7 8 9], hdr), 'PRESTUS:transform:onePoint');
        end

        function test_rejects_column_vector_affine(tc)
            % Negative control: a column-vector (transposed) affine is refused,
            % not silently applied the wrong way round.
            T = [1 0 0 10; 0 1 0 20; 0 0 1 30; 0 0 0 1];
            tc.verifyError(@() ras_to_voxel([0 0 0], T), 'PRESTUS:transform:badAffine');
        end

    end

end

function T = oblique_affine()
% 10 and 20 degree rotations about z and x, 0.8 mm voxels, an offset
    a = deg2rad(10); b = deg2rad(20);
    Rz = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
    Rx = [1 0 0; 0 cos(b) -sin(b); 0 sin(b) cos(b)];
    M = Rz * Rx * diag([0.8 0.8 0.8]);
    T = [M' zeros(3, 1); -80 -110 -60 1];
end
