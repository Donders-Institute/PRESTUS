classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_ras_round_trip < matlab.unittest.TestCase
% TEST_RAS_ROUND_TRIP  RAS <-> voxel against the NIfTI header itself (audit item I1).
%
% Writes a NIfTI with an oblique, non-isotropic sform, reads it back with
% niftiinfo and takes the ground truth from the raw header rows
% (srow_x/y/z), which by the NIfTI-1 standard map 0-based voxel indices to
% world mm. RAS_TO_GRID must return the 1-based MATLAB index of that voxel,
% GRID_TO_RAS its inverse, both must round to the nearest voxel centre, and
% GRID_TO_RAS_AFFINE (identity planning transform) must agree with
% GRID_TO_RAS. test_transform checks the same convention with a bright voxel
% read through niftiread; this test is independent of MATLAB's Transform
% field.
%
%   results = runtests('tests/test_ras_round_trip.m');

    properties
        Hdr     % niftiinfo of the written file
        S       % 3x4 raw sform rows (0-based voxel -> mm)
        Dims = [30 40 25]
    end

    methods (TestClassSetup)
        function write_nifti(tc)
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            Rx = rotx_deg(17); Rz = rotz_deg(-11); Ry = roty_deg(6);
            A = eye(4);
            A(1:3, 1:3) = Rz * Ry * Rx * diag([0.9 1.2 2.0]);
            A(1:3, 4) = [-92.5; 118.25; -60.75];
            img = zeros(tc.Dims, 'single');
            f = fullfile(tmp, 'oblique.nii');
            niftiwrite(img, f);
            info = niftiinfo(f);
            info.Transform = affine3d(A');
            info.TransformName = 'Sform';
            niftiwrite(img, f, info);
            tc.Hdr = niftiinfo(f);
            r = tc.Hdr.raw;
            tc.S = double([r.srow_x; r.srow_y; r.srow_z]);
            tc.assertGreaterThan(double(r.sform_code), 0, 'sform not written');
            tc.assertEqual(tc.S, A(1:3, :), 'AbsTol', 1e-4, 'sform stored differently from written');
        end
    end

    methods (Test, TestTags = {'ras_round_trip'})

        function test_ras_to_grid_returns_1based_index(tc)
            for v = test_voxels(tc.Dims)'
                world = tc.S * [v - 1; 1];                       % NIfTI: 0-based
                tc.verifyEqual(ras_to_grid(world', tc.Hdr), v', sprintf('voxel %s', mat2str(v')));
            end
        end

        function test_grid_to_ras_matches_header(tc)
            for v = test_voxels(tc.Dims)'
                tc.verifyEqual(grid_to_ras(v', tc.Hdr), (tc.S * [v - 1; 1])', 'AbsTol', 1e-6, ...
                    sprintf('voxel %s', mat2str(v')));
            end
        end

        function test_round_trip_is_identity(tc)
            for v = test_voxels(tc.Dims)'
                tc.verifyEqual(ras_to_grid(grid_to_ras(v', tc.Hdr), tc.Hdr), v');
            end
        end

        function test_rounds_to_nearest_voxel_centre(tc)
            % Offsets below half a voxel stay, above half a voxel move one index
            v = [12 17 9];
            c = grid_to_ras(v, tc.Hdr);
            step_i = tc.S(:, 1)';                               % world step of +1 in i
            step_k = tc.S(:, 3)';
            tc.verifyEqual(ras_to_grid(c + 0.4 * step_i - 0.4 * step_k, tc.Hdr), v);
            tc.verifyEqual(ras_to_grid(c + 0.6 * step_i, tc.Hdr), v + [1 0 0]);
            tc.verifyEqual(ras_to_grid(c - 0.6 * step_k, tc.Hdr), v - [0 0 1]);
        end

        function test_column_and_row_input_agree(tc)
            w = grid_to_ras([5 6 7], tc.Hdr);
            tc.verifyEqual(ras_to_grid(w', tc.Hdr), ras_to_grid(w, tc.Hdr));
        end

        function test_grid_to_ras_affine_agrees(tc)
            planimg = struct('inv_transf', eye(4), 't1_header', tc.Hdr);
            T = grid_to_ras_affine(planimg);
            for v = test_voxels(tc.Dims)'
                w = [v' 1] * T;
                tc.verifyEqual(w(1:3), grid_to_ras(v', tc.Hdr), 'AbsTol', 1e-6);
            end
        end

    end

end

% =========================================================================

function V = test_voxels(dims)
% Corners and interior points, 1-based.
    V = [1 1 1; dims; dims(1) 1 1; 1 dims(2) 1; 1 1 dims(3); 7 13 21; 15 20 12];
end

function R = rotx_deg(a)
    R = [1 0 0; 0 cosd(a) -sind(a); 0 sind(a) cosd(a)];
end

function R = roty_deg(a)
    R = [cosd(a) 0 sind(a); 0 1 0; -sind(a) 0 cosd(a)];
end

function R = rotz_deg(a)
    R = [cosd(a) -sind(a) 0; sind(a) cosd(a) 0; 0 0 1];
end
