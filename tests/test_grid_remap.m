classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_grid_remap < matlab.unittest.TestCase
% TEST_GRID_REMAP  Mapping fields between simulation grids through world space:
% GRID_WORLD_AFFINE, GRID_REMAP_AFFINE and REMAP_ACOUSTIC_FIELD, on synthetic
% grids with known answers (identity, FOV box offset, world shift, resolution
% change, partial overlap), every quantity, and negative controls.

    methods (Test)

        function test_identity_keeps_the_field(tc)
            A = grid_affine(1, [0 0 0]);
            [M, info] = grid_remap_affine(A, A);
            tc.verifyEqual(M, eye(4), 'AbsTol', 1e-12);
            tc.verifyTrue(info.rigid);
            tc.verifyEqual(info.shift_vox, [0 0 0], 'AbsTol', 1e-12);
            v = rand(6, 7, 8, 'single') + 1;
            [out, rinfo] = remap_acoustic_field(v, M, size(v), 'pressure');
            tc.verifyEqual(out, v, 'AbsTol', 1e-5);
            tc.verifyClass(out, 'single');
            tc.verifyEqual(rinfo.coverage, 1);
        end

        function test_fov_box_offset(tc)
            % a box stored from full-grid voxel [11 3 1] onto the full grid
            A = grid_affine(1, [0 0 0]);
            [M, info] = grid_remap_affine(A, A, [11 3 1], [1 1 1]);
            tc.verifyEqual(info.shift_vox, [10 2 0], 'AbsTol', 1e-12);
            box = reshape(1:4*5*6, [4 5 6]);
            out = remap_acoustic_field(box, M, [20 10 6], 'linear');
            tc.verifyEqual(out(11:14, 3:7, :), box, 'AbsTol', 1e-5);
            out(11:14, 3:7, :) = 0;
            tc.verifyEqual(nnz(abs(out) > 1e-5), 0);      % zero outside the box
        end

        function test_world_shift_between_grids(tc)
            % destination grid starts 2 mm further along x: source voxel 5 -> destination voxel 3
            [M, info] = grid_remap_affine(grid_affine(1, [0 0 0]), grid_affine(1, [2 0 0]));
            tc.verifyEqual(info.shift_vox, [-2 0 0], 'AbsTol', 1e-12);
            v = zeros(10, 4, 4); v(5, 2, 2) = 3;
            out = remap_acoustic_field(v, M, [10 4 4], 'pressure');
            tc.verifyEqual(out(3, 2, 2), 3, 'AbsTol', 1e-5);
            tc.verifyEqual(nnz(out > 1e-5), 1);
        end

        function test_resolution_change_keeps_intensity(tc)
            % 1 mm source onto a 0.5 mm grid covering the same extent: scale 2,
            % rigid; a smooth field keeps its integral of p^2
            A_src = grid_affine(1.0, [0 0 0]);
            A_dst = grid_affine(0.5, [0 0 0]);
            [M, info] = grid_remap_affine(A_src, A_dst);
            tc.verifyEqual(info.scale, [2 2 2], 'AbsTol', 1e-12);
            tc.verifyTrue(info.rigid);
            [x, y, z] = ndgrid(1:24, 1:24, 1:24);
            v = exp(-((x - 12.5).^2 + (y - 12.5).^2 + (z - 12.5).^2) / (2 * 3^2));
            out = remap_acoustic_field(v, M, [47 47 47], 'pressure');
            e_src = sum(v(:).^2) * 1.0^3;
            e_dst = sum(out(:).^2) * 0.5^3;
            tc.verifyEqual(e_dst, e_src, 'RelTol', 0.02);
            tc.verifyEqual(max(out(:)), max(v(:)), 'RelTol', 0.05);
        end

        function test_pressure_is_resampled_as_p_squared(tc)
            % between voxels of 1 and 3, p^2 interpolation gives sqrt(5) and
            % linear interpolation 2 (negative control for the quantity rule)
            v = ones(2, 3, 3); v(2, :, :) = 3;
            M = eye(4); M(1, 1) = 2; M(4, 1) = -1;   % source x 1 -> 1, 2 -> 3 (midpoint at 2)
            p  = remap_acoustic_field(v, M, [3 3 3], 'pressure');
            li = remap_acoustic_field(v, M, [3 3 3], 'linear');
            tc.verifyEqual(p(2, 2, 2), sqrt(5), 'AbsTol', 1e-5);
            tc.verifyEqual(li(2, 2, 2), 2, 'AbsTol', 1e-5);
        end

        function test_every_quantity(tc)
            M = grid_remap_affine(grid_affine(1, [0 0 0]), grid_affine(1, [0 0 0]));
            v = -rand(5, 5, 5) - 0.5;
            n = remap_acoustic_field(v, M, [5 5 5], 'negative');
            tc.verifyEqual(n, v, 'AbsTol', 1e-5);                       % sign kept
            c = complex(rand(5, 5, 5), rand(5, 5, 5));
            tc.verifyEqual(remap_acoustic_field(c, M, [5 5 5], 'complex'), c, 'AbsTol', 1e-5);
            m = rand(5, 5, 5) > 0.5;
            mo = remap_acoustic_field(m, M, [5 5 5], 'mask');
            tc.verifyClass(mo, 'logical');
            tc.verifyEqual(mo, m);
            tc.verifyError(@() remap_acoustic_field(v, M, [5 5 5], 'nope'), 'remap_acoustic_field:quantity');
        end

        function test_partial_overlap_coverage(tc)
            % destination grid shifted by half its x extent: half of it is covered
            [M, ~] = grid_remap_affine(grid_affine(1, [0 0 0]), grid_affine(1, [5 0 0]));
            [~, info] = remap_acoustic_field(ones(10, 4, 4), M, [10 4 4], 'pressure');
            tc.verifyEqual(info.coverage, 0.5, 'AbsTol', 1e-12);
        end

        function test_shear_is_not_rigid(tc)
            A = grid_affine(1, [0 0 0]);
            B = A; B(2, 1) = 0.3;                                         % shear
            [~, info] = grid_remap_affine(A, B);
            tc.verifyFalse(info.rigid);
        end

        function test_world_affine_from_planimg(tc)
            T1 = grid_affine(1.0, [-90 -126 -72]);
            pl = struct('inv_transf', grid_affine(0.5, [3 4 5]), 't1_header', struct('Transform', struct('T', T1)));
            tc.verifyEqual(grid_world_affine(pl), pl.inv_transf * T1, 'AbsTol', 1e-12);
            pl.t1_header = [];
            tc.verifyEmpty(grid_world_affine(pl));                         % phantom: no T1 header
            tc.verifyEmpty(grid_world_affine(struct('parameters', struct())));
            tc.verifyEmpty(grid_world_affine([]));
        end
    end
end

function A = grid_affine(dx, origin)
% voxel -> world, row-vector convention: [w 1] = [v 1] * A
A = diag([dx dx dx 1]);
A(4, 1:3) = origin;
end
