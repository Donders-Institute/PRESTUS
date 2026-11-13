classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_transform < matlab.unittest.TestCase
% TEST_TRANSFORM  Unit tests for functions/transform/
%
%   Run with:   results = runtests('tests/test_transform.m');

    methods (Test, TestTags = {'ras_to_grid'})

        % Convention: PRESTUS voxels are 1-based MATLAB indices into the
        % niftiread array; niftiinfo's Transform.T maps 0-based indices.

        function test_ras_to_grid_identity_affine(tc)
            % Identity affine: RAS == 0-based index, so the 1-based index is RAS + 1
            hdr.Transform.T = eye(4);
            grid = ras_to_grid([5, 6, 7], hdr);
            tc.verifyEqual(grid, [6, 7, 8]);
        end

        function test_ras_to_grid_translation(tc)
            % Affine = shift by [10,20,30]
            T = eye(4);
            T(1:3, 4) = [10; 20; 30];
            hdr.Transform.T = T';  % niftiinfo stores T transposed
            grid = ras_to_grid([15, 26, 37], hdr);
            tc.verifyEqual(grid, [6, 7, 8]);
        end

        function test_ras_to_grid_column_input(tc)
            hdr.Transform.T = eye(4);
            grid = ras_to_grid([1; 2; 3], hdr);
            tc.verifyEqual(grid, [2, 3, 4]);
        end

        function test_grid_to_ras_inverts_ras_to_grid(tc)
            A = [-1.1 0.05 0 90.3; 0.02 1.2 0.1 -126.7; 0 -0.08 0.9 -72.2; 0 0 0 1];
            hdr.Transform.T = A';
            vox = [11 21 31; 3 4 5];
            ras = grid_to_ras(vox, hdr);
            tc.verifyEqual(ras(1,:), (A(1:3,:) * [10; 20; 30; 1])', 'AbsTol', 1e-9);
            tc.verifyEqual(ras_to_grid(ras(2,:), hdr), [3 4 5]);
            tc.verifyEqual(grid_to_ras([11; 21; 31], hdr), ras(1,:), 'AbsTol', 1e-9);
        end

        function test_convention_matches_niftiread(tc)
            % A single bright voxel: its RAS centre must map to its array index
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            A = [-1.1 0.05 0 90.3; 0.02 1.2 0.1 -126.7; 0 -0.08 0.9 -72.2; 0 0 0 1];
            img = zeros(40, 50, 60, 'single');
            img(11, 21, 31) = 1;
            f = fullfile(tmp, 'dot.nii');
            niftiwrite(img, f);
            info = niftiinfo(f);
            info.Transform = affine3d(A');
            info.TransformName = 'Sform';
            niftiwrite(img, f, info);
            hdr = niftiinfo(f);
            vol = niftiread(hdr);
            ras = grid_to_ras([11 21 31], hdr);
            g = ras_to_grid(ras, hdr);
            tc.verifyEqual(g, [11 21 31]);
            tc.verifyEqual(vol(g(1), g(2), g(3)), single(1));
            tc.verifyEqual(transform_coordinates(struct(), ras, 'ras_plus', 'grid', hdr), [11 21 31]);
            tc.verifyEqual(transform_coordinates(struct(), [11 21 31], 'grid', 'ras_plus', hdr), ras, ...
                'AbsTol', 1e-9);
        end

        function test_nifti_affine_from_1based(tc)
            % [ras 1] = [vox1 1] * T1  ==>  [ras 1] = [vox0 1] * nifti_affine_from_1based(T1)
            T1 = [2 0 0 0; 0 3 0 0; 0 0 4 0; 5 6 7 1];
            Tn = nifti_affine_from_1based(T1);
            tc.verifyEqual([2 3 4 1] * Tn, [3 4 5 1] * T1, 'AbsTol', 1e-12);
        end

        function test_ensure_ras_plus_keeps_world_position(tc)
            % Reorienting (flip + axis permutation) must not move a voxel in RAS
            A = [0 0 -2 50; -1 0 0 30; 0 1.5 0 -40; 0 0 0 1];   % dims: -Y, +Z, -X
            hdr.Transform = affine3d(A');
            hdr.PixelDimensions = [1 1.5 2];
            img = zeros(10, 12, 14, 'single');
            img(3, 5, 7) = 1;
            ras = grid_to_ras([3 5 7], hdr);
            [img2, hdr2] = ensure_ras_plus(img, hdr);
            [i, j, k] = ind2sub(size(img2), find(img2 == 1));
            tc.verifyEqual(grid_to_ras([i j k], hdr2), ras, 'AbsTol', 1e-9);
            tc.verifyEqual(ras_to_grid(ras, hdr2), [i j k]);
        end

        function test_ensure_ras_plus_vox_map(tc)
            % vox_map carries input-grid voxel indices to the reoriented grid
            A = [0 0 -2 50; -1 0 0 30; 0 1.5 0 -40; 0 0 0 1];   % dims: -Y, +Z, -X
            hdr.Transform = affine3d(A');
            hdr.PixelDimensions = [1 1.5 2];
            img = zeros(10, 12, 14, 'single');
            img(3, 5, 7) = 1;
            [img2, ~, vox_map] = ensure_ras_plus(img, hdr);
            [i, j, k] = ind2sub(size(img2), find(img2 == 1));
            tc.verifyEqual(affine_apply_pts([3 5 7], vox_map), [i j k], 'AbsTol', 1e-9);

            % already RAS+: identity
            hdr.Transform = affine3d(eye(4));
            hdr.PixelDimensions = [1 1 1];
            [~, ~, vox_map] = ensure_ras_plus(img, hdr);
            tc.verifyEqual(vox_map, eye(4));
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


    % convert_final_to_MNI_simnibs
    methods (Test, TestTags = {'mni'})

        function test_native_mni_warp_keeps_voxel_position(tc)
            tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            % Identity warp on a 20^3 MNI grid (origin -10 mm); input grid
            % with origin -5 mm. A voxel at world (2,3,4) mm must land on
            % MNI voxel (13,14,15), 1-based.
            m2m = fullfile(tmp, 'm2m_sub-001');
            mkdir(fullfile(m2m, 'toMNI'));
            n = 20;
            [i, j, k] = ndgrid(0:n-1, 0:n-1, 0:n-1);
            warp = single(cat(4, i - 10, j - 10, k - 10));
            write_nii_1mm(fullfile(m2m, 'toMNI', 'MNI2Conform_nonl'), warp, -10);
            write_nii_1mm(fullfile(m2m, 'toMNI', 'final_tissues_MNI'), zeros(n, n, n, 'uint8'), -10);
            vol = zeros(n, n, n, 'single');
            vol(8, 9, 10) = 1;
            in = write_nii_1mm(fullfile(tmp, 'in'), vol, -5);
            out = fullfile(tmp, 'out_MNI.nii.gz');
            convert_final_to_MNI_simnibs(in, m2m, out, struct(), ...
                'method', 'native', 'interpolation_order', 0);
            res = niftiread(out);
            [~, idx] = max(res(:));
            [a, b, c] = ind2sub(size(res), idx);
            tc.verifyEqual([a b c], [13 14 15]);
        end

    end
end

function f = write_nii_1mm(base, vol, origin_mm)
% 1 mm isotropic NIfTI whose 0-based voxel [0 0 0] sits at origin_mm on each axis.
    niftiwrite(vol, base);
    info = niftiinfo([base '.nii']);
    T = eye(4);
    T(4, 1:3) = origin_mm;
    info.Transform = affine3d(T);
    info.TransformName = 'Sform';   % 'None' would drop the affine on write
    niftiwrite(vol, base, info, 'Compressed', true);
    delete([base '.nii']);
    f = [base '.nii.gz'];
end
