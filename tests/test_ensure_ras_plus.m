classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_ensure_ras_plus < matlab.unittest.TestCase
% TEST_ENSURE_RAS_PLUS  Reordering a volume to RAS+ keeps every voxel at its
% world position and value, also for an oblique affine; an image without
% rotation gets a positive diagonal affine. Pure logic: no data.
%
%   results = runtests('tests/test_ensure_ras_plus.m');

    properties (TestParameter)
        % voxel axes as world directions (columns, before rotation) and rotation about z [deg]
        layout = struct( ...
            'las_oblique',      struct('P', diag([-1 1 1]),        'rot', 10), ...
            'permuted_oblique', struct('P', [0 0 -1; 1 0 0; 0 1 0], 'rot', -12), ...
            'las_straight',     struct('P', diag([-1 1 1]),        'rot', 0), ...
            'permuted_straight',struct('P', [0 0 -1; 1 0 0; 0 1 0], 'rot', 0));
    end

    methods (Test)
        function test_world_position_and_value_kept(tc, layout)
            [img, hdr] = make_volume(layout.P, layout.rot);
            [out, hdr_out, vox_map] = ensure_ras_plus(img, hdr);
            [i, j, k] = ndgrid(1:size(img, 1), 1:size(img, 2), 1:size(img, 3));
            v_in  = [i(:), j(:), k(:)];
            v_out = round(affine_apply_pts(v_in, vox_map));
            w_in  = [v_in - 1,  ones(size(v_in, 1), 1)] * hdr.Transform.T;
            w_out = [v_out - 1, ones(size(v_out, 1), 1)] * hdr_out.Transform.T;
            tc.verifyEqual(w_out, w_in, 'AbsTol', 1e-9, 'every voxel keeps its world position');
            tc.verifyEqual(out(sub2ind(size(out), v_out(:, 1), v_out(:, 2), v_out(:, 3))), img(:));
            tc.verifyEqual(hdr_out.ImageSize, size(out));
        end

        function test_output_axes_point_right_anterior_superior(tc, layout)
            [img, hdr] = make_volume(layout.P, layout.rot);
            [~, hdr_out] = ensure_ras_plus(img, hdr);
            R = hdr_out.Transform.T(1:3, 1:3)';
            tc.verifyTrue(all(diag(R) > 0));
            [~, dominant] = max(abs(R), [], 1);
            tc.verifyEqual(dominant, [1 2 3]);
            tc.verifyEqual(sqrt(sum(R.^2, 1)), hdr_out.PixelDimensions, 'AbsTol', 1e-12);
            if layout.rot == 0
                tc.verifyEqual(R, diag(hdr_out.PixelDimensions), 'AbsTol', 1e-12);
            end
        end

        function test_ras_plus_input_is_unchanged(tc)
            % Negative control: an image already RAS+ (oblique or not) is returned as is
            [img, hdr] = make_volume(eye(3), 8);
            [out, hdr_out, vox_map] = ensure_ras_plus(img, hdr);
            tc.verifyEqual(out, img);
            tc.verifyEqual(hdr_out.Transform.T, hdr.Transform.T);
            tc.verifyEqual(vox_map, eye(4));
        end
    end
end

function [img, hdr] = make_volume(P, rot_deg)
% A [5 6 7] volume with distinct values, anisotropic voxels, and the
% column-vector affine Rz(rot) * P * diag(pix) + t
    sz  = [5 6 7];
    pix = [1.0 1.2 0.9];
    img = reshape(1:prod(sz), sz);
    c = cosd(rot_deg); s = sind(rot_deg);
    Rz = [c -s 0; s c 0; 0 0 1];
    A = eye(4);
    A(1:3, 1:3) = Rz * P * diag(pix);
    A(1:3, 4)   = [-40; 25; -10];
    hdr = struct('Transform', affine3d(A'), 'PixelDimensions', pix, 'ImageSize', sz);
end
