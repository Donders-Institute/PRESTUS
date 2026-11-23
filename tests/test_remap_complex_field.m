classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_remap_complex_field < matlab.unittest.TestCase
% TEST_REMAP_COMPLEX_FIELD  Complex fields resampled between grids keep their phase
% (remap_complex_field): a plane wave along the beam axis is exact after removing
% the carrier, while plain interpolation loses amplitude between the samples.
% Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_remap_complex_field.m');

    methods (Test)

        function test_plane_wave_rotated_grid(tc)
            % 6 points per wavelength, output grid rotated by 30 deg and shifted by
            % half a voxel: demodulated resampling is exact in the interior
            ppw = 6; k = 2 * pi / ppw;
            sz = [40 40 40];
            n = [0.3 0.2 1]; n = n / norm(n);
            o = [20 20 1];
            [x, y, z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
            P = 1.7 * exp(-1i * k * ((x - o(1)) * n(1) + (y - o(2)) * n(2) + (z - o(3)) * n(3)) + 0.4i);
            M = rot_shift(30, [0.5 0.5 0.5], sz / 2);
            [Pd, info] = remap_complex_field(P, M, sz, KPerVoxel=k, Axis=n, Origin=o);
            Pp = remap_complex_field(P, M, sz, Demodulate=false);
            tc.verifyTrue(info.demodulated);
            % analytic field at each output voxel's source position
            Minv = inv(M);
            [xo, yo, zo] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
            v = [xo(:), yo(:), zo(:), ones(numel(xo), 1)] * Minv;
            truth = reshape(1.7 * exp(-1i * k * ((v(:, 1:3) - o) * n(:)) + 0.4i), sz);
            inside = reshape(all(v(:, 1:3) >= 3 & v(:, 1:3) <= sz - 2, 2), sz);
            err_d = max(abs(Pd(inside) - truth(inside)));
            err_p = max(abs(Pp(inside) - truth(inside)));
            tc.verifyLessThan(err_d, 1e-6);   % AFFINE_RESAMPLE_3D interpolates in single precision
            tc.verifyGreaterThan(err_p, 0.05 * 1.7, 'plain interpolation loses amplitude at 6 PPW');
        end

        function test_identity_is_unchanged(tc)
            P = complex(randn(8, 9, 10), randn(8, 9, 10));
            out = remap_complex_field(P, eye(4), [8 9 10], KPerVoxel=0.5, Axis=[1 0 0], Origin=[1 1 1]);
            tc.verifyEqual(out, P, 'AbsTol', 1e-6);
        end

        function test_coherent_sum_after_remap(tc)
            % Two opposed plane waves on a shifted grid: the standing-wave pattern
            % of the sum survives when each is remapped along its own axis
            ppw = 6; k = 2 * pi / ppw; sz = [12 12 48];
            [x, y, z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
            P1 = exp(-1i * k * z);
            P2 = exp(1i * k * z);
            M = eye(4); M(4, 1:3) = [0 0 0.5];
            R1 = remap_complex_field(P1, M, sz, KPerVoxel=k, Axis=[0 0 1], Origin=[0 0 0]);
            R2 = remap_complex_field(P2, M, sz, KPerVoxel=k, Axis=[0 0 -1], Origin=[0 0 0]);
            truth = 2 * cos(k * (z - 0.5));
            in = z >= 3 & z <= sz(3) - 2;
            tc.verifyLessThan(max(abs(abs(R1(in) + R2(in)) - abs(truth(in)))), 1e-6);
            clear x y
        end

        function test_wavenumber_required(tc)
            tc.verifyError(@() remap_complex_field(ones(4, 4, 4), eye(4), [4 4 4]), ...
                'remap_complex_field:wavenumber');
        end

    end
end

function M = rot_shift(deg, shift, centre)
% Rotation about z through centre, then a shift (row-vector convention)
a = deg2rad(deg);
R = [cos(a) sin(a) 0; -sin(a) cos(a) 0; 0 0 1];
T1 = eye(4); T1(4, 1:3) = -centre;
T2 = eye(4); T2(1:3, 1:3) = R;
T3 = eye(4); T3(4, 1:3) = centre + shift;
M = T1 * T2 * T3;
end
