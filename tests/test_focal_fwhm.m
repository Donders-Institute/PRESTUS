classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_focal_fwhm < matlab.unittest.TestCase
% TEST_FOCAL_FWHM  Unit tests for COMPUTE_FOCAL_FWHM (−3 dB axial and lateral
% focal widths) on Gaussian foci with known widths.
%
%   Run with:   results = runtests('tests/test_focal_fwhm.m');

    methods (Test)

        function test_axis_aligned_gaussian(tc)
            % FWHM of exp(-x²/(2σ²)) is 2σ·sqrt(2 ln 2)
            sig = [3, 3, 8];
            I = gaussian3([41 41 81], [21 21 41], sig, eye(3));
            fw = compute_focal_fwhm(I, [21 21 41], [0 0 1], 0.5);
            k = 2 * sqrt(2 * log(2)) * 0.5;
            tc.verifyEqual(fw.axial_mm, sig(3) * k, 'RelTol', 0.02);
            tc.verifyEqual(fw.lateral_min_mm, sig(1) * k, 'RelTol', 0.02);
            tc.verifyEqual(fw.lateral_max_mm, sig(1) * k, 'RelTol', 0.02);
        end

        function test_anisotropic_lateral(tc)
            sig = [2, 4, 8];
            I = gaussian3([41 41 81], [21 21 41], sig, eye(3));
            fw = compute_focal_fwhm(I, [21 21 41], [0 0 -1], 1);
            k = 2 * sqrt(2 * log(2));
            tc.verifyEqual(fw.lateral_min_mm, sig(1) * k, 'RelTol', 0.02);
            tc.verifyEqual(fw.lateral_max_mm, sig(2) * k, 'RelTol', 0.02);
            tc.verifyEqual(fw.axial_mm, sig(3) * k, 'RelTol', 0.02);
        end

        function test_oblique_axis(tc)
            % Elongated focus along a tilted beam axis
            u = [1 0 1] / sqrt(2);
            R = [cross([0 1 0], u); 0 1 0; u];      % rows: local x, y, beam axis
            sig = [3, 3, 9];
            I = gaussian3([61 61 61], [31 31 31], sig, R);
            fw = compute_focal_fwhm(I, [31 31 31], u * 10, 1);
            k = 2 * sqrt(2 * log(2));
            tc.verifyEqual(fw.axial_mm, sig(3) * k, 'RelTol', 0.02);
            tc.verifyEqual(fw.lateral_min_mm, sig(1) * k, 'RelTol', 0.03);
            tc.verifyEqual(fw.lateral_max_mm, sig(1) * k, 'RelTol', 0.03);
        end

        function test_truncated_lobe_is_nan(tc)
            % Axial lobe running into the grid edge: no crossing on one side
            I = gaussian3([31 31 31], [16 16 2], [3 3 20], eye(3));
            fw = compute_focal_fwhm(I, [16 16 2], [0 0 1], 1);
            tc.verifyTrue(isnan(fw.axial_mm));
            tc.verifyTrue(isfinite(fw.lateral_min_mm));
        end

        function test_2d_grid(tc)
            [x, z] = ndgrid(1:41, 1:81);
            I = exp(-(x - 21).^2 / (2 * 3^2) - (z - 41).^2 / (2 * 8^2));
            fw = compute_focal_fwhm(I, [21 41], [0 1], 1);
            k = 2 * sqrt(2 * log(2));
            tc.verifyEqual(fw.axial_mm, 8 * k, 'RelTol', 0.02);
            tc.verifyEqual(fw.lateral_min_mm, 3 * k, 'RelTol', 0.02);
            tc.verifyEqual(fw.lateral_max_mm, fw.lateral_min_mm);
        end

    end
end

function I = gaussian3(sz, c, sig, R)
% Gaussian intensity centred at c with widths sig along the rows of R.
    [x, y, z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
    d = [x(:) - c(1), y(:) - c(2), z(:) - c(3)] * R';
    I = reshape(exp(-sum(d.^2 ./ (2 * sig.^2), 2)), sz);
end
