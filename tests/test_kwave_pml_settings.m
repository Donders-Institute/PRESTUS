classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_kwave_pml_settings < matlab.unittest.TestCase
% TEST_KWAVE_PML_SETTINGS  KWAVE_PML_SETTINGS, the PML resolution shared by
% ACOUSTIC_WRAPPER and REVERSE_POINT_SOURCE_FIELD: numeric sizes pass through,
% 'auto' honours grid.pml_size_range (and WSWA on axisymmetric grids), and
% grid.pml_alpha is returned.
%
%   Run with:   results = runtests('tests/test_kwave_pml_settings.m');

    methods (Test)

        function test_numeric_size_passes_through(tc)
            p.grid = struct('pml_size', 12);
            [pml, alpha] = kwave_pml_settings(p, kWaveGrid(64, 1e-3, 64, 1e-3));
            tc.verifyEqual(pml, 12);
            tc.verifyEmpty(alpha);
        end

        function test_auto_uses_range_and_alpha(tc)
            kgrid = kWaveGrid(70, 1e-3, 70, 1e-3, 90, 1e-3);
            p.grid = struct('pml_size', "auto");
            tc.verifyEqual(kwave_pml_settings(p, kgrid), getOptimalPMLSize(kgrid));
            p.grid.pml_size_range = [16 24];
            p.grid.pml_alpha = 1.5;
            [pml, alpha] = kwave_pml_settings(p, kgrid);
            tc.verifyEqual(pml, getOptimalPMLSize(kgrid, [16 24]));
            tc.verifyTrue(all(pml >= 16 & pml <= 24));
            tc.verifyEqual(alpha, 1.5);
        end

        function test_axisymmetric_uses_wswa(tc)
            kgrid = kWaveGrid(90, 1e-3, 50, 1e-3);
            p.grid = struct('pml_size', 'auto', 'axisymmetric', 1);
            tc.verifyEqual(kwave_pml_settings(p, kgrid), getOptimalPMLSize(kgrid, [], 'WSWA'));
            p.grid.axisymmetric = 0;   % 2D Cartesian: no WSWA
            tc.verifyEqual(kwave_pml_settings(p, kgrid), getOptimalPMLSize(kgrid, []));
        end

    end
end
