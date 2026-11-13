classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_water_baseline_box < matlab.unittest.TestCase
% TEST_WATER_BASELINE_BOX  The free-water baseline on the beam box equals the full-grid one
%
% 2D k-Wave (MATLAB, CPU, seconds): a focused arc in water on a grid much
% wider than the beam. WATER_BASELINE with grid.water_baseline_grid 'beam'
% (WATER_BASELINE_BOX) and 'full' must give the same free-water ISPPA; the box
% must keep every source point and be clearly smaller than the grid.
%
%   Run with:   results = runtests('tests/test_water_baseline_box.m');

    properties
        P; Kgrid; Source; Sensor
    end

    methods (TestMethodSetup)
        function setup_case(tc)
            res = 0.5; dx = res * 1e-3; f = 500e3; c = 1500;
            Nx = 288; Ny = 448;                         % beam along x, wide in y
            kgrid = kWaveGrid(Nx, dx, Ny, dx);
            kgrid.makeTime(c, 0.3, 2 * Nx * dx / c);
            R = 128;                                    % 64 mm curvature radius [voxels]
            apex  = [12, Ny / 2];
            focus = [12 + R, Ny / 2];
            arc = makeArc([Nx, Ny], apex, R, 2 * 60 + 1, focus);
            source = struct('p_mask', arc);
            n = nnz(arc);
            source.p = createCWSignals(kgrid.t_array, f, 1e5 * ones(n, 1), zeros(n, 1));
            sensor = struct('mask', true(Nx, Ny), 'record', {{'p_max_all'}});

            p = struct();
            p.grid = struct('resolution_mm', res, 'pml_size', 20, 'axisymmetric', 0, 'dims', [Nx, Ny]);
            p.simulation = struct('code_type', 'matlab_cpu', 'precision', 'single', 'interactive', 0, ...
                'debug', 0, 'medium', 'water');
            p.io = struct('output_affix', '_wbbox', 'dir_cache', tempdir);
            p.medium_properties.water = struct('sound_speed', c, 'density', 1000, 'alpha_coeff', 0, 'alpha_power', 1.43);
            p.transducer = struct('type', 'annular', 'trans_pos', apex, 'focus_pos', focus, ...
                'annular', struct('curv_radius_mm', R * res, 'elem_p_amp', 1e5));
            tc.P = p; tc.Kgrid = kgrid; tc.Source = source; tc.Sensor = sensor;
        end
    end

    methods (Test, TestTags = {'water_baseline'})

        function test_box_keeps_source_and_is_smaller(tc)
            [kb, sb, ~, off] = water_baseline_box(tc.P, tc.Kgrid, tc.Source, tc.Sensor);
            tc.verifyEqual(nnz(sb.p_mask), nnz(tc.Source.p_mask));
            tc.verifyLessThan(kb.Nx * kb.Ny, 0.8 * tc.Kgrid.Nx * tc.Kgrid.Ny);
            tc.verifyEqual(kb.Nt, tc.Kgrid.Nt);
            tc.verifyGreaterThanOrEqual(off, [0 0]);
        end

        function test_same_isppa_as_full_grid(tc)
            pb = tc.P; pb.grid.water_baseline_grid = 'beam';
            pf = tc.P; pf.grid.water_baseline_grid = 'full';
            i_box  = water_baseline(pb, tc.Kgrid, tc.Source, tc.Sensor);
            i_full = water_baseline(pf, tc.Kgrid, tc.Source, tc.Sensor);
            tc.verifyEqual(i_box, i_full, 'RelTol', 0.01);
        end
    end
end
