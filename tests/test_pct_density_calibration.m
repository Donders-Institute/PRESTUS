classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_pct_density_calibration < matlab.unittest.TestCase
% TEST_PCT_DENSITY_CALIBRATION  pct.density_calibration for mapping_density 'k-plan'.
%
% Empty selects k-Plan's default table (results unchanged from the former
% hard-coded table); a shipped table is selected by name, any CSV by path;
% bad selections and tables error; treeby without its 80 kVp table warns.
%
% Run with:   results = runtests('tests/test_pct_density_calibration.m');

    properties
        TmpDir
    end

    methods (TestMethodSetup)
        function setup(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
            tc.addTeardown(@() rmdir(tc.TmpDir, 's'));
        end
    end

    methods (Test)
        function test_default_is_kplan_table(tc)
            [hu, rho, src] = medium_pct_density_calibration(params(''));
            tc.verifyEqual(hu, [-990 60 1000 1950]);
            tc.verifyEqual(rho, [1.2 1060 1530 2150]);
            tc.verifySubstring(src, 'kplan_default.csv');
            % same densities as the former hard-coded table
            x = [-500 200 1000 1500 2500];
            d = run_kplan(params(''), x);
            % the mapping itself; the water floor is a limit (pct.limits.density_min, PCT_MAP_PROPERTIES)
            tc.verifyEqual(d, fit_pairwiselinear(x, [-990 60 1000 1950], [1.2 1060 1530 2150]), 'AbsTol', 1e-9);
        end

        function test_shipped_table_by_name(tc)
            [hu, rho] = medium_pct_density_calibration(params('ucl_lowdose_80kvp'));
            tc.verifyEqual(rho([1 end]), [1.2 2150], 'AbsTol', 1e-3);
            tc.verifyEqual(hu(end), 2614.7998, 'AbsTol', 1e-3);
            % the same density needs more HU at 80 kVp than in the default table
            tc.verifyLessThan(run_kplan(params('ucl_lowdose_80kvp'), 1500), run_kplan(params(''), 1500));
            [~, rho] = medium_pct_density_calibration(params('cph2025'));
            tc.verifyEqual(rho(end), 2150, 'AbsTol', 1e-3);
        end

        function test_csv_path_unsorted_rows(tc)
            f = fullfile(tc.TmpDir, 'my_ct.csv');
            writelines(["# HU, density", "1000, 1600", "0, 1000", "2000, 2200"], f);
            [hu, rho, src] = medium_pct_density_calibration(params(f));
            tc.verifyEqual(hu, [0 1000 2000]); tc.verifyEqual(rho, [1000 1600 2200]);
            tc.verifyEqual(src, f);
            tc.verifyEqual(run_kplan(params(f), 500), 1300, 'AbsTol', 1e-9);
        end

        function test_treeby_pairing_warns(tc)
            p = params(''); p.pct.skull_mapping = 'treeby';
            tc.verifyWarning(@() run_kplan(p, 1000), 'prestus:pct:calibrationMismatch');
            p.pct.density_calibration = 'ucl_lowdose_80kvp';
            tc.verifyWarningFree(@() run_kplan(p, 1000));
            p.pct.skull_mapping = 'miscouridou'; p.pct.density_calibration = '';
            tc.verifyWarningFree(@() run_kplan(p, 1000));
        end

        function test_errors(tc)
            tc.verifyError(@() medium_pct_density_calibration(params('no_such_table')), 'prestus:pct:densityCalibration');
            f = fullfile(tc.TmpDir, 'bad.csv');
            writelines(["0, 1000", "1000, 900"], f);              % density decreasing
            tc.verifyError(@() medium_pct_density_calibration(params(f)), 'prestus:pct:densityCalibration');
        end
    end
end

function p = params(sel)
p = struct();
p.pct.density_calibration = sel;
p.simulation.debug = 0;
p.medium_properties.water.density = 1000;
end

function d = run_kplan(p, hu)
d = medium_pct_density(p, zeros(size(hu)), hu, 1:numel(hu), 'k-plan');
end
