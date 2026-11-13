classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_thermal_dose < matlab.unittest.TestCase
% TEST_THERMAL_DOSE  ISO CEM43 increments against analytic values (Sapareto & Dewey 1984).
%
%   CEM43 = t[min] * R^(43 - T), R = 0.5 at/above 43 °C, 0.25 for 39–43 °C,
%   0 below 39 °C; once 43 °C was reached R stays 0.5; >= 57 °C is Inf.
%
%   results = runtests('tests/test_thermal_dose.m');


    methods (Test, TestTags = {'thermal_dose'})

        function test_one_minute_at_43_is_one_minute(tc)
            tc.verifyEqual(cem43_iso_increment(43, 43, 60), 1, 'AbsTol', 1e-12);
        end

        function test_each_degree_above_43_doubles(tc)
            tc.verifyEqual(cem43_iso_increment(44, 44, 60), 2, 'AbsTol', 1e-12);
            tc.verifyEqual(cem43_iso_increment(46, 46, 60), 8, 'AbsTol', 1e-12);
        end

        function test_each_degree_below_43_quarters(tc)
            tc.verifyEqual(cem43_iso_increment(42, 42, 60), 0.25, 'AbsTol', 1e-12);
            tc.verifyEqual(cem43_iso_increment(40, 40, 60), 0.25^3, 'AbsTol', 1e-12);
        end

        function test_below_39_contributes_nothing(tc)
            tc.verifyEqual(cem43_iso_increment(37, 37, 60), 0);
            tc.verifyEqual(cem43_iso_increment(38.9, 38.9, 600), 0);
        end

        function test_after_reaching_43_rate_stays_half(tc)
            % Cooled to 42 °C after having reached 44 °C
            tc.verifyEqual(cem43_iso_increment(42, 44, 60), 0.5, 'AbsTol', 1e-12);
        end

        function test_scales_with_step_duration(tc)
            tc.verifyEqual(cem43_iso_increment(43, 43, 30), 0.5, 'AbsTol', 1e-12);
            tc.verifyEqual(cem43_iso_increment(43, 43, 0), 0);
        end

        function test_irreversible_damage_is_inf(tc)
            tc.verifyEqual(cem43_iso_increment(57, 57, 1), Inf);
            tc.verifyEqual(cem43_iso_increment(40, 58, 1), Inf);   % reached 57 °C before
        end

        function test_elementwise_on_volumes(tc)
            T = cat(3, [37 42; 43 44], [57 39; 41 45]);
            inc = cem43_iso_increment(T, T, 60);
            tc.verifySize(inc, size(T));
            tc.verifyEqual(inc(:, :, 1), [0 0.25; 1 2], 'AbsTol', 1e-12);
            tc.verifyEqual(inc(1, 1, 2), Inf);
        end

    end

end
