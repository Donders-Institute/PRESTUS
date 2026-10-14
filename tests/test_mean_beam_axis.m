classdef test_mean_beam_axis < matlab.unittest.TestCase
% TEST_MEAN_BEAM_AXIS  Unit tests for MEAN_BEAM_AXIS (grid.mode = 'focal_axis'
% orientation for several transducers).
%
%   Run with:   results = runtests('tests/test_mean_beam_axis.m');

    methods (Test, TestTags = {'grid'})

        function test_single_beam(tc)
            u = [1 2 2] / 3;
            tc.verifyEqual(mean_beam_axis(u), u, 'AbsTol', 1e-12);
        end

        function test_converging_beams_bisector(tc)
            % same side, 60 deg apart: bisector, 30 deg from each
            u = [sind(30) 0 cosd(30); -sind(30) 0 cosd(30)];
            [a, n] = mean_beam_axis(u);
            tc.verifyEqual(a, [0 0 1], 'AbsTol', 1e-12);
            tc.verifyEqual(n, cosd(30), 'AbsTol', 1e-12);
        end

        function test_opposed_beams_do_not_cancel(tc)
            % exactly opposed: the shared axis, oriented along beam 1
            u = [0 0 1; 0 0 -1];
            [a, n] = mean_beam_axis(u);
            tc.verifyEqual(a, [0 0 1], 'AbsTol', 1e-12);
            tc.verifyEqual(n, 1, 'AbsTol', 1e-12);
        end

        function test_opposed_tilted_beams_split_the_angle(tc)
            % beam 2 is 40 deg off antiparallel: each beam ends 20 deg from the axis
            u = [0 0 1; sind(40) 0 -cosd(40)];
            a = mean_beam_axis(u);
            tc.verifyEqual(acosd(a * u(1,:)'), 20, 'AbsTol', 1e-9);
            tc.verifyEqual(acosd(abs(a * u(2,:)')), 20, 'AbsTol', 1e-9);
            tc.verifyGreaterThan(a(3), 0);   % points along beam 1
        end

        function test_two_beams_never_fall_back(tc)
            % perpendicular beams are the worst case for two: norm cos(45 deg)
            u = [0 0 1; 1 0 0];
            [~, n] = mean_beam_axis(u);
            tc.verifyEqual(n, cosd(45), 'AbsTol', 1e-12);
            tc.verifyWarningFree(@() mean_beam_axis(u));
        end

        function test_spread_beams_fall_back_to_beam_1(tc)
            % three mutually perpendicular beams: norm 1/sqrt(3) < 0.6
            u = eye(3);
            a = tc.verifyWarning(@() mean_beam_axis(u, 0.6), 'PRESTUS:gridAlignment:ambiguousFocalAxis');
            tc.verifyEqual(a, [1 0 0], 'AbsTol', 1e-12);
        end

        function test_unnormalised_input(tc)
            a = mean_beam_axis([0 0 5; 0 0 -2]);
            tc.verifyEqual(a, [0 0 1], 'AbsTol', 1e-12);
        end

    end
end
