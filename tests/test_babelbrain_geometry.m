classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_babelbrain_geometry < matlab.unittest.TestCase
% TEST_BABELBRAIN_GEOMETRY  Unit tests for BABELBRAIN_CHECK_GEOMETRY: the
% BabelBrain backend takes one transducer, with its beam along the grid z axis.
%
%   Run with:   results = runtests('tests/test_babelbrain_geometry.m');

    methods (Test, TestTags = {'babelbrain'})

        function test_single_transducer_on_axis_passes(tc)
            p = params(1);
            tc.verifyWarningFree(@() babelbrain_check_geometry(p));
            tc.verifyWarningFree(@() babelbrain_check_geometry(p, [40 40 2], [40 40 80]));
            tc.verifyWarningFree(@() babelbrain_check_geometry(p, [40 40 80], [40 40 2]));   % beam -z
        end

        function test_two_transducers_error(tc)
            tc.verifyError(@() babelbrain_check_geometry(params(2)), ...
                'PRESTUS:babelbrain:multiTransducer');
        end

        function test_oblique_beam_errors(tc)
            % 20 deg off z, as each beam of an opposed pair on a mean-axis grid
            tp = [40 40 2]; fp = tp + 78 * [sind(20) 0 cosd(20)];
            tc.verifyError(@() babelbrain_check_geometry(params(1), tp, fp), ...
                'PRESTUS:babelbrain:obliqueBeam');
        end

        function test_small_tilt_passes(tc)
            % voxel rounding of a beam on z: 1 voxel lateral over 78 (0.7 deg)
            tc.verifyWarningFree(@() babelbrain_check_geometry(params(1), [40 40 2], [41 40 80]));
        end

    end
end

function p = params(n)
p = struct();
p.grid.mode = 'focal_axis';
p.transducer = repmat(struct('trans_pos', [1 1 1], 'focus_pos', [1 1 2]), 1, n);
end
