classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_thermal_plane_frames < matlab.unittest.TestCase
% TEST_THERMAL_PLANE_FRAMES  Focal-plane frame selection under thermal.max_focal_plane_gb.
%
% Every combination of protocol length, frame size, and budget: the first and
% last timepoints are always kept, kept timepoints are evenly spaced, the
% kept frames fit the budget (up to the final frame), and budget 0 or a
% budget above the full series keeps every timepoint. A negative budget is
% rejected (negative control). No k-Wave run, no subject data.
%
%   results = runtests('tests/test_thermal_plane_frames.m');

    properties (TestParameter)
        n_time = struct('single_seg', 2, 'short', 101, 'mouse', 6521);
        n_vox  = struct('small', 1e3, 'grid_2d', 792 * 792);
        budget = struct('none', 0, 'tiny', 1e-3, 'default', 2, 'large', 1e4);
    end

    methods (Test)

        function test_selection(tc, n_time, n_vox, budget)
            [keep, stride] = thermal_plane_frames(n_time, n_vox, budget);

            tc.verifySize(keep, [1 n_time]);
            tc.verifyClass(keep, 'logical');
            tc.verifyTrue(keep(1) && keep(end), 'Initial and final state must be kept.');

            frame_bytes = 3 * 8 * n_vox;
            if budget == 0 || n_time * frame_bytes <= budget * 1e9
                tc.verifyEqual(stride, 1);
                tc.verifyTrue(all(keep), 'Within budget, every timepoint is kept.');
            else
                tc.verifyGreaterThan(stride, 1);
                gaps = diff(find(keep(1:end-1)));
                tc.verifyTrue(all(gaps == stride), 'Kept timepoints are every stride-th.');
                tc.verifyLessThanOrEqual((nnz(keep) - 1) * frame_bytes, budget * 1e9 + frame_bytes);
            end
        end

        function test_negative_budget_rejected(tc)
            tc.verifyError(@() thermal_plane_frames(10, 10, -1), 'MATLAB:validators:mustBeNonnegative');
        end

    end
end
