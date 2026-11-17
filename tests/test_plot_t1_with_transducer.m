classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_plot_t1_with_transducer < matlab.unittest.TestCase
% TEST_PLOT_T1_WITH_TRANSDUCER  Focus box drawing on the T1 slice, inside and at the edge.
%
% A focus whose 5x5 box leaves the slice raises the 'Boxes may be out of
% bounds' warning instead of being dropped silently; a central focus draws
% the box without that warning. Needs k-Wave's makeBowl on the path.
%
%   results = runtests('tests/test_plot_t1_with_transducer.m');

    methods (Test, TestTags = {'plot'})

        function test_edge_focus_warns(tc)
            % `catch warn(...)` parsed warn as the exception variable: no warning
            [t1, p] = tc.inputs();
            lastwarn('');
            res = plot_t1_with_transducer(t1, 1, [20 20 30], [1 20 5], p);
            msg = lastwarn();
            tc.verifySubstring(msg, 'Boxes may be out of bounds');
            tc.verifySize(res, [40 40 3]);
        end

        function test_central_focus_draws_box(tc)
            % negative control: box inside the slice, no warning
            [t1, p] = tc.inputs();
            lastwarn('');
            res = plot_t1_with_transducer(t1, 1, [20 20 35], [20 20 15], p);
            tc.verifyEmpty(strfind(lastwarn(), 'Boxes may be out of bounds'));
            tc.verifyEqual(res(18:22, 13:17, 2), ones(5), 'Focus box not drawn.');
        end
    end

    methods (Access = private)
        function [t1, p] = inputs(tc)
            tc.assumeTrue(exist('makeBowl', 'file') > 0, 'k-Wave (makeBowl) not on the path.');
            t1 = rand(40, 40, 40);
            p = make_minimal_parameters(tempname());
            tc.addTeardown(@() rmdir(p.path.sim, 's'));
        end
    end
end
