classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_outputs_follow_recompute < matlab.unittest.TestCase
% TEST_OUTPUTS_FOLLOW_RECOMPUTE  A recomputed stage replaces its existing outputs.
%
% Regression for io.overwrite_files 'never': a re-run without caches
% recomputed the fields and rewrote the CSV but kept the earlier NIfTIs.
% Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_outputs_follow_recompute.m');


    methods (Static, Access = private)
        function p = params(mode, recomputed)
            p = struct();
            p.io.overwrite_files = mode;
            if ~isempty(recomputed)
                p.state.recomputed = struct('acoustic', recomputed);
            end
        end
    end

    methods (Test)
        function test_recomputed_stage_overrides_never(tc)
            p = outputs_follow_recompute(tc.params('never', true), 'acoustic');
            tc.verifyEqual(p.io.overwrite_files, 'always');
        end

        function test_cached_stage_keeps_never(tc)
            p = outputs_follow_recompute(tc.params('never', false), 'acoustic');
            tc.verifyEqual(p.io.overwrite_files, 'never');
        end

        function test_unknown_stage_state_keeps_never(tc)
            p = outputs_follow_recompute(tc.params('never', []), 'acoustic');
            tc.verifyEqual(p.io.overwrite_files, 'never');
            p = outputs_follow_recompute(tc.params('never', true), 'thermal');
            tc.verifyEqual(p.io.overwrite_files, 'never');
        end

        function test_force_untouched(tc)
            p = outputs_follow_recompute(tc.params('force', true), 'acoustic');
            tc.verifyEqual(p.io.overwrite_files, 'force');
        end
    end
end
