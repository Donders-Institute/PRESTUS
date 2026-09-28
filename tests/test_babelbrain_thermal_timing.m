classdef test_babelbrain_thermal_timing < matlab.unittest.TestCase
% TEST_BABELBRAIN_THERMAL_TIMING  Unit tests for the PRESTUS -> BabelBrain thermal timing map.
%
% BabelBrain runs Repetitions x (DurationUS on, DurationOff off) and then, with
% a single sonication group, one PauseBetweenGroupedSonications off period.
% The mapping must reproduce the native PRESTUS protocol, including the
% post_ptri_dur cooldown after ptrd.
%
%   Run with:   results = runtests('tests/test_babelbrain_thermal_timing.m');

    methods (Test, TestTags = {'unit', 'quick', 'nogpu'})

        function test_two_train_protocol_with_cooldown(tc)
            t = tc.base_timing();
            bt = babelbrain_thermal_timing(t);
            tc.verifyEqual(bt.duty_cycle,     0.3,  'AbsTol', 1e-12);
            tc.verifyEqual(bt.prf_hz,         5,    'AbsTol', 1e-12);
            tc.verifyEqual(bt.duration_s,     10);
            tc.verifyEqual(bt.duration_off_s, 20);
            tc.verifyEqual(bt.repetitions,    2);
            tc.verifyEqual(bt.pause_after_s,  60);
            % BabelBrain total simulated time equals the native protocol span.
            total = bt.repetitions * (bt.duration_s + bt.duration_off_s) + bt.pause_after_s;
            tc.verifyEqual(total, t.ptrd + t.post_ptri_dur, 'AbsTol', 1e-12);
        end

        function test_no_cooldown_when_unset_or_nonpositive(tc)
            t = rmfield(tc.base_timing(), 'post_ptri_dur');
            tc.verifyEqual(babelbrain_thermal_timing(t).pause_after_s, 0);
            for v = {0, NaN, [], -5}
                t.post_ptri_dur = v{1};
                tc.verifyEqual(babelbrain_thermal_timing(t).pause_after_s, 0);
            end
        end

        function test_single_continuous_train(tc)
            t = tc.base_timing();
            t.ptri = t.ptd;  t.ptrd = t.ptd;  t.n_ptri_reps = 1;
            bt = babelbrain_thermal_timing(t);
            tc.verifyEqual(bt.duration_off_s, 0);
            tc.verifyEqual(bt.repetitions,    1);
            tc.verifyEqual(bt.pause_after_s,  60);
        end

        function test_invalid_timing_errors(tc)
            t = rmfield(tc.base_timing(), 'ptd');
            tc.verifyError(@() babelbrain_thermal_timing(t), 'babelbrain_thermal_timing:missingPtd');
            t = tc.base_timing();  t.pd = 2 * t.pri;
            tc.verifyError(@() babelbrain_thermal_timing(t), 'babelbrain_thermal_timing:invalidDutyCycle');
        end

    end

    methods (Static)
        function t = base_timing()
            t = struct('pd', 0.06, 'pri', 0.2, 'ptd', 10, 'ptri', 30, 'ptrd', 60, ...
                       'n_ptri_reps', 2, 'post_ptri_dur', 60);
        end
    end
end
