classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_cw_ramp_periods < matlab.unittest.TestCase
% TEST_CW_RAMP_PERIODS  Start ramp of CW sources (CW_RAMP_PERIODS): k-Wave's 4
% periods for runs of 8 periods or more, half the run below that, so that
% CREATECWSIGNALS works for short runs (codetypes set C: 464 steps,
% 2026-11-17).

    methods (Test)

        function test_long_run_keeps_default(tc)
            t = (0:9999) * 4e-8;                                 % 400 us at 300 kHz: 120 periods
            tc.verifyEqual(cw_ramp_periods(t, 300e3), 4);
            t = (0:ceil(8 / 300e3 / 4e-8)) * 4e-8;               % just over 8 periods
            tc.verifyEqual(cw_ramp_periods(t, 300e3), 4);
        end

        function test_short_run_halves(tc)
            dt = 2.67e-8;
            t = (0:463) * dt;                                    % 464 steps, as set C
            n = cw_ramp_periods(t, 300e3);
            tc.verifyEqual(n, numel(t) * dt * 300e3 / 2, 'RelTol', 1e-12);
            tc.verifyLessThan(n, 4);
            s = createCWSignals(t, 300e3, [1; 2], [0; 0.5], n);  % no index error
            tc.verifySize(s, [2 numel(t)]);
        end

        function test_default_ramp_fails_on_short_run(tc)
            % negative control: the case the helper exists for
            t = (0:463) * 2.67e-8;
            tc.verifyError(@() createCWSignals(t, 300e3, 1, 0), ?MException);
        end

        function test_too_short_for_a_ramp(tc)
            t = (0:1) * 1e-8;                                    % ramp would span < 2 steps: none
            tc.verifyEqual(cw_ramp_periods(t, 300e3), 0);
            tc.verifyEqual(cw_ramp_periods(0, 300e3), 0);
        end
    end
end
