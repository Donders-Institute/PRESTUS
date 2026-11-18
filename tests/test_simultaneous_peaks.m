classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_simultaneous_peaks < matlab.unittest.TestCase
% TEST_SIMULTANEOUS_PEAKS  Peak pressures of tones firing together (SIMULTANEOUS_PEAKS)
% against closed forms: one tone, two tones in phase, two tones in antiphase
% at the peak of one, a pulse shorter than the common period, and the upper
% bound when the common period needs too many samples.

    methods (Test)

        function test_single_tone_is_its_amplitude(tc)
            P = {[2; 5 * exp(1i * 0.7)]};
            [pp, pn, info] = simultaneous_peaks(P, 300e3, NaN);
            tc.verifyEqual(pp, [2; 5], 'RelTol', 2e-3);
            tc.verifyEqual(pn, -[2; 5], 'RelTol', 2e-3);
            tc.verifyEqual(info.method, 'phase_sampled');
            tc.verifyEqual(info.common_period_s, 1 / 300e3, 'RelTol', 1e-12);
        end

        function test_two_tones_in_phase_add(tc)
            % cos(wt) + cos(2wt) peaks at t = 0 with 2; its minimum is -9/8
            [pp, pn, info] = simultaneous_peaks({1, 1}, [300e3 600e3], NaN);
            tc.verifyEqual(pp, 2, 'RelTol', 2e-3);
            tc.verifyEqual(pn, -9/8, 'RelTol', 2e-3);
            tc.verifyEqual(info.common_period_s, 1 / 300e3, 'RelTol', 1e-12);
        end

        function test_sum_of_amplitudes_is_not_reached_out_of_phase(tc)
            % negative control: cos(wt) - cos(3wt) = 4c - 4c^3 never reaches 2 (max 8/(3*sqrt(3)) ~ 1.54)
            [pp, ~] = simultaneous_peaks({1, -1}, [300e3 900e3], NaN);
            tc.verifyLessThan(pp, 1.6);
            tc.verifyGreaterThan(pp, 1.5);
        end

        function test_short_pulse_limits_the_window(tc)
            % 300 + 450 kHz: common period 6.67 us; a 1 us pulse is evaluated over 1 us only
            [~, ~, info] = simultaneous_peaks({1, 1}, [300e3 450e3], 1e-6);
            tc.verifyEqual(info.evaluated_s, 1e-6, 'RelTol', 1e-12);
            [~, ~, info] = simultaneous_peaks({1, 1}, [300e3 450e3], 1);
            tc.verifyEqual(info.evaluated_s, 1 / 150e3, 'RelTol', 1e-9);
        end

        function test_upper_bound_for_long_common_period(tc)
            % 300000 Hz + 300001 Hz: gcd 1 Hz, common period 1 s -> too many samples
            [pp, pn, info] = simultaneous_peaks({2, 3}, [300e3 300001], NaN);
            tc.verifyEqual(info.method, 'upper_bound');
            tc.verifyEqual(pp, 5);
            tc.verifyEqual(pn, -5);
        end
    end
end
