classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_thermal_overlap_terms < matlab.unittest.TestCase
% TEST_THERMAL_OVERLAP_TERMS  Interference of partly overlapping pulses
% (thermal_overlap_terms): which segments carry the interference term, and its
% phase for a shared clock and for per-pulse bursts. Pure logic.
%
%   results = runtests('tests/test_thermal_overlap_terms.m');

    properties
        p   % two transducers, 20 ms pulses every 200 ms for 0.4 s; T2 later by d
    end

    methods (TestMethodSetup)
        function setup(tc)
            p = struct();
            p.timing = struct('pd', 0.02, 'pri', 0.2, 'ptd', 0.4, 'ptri', 0.4, 'ptrd', 0.4, ...
                'pt_timestep', 0.01, 'post_pt_timestep', 0.1, 'post_ptri_dur', 0, ...
                'equal_step_duration', 0, 'async_schedule', 'custom');
            p.thermal = struct();
            p.modules.run_heating_sims = 1;
            d = 0.01 + 1 / (4 * 5e5);          % 10 ms plus a quarter period at 500 kHz
            p.transducer = struct('freq_hz', {5e5, 5e5}, 'timing_offset_s', {0, d}, ...
                                  'timing', {struct(), struct()});
            tc.p = p;
        end
    end

    methods (Test)

        function test_overlap_segments_and_clock_phase(tc)
            s = schedule(tc.p);
            ov = thermal_overlap_terms(s, [5e5 5e5], [1 2], 'clock');
            tc.verifyEqual(ov.pairs, [1 2]);
            % both on for 10 ms minus a quarter period in each of the two PRIs
            tc.verifyEqual(sum(s.dur(ov.active)), 2 * (0.02 - 0.01 - 0.5e-6), 'AbsTol', 1e-9);
            tc.verifyEqual(ov.cos_w(ov.active), ones(nnz(ov.active), 1));
            tc.verifyEqual(ov.thetas, {0});
            tc.verifyEqual(ov.cos_w(~ov.active), zeros(nnz(~ov.active), 1));
        end

        function test_burst_phase_from_onsets(tc)
            % theta = -2*pi*f*d: 10 ms is a whole number of periods, the quarter period is -pi/2
            s = schedule(tc.p);
            ov = thermal_overlap_terms(s, [5e5 5e5], [1 2], 'burst');
            tc.verifyEqual(ov.thetas{1}, -pi / 2, 'AbsTol', 1e-6);
            tc.verifyEqual(ov.sin_w(ov.active), -ones(nnz(ov.active), 1), 'AbsTol', 1e-6);
            tc.verifyEqual(ov.cos_w(ov.active), zeros(nnz(ov.active), 1), 'AbsTol', 1e-6);
        end

        function test_no_pair_across_frequencies_or_within_a_group(tc)
            % negative controls: different frequencies average out; one group is summed when combined
            s = schedule(tc.p);
            tc.verifyEmpty(thermal_overlap_terms(s, [5e5 3e5], [1 2], 'clock').pairs);
            tc.verifyEmpty(thermal_overlap_terms(s, [5e5 5e5], [1 1], 'clock').pairs);
        end

        function test_no_pair_without_overlap(tc)
            p = tc.p;
            p.transducer(2).timing_offset_s = 0.05;   % after T1's pulse ends
            s = schedule(p);
            tc.verifyEmpty(thermal_overlap_terms(s, [5e5 5e5], [1 2], 'clock').pairs);
        end

        function test_train_averaged_is_refused(tc)
            s = schedule(tc.p);
            s.pulse_model = 'train_averaged';
            tc.verifyError(@() thermal_overlap_terms(s, [5e5 5e5], [1 2], 'clock'), ...
                'thermal_overlap_terms:trainAveraged');
        end

    end
end

function s = schedule(p)
% thermal_schedule refuses overlap at one frequency; build the windows with the
% frequencies told apart, the pulse windows do not depend on them
w = warning('off', 'thermal_schedule:overlap');
c = onCleanup(@() warning(w));
s = thermal_schedule(p, 2, SourceTx=[1 2], SourceFreq=[5e5 5e5 + 1]);
end
