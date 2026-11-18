classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_firing_groups < matlab.unittest.TestCase
% TEST_FIRING_GROUPS  Unit tests for firing_groups(): which transducers fire
% together in an async run, and the heat sources (one per group and frequency).
%
%   Run with:   results = runtests('tests/test_firing_groups.m');

    properties
        base  % async 'custom' run with three transducers at 500 kHz
    end

    methods (TestMethodSetup)
        function setup(tc)
            p = struct();
            p.simulation.transducer_coupling = 'async';
            p.timing = struct('pd', 0.02, 'pri', 0.1, 'ptd', 1, 'ptri', 2, 'ptrd', 10, ...
                              'async_schedule', 'custom');
            p.transducer = struct('freq_hz', {5e5, 5e5, 5e5}, 'timing_offset_s', {0, 0, 20}, ...
                                  'timing', {struct(), struct(), struct()});
            tc.base = p;
        end
    end

    methods (Test)

        function test_same_start_and_timing_is_one_group(tc)
            fg = firing_groups(tc.base);
            tc.verifyEqual(fg.group, [1 1 2]);
            tc.verifyEqual(fg.coherent, [true false]);
            tc.verifyEqual(fg.members, {[1 2], 3});
            tc.verifyEqual(fg.source_tx, [1 3]);
            tc.verifyEqual(fg.source_freq, [5e5 5e5]);
            tc.verifyEqual(fg.source_group, [1 2]);
        end

        function test_different_timing_is_not_a_group(tc)
            % Same start, but transducer 2 sonicates longer: overlapping, not together
            p = tc.base;
            p.transducer(2).timing = struct('ptrd', 4);
            fg = firing_groups(p);
            tc.verifyEqual(fg.group, [1 2 3]);
            tc.verifyFalse(any(fg.coherent));
        end

        function test_override_equal_to_block_is_a_group(tc)
            % An override that repeats the timing block keeps the transducer in the group
            p = tc.base;
            p.transducer(2).timing = struct('ptrd', 10);
            fg = firing_groups(p);
            tc.verifyEqual(fg.group, [1 1 2]);
        end

        function test_different_start_is_not_a_group(tc)
            p = tc.base;
            p.transducer(2).timing_offset_s = 10;
            fg = firing_groups(p);
            tc.verifyEqual(fg.group, [1 2 3]);
        end

        function test_non_custom_schedule_has_no_groups(tc)
            % 'pulses' and 'trains' interleave the transducers by construction
            for mode = {'pulses', 'trains'}
                p = tc.base;
                p.timing.async_schedule = mode{1};
                fg = firing_groups(p);
                tc.verifyEqual(fg.group, [1 2 3], mode{1});
                tc.verifyEqual(fg.source_tx, [1 2 3], mode{1});
            end
        end

        function test_coherent_coupling_has_no_groups_at_one_frequency(tc)
            % One coherent k-Wave run, not the async pipeline: firing_groups is not used,
            % but must not invent groups
            p = tc.base;
            p.simulation.transducer_coupling = 'coherent';
            fg = firing_groups(p);
            tc.verifyEqual(fg.group, [1 2 3]);
        end

        function test_simultaneous_multifreq_is_one_group(tc)
            p = tc.base;
            p.simulation.transducer_coupling = 'coherent';
            p.transducer(3).freq_hz = 3e5;
            fg = firing_groups(p);
            tc.verifyEqual(fg.group, [1 1 1]);
            tc.verifyEqual(fg.source_freq, [3e5 5e5]);   % one source per frequency, sorted
            tc.verifyEqual(fg.source_tx, [1 1]);
        end

        function test_group_at_two_frequencies_gives_two_sources(tc)
            p = tc.base;
            p.transducer(2).freq_hz = 2.5e5;
            fg = firing_groups(p);
            tc.verifyEqual(fg.group, [1 1 2]);
            tc.verifyEqual(fg.source_tx, [1 1 3]);
            tc.verifyEqual(fg.source_freq, [2.5e5 5e5 5e5]);
            tc.verifyEqual(fg.source_group, [1 1 2]);
        end

        function test_single_transducer(tc)
            p = tc.base;
            p.transducer = p.transducer(1);
            fg = firing_groups(p);
            tc.verifyEqual(fg.group, 1);
            tc.verifyFalse(fg.coherent);
        end

        function test_sources_feed_thermal_schedule(tc)
            % The group's windows are used once; without the map the overlap errors
            p = tc.base;
            p.timing.pt_timestep = 0.02; p.timing.post_pt_timestep = 1;
            p.timing.post_ptri_dur = 0; p.timing.equal_step_duration = 0;
            p.thermal = struct();
            p.modules.run_heating_sims = 1;
            fg = firing_groups(p);
            s = thermal_schedule(p, numel(fg.source_tx), SourceTx=fg.source_tx, SourceFreq=fg.source_freq);
            tc.verifyEqual([s.source.offset_s], [0 20]);
            tc.verifyError(@() thermal_schedule(p, 3), 'thermal_schedule:coherentOverlap');
        end

    end

end
