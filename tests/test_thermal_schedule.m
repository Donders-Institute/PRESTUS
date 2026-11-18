classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_thermal_schedule < matlab.unittest.TestCase
% TEST_THERMAL_SCHEDULE  Unit tests for thermal_schedule(): the async 'custom'
% schedule as a sequence of stances (per-transducer ptrd and start delay).
%
%   Run with:   results = runtests('tests/test_thermal_schedule.m');

    properties
        base  % parameters with two transducers and a 10 s protocol
    end

    methods (TestMethodSetup)
        function setup(tc)
            p = struct();
            p.timing.pd               = 0.02;   % 20 ms pulse
            p.timing.pri              = 0.1;    % 100 ms PRI
            p.timing.ptd              = 1.0;    % 1 s pulse train
            p.timing.ptri             = 2.0;    % 2 s PTRI
            p.timing.ptrd             = 10.0;   % 10 s, 5 trains
            p.timing.pt_timestep      = 0.02;
            p.timing.post_pt_timestep = 1.0;
            p.timing.post_ptri_dur    = 0.0;
            p.timing.equal_step_duration = 0;
            p.timing.async_schedule   = 'custom';
            p.modules.run_heating_sims = 1;
            p.thermal = struct();
            tx = struct('freq_hz', {5e5, 5e5}, 'timing_offset_s', {0, 15}, ...
                        'timing', {struct(), struct()});
            p.transducer = tx;
            tc.base = p;
        end
    end

    methods (Test)

        function test_single_source_ends_at_ptrd(tc)
            p = tc.base;
            p.transducer = p.transducer(1);
            s = thermal_schedule(p);
            tc.verifyEqual(s.t_sonication, 10, 'AbsTol', 1e-9);
        end

        function test_ptrd_override_sets_stance_duration(tc)
            p = tc.base;
            p.transducer(2).timing = struct('ptrd', 4);
            s = thermal_schedule(p, 2);
            tc.verifyEqual(s.source(1).n_ptri_reps, 5);
            tc.verifyEqual(s.source(2).n_ptri_reps, 2);
            tc.verifyEqual(s.source(2).on_time_s, 0.02 * 10 * 2, 'AbsTol', 1e-12);
            tc.verifyEqual(s.t_sonication, 15 + 4, 'AbsTol', 1e-9);
        end

        function test_custom_ends_at_offset_plus_ptrd(tc)
            % The last train of transducer 2 ends at 15 + 4*2 + 1 = 24 s; its
            % protocol spans 15 + 10 = 25 s, as a single run spans its PTRD.
            s = thermal_schedule(tc.base, 2);
            tc.verifyEqual(s.t_sonication, 25, 'AbsTol', 1e-9);
            tc.verifyEqual(s.t_end, 25, 'AbsTol', 1e-9);
        end

        function test_pause_between_stances_is_coarse(tc)
            s = thermal_schedule(tc.base, 2);
            in_pause = s.t0 >= 10 - 1e-9 & s.t0 < 15 - 1e-9;
            tc.verifyTrue(any(in_pause));
            tc.verifyFalse(any(s.fine(in_pause)));
            tc.verifyFalse(any(any(s.active(in_pause, :))));
            tc.verifyEqual(s.dur(in_pause), ones(nnz(in_pause), 1), 'AbsTol', 1e-9);
        end

        function test_ptrd_override_needs_custom(tc)
            p = tc.base;
            p.timing.async_schedule = 'trains';
            p.transducer(2).timing_offset_s = 0;
            p.transducer(2).timing = struct('ptrd', 4);
            tc.verifyError(@() thermal_schedule(p, 2), 'thermal_schedule:overrideNeedsCustom');
        end

        function test_ptrd_override_must_hold_whole_trains(tc)
            p = tc.base;
            p.transducer(2).timing = struct('ptrd', 3);    % 3 s / 2 s PTRI
            tc.verifyError(@() thermal_schedule(p, 2), ?MException);
        end

        function test_overlapping_stances_same_frequency_error(tc)
            p = tc.base;
            p.transducer(2).timing_offset_s = 0;
            tc.verifyError(@() thermal_schedule(p, 2), 'thermal_schedule:coherentOverlap');
        end

        function test_group_sources_share_windows(tc)
            % Two fields of one firing group (two frequencies) follow transducer 1's
            % windows together: no overlap error or warning; transducer 2 keeps its own.
            p = tc.base;
            p.transducer(1).freq_hz = 3e5;
            s = tc.verifyWarningFree(@() thermal_schedule(p, 3, SourceTx=[1 1 2], ...
                                                          SourceFreq=[3e5 5e5 5e5]));
            tc.verifyEqual([s.source.offset_s], [0 0 15]);
            on = s.active(:, 1) | s.active(:, 2);
            tc.verifyEqual(s.active(on, 1), s.active(on, 2));
        end

        function test_group_follows_its_transducer_timing(tc)
            % A single source of a group still reads its transducer's offset and override
            p = tc.base;
            p.transducer(2).timing = struct('ptrd', 4);
            s = thermal_schedule(p, 1, SourceTx=2, SourceFreq=5e5);
            tc.verifyEqual(s.mode, 'custom');
            tc.verifyEqual(s.source.offset_s, 15);
            tc.verifyEqual(s.source.n_ptri_reps, 2);
        end

        function test_groups_overlapping_at_same_frequency_error(tc)
            % Negative control: sources of different groups at one frequency still error
            p = tc.base;
            p.transducer(2).timing_offset_s = 0;
            tc.verifyError(@() thermal_schedule(p, 2, SourceTx=[1 2], SourceFreq=[5e5 5e5]), ...
                'thermal_schedule:coherentOverlap');
        end

        function test_bad_source_map_errors(tc)
            tc.verifyError(@() thermal_schedule(tc.base, 2, SourceTx=1), 'thermal_schedule:badSourceTx');
            tc.verifyError(@() thermal_schedule(tc.base, 2, SourceTx=[1 2], SourceFreq=5e5), ...
                'thermal_schedule:badSourceFreq');
        end

    end

end
