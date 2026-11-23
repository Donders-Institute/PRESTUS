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

        function test_ramp_factor_per_source_pd(tc)
            % custom schedule, transducer 2 with its own PD: each source gets the ramp
            % factor of its own PD; without a ramp all factors are 1; a ramp longer
            % than half of one source's PD stops the run
            p = tc.base;
            p.transducer(2).timing = struct('pd', 0.05);
            p.timing.pulse_ramp_dur = 0.004;
            s = thermal_schedule(p, 2);
            pt = thermal_parameters(p, true);
            eta = thermal_ramp_factors(pt, s);
            tc.verifyEqual(eta, [compute_ramp_factor(0.02, 0.004, 'linear'), ...
                                 compute_ramp_factor(0.05, 0.004, 'linear')], 'AbsTol', 1e-12);
            tc.verifyNotEqual(eta(1), eta(2));
            p0 = tc.base;
            tc.verifyEqual(thermal_ramp_factors(thermal_parameters(p0, true), thermal_schedule(p0, 2)), [1 1]);
            p.timing.pulse_ramp_dur = 0.009;              % fits the global 20 ms PD
            p.transducer(2).timing = struct('pd', 0.015); % but not 2 x 9 ms in 15 ms
            tc.verifyError(@() thermal_ramp_factors(thermal_parameters(p, true), thermal_schedule(p, 2)), ...
                'PRESTUS:thermal:ramp');
        end

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

        function test_overlapping_stances_same_frequency_interfere(tc)
            % Same frequency, both on: the pair interferes in every segment where both
            % are on, with the clock phase by default and the onset delay with 'burst'
            p = tc.base;
            p.transducer(2).timing_offset_s = 0.01 + 1 / (4 * 5e5);   % 10 ms + a quarter period
            s = thermal_schedule(p, 2);
            tc.verifyEqual(s.overlap.pairs, [1 2]);
            both = s.active(:, 1) & s.active(:, 2);
            tc.verifyEqual(s.overlap.active(:, 1), both);
            tc.verifyEqual(s.overlap.thetas, {0});
            p.timing.overlap_phase = 'burst';
            s = thermal_schedule(p, 2);
            tc.verifyEqual(s.overlap.thetas{1}, -pi / 2, 'AbsTol', 1e-6);
            p.timing.overlap_phase = 'random';
            tc.verifyError(@() thermal_schedule(p, 2), 'thermal_schedule:overlapPhase');
        end

        function test_overlap_needs_known_frequencies(tc)
            p = tc.base;
            p.transducer(2).timing_offset_s = 0;
            p.transducer(2).freq_hz = [];
            tc.verifyError(@() thermal_schedule(p, 2), 'thermal_schedule:overlapUnknownFrequency');
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

        function test_groups_overlapping_at_same_frequency_interfere(tc)
            % Sources of different groups at one frequency form an interfering pair;
            % sources of one group never do (negative control)
            p = tc.base;
            p.transducer(2).timing_offset_s = 0;
            s = thermal_schedule(p, 2, SourceTx=[1 2], SourceFreq=[5e5 5e5]);
            tc.verifyEqual(s.overlap.pairs, [1 2]);
            s = thermal_schedule(p, 2, SourceTx=[1 1], SourceFreq=[3e5 5e5]);
            tc.verifyEmpty(s.overlap.pairs);
        end

        function test_train_averaged_keeps_heat_and_drops_pulses(tc)
            % One source, 5 trains of 10 pulses: 'train_averaged' heats through each
            % train with weight PD/PRI, so heat x time equals the resolved pulses',
            % with one window per train instead of one per pulse
            p = tc.base;
            p.transducer = p.transducer(1);
            a = thermal_schedule(p);                      % default
            p.thermal.pulse_model = 'pulse_resolved';
            r = thermal_schedule(p);
            tc.verifyEqual(r.pulse_model, 'pulse_resolved');
            tc.verifyEqual(a.pulse_model, 'train_averaged');
            tc.verifyEqual(sum(a.weight .* a.dur), sum(r.weight .* r.dur), 'RelTol', 1e-12);
            tc.verifyEqual(sum(r.weight .* r.dur), 5 * 10 * 0.02, 'RelTol', 1e-12);
            tc.verifyEqual(unique(a.weight(a.weight > 0)), 0.02 / 0.1, 'AbsTol', 1e-12);
            tc.verifyEqual(sum(a.dur(a.active)), 5 * 1.0, 'AbsTol', 1e-9);      % heats through each 1 s train
            tc.verifyLessThan(numel(a.t0), numel(r.t0));
            tc.verifyEqual(a.t_end, r.t_end, 'AbsTol', 1e-9);
        end

        function test_interleaved_async_resolves_pulses(tc)
            % 'custom' with transducers in turn cannot average within a train: the
            % default falls back to resolved pulses; an unknown model still errors
            p = tc.base;
            s = thermal_schedule(p, 2);
            tc.verifyEqual(s.pulse_model, 'pulse_resolved');
            tc.verifyEqual(unique(s.weight(s.weight > 0)), 1);
            p.thermal.pulse_model = 'smooth';
            tc.verifyError(@() thermal_schedule(p, 2), 'thermal_schedule:pulseModel');
        end

        function test_bad_source_map_errors(tc)
            tc.verifyError(@() thermal_schedule(tc.base, 2, SourceTx=1), 'thermal_schedule:badSourceTx');
            tc.verifyError(@() thermal_schedule(tc.base, 2, SourceTx=[1 2], SourceFreq=5e5), ...
                'thermal_schedule:badSourceFreq');
        end

    end

end
