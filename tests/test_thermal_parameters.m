classdef test_thermal_parameters < matlab.unittest.TestCase
% TEST_THERMAL_PARAMETERS  Unit tests for thermal_parameters()
%
%   Run with:   results = runtests('tests/test_thermal_parameters.m');

    methods (TestMethodSetup)
        function setup(tc)
            % Minimal valid timing struct used as base for all tests
            tc.base = struct();
            tc.base.timing.pd              = 0.02;    % 20 ms pulse
            tc.base.timing.pri             = 0.1;     % 100 ms PRI  → DC = 0.2
            tc.base.timing.ptd             = 1.0;     % 1 s pulse train
            tc.base.timing.pt_timestep     = 0.02;
            tc.base.timing.ptri            = 2.0;     % 2 s PTRI
            tc.base.timing.ptrd            = 10.0;    % 10 s total
            tc.base.timing.post_pt_timestep = 1.0;
            tc.base.timing.post_ptri_dur   = 0.0;
            tc.base.timing.equal_step_duration = 0;
            tc.base.modules.run_heating_sims = 1;
            tc.base.thermal = struct();  % thermal_parameters also reads parameters.thermal
        end
    end

    properties
        base  % base parameter struct
    end

    methods (Test, TestTags = {'duty_cycle'})

        function test_duty_cycle(tc)
            p = thermal_parameters(tc.base, true);
            tc.verifyEqual(p.dc, tc.base.timing.pd / tc.base.timing.pri, 'AbsTol', 1e-10);
        end

        function test_prf(tc)
            p = thermal_parameters(tc.base, true);
            tc.verifyEqual(p.prf, 1 / tc.base.timing.pri, 'AbsTol', 1e-10);
        end

    end

    methods (Test, TestTags = {'pulse_counts'})

        function test_pulses_per_train(tc)
            p = thermal_parameters(tc.base, true);
            expected = round(tc.base.timing.ptd / tc.base.timing.pri);
            tc.verifyEqual(p.n_pulses_per_pt, expected);
        end

        function test_ptri_repetitions(tc)
            p = thermal_parameters(tc.base, true);
            expected_reps = round(tc.base.timing.ptrd / tc.base.timing.ptri);
            tc.verifyEqual(p.n_ptri_reps, expected_reps);
        end

    end

    methods (Test, TestTags = {'step_counts'})

        function test_on_off_steps_sum_to_pri(tc)
            % ON + OFF step durations should sum to PRI within one pulse
            p = thermal_parameters(tc.base, true);
            on_dur  = p.pt_on_steps_n  * p.pt_on_steps_dur;
            off_dur = p.pt_off_steps_n * p.pt_off_steps_dur;
            tc.verifyEqual(on_dur + off_dur, tc.base.timing.pri, 'AbsTol', 1e-9);
        end

    end

    methods (Test, TestTags = {'validation'})

        function test_non_integer_pulses_per_train_errors(tc)
            p = tc.base;
            p.timing.ptd = 0.15;  % 0.15 / 0.1 = 1.5 → not integer
            tc.verifyError(@() thermal_parameters(p, true), '');
        end

    end

end
