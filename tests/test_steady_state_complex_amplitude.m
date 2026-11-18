classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_steady_state_complex_amplitude < matlab.unittest.TestCase
% TEST_STEADY_STATE_COMPLEX_AMPLITUDE  Complex amplitude extraction from a recorded window
%
%   Run with:   results = runtests('tests/test_steady_state_complex_amplitude.m');

    methods (Test)

        function test_exact_for_steady_tone(tc)
            % Whole and fractional points per period, 1-3 recorded periods
            dt = 40.1606e-9; n_sim = 6268;
            rng(1); n = 500;
            a = 1e5 * (0.5 + rand(n, 1)); phi = 2 * pi * rand(n, 1) - pi;
            truth = a .* exp(1i * phi);
            for f = [250e3, 300e3, 500e3, 650e3]
                for periods = 1:3
                    n_rec = round(periods / (f * dt));
                    rsi = n_sim - n_rec + 1;
                    t = ((rsi - 1) + (0:n_rec-1)) * dt;
                    p = single(a .* cos(2 * pi * f * t + phi));
                    pc = steady_state_complex_amplitude(p, dt, f, rsi);
                    tc.verifySize(pc, [n 1]);
                    tc.verifyLessThan(max(abs(pc - truth) ./ abs(truth)), 1e-5, ...
                        sprintf('%g kHz, %d periods', f / 1e3, periods));
                end
            end
        end

        function test_phase_referred_to_simulation_start(tc)
            % Same field recorded from two start indices gives the same p_complex
            dt = 50e-9; f = 400e3; a = 2; phi = 0.7;
            n_rec = 120;
            p1 = a * cos(2 * pi * f * ((100 - 1) + (0:n_rec-1)) * dt + phi);
            p2 = a * cos(2 * pi * f * ((737 - 1) + (0:n_rec-1)) * dt + phi);
            pc1 = steady_state_complex_amplitude(p1, dt, f, 100);
            pc2 = steady_state_complex_amplitude(p2, dt, f, 737);
            tc.verifyEqual(pc1, a * exp(1i * phi), 'AbsTol', 1e-9);
            tc.verifyEqual(pc2, pc1, 'AbsTol', 1e-9);
        end

        function test_quarter_period_window(tc)
            % The window acoustic_simulation records for one frequency: exact for a steady tone
            dt = 23.9017e-9; f = 300e3; n_sim = 10365;
            ppp = 1 / (f * dt);
            n_rec = max(4, ceil(ppp / 4) + 1);
            rsi = n_sim - n_rec + 1;
            rng(3); n = 200;
            a = 1e5 * (0.5 + rand(n, 1)); phi = 2 * pi * rand(n, 1) - pi;
            t = ((rsi - 1) + (0:n_rec-1)) * dt;
            p = single(a .* cos(2 * pi * f * t + phi));
            pc = steady_state_complex_amplitude(p, dt, f, rsi);
            tc.verifyLessThan(max(abs(pc - a .* exp(1i * phi)) ./ a), 1e-4);
        end

        function test_two_tones_jointly_resolved(tc)
            % One fit per frequency: a second tone leaks into the first unless the
            % window resolves them; with ~4 periods of the difference frequency it stays small
            dt = 40e-9; f = [300e3, 450e3]; n_t = round(4 / ((f(2) - f(1)) * dt));
            t = (0:n_t-1) * dt;
            p = 1.0 * cos(2 * pi * f(1) * t + 0.3) + 0.5 * cos(2 * pi * f(2) * t - 1.1);
            pc = steady_state_complex_amplitude(p, dt, f, 1);
            tc.verifyEqual(pc, [1.0 * exp(1i * 0.3), 0.5 * exp(-1i * 1.1)], 'AbsTol', 0.05);
        end

        function test_blocks_match_single_pass(tc)
            % More points than one ~200 MB block holds: block edges change nothing
            dt = 40e-9; f = 300e3; n_t = 30;
            n = floor(2e8 / (8 * n_t)) + 17;
            p = single(repmat(cos(2 * pi * f * (0:n_t-1) * dt), n, 1));
            pc = steady_state_complex_amplitude(p, dt, f, 1);
            tc.verifyLessThan(max(abs(pc - 1)), 1e-5);
        end

    end
end
