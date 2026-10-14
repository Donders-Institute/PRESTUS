classdef test_coherent_drive_calibration < matlab.unittest.TestCase
% TEST_COHERENT_DRIVE_CALIBRATION  Unit tests for COHERENT_DRIVE_CALIBRATION:
% each transducer's drive is set to its own ISPPA target before a coherent run.
% The free-water simulation is replaced by a model: ISPPA ∝ amplitude² of the
% one transducer left firing.
%
%   Run with:   results = runtests('tests/test_coherent_drive_calibration.m');

    methods (Test, TestTags = {'isppa'})

        function test_simulated_baselines_reach_targets(tc)
            p = params([20 40], [NaN NaN]);
            [p2, cal] = coherent_drive_calibration(p, @model_baseline);
            for k = 1:2
                tc.verifyEqual(model_isppa(p2, k), cal(k).target_isppa_wcm2, 'RelTol', 1e-12);
                tc.verifyEqual(cal(k).baseline_source, 'simulation');
            end
        end

        function test_baseline_is_of_one_transducer(tc)
            % the simulated source has every other transducer at zero amplitude
            p = params([20 40], [NaN NaN]);
            seen = {};
            function v = spy(q)
                seen{end+1} = arrayfun(@(t) max(t.annular.elem_p_amp), q.transducer);
                v = model_baseline(q);
            end
            coherent_drive_calibration(p, @spy);
            tc.verifyEqual(seen{1}(2), 0); tc.verifyGreaterThan(seen{1}(1), 0);
            tc.verifyEqual(seen{2}(1), 0); tc.verifyGreaterThan(seen{2}(2), 0);
        end

        function test_matching_calibration_needs_no_simulation(tc)
            % T1 calibrated at 30 W/cm² with target 30: no simulation, no rescale
            p = params([30 40], [30 NaN]);
            n = 0;
            function v = count(q), n = n + 1; v = model_baseline(q); end
            [p2, cal] = coherent_drive_calibration(p, @count);
            tc.verifyEqual(n, 1);   % only T2
            tc.verifyEqual(cal(1).baseline_source, 'calibration');
            tc.verifyEqual(cal(1).amplitude_scale, 1);
            tc.verifyEqual(p2.transducer(1).annular.elem_p_amp, p.transducer(1).annular.elem_p_amp);
        end

        function test_no_target_keeps_drive(tc)
            p = params([NaN 40], [NaN NaN]);
            [p2, cal] = coherent_drive_calibration(p, @model_baseline);
            tc.verifyEqual(p2.transducer(1).annular.elem_p_amp, p.transducer(1).annular.elem_p_amp);
            tc.verifyEqual(cal(1).amplitude_scale, 1);
        end

        function test_vector_target_errors(tc)
            p = params([20 40], [NaN NaN]);
            p.transducer(2).target_isppa_wcm2 = [10 20];
            tc.verifyError(@() coherent_drive_calibration(p, @model_baseline), 'prestus:isppa:coherentSweep');
        end

    end
end

function p = params(targets, cal)
% Two annular transducers, 3 rings; model baseline 10 W/cm² per (100 kPa)²
p = struct();
p.modules.run_water_baseline = 'auto';
for k = 1:2
    t = struct('type', 'annular', 'target_isppa_wcm2', targets(k));
    t.annular = struct('elem_p_amp', 1e5 * [1 0.9 0.8] * k);
    t.tx_calibration = struct('isppa_wcm2', NaN);
    if isfinite(cal(k))
        t.tx_calibration.isppa_wcm2 = cal(k);
        t.annular.elem_p_amp = 1e5 * [1 0.9 0.8] * sqrt(cal(k) / 10);   % drive that gives cal(k)
    end
    p.transducer(k) = t;
end
end

function v = model_isppa(p, k)
v = 10 * (max(p.transducer(k).annular.elem_p_amp) / 1e5)^2;
end

function v = model_baseline(p)
% ISPPA of the only transducer firing
on = find(arrayfun(@(t) max(t.annular.elem_p_amp) > 0, p.transducer));
assert(isscalar(on), 'exactly one transducer should fire');
v = model_isppa(p, on);
end
