classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_isppa_policy < matlab.unittest.TestCase
% TEST_ISPPA_POLICY  When the free-water baseline runs, and when a target ISPPA can be met.
%
% Covers NEEDS_WATER_BASELINE (modules.run_water_baseline = 1 | 'auto' | 0),
% ISPPA_BASELINE (baseline lookup incl. legacy caches), ISPPA_SCALE_FACTOR and
% APPLY_ISPPA_SCALING. A target ISPPA on an uncalibrated transducer is never
% reported at the raw drive level.
% Pure logic: no k-Wave, no data.
%
%   results = runtests('tests/test_isppa_policy.m');


    methods (Static, Access = private)
        function p = params(mode, target)
            p = struct();
            p.modules.run_water_baseline = mode;
            p.transducer = struct('target_isppa_wcm2', target);
        end
        function prov = cal(isppa)
            prov = struct('tx_calibration', struct('isppa_wcm2', isppa, 'isppa_source', 'library'));
        end
    end

    %% needs_water_baseline
    methods (Test, TestTags = {'isppa_policy'})

        function test_requested_always_runs(tc)
            [run, why] = needs_water_baseline(tc.params(1, NaN), struct());
            tc.verifyTrue(run);
            tc.verifyEqual(why, 'requested');
        end

        function test_auto_no_target_does_not_run(tc)
            tc.verifyFalse(needs_water_baseline(tc.params('auto', NaN), struct()));
        end

        function test_auto_target_without_calibration_runs(tc)
            [run, why] = needs_water_baseline(tc.params('auto', 30), struct());
            tc.verifyTrue(run);
            tc.verifyEqual(why, 'no_calibration');
        end

        function test_auto_string_type_accepted(tc)
            % yaml.loadFile may return a string object
            tc.verifyTrue(needs_water_baseline(tc.params("auto", 30), struct()));
        end

        function test_auto_target_matching_calibration_skips(tc)
            tc.verifyFalse(needs_water_baseline(tc.params('auto', 30), tc.cal(30.1)));
        end

        function test_target_differing_from_calibration_runs_even_when_off(tc)
            [run, why] = needs_water_baseline(tc.params(0, 20), tc.cal(30));
            tc.verifyTrue(run);
            tc.verifyEqual(why, 'target_differs');
        end

        function test_off_target_without_calibration_does_not_run(tc)
            tc.verifyFalse(needs_water_baseline(tc.params(0, 30), struct()));
        end

        function test_missing_flag_defaults_to_auto(tc)
            p = struct('transducer', struct('target_isppa_wcm2', 30));
            tc.verifyTrue(needs_water_baseline(p, struct()));
        end

        function test_invalid_mode_errors(tc)
            tc.verifyError(@() needs_water_baseline(tc.params('sometimes', 30), struct()), ...
                'prestus:invalidRunWaterBaseline');
        end

    end

    %% isppa_baseline / isppa_scale_factor
    methods (Test, TestTags = {'isppa_policy'})

        function test_baseline_prefers_freefield(tc)
            prov = tc.cal(30);
            prov.freefield_isppa_wcm2 = 25;
            prov.freefield_source = 'kwave_simulation';
            [b, src] = isppa_baseline(prov);
            tc.verifyEqual(b, 25);
            tc.verifyEqual(src, 'kwave_simulation');
        end

        function test_legacy_freefield_without_calibration_is_simulated(tc)
            % Caches from before freefield_source: only a simulation wrote it
            [b, src] = isppa_baseline(struct('freefield_isppa_wcm2', 12));
            tc.verifyEqual(b, 12);
            tc.verifyEqual(src, 'kwave_simulation');
        end

        function test_calibration_only_baseline(tc)
            [b, src] = isppa_baseline(tc.cal(30));
            tc.verifyEqual(b, 30);
            tc.verifyEqual(src, 'library');
        end

        function test_no_baseline_is_nan(tc)
            [b, src] = isppa_baseline(struct());
            tc.verifyTrue(isnan(b));
            tc.verifyEqual(src, '');
        end

        function test_scale_factor_from_simulated_baseline(tc)
            prov = struct('freefield_isppa_wcm2', 10, 'freefield_source', 'jwave_water_sim');
            tc.verifyEqual(isppa_scale_factor(40, prov), 2, 'AbsTol', 1e-12);
        end

        function test_no_target_is_unity(tc)
            tc.verifyEqual(isppa_scale_factor(NaN, struct()), 1);
            tc.verifyEqual(isppa_scale_factor([], struct()), 1);
        end

        function test_target_without_baseline_errors(tc)
            % A target that needs scaling without a baseline errors (no unscaled field)
            tc.verifyError(@() isppa_scale_factor(30, struct()), 'prestus:isppa:noBaseline');
        end

        function test_target_matching_calibration_needs_no_rescale(tc)
            tc.verifyEqual(isppa_scale_factor(30, tc.cal(30)), 1);
        end

        function test_target_differing_from_unsimulated_baseline_errors(tc)
            tc.verifyError(@() isppa_scale_factor(20, tc.cal(30)), ...
                'prestus:isppa:unsimulatedBaseline');
        end

        function test_calibration_baseline_allowed_when_requested(tc)
            s = isppa_scale_factor(20, tc.cal(80), AllowCalibrationBaseline=true);
            tc.verifyEqual(s, 0.5, 'AbsTol', 1e-12);
        end

        function test_large_scale_warns(tc)
            prov = struct('freefield_isppa_wcm2', 1, 'freefield_source', 'kwave_simulation');
            tc.verifyWarning(@() isppa_scale_factor(100, prov), 'prestus:isppa:largeScale');
        end

    end

    %% apply_isppa_scaling
    methods (Test, TestTags = {'isppa_policy'})

        function test_apply_scales_all_pressure_fields(tc)
            sd = struct('p_max_all', ones(3), 'p_min_all', -ones(3), 'p_final', ones(3));
            prov = struct('freefield_isppa_wcm2', 10, 'freefield_source', 'kwave_simulation');
            out = apply_isppa_scaling(sd, prov, tc.params('auto', 40));
            tc.verifyEqual(out.p_max_all, 2 * ones(3), 'AbsTol', 1e-12);
            tc.verifyEqual(out.p_min_all, -2 * ones(3), 'AbsTol', 1e-12);
            tc.verifyEqual(out.p_final, 2 * ones(3), 'AbsTol', 1e-12);
            tc.verifyEqual(out.isppa_scale_p, 2, 'AbsTol', 1e-12);
        end

        function test_apply_errors_for_unmet_target(tc)
            sd = struct('p_max_all', ones(3));
            tc.verifyError(@() apply_isppa_scaling(sd, struct(), tc.params(0, 30)), ...
                'prestus:isppa:noBaseline');
        end

        function test_apply_skips_combined_caches(tc)
            sd = struct('p_max_all', ones(3));
            prov = struct('combination_mode', 'async_dual_field');
            out = apply_isppa_scaling(sd, prov, tc.params(0, 30));
            tc.verifyEqual(out.p_max_all, ones(3));
        end

    end

end
