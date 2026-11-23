classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_legacy_bridge < matlab.unittest.TestCase
% TEST_LEGACY_BRIDGE  Study scripts written for earlier PRESTUS versions keep running, with warnings.
%
% Covers LEGACY_BRIDGE_ENTRY (pre-v0.6 call form), LEGACY_BRIDGE_STRUCT (old
% key names set on the parameters struct after LOAD_PARAMETERS, incl.
% transducer arrays and sequential configs), LEGACY_BRIDGE_LAYER (top-level
% expected_focal_distance_mm made redundant by the transducers), the
% deprecation entries of the v0.5 development builds, and the v0.5.0 entry
% functions. Negative controls: a current config passes unchanged and without
% warnings. No simulation runs.
%
%   results = runtests('tests/test_legacy_bridge.m');

    methods (Test)

        %% Old names set on the struct

        function test_struct_old_names_are_converted(tc)
            p = struct('subject_id', 1, 'code_type', 'cpp_gpu', ...
                       'simulation', struct('code_type', 'matlab_cpu', 'medium', 'water'), ...
                       'results_filename_affix', '_x');
            tc.verifyWarning(@() legacy_bridge_struct(p), 'prestus:renamedConfigKey');
            q = tc.quiet(@() legacy_bridge_struct(p));
            tc.verifyEqual(q.simulation.code_type, 'cpp_gpu');    % the script's value wins
            tc.verifyEqual(q.simulation.medium, 'water');        % siblings kept
            tc.verifyEqual(q.io.output_affix, '_x');
            tc.verifyFalse(isfield(q, 'code_type'));
            tc.verifyFalse(isfield(q, 'results_filename_affix'));
        end

        function test_struct_dev_build_names_are_converted(tc)
            p = struct('hpc_timelimit', '05:00:00', 'hpc_memorylimit', 40, 'hpc_wait_for_job', true, ...
                       'run_nifti_creation', 0, 'use_pseudoCT', 1, 'smooth_fwhm_mm', 2);
            q = tc.quiet(@() legacy_bridge_struct(p));
            tc.verifyEqual(q.hpc.timelimit, '05:00:00');
            tc.verifyEqual(q.hpc.memorylimit, 40);
            tc.verifyTrue(logical(q.hpc.wait_for_job));
            tc.verifyEqual(q.modules.run_nifti_creation, 0);
            tc.verifyEqual(q.pct.enabled, 1);
            tc.verifyEqual(q.headmodel.smooth_fwhm_mm, 2);
        end

        function test_struct_transducer_fields_per_element(tc)
            t1 = struct('type', 'annular', 'set_intensity_w_per_cm2', 80, 'focal_distance_ep', 40);
            t2 = struct('type', 'annular', 'set_intensity_w_per_cm2', 30, 'focal_distance_ep', 60);
            p = struct('transducer', [t1, t2]);
            q = tc.quiet(@() legacy_bridge_struct(p));
            tc.verifySize(q.transducer, [1 2]);
            tc.verifyEqual(q.transducer(1).tx_calibration.isppa_wcm2, 80);
            tc.verifyEqual(q.transducer(2).tx_calibration.isppa_wcm2, 30);
            tc.verifyEqual(q.transducer(2).focal_distance_ep, 60);
            tc.verifyEqual(q.transducer(1).type, 'annular');
            tc.verifyFalse(isfield(q.transducer, 'set_intensity_w_per_cm2'));
        end

        function test_struct_transducer_v050_fields_keep_given_type(tc)
            p = struct('transducer', struct('type', 'annular', 'n_elements', 4, 'focal_distance_ep', 50));
            q = tc.quiet(@() legacy_bridge_struct(p));
            tc.verifyEqual(q.transducer.annular.elem_n, 4);
            tc.verifyEqual(q.transducer.type, 'annular');
        end

        function test_struct_unconvertible_key_stops(tc)
            p = struct('transducer', struct('expected_focal_distance_mm', 60));
            tc.verifyError(@() tc.quiet(@() legacy_bridge_struct(p)), 'prestus:legacyConfigKey');
        end

        %% Negative controls

        function test_current_struct_is_unchanged_and_silent(tc)
            p = struct('subject_id', 3, 'simulation', struct('code_type', 'cpp_gpu'), ...
                       'transducer', struct('type', 'annular', 'focal_distance_ep', 60));
            tc.verifyWarningFree(@() legacy_bridge_struct(p));
            tc.verifyEqual(legacy_bridge_struct(p), p);
        end

        function test_loaded_default_config_is_unchanged_and_silent(tc)
            % Every field LOAD_PARAMETERS returns must pass untouched: the bridge
            % must never take a current field for an old one.
            ws = warning('off', 'all');
            restore = onCleanup(@() warning(ws));
            evalc('P = load_parameters();');
            clear restore
            P.subject_id = 1;
            tc.verifyWarningFree(@() legacy_bridge_struct(P));
            tc.verifyEqual(legacy_bridge_struct(P), P);
            tc.verifyWarningFree(@() legacy_bridge_entry('prestus_pipeline_start', P, struct(), []));
        end

        %% Call forms

        function test_pre_v06_call_form(tc)
            p = struct('simulation', struct('code_type', 'cpp_gpu'));
            tc.verifyWarning(@() legacy_bridge_entry('prestus_pipeline_start', 7, p, []), ...
                'prestus:legacyCallForm');
            [q, o] = tc.quiet2(@() legacy_bridge_entry('prestus_pipeline_start', 7, p, []));
            tc.verifyEqual(q.subject_id, 7);
            tc.verifyEqual(o, struct());
            seq = struct('sequential_configs', struct('config_2', struct('code_type', 'cpp_gpu')));
            [~, o] = tc.quiet2(@() legacy_bridge_entry('prestus_pipeline_start', 7, p, seq));
            tc.verifyEqual(o.sequential_configs.config_2.simulation.code_type, 'cpp_gpu');
        end

        function test_current_call_form_passes(tc)
            p = struct('subject_id', 2);
            [q, o] = legacy_bridge_entry('prestus_pipeline_start', p, struct(), []);
            tc.verifyEqual(q, p);
            tc.verifyEqual(o, struct());
        end

        function test_bad_inputs_stop(tc)
            tc.verifyError(@() legacy_bridge_entry('f', struct(), struct(), struct('a', 1)), 'prestus:tooManyInputs');
            tc.verifyError(@() legacy_bridge_entry('f', 'abc', struct(), []), 'prestus:invalidParameters');
            tc.verifyError(@() legacy_bridge_entry('f', struct(), 5, []), 'prestus:invalidOptions');
        end

        function test_pipeline_entry_points_take_the_old_form(tc)
            % The argument blocks accept a subject_id first (validation is in the bridge)
            for f = {'prestus_pipeline_start', 'prestus_pipeline'}
                tc.verifyEqual(nargin(f{1}), 3, f{1});
            end
            tc.verifyError(@() prestus_pipeline_start('abc'), 'prestus:invalidParameters');
        end

        %% Redundant expected_focal_distance_mm

        function test_expected_focal_distance_redundant_with_transducers(tc)
            c = struct('expected_focal_distance_mm', 64, ...
                       'transducer', [struct('focal_distance_ep', 85), struct('focal_distance_ep', 85)]);
            tc.verifyWarning(@() upgrade_config_keys(c, 'x'), 'prestus:legacyRedundantKey');
            u = tc.quiet(@() upgrade_config_keys(c, 'x'));
            tc.verifyFalse(isfield(u, 'expected_focal_distance_mm'));
            tc.verifyEqual([u.transducer.focal_distance_ep], [85 85]);
        end

        function test_expected_focal_distance_without_transducer_distance_stops(tc)
            c = struct('expected_focal_distance_mm', 64, ...
                       'transducer', [struct('focal_distance_ep', 85), struct('focal_distance_ep', [])]);
            tc.verifyError(@() tc.quiet(@() upgrade_config_keys(c, 'x')), 'prestus:legacyConfigKey');
        end

        %% v0.5.0 entry functions

        function test_v050_entry_functions_exist(tc)
            for f = {'single_subject_pipeline', 'single_subject_pipeline_with_slurm', 'single_subject_pipeline_with_qsub'}
                tc.verifyNotEmpty(which(f{1}), f{1});
            end
        end

        function test_v050_adopted_heatmap_stops_with_instructions(tc)
            tc.verifyError(@() single_subject_pipeline(1, struct(), adopted_heatmap=ones(2,2,2)), ...
                'prestus:legacyFunction');
        end
    end

    methods (Access = private)
        function out = quiet(~, f)
            ws = warning('off', 'all');
            c = onCleanup(@() warning(ws));
            out = f();
        end

        function [a, b] = quiet2(~, f)
            ws = warning('off', 'all');
            c = onCleanup(@() warning(ws));
            [a, b] = f();
        end
    end
end
