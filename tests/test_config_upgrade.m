classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_config_upgrade < matlab.unittest.TestCase
% TEST_CONFIG_UPGRADE  Conversion of config keys from PRESTUS v0.5.0 onwards.
%
% Covers UPGRADE_CONFIG_KEYS and CONFIG_DEPRECATIONS: the default configs of
% v0.5.0, v0.6.0 and v0.6.1 (tests/fixtures/config_history) convert without
% unknown keys left over; renames, exact value conversions, removed keys,
% errors, conflicts, transducer arrays and key/value overrides; and
% LOAD_PARAMETERS applying the conversion to a study config. Writes only to
% temporary folders.
%
%   results = runtests('tests/test_config_upgrade.m');

    properties
        Repo
        Fixtures
        Dir
    end

    methods (TestClassSetup)
        function add_paths(tc)
            tc.Repo = fileparts(fileparts(mfilename('fullpath')));
            tc.Fixtures = fullfile(tc.Repo, 'tests', 'fixtures', 'config_history');
        end
    end

    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tempname;
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
        end
    end

    methods (Test, TestTags = {'config_upgrade'})

        %% Released default configs

        function test_v060_v061_defaults_leave_no_unknown_keys(tc)
            for v = {'v060', 'v061'}
                c = tc.release_config(v{1});
                u = tc.upgrade_quiet(c);
                tc.verifyEmpty(tc.unknown(u), v{1});
            end
        end

        function test_v050_default_stops_only_at_smooth_window(tc)
            c = tc.release_config('v050');
            msg = '';
            try
                tc.upgrade_quiet(c);
            catch err
                msg = err.message;
            end
            tc.verifySubstring(msg, 'smooth_window');
            tc.verifyEqual(numel(regexp(msg, '\n  ')), 1, 'only smooth_window cannot be converted');
            c = rmfield(c, 'smooth_window');
            u = tc.upgrade_quiet(c);
            tc.verifyEmpty(tc.unknown(u));
            tc.verifyEqual(char(u.simulation.medium), 'layered');
            tc.verifyEqual(u.grid.resolution_mm, 0.5);
            tc.verifyTrue(isfield(u.medium_properties.water, 'sound_speed'));
            tc.verifyTrue(isfield(u.timing, 'pd'));
            % kosciessa preset, with yakuub now called mueller
            tc.verifyEqual({u.pct.mapping_density, u.pct.mapping_soundspeed, u.pct.mapping_attenuation}, ...
                {'k-plan', 'marsac', 'mueller'});
            % v0.5.0 weights 1 = uniform -> opt_weights 0
            tc.verifyEqual(u.calibration.opt_weights, 0);
            tc.verifyFalse(isfield(u, 'seg_labels'));
            tc.verifyFalse(isfield(u, 'subject_subfolder'));
        end

        function test_renamed_keys_warn(tc)
            tc.verifyWarning(@() upgrade_config_keys(struct('grid_step_mm', 1), 'x'), ...
                'prestus:renamedConfigKey');
        end

        function test_current_config_is_unchanged(tc)
            c = struct('grid', struct('resolution_mm', 1), 'simulation', struct('debug', 0));
            tc.verifyWarningFree(@() upgrade_config_keys(c, 'x'));
            tc.verifyEqual(upgrade_config_keys(c, 'x'), c);
        end

        %% Renames and conflicts

        function test_group_rename(tc)
            u = tc.upgrade_quiet(struct('medium', struct('brain', struct('density', 1040))));
            tc.verifyEqual(u.medium_properties.brain.density, 1040);
            tc.verifyFalse(isfield(u, 'medium'));
        end

        function test_old_and_new_key_with_different_values_stop(tc)
            c = struct('grid_step_mm', 1, 'grid', struct('resolution_mm', 0.5));
            tc.verifyError(@() tc.upgrade_quiet(c), 'prestus:legacyConfigKey');
        end

        function test_old_and_new_key_with_same_value_merge(tc)
            u = tc.upgrade_quiet(struct('grid_step_mm', 1, 'grid', struct('resolution_mm', 1)));
            tc.verifyEqual(u.grid.resolution_mm, 1);
            tc.verifyFalse(isfield(u, 'grid_step_mm'));
        end

        function test_alias_yields_to_set_key(tc)
            u = tc.upgrade_quiet(struct('hpc_gpu', 'gpu:1', 'hcp_gpu', 'gpu:2'));
            tc.verifyEqual(u.hpc.gpu, 'gpu:1');
        end

        %% Transducers

        function test_v050_transducer_becomes_annular(tc)
            t = struct('n_elements', 2, 'Elements_ID_mm', [0 30], 'Elements_OD_mm', [29 40], ...
                'curv_radius_mm', 63, 'dist_to_plane_mm', 52, 'source_amp', 9e4, ...
                'source_phase_deg', [0 10], 'source_freq_hz', 25e4);
            u = tc.upgrade_quiet(struct('transducer', t));
            a = u.transducer.annular;
            tc.verifyEqual([a.elem_n, a.curv_radius_mm, a.dist_geom_ep_mm, a.elem_p_amp], [2 63 52 9e4]);
            tc.verifyEqual(a.elem_phase_deg, [0 10]);
            tc.verifyEqual(u.transducer.freq_hz, 25e4);
            tc.verifyEqual(u.transducer.type, 'annular');
        end

        function test_every_transducer_of_an_array_is_converted(tc)
            t = struct('source_amp', {1, 2}, 'source_freq_hz', {5e5, 2.5e5});
            u = tc.upgrade_quiet(struct('transducer', t));
            tc.verifyEqual(numel(u.transducer), 2);
            tc.verifyEqual([u.transducer(1).annular.elem_p_amp, u.transducer(2).annular.elem_p_amp], [1 2]);
            tc.verifyEqual([u.transducer.freq_hz], [5e5 2.5e5]);
        end

        function test_calibration_file_intensity_is_the_calibration_isppa(tc)
            u = tc.upgrade_quiet(struct('transducer', struct('type', 'annular', 'set_intensity_w_per_cm2', 60)));
            tc.verifyEqual(u.transducer.tx_calibration.isppa_wcm2, 60);
            tc.verifyFalse(isfield(u.transducer, 'target_isppa_wcm2'));
        end

        function test_expected_focal_distance_stops(tc)
            tc.verifyError(@() tc.upgrade_quiet(struct('expected_focal_distance_mm', 60)), ...
                'prestus:legacyConfigKey');
            tc.verifyError(@() tc.upgrade_quiet(struct('transducer', ...
                struct('expected_focal_distance_mm', 60))), 'prestus:legacyConfigKey');
        end

        %% Value conversions

        function test_output_location(tc)
            u = tc.upgrade_quiet(struct('data_path', '/data/study', 'output_location', 'sims'));
            tc.verifyEqual(u.path.sim, fullfile('/data/study', 'sims'));
            u = tc.upgrade_quiet(struct('output_location', '/abs/out'));
            tc.verifyEqual(u.path.sim, '/abs/out');
            tc.verifyError(@() tc.upgrade_quiet(struct('output_location', 'sims')), 'prestus:legacyConfigKey');
        end

        function test_pct_presets_and_custom_mapping(tc)
            u = tc.upgrade_quiet(struct('pseudoCT_variant', 'carpino'));
            tc.verifyEqual({u.pct.mapping_density, u.pct.mapping_soundspeed, u.pct.mapping_attenuation}, ...
                {'k-wave', 'k-plan', 'mueller'});
            m = struct('density', 'marsac', 'soundspeed', 'k-plan', 'attenuation', 'yakuub');
            u = tc.upgrade_quiet(struct('pseudoCT_variant', 'custom', 'pct_skullmapping', m));
            tc.verifyEqual({u.pct.mapping_density, u.pct.mapping_soundspeed, u.pct.mapping_attenuation}, ...
                {'marsac', 'k-plan', 'mueller'});
            tc.verifyError(@() tc.upgrade_quiet(struct('pseudoCT_variant', 'marquet')), 'prestus:legacyConfigKey');
            tc.verifyError(@() tc.upgrade_quiet(struct('pct_skullmapping', m)), 'prestus:legacyConfigKey');
        end

        function test_calibration_weights(tc)
            u = tc.upgrade_quiet(struct('calibration', struct('weights', 2, 'optmethod', 'FEXminimize')));
            tc.verifyEqual(u.calibration.opt_weights, 0);
            tc.verifyEqual(u.calibration.opt_method, 'FEXminimize');
            tc.verifyError(@() tc.upgrade_quiet(struct('calibration', struct('weights', 0))), ...
                'prestus:legacyConfigKey');
        end

        %% Removed keys

        function test_removed_key_with_supported_value_is_dropped(tc)
            u = tc.upgrade_quiet(struct('subject_subfolder', 1, 'segmentation_software', 'charm'));
            tc.verifyEmpty(fieldnames(u));
        end

        function test_removed_key_with_unsupported_value_stops(tc)
            tc.verifyError(@() tc.upgrade_quiet(struct('subject_subfolder', 0)), 'prestus:legacyConfigKey');
            labels = charm_seg_labels();
            labels.csf = 99;
            tc.verifyError(@() tc.upgrade_quiet(struct('seg_labels', labels)), 'prestus:legacyConfigKey');
        end

        %% Overrides and load_parameters

        function test_overrides(tc)
            ov = tc.upgrade_quiet({'grid_step_mm', 1; 'transducer(2).source_amp', 5; 'medium.brain.density', 1});
            tc.verifyEqual(ov(:, 1)', {'grid.resolution_mm', 'transducer(2).annular.elem_p_amp', ...
                'medium_properties.brain.density'});
            tc.verifyError(@() tc.upgrade_quiet({'smooth_window', 3}), 'prestus:legacyConfigKey');
        end

        function test_load_parameters_converts_a_study_config(tc)
            write_text(fullfile(tc.Dir, 'config_old.yaml'), { ...
                'simulation_medium: water', 'grid_step_mm: 1', 'interactive: 0', ...
                'run_heating_sims: 0', 'subject_id: 3'});
            tc.verifyWarning(@() tc.load_quiet_except('config_old.yaml', 'prestus:renamedConfigKey'), ...
                'prestus:renamedConfigKey');
            p = tc.load_quiet('config_old.yaml');
            tc.verifyEqual(char(p.simulation.medium), 'water');
            tc.verifyEqual(p.grid.resolution_mm, 1);
            tc.verifyEqual(p.modules.run_heating_sims, 0);
        end
    end

    methods (Access = private)
        function c = release_config(tc, v)
            c = yaml.loadFile(fullfile(tc.Fixtures, ['config_default_' v '.yaml']), 'ConvertToArray', true);
            cal = yaml.loadFile(fullfile(tc.Fixtures, ['config_calibration_' v '.yaml']), 'ConvertToArray', true);
            if isfield(cal, 'calibration'), cal = cal.calibration; end
            c.calibration = cal;
        end

        function unk = unknown(tc, cfg)
            known = yaml.loadFile(fullfile(tc.Repo, 'config', 'config_default.yaml'), 'ConvertToArray', true);
            known.calibration = yaml.loadFile(fullfile(tc.Repo, 'config', 'config_calibration.yaml'), ...
                'ConvertToArray', true);
            unk = config_unknown_keys(cfg, known, config_known_extra_keys(tc.Repo));
        end

        function out = upgrade_quiet(~, cfg)
            s = warning('off', 'prestus:renamedConfigKey');
            restore = onCleanup(@() warning(s));
            out = upgrade_config_keys(cfg, 'test');
        end

        function p = load_quiet(tc, file)
            s = warning('off', 'all');
            restore = onCleanup(@() warning(s));
            p = load_parameters(file, tc.Dir);
        end

        function p = load_quiet_except(tc, file, id)
            state = warning('query', id);
            s = warning('off', 'all');
            warning(state.state, id);
            restore = onCleanup(@() warning(s));
            p = load_parameters(file, tc.Dir);
        end
    end
end

function write_text(path, lines)
    fid = fopen(path, 'w');
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end
