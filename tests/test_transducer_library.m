classdef test_transducer_library < matlab.unittest.TestCase
% TEST_TRANSDUCER_LIBRARY  Unit tests for transducer serial resolution.
%
% Tests resolve_transducer_from_serial, the serial→geometry→library
% pipeline path in load_transducer_parameters, and combo-name logic.
% Does not require k-Wave or SimNIBS; uses the real equipment YAMLs.
%
%   Run with:   results = runtests('tests/test_transducer_library.m');
%   Run tag:    results = runtests('tests/test_transducer_library.m', 'Tag', 'geometry');

    properties
        equip_param      % loaded once per test class
        library_path     % path to config/calibration/
        repo_root        % PRESTUS root
        tmp_lib_dir      % temp dir for synthetic calibration YAMLs
    end

    methods (TestClassSetup)
        function setup_paths(tc)
            test_dir       = fileparts(mfilename('fullpath'));
            tc.repo_root   = fileparts(test_dir);
            addpath(genpath(fullfile(tc.repo_root, 'functions')));
            addpath(genpath(fullfile(tc.repo_root, 'external')));
            tc.equip_param  = load_equipment_config();
            tc.library_path = fullfile(tc.repo_root, 'config', 'calibration');
        end
    end

    methods (TestClassTeardown)
        function cleanup_tmp(tc)
            if ~isempty(tc.tmp_lib_dir) && isfolder(tc.tmp_lib_dir)
                rmdir(tc.tmp_lib_dir, 's');
            end
        end
    end

    % ------------------------------------------------------------------ %
    %% load_equipment_config
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'equipment'})

        function test_equipment_config_loads(tc)
            eq = tc.equip_param;
            tc.verifyClass(eq, 'struct');
            tc.verifyTrue(isfield(eq, 'trans'),  'Missing eq.trans');
            tc.verifyTrue(isfield(eq, 'combos'), 'Missing eq.combos');
        end

        function test_known_serials_present(tc)
            serials = fieldnames(tc.equip_param.trans);
            tc.verifyFalse(isempty(serials), 'No transducers found in equipment config');
        end

        function test_transducer_has_geometry(tc)
            serials = fieldnames(tc.equip_param.trans);
            serial  = serials{1};
            tran    = tc.equip_param.trans.(serial);
            tc.verifyTrue(isfield(tran, 'transducer'), ...
                sprintf('%s missing .transducer geometry field', serial));
        end

    end

    % ------------------------------------------------------------------ %
    %% resolve_transducer_from_serial — geometry only
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'geometry'})

        function test_geometry_populated_from_serial(tc)
            serial = tc.first_annular_serial();
            tr     = tc.minimal_tr(serial);
            out    = resolve_transducer_from_serial(tr, tc.equip_param, tc.library_path, 1);
            tc.verifyTrue(isfield(out, 'type'), 'type not populated');
            tc.verifyTrue(isfield(out, 'freq_hz'), 'freq_hz not populated');
            tc.verifyTrue(isfield(out, 'annular'), 'annular not populated');
            tc.verifyTrue(isfield(out.annular, 'curv_radius_mm'), 'curv_radius_mm missing');
        end

        function test_inline_fields_override_equipment(tc)
            serial   = tc.first_annular_serial();
            tr       = tc.minimal_tr(serial);
            tr.freq_hz = 123456;    % intentionally wrong — should survive merge
            out = resolve_transducer_from_serial(tr, tc.equip_param, tc.library_path, 1);
            tc.verifyEqual(out.freq_hz, 123456, ...
                'Inline freq_hz should not be overwritten by equipment YAML');
        end

        function test_unknown_serial_errors(tc)
            tr.serial = 'NONEXISTENT_DEVICE_XYZ';
            tc.verifyError( ...
                @() resolve_transducer_from_serial(tr, tc.equip_param, tc.library_path, 1), ...
                ?MException);
        end

    end

    % ------------------------------------------------------------------ %
    %% combo-name resolution
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'combo'})

        function test_generic_combo_name(tc)
            % No ds_serial → key is just the serial
            serial = tc.first_annular_serial();
            tr     = tc.minimal_tr(serial);
            name   = tc.get_combo_name(tr, serial);
            tc.verifyEqual(name, serial);
        end

        function test_ds_specific_combo_name(tc)
            serial        = tc.first_annular_serial();
            tr            = tc.minimal_tr(serial);
            tr.combo.ds_serial = 'MY_DS';
            name = tc.get_combo_name(tr, serial);
            tc.verifyEqual(name, [serial '_MY_DS']);
        end

    end

    % ------------------------------------------------------------------ %
    %% library lookup with a synthetic YAML
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'library'})

        function test_phases_resolved_from_library(tc)
            serial = tc.first_annular_serial();
            tran   = tc.equip_param.trans.(serial);
            n_elem = tran.transducer.annular.elem_n;

            lib_dir = tc.get_tmp_lib();
            tc.write_synthetic_library(lib_dir, serial, n_elem, 60, 32);

            tr = tc.minimal_tr_with_targets(serial, 60, 32);
            out = resolve_transducer_from_serial(tr, tc.equip_param, lib_dir, 1);

            tc.verifyTrue(isfield(out.annular, 'elem_phase_deg'), ...
                'elem_phase_deg not resolved');
            tc.verifyEqual(numel(out.annular.elem_phase_deg), n_elem, ...
                'Wrong number of phase elements');
            tc.verifyTrue(isfield(out.annular, 'elem_p_amp'), 'elem_p_amp not resolved');
        end

        function test_inline_phases_not_overwritten_by_library(tc)
            serial = tc.first_annular_serial();
            tran   = tc.equip_param.trans.(serial);
            n_elem = tran.transducer.annular.elem_n;

            lib_dir = tc.get_tmp_lib();
            tc.write_synthetic_library(lib_dir, serial, n_elem, 60, 32);

            tr = tc.minimal_tr_with_targets(serial, 60, 32);
            tr.annular.elem_phase_deg = zeros(1, n_elem);  % inline override
            tr.annular.elem_p_amp       = ones(1, n_elem) * 99999;

            out = resolve_transducer_from_serial(tr, tc.equip_param, lib_dir, 1);

            tc.verifyEqual(out.annular.elem_phase_deg, zeros(1, n_elem), ...
                'Inline phases were overwritten');
            tc.verifyEqual(out.annular.elem_p_amp, ones(1, n_elem)*99999, ...
                'Inline amplitude was overwritten');
        end

        function test_missing_library_errors_with_message(tc)
            serial  = tc.first_annular_serial();
            empty_dir = tc.get_tmp_lib();  % tmp dir with no YAMLs
            tr = tc.minimal_tr_with_targets(serial, 60, 32);
            tc.verifyError( ...
                @() resolve_transducer_from_serial(tr, tc.equip_param, empty_dir, 1), ...
                ?MException);
        end

        function test_no_lookup_when_targets_absent(tc)
            % Without focal_distance_ep or target_isppa_wcm2, resolution
            % should silently skip the library step (no error).
            serial  = tc.first_annular_serial();
            empty_dir = tc.get_tmp_lib();
            tr = tc.minimal_tr(serial);  % no targets set
            tc.verifyWarningFree( ...
                @() resolve_transducer_from_serial(tr, tc.equip_param, empty_dir, 1));
        end

    end

    % ------------------------------------------------------------------ %
    %% static offset + optional depth-specific refinement
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'library'})

        function test_static_offset_added_to_geometric(tc)
            [serial, tran, n] = tc.first_annular();
            corr = mod((1:n) * 7, 360);
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, corr, struct());
            out = load_transducer_from_library(serial, 40, 30, tc.equip_param, lib_dir);
            tc.verifyEqual(tc.wrap(out.transducer.annular.elem_phase_deg - tc.geo(tran, 40)), ...
                tc.wrap(corr), 'AbsTol', 1e-9);
            tc.verifyEqual(char(out.transducer.tx_calibration.phase_model), 'static');
        end

        function test_refinement_exact_at_knots_and_blends_outside(tc)
            [serial, tran, n] = tc.first_annular();
            corr = zeros(1, n);
            ref.depths_ep_mm = [40, 50];
            ref.elem_phase_correction_deg = {repmat(20, 1, n), repmat(40, 1, n)};
            ref.blend_mm = 5;
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, corr, ref);
            off = @(d) tc.wrap(tc.phases_at(serial, d, lib_dir) - tc.geo(tran, d));
            tc.verifyEqual(off(40), repmat(20, 1, n), 'AbsTol', 1e-9);   % knot
            tc.verifyEqual(off(45), repmat(30, 1, n), 'AbsTol', 1e-9);   % between knots
            tc.verifyEqual(off(52.5), repmat(20, 1, n), 'AbsTol', 1e-6); % half-way fade: angle(0.5+0.5e^{i40°}) = 20°
            tc.verifyEqual(off(60), zeros(1, n), 'AbsTol', 1e-9);        % beyond blend → static
        end

        function test_refinement_drive_and_gains_used(tc)
            % Per-depth drives and gains: exact at knots, interpolated between,
            % blended back to the global amplitude outside, scaled as sqrt(I)
            [serial, ~, n] = tc.first_annular();
            g1 = ones(1, n); g1(1) = 0.5; g1 = g1 / mean(g1);
            ref.depths_ep_mm = [40, 50];
            ref.elem_phase_correction_deg = {zeros(1, n), zeros(1, n)};
            ref.elem_p_amp = [80000, 120000];
            ref.isppa_wcm2 = 30;
            ref.elem_amp_weight = {g1, ones(1, n)};
            ref.blend_mm = 5;
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, zeros(1, n), ref);   % global i30 = 100000
            amp = @(d, I) load_transducer_from_library(serial, d, I, tc.equip_param, lib_dir);
            o40 = amp(40, 30);
            tc.verifyEqual(o40.transducer.annular.elem_p_amp, round(80000 * g1), 'AbsTol', 1);
            tc.verifyEqual(char(o40.transducer.tx_calibration.amplitude_source), 'refinement');
            g45 = (g1 + 1) / 2; g45 = g45 / mean(g45);
            tc.verifyEqual(amp(45, 30).transducer.annular.elem_p_amp, round(100000 * g45), 'AbsTol', 1);
            tc.verifyEqual(amp(60, 30).transducer.annular.elem_p_amp, repmat(100000, 1, n), 'AbsTol', 1);
            tc.verifyEqual(amp(50, 7.5).transducer.annular.elem_p_amp, repmat(60000, 1, n), 'AbsTol', 1);
        end

        function test_unverified_amplitude_warns(tc)
            % Analytical-only drive (no refinement.amp_source, or 'analytical'):
            % amp_kwave_verified is false and the loader warns.
            [serial, ~, n] = tc.first_annular();
            ref.depths_ep_mm = 40;
            ref.elem_phase_correction_deg = {zeros(1, n)};
            ref.elem_p_amp = 80000;
            ref.isppa_wcm2 = 30;
            ref.amp_source = {'analytical'};
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, zeros(1, n), ref);
            tc.verifyWarning(@() load_transducer_from_library(serial, 40, 30, tc.equip_param, lib_dir), ...
                'load_transducer_from_library:unverifiedAmplitude');
            w = warning('off', 'load_transducer_from_library:unverifiedAmplitude');
            out = load_transducer_from_library(serial, 40, 30, tc.equip_param, lib_dir);
            warning(w);
            tc.verifyFalse(out.transducer.tx_calibration.amp_kwave_verified);
        end

        function test_kwave_verified_amplitude_silent(tc)
            % A refinement whose amp_source is 'kwave', resolved exactly at its
            % own depth (no blend with the analytical global model): verified,
            % no warning.
            [serial, ~, n] = tc.first_annular();
            ref.depths_ep_mm = 40;
            ref.elem_phase_correction_deg = {zeros(1, n)};
            ref.elem_p_amp = 80000;
            ref.isppa_wcm2 = 30;
            ref.amp_source = {'kwave'};
            ref.kwave_r2 = 0.98;
            ref.kwave_focal_err_mm = 0.3;
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, zeros(1, n), ref);
            tc.verifyWarningFree(@() load_transducer_from_library(serial, 40, 30, tc.equip_param, lib_dir));
            out = load_transducer_from_library(serial, 40, 30, tc.equip_param, lib_dir);
            tc.verifyTrue(out.transducer.tx_calibration.amp_kwave_verified);
            tc.verifyEqual(out.transducer.tx_calibration.amp_kwave_r2, 0.98);
        end

        function test_blend_with_kwave_refinement_stays_unverified(tc)
            % Between the refined knot and the blend edge, the drive is a mix of
            % the (verified) refinement and the (always-analytical) global
            % model — amp_kwave_verified must be false there even though the
            % nearest refinement is 'kwave', since the blended amplitude is not.
            [serial, ~, n] = tc.first_annular();
            ref.depths_ep_mm = 40;
            ref.elem_phase_correction_deg = {zeros(1, n)};
            ref.elem_p_amp = 80000;
            ref.isppa_wcm2 = 30;
            ref.amp_source = {'kwave'};
            ref.blend_mm = 5;
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, zeros(1, n), ref);   % global i30 = 100000
            out_knot  = load_transducer_from_library(serial, 40, 30, tc.equip_param, lib_dir);
            out_blend = load_transducer_from_library(serial, 42, 30, tc.equip_param, lib_dir);
            tc.verifyTrue(out_knot.transducer.tx_calibration.amp_kwave_verified);
            tc.verifyFalse(out_blend.transducer.tx_calibration.amp_kwave_verified);
        end

        function test_model_at_depth_matches_loader(tc)
            % LIBRARY_MODEL_AT_DEPTH on the in-memory entry = the file-based loader
            [serial, tran, n] = tc.first_annular();
            ref.depths_ep_mm = [40, 50];
            ref.elem_phase_correction_deg = {repmat(20, 1, n), repmat(40, 1, n)};
            ref.elem_p_amp = [80000, 120000];
            ref.isppa_wcm2 = 30;
            ref.blend_mm = 5;
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, mod((1:n) * 7, 360), ref);
            lib = yaml.loadFile(fullfile(lib_dir, [serial '.yaml']), 'ConvertToArray', true);
            for d = [35 40 45 52.5 60]
                out = load_transducer_from_library(serial, d, 25, tc.equip_param, lib_dir);
                m   = library_model_at_depth(lib, tran, d, 25);
                tc.verifyEqual(m.elem_phase_deg, out.transducer.annular.elem_phase_deg, 'AbsTol', 1e-12);
                tc.verifyEqual(m.elem_p_amp, out.transducer.annular.elem_p_amp);
                tc.verifyEqual(m.phase_model, char(out.transducer.tx_calibration.phase_model));
            end
        end

        function test_evaluate_entry_at_unmeasured_depths(tc)
            % Runtime-path prediction: refined at the knot, blended just outside,
            % static beyond blend_mm; metrics only where a profile was measured
            [serial, tran, n] = tc.first_annular();
            d0 = tran.transducer.annular.dist_geom_ep_mm;
            lib.meta.tran_serial = serial;
            lib.global_model.elem_phase_correction_deg = zeros(1, n);
            lib.global_model.fit_z_shift_mm = 2;            % diagnostic: must not be applied
            lib.refinement = struct('depths_ep_mm', d0, ...
                'elem_phase_correction_deg', {{mod((1:n) * 30, 360)}}, 'blend_mm', 5);
            z = (0:0.25:1.5 * (d0 + 10))';
            prof = struct('depths_ep_mm', [d0, d0 + 2.5, d0 + 10], 'z_ep_mm', z, ...
                'I', [exp(-((z - d0) / 5).^2), nan(numel(z), 2)]);
            ev = evaluate_library_entry(tran, lib, prof, struct('NearFieldSkipMm', 0));
            tc.verifyEqual({ev.settings.phase_model}, {'refined', 'blend', 'static'});
            tc.verifyEqual(ev.measured, [true false false]);
            tc.verifyTrue(isfinite(ev.metrics.r2(1)));
            tc.verifyTrue(all(isnan(ev.metrics.r2(2:3))));
            tc.verifyTrue(all(isnan(ev.I_meas(:, 2:3)), 'all'));
            % Static depth equals geometric steering (zero offset, no z shift)
            geo = fit_transducer_model(tran, struct('depths_ep_mm', d0 + 10, 'z_ep_mm', z, ...
                'I', ones(numel(z), 1)), struct('EvaluateOnly', true, 'NearFieldSkipMm', 0));
            tc.verifyEqual(ev.I_model(:, 3), geo.I_model, 'AbsTol', 1e-12);
            tc.verifyEqual(ev.metrics.model_focus_ep_mm(3), ...
                compute_flhm_center(geo.I_model, geo.z_ep_mm), 'AbsTol', 1e-12);
        end

        function test_interpolated_profile_lands_at_desired_depth(tc)
            % Focus-aligned blend: the interpolated focus is at the desired
            % depth (the legacy EXTRACT_REAL_INTENSITY_PROFILE put it at a
            % fixed position for every desired depth)
            z = (0:0.25:120)';
            g = @(c, w) exp(-((z - c) / w).^2);
            meas = struct('depths_ep_mm', [30 50], 'z_ep_mm', z, 'I', 20 * [g(30, 6), g(50, 8)]);
            for align = {'label', 'peak'}
                p = interpolate_axial_profiles(meas, [30 35 40 45 50], struct('Align', align{1}));
                [~, k] = max(p.I);
                tc.verifyEqual(z(k)', [30 35 40 45 50], 'AbsTol', 0.25, align{1});
            end
            tc.verifyEqual(p.profile_type, {'measured', 'interpolated', 'interpolated', 'interpolated', 'measured'});
            tc.verifyEqual(p.I(:, 1), meas.I(:, 1));                    % measured column passed through
            tc.verifyEqual(p.weights(3, :), [0.5 0.5], 'AbsTol', 1e-12);
            tc.verifyEqual(p.source_depths(2, :), [30 50]);
            % FWHM blends between the two widths
            fw = @(y) nnz(y >= max(y) / 2) * 0.25;
            tc.verifyTrue(fw(p.I(:, 3)) > fw(meas.I(:, 1)) && fw(p.I(:, 3)) < fw(meas.I(:, 2)));
            % Legacy wrapper returns the same profile on the measured grid
            par.calibration = struct('skip_front_peak_mm', 10, 'path_output_profiles', tc.fresh_dir());
            leg = extract_real_intensity_profile(par, [30 50], 40, meas.I, 'test', z);
            tc.verifyEqual(leg(:), p.I(:, 3), 'AbsTol', 1e-12);
            tc.verifyNotEmpty(dir(fullfile(par.calibration.path_output_profiles, 'Profile_interpolated_at_*.csv')));
        end

        function test_profile_extrapolation_is_opt_in(tc)
            z = (0:0.5:120)';
            meas = struct('depths_ep_mm', [30 50], 'z_ep_mm', z, 'I', [exp(-((z - 30) / 6).^2), exp(-((z - 50) / 8).^2)]);
            p0 = interpolate_axial_profiles(meas, [25 55]);
            tc.verifyEqual(p0.profile_type, {'unavailable', 'unavailable'});
            tc.verifyTrue(all(isnan(p0.I), 'all'));
            p = interpolate_axial_profiles(meas, [25 55 70], struct('MaxExtrapolationMm', 5));
            tc.verifyEqual(p.profile_type, {'extrapolated', 'extrapolated', 'unavailable'});
            [~, k] = max(p.I(:, 1:2));
            tc.verifyEqual(z(k)', [25 55], 'AbsTol', 0.5);
            tc.verifyEqual(p.shift_mm(1:2, 1)', [-5 5]);
            par.calibration = struct('skip_front_peak_mm', 10, 'path_output_profiles', tc.fresh_dir());
            tc.verifyError(@() extract_real_intensity_profile(par, [30 50], 25, meas.I, 'test', z), ?MException);
        end

        function test_write_axial_profiles_round_trip(tc)
            z = (0:0.5:100)';
            meas = struct('depths_ep_mm', [30 50], 'z_ep_mm', z, 'I', 10 * [exp(-((z - 30) / 6).^2), exp(-((z - 50) / 8).^2)]);
            p = interpolate_axial_profiles(meas, [30 42.5 60]);           % 60: unavailable, not written
            f = fullfile(tc.fresh_dir(), 'profiles.csv');
            write_axial_profiles(p, f, struct('combo', 'X'));
            back = read_axial_profiles(f);
            tc.verifyEqual(back.depths_ep_mm, [30 42.5]);
            tc.verifyEqual(back.z_ep_mm, z);
            tc.verifyEqual(back.I, p.I(:, 1:2), 'AbsTol', 1e-8);
            sc = yaml.loadFile(strrep(f, '.csv', '.yaml'), 'ConvertToArray', true);
            tc.verifyEqual(string(sc.profile_type), ["measured", "interpolated"]);
            dv = sc.derived; if iscell(dv); dv = dv{1}; end
            tc.verifyEqual(double(dv.source_depths_ep_mm), [30 50]);
        end

        function test_refinement_records_profile_type(tc)
            [~, ~, n] = tc.first_annular();
            f = struct('fit_depths_ep_mm', 42.5, 'elem_phase_correction_deg', zeros(1, n), ...
                'drive_velocity', 0.05, 'measured_peak_wcm2', 30, 'z_shift_mm', 0, ...
                'free_params', {{'phase'}}, 'forward_model', 'rayleigh', 'profile_type', 'interpolated');
            lib = library_refinement_from_fits(struct(), {f}, 30);
            tc.verifyEqual(lib.refinement.profile_type, {'interpolated'});
            lib = library_add_refinements(lib, 50, {zeros(1, n)});        % no type recorded for this depth
            tc.verifyEqual(lib.refinement.profile_type, {'interpolated', ''});
        end

        function test_settings_table_matches_runtime(tc)
            [serial, tran, n] = tc.first_annular();
            ref.depths_ep_mm = [40, 50];
            ref.elem_phase_correction_deg = {repmat(20, 1, n), repmat(40, 1, n)};
            ref.elem_p_amp = [80000, 120000]; ref.isppa_wcm2 = 30; ref.blend_mm = 5;
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, zeros(1, n), ref);
            lib = yaml.loadFile(fullfile(lib_dir, [serial '.yaml']), 'ConvertToArray', true);
            csv = fullfile(lib_dir, 'settings.csv');
            T = library_settings_table(tran, lib, [35 45 52.5], struct('IsppaWcm2', 20, 'CsvPath', csv));
            tc.verifyEqual(T.phase_model', ["static", "refined", "blend"]);
            for k = 1:3
                out = load_transducer_from_library(serial, T.depth_ep_mm(k), 20, tc.equip_param, lib_dir);
                tc.verifyEqual(T{k, compose('elem_phase_deg_%d', 1:n)}, out.transducer.annular.elem_phase_deg, 'AbsTol', 1e-9);
                tc.verifyEqual(T{k, compose('elem_p_amp_pa_%d', 1:n)}, out.transducer.annular.elem_p_amp);
            end
            tc.verifyEqual(height(readtable(csv)), 3);
        end

        function test_refine_at_unmeasured_depth_exports_profiles(tc)
            % Headless: a requested depth between measured ones is calibrated
            % on its interpolated profile; profiles and settings are exported
            [serial, tran, n] = tc.first_annular();
            d0 = tran.transducer.annular.dist_geom_ep_mm;
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, zeros(1, n), struct());
            z = (0:0.25:1.6 * d0)';
            dm = d0 + [-10 0 10];
            I = zeros(numel(z), 3);
            for k = 1:3
                m = fit_transducer_model(tran, struct('depths_ep_mm', dm(k), 'z_ep_mm', z, 'I', ones(numel(z), 1)), ...
                    struct('EvaluateOnly', true, 'NearFieldSkipMm', 0));
                I(:, k) = 30 * interp1(m.z_ep_mm, m.I_model, z, 'linear', 0);
            end
            csv = fullfile(lib_dir, 'meas.csv');
            write_axial_profiles(struct('depths_ep_mm', dm, 'z_ep_mm', z, 'I', I), csv);
            [out_file, fits, rep] = refine_transducer_library(serial, csv, lib_dir, struct( ...
                'Depths', [d0 - 10, d0 - 5], 'FitOptions', struct('NStarts', 2), ...
                'EquipParam', tc.equip_param, 'ExportSettings', true, 'ShippedLibrary', lib_dir));
            tc.verifyEqual(numel(fits), 2);
            tc.verifyEqual(rep.profile_type', ["measured", "interpolated"]);
            lib = yaml.loadFile(out_file, 'ConvertToArray', true);
            tc.verifyEqual(string(lib.refinement.profile_type), ["measured", "interpolated"]);
            tc.verifyTrue(isfile(lib.refinement.calibration_profiles_csv));
            back = read_axial_profiles(lib.refinement.calibration_profiles_csv);
            tc.verifyEqual(back.depths_ep_mm, round([d0 - 10, d0 - 5], 2));
            tc.verifyTrue(isfile(fullfile(lib_dir, [serial '_settings.csv'])));
        end

        function test_add_refinements_replaces_depth_and_kwave_result(tc)
            % A refit depth replaces its drive and drops its k-Wave check
            [~, ~, n] = tc.first_annular();
            amp = @(a) struct('elem_p_amp', a, 'isppa_wcm2', 30);
            lib = library_add_refinements(struct(), [40 50], {zeros(1, n), zeros(1, n)}, [], '', '', amp([1e5 2e5]));
            lib.refinement.kwave_isppa_wcm2 = [25 35];
            lib.refinement.kwave_r2 = [0.9 0.8]; lib.refinement.kwave_focal_err_mm = [0.1 0.2];
            lib = library_add_refinements(lib, 50, {repmat(10, 1, n)}, [], '', '', ...
                struct('elem_p_amp', 3e5, 'isppa_wcm2', 120));    % other ISPPA: rescaled to 30
            r = lib.refinement;
            tc.verifyEqual(r.depths_ep_mm, [40 50]);
            tc.verifyEqual(r.elem_p_amp, [1e5 1.5e5], 'AbsTol', 1);
            tc.verifyEqual(r.kwave_isppa_wcm2(1), 25);
            tc.verifyTrue(isnan(r.kwave_isppa_wcm2(2)));
            tc.verifyEqual(string(r.amp_source), ["analytical", "analytical"]);
        end

        function test_refinement_from_fits_stores_drive_gains_and_zshift(tc)
            [~, ~, n] = tc.first_annular();
            g = linspace(0.8, 1.2, n);
            f = struct('fit_depths_ep_mm', 45, 'elem_phase_correction_deg', zeros(1, n), ...
                'elem_amp_weight', g, 'z_shift_mm', 1.5, 'free_params', {{'phase', 'amp', 'z_shift'}}, ...
                'forward_model', 'rayleigh', 'drive_velocity', 0.05, 'measured_peak_wcm2', 120);
            lib = library_refinement_from_fits(struct(), {f}, 30);
            r = lib.refinement;
            tc.verifyEqual(r.elem_p_amp, round(0.05 * 994 * 1500 * sqrt(30 / 120)), 'AbsTol', 1);
            tc.verifyEqual(r.elem_amp_weight{1}, g, 'AbsTol', 1e-12);
            tc.verifyEqual(r.fit_z_shift_mm, 1.5);
        end

        function test_correction_length_mismatch_errors(tc)
            [serial, ~, n] = tc.first_annular();
            lib_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, zeros(1, n + 1), struct());
            tc.verifyError(@() load_transducer_from_library(serial, 40, 30, tc.equip_param, lib_dir), ...
                'load_transducer_from_library:correctionLength');
        end

        function test_user_library_path_takes_precedence(tc)
            [serial, tran, n] = tc.first_annular();
            shipped = tc.fresh_dir(); user = tc.fresh_dir();
            tc.write_static_library(shipped, serial, zeros(1, n), struct());
            tc.write_static_library(user,    serial, repmat(90, 1, n), struct());
            tr = tc.minimal_tr_with_targets(serial, 40, 30);
            tr.library_path = user;
            out = resolve_transducer_from_serial(tr, tc.equip_param, shipped, 1);
            tc.verifyEqual(tc.wrap(out.annular.elem_phase_deg - tc.geo(tran, 40)), ...
                repmat(90, 1, n), 'AbsTol', 1e-9);
        end

        function test_ingested_refinement_not_double_corrected(tc)
            % A single-depth calibration run stores absolute phases; after
            % update_transducer_library the loader must reproduce them exactly
            % at that depth even though a static offset is present.
            [serial, tran, n] = tc.first_annular();
            lib_dir = tc.fresh_dir(); run_dir = tc.fresh_dir();
            tc.write_static_library(lib_dir, serial, repmat(15, 1, n), struct());
            abs_deg = mod(tc.geo(tran, 45) + (1:n) * 11, 360);
            d.meta.tran_serial = serial;
            d.transducer.annular.elem_phase_deg = abs_deg;
            d.transducer.annular.elem_p_amp     = repmat(1e5, 1, n);
            yaml.dumpFile(fullfile(run_dir, [serial '-F45mm-I30wpercm2.yaml']), d);
            update_transducer_library(serial, run_dir, lib_dir);
            lib = yaml.loadFile(fullfile(lib_dir, [serial '.yaml']), 'ConvertToArray', true);
            tc.verifyEqual(lib.global_model.elem_phase_correction_deg(:)', repmat(15, 1, n), ...
                'static offset must be left untouched');
            tc.verifyEqual(tc.wrap(tc.phases_at(serial, 45, lib_dir) - abs_deg), zeros(1, n), 'AbsTol', 1e-6);
            tc.verifyEqual(tc.wrap(tc.phases_at(serial, 70, lib_dir) - tc.geo(tran, 70)), ...
                repmat(15, 1, n), 'AbsTol', 1e-9, 'far from the refinement the static offset applies');
        end

    end

    % ------------------------------------------------------------------ %
    %% load_transducer_parameters integration
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'integration'})

        function test_serial_triggers_geometry_in_pipeline(tc)
            % Verify that resolve_transducer_from_serial (the function wired
            % into load_transducer_parameters) populates type and freq_hz from
            % the equipment YAML when only a serial is given.
            serial = tc.first_annular_serial();
            tran   = tc.equip_param.trans.(serial);
            n      = tran.transducer.annular.elem_n;

            tr.serial              = serial;
            tr.focal_distance_ep   = 60;
            tr.target_isppa_wcm2   = 5;
            tr.annular.elem_phase_deg = zeros(1, n);
            tr.annular.elem_p_amp       = ones(1, n) * 100000;

            out = resolve_transducer_from_serial(tr, tc.equip_param, tc.library_path, 1);

            tc.verifyEqual(out.freq_hz, tran.transducer.freq_hz, ...
                'freq_hz should be populated from equipment YAML');
            tc.verifyEqual(char(out.type), 'annular', ...
                'type should be set to annular from equipment YAML');
        end

    end

    % ------------------------------------------------------------------ %
    %% Private helpers
    % ------------------------------------------------------------------ %
    methods (Access = private)

        function [serial, tran, n] = first_annular(tc)
            serial = tc.first_annular_serial();
            tran   = tc.equip_param.trans.(serial);
            tran.n_elem = tran.transducer.annular.elem_n;
            n      = tran.n_elem;
        end

        function d = fresh_dir(~)
            d = tempname(); mkdir(d);
        end

        function g = geo(~, tran, depth)
            p.medium_properties.water.sound_speed = 1500;
            p.transducer = tran.transducer; p.transducer.n_elem = tran.transducer.annular.elem_n;
            g = reshape(set_real_phases(tran, depth, p), 1, []);
        end

        function ph = phases_at(tc, serial, depth, lib_dir)
            out = load_transducer_from_library(serial, depth, 30, tc.equip_param, lib_dir);
            ph  = out.transducer.annular.elem_phase_deg(:)';
        end

        function w = wrap(~, deg)
            % map to [0, 360) with round-off at 360 folded to 0
            w = mod(deg(:)', 360);
            w(w > 360 - 1e-6) = 0;
        end

        function write_static_library(~, lib_dir, serial, corr, ref)
            data.meta.tran_serial = serial;
            data.meta.calibration_type = 'test';
            data.global_model.elem_phase_correction_deg = corr;
            data.global_model.amplitude_scaling.i30 = 100000;
            if ~isempty(fieldnames(ref)); data.refinement = ref; end
            yaml.dumpFile(fullfile(lib_dir, [serial '.yaml']), data);
        end

        function serial = first_annular_serial(tc)
            serials = fieldnames(tc.equip_param.trans);
            serial  = '';
            for i = 1:numel(serials)
                tran = tc.equip_param.trans.(serials{i});
                if isfield(tran, 'transducer') && isfield(tran.transducer, 'annular')
                    serial = serials{i};
                    return;
                end
            end
            error('test_transducer_library: no annular transducer found in equipment config');
        end

        function tr = minimal_tr(~, serial)
            tr.serial = serial;
        end

        function tr = minimal_tr_with_targets(~, serial, focal_ep, isppa)
            tr.serial              = serial;
            tr.focal_distance_ep   = focal_ep;
            tr.target_isppa_wcm2   = isppa;
        end

        function d = get_tmp_lib(tc)
            if isempty(tc.tmp_lib_dir)
                tc.tmp_lib_dir = tempname();
                mkdir(tc.tmp_lib_dir);
            end
            d = tc.tmp_lib_dir;
        end

        function name = get_combo_name(~, tr, serial)
            % Mirror the logic in resolve_combo_name without calling it directly.
            if isfield(tr, 'combo') && isstruct(tr.combo) && ...
                    isfield(tr.combo, 'ds_serial') && ~isempty(tr.combo.ds_serial)
                name = [serial '_' char(tr.combo.ds_serial)];
            else
                name = serial;
            end
        end

        function write_synthetic_library(~, lib_dir, serial, n_elem, focal_mm, isppa)
            % Write a minimal library YAML to lib_dir/{serial}.yaml so that
            % load_transducer_from_library can parse it.
            fkey = sprintf('f%s', strrep(num2str(focal_mm), '.', 'p'));
            ikey = sprintf('i%s', strrep(num2str(isppa),   '.', 'p'));

            data.meta.tran_serial         = serial;
            data.meta.created_by          = 'test_transducer_library';
            data.calibration.focal_depths.(fkey).precession_mode   = 'linear';
            data.calibration.focal_depths.(fkey).phase_start_deg   = 0;
            data.calibration.focal_depths.(fkey).phase_step_deg    = 10;
            data.calibration.focal_depths.(fkey).amplitude_scaling.(ikey) = 100000;

            yaml_path = fullfile(lib_dir, [serial '.yaml']);
            yaml.dumpFile(yaml_path, data);
        end

    end

end
