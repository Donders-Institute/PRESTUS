classdef test_transducer_model < matlab.unittest.TestCase
% TEST_TRANSDUCER_MODEL  Unit tests for the static transducer model and its fit.
%
% Tags:
%   'model' — forward-model parity and library fields (fast, no optimisation)
%   'fit'   — synthetic parameter recovery with FIT_TRANSDUCER_MODEL (runs lsqnonlin)
%
%   results = runtests('tests/test_transducer_model.m', 'Tag', 'model');
%   results = runtests('tests/test_transducer_model.m', 'Tag', 'fit');

    properties
        equip_param
    end

    methods (TestClassSetup)
        function setup_paths(tc)
            repo_root = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(repo_root, 'functions')));
            addpath(genpath(fullfile(repo_root, 'external')));
            tc.equip_param = load_equipment_config();
        end
    end

    methods (Test, TestTags = {'model'})

        function test_transfer_matrix_matches_rayleigh(tc)
            tran = tc.equip_param.trans.CTX_250_001;
            n    = tran.transducer.annular.elem_n;
            z    = (10:0.5:90)';
            ph   = linspace(0, 2*pi, n + 1); ph = ph(1:n);
            p.transducer = tran.transducer;
            p.medium_properties.water.sound_speed = 1500;
            p.medium_properties.water.density     = 994;
            I_ref = rayleigh_axial_intensity(ph, 0.05, p, z);
            H     = axial_transfer_matrix(tran, z, 1500, 994, 'rayleigh');
            I     = abs(H * (0.05 * exp(1i * ph(:)))).^2 / (2 * 994 * 1500) * 1e-4;
            tc.verifyEqual(I, I_ref, 'RelTol', 1e-10);
        end

        function test_oneil_transfer_matrix_matches_kwave(tc)
            % per-ring columns must be complex so that H*w reproduces focusedAnnulusONeil
            tran = tc.equip_param.trans.CTX_250_001;
            ann  = tran.transducer.annular;
            n    = ann.elem_n;
            z    = (10:0.5:90)';
            ph   = linspace(0.3, 5.9, n);
            ref  = focusedAnnulusONeil(ann.curv_radius_mm * 1e-3, [ann.elem_id_mm(:)'; ann.elem_od_mm(:)'] * 1e-3, ...
                0.05 * ones(1, n), ph, tran.transducer.freq_hz, 1500, 994, z * 1e-3);
            H    = axial_transfer_matrix(tran, z, 1500, 994, 'oneil');
            tc.verifyEqual(abs(H * (0.05 * exp(1i * ph(:)))), ref(:), 'RelTol', 1e-10);
        end

        function test_oneil_and_rayleigh_place_focus_alike(tc)
            % Both forward models take distances from the bowl apex, so for the
            % same geometric drive they must put the axial peak and FLHM centre
            % at the same depth (their amplitudes differ; positions do not).
            for serial = {'CTX_250_001', 'CTX_500_026', 'IS_PCD15473_01001'}
                tran = tc.equip_param.trans.(serial{1});
                ann  = tran.transducer.annular;
                tran.n_elem = ann.elem_n;
                bowl_off = ann.curv_radius_mm - ann.dist_geom_ep_mm;
                for d = ann.dist_geom_ep_mm * [0.7 1 1.2]
                    z_ep = (10:0.02:2*d)';
                    ph   = compute_phases(1500, tran, d);
                    w    = 0.05 * exp(1i * ph(:) * pi/180);
                    I    = @(m) abs(axial_transfer_matrix(tran, z_ep + bowl_off, 1500, 994, m) * w).^2;
                    I_r  = I('rayleigh');
                    I_o  = I('oneil');
                    [~, k_r] = max(I_r);
                    [~, k_o] = max(I_o);
                    msg = sprintf('%s at %.1f mm', serial{1}, d);
                    tc.verifyEqual(z_ep(k_o), z_ep(k_r), 'AbsTol', 0.1, ['peak, ' msg]);
                    tc.verifyEqual(compute_flhm_center(I_o / max(I_o), z_ep), ...
                        compute_flhm_center(I_r / max(I_r), z_ep), 'AbsTol', 0.1, ['FLHM centre, ' msg]);
                end
            end
        end

        function test_steering_paths_match_compute_phases(tc)
            for serial = {'CTX_250_001', 'IS_PCD15473_01001'}
                tran = tc.equip_param.trans.(serial{1});
                tran.n_elem = tran.transducer.annular.elem_n;
                for d = [20 45 70]
                    ref = compute_phases(1480, tran, d);
                    ph  = mod(360 * tran.transducer.freq_hz * steering_path_lengths(tran, d) / 1480, 360);
                    tc.verifyEqual(angle(exp(1i * (ph - ref) * pi/180)), zeros(size(ref)), 'AbsTol', 1e-9);
                end
            end
        end

        function test_evaluate_only_does_not_fit(tc)
            tran = tc.equip_param.trans.CTX_250_001;
            prof = tc.synthetic_profiles(tran, struct(), [20 35 50]);
            m = fit_transducer_model(tran, prof, struct('EvaluateOnly', true));
            tc.verifyTrue(isnan(m.exitflag));
            tc.verifyEqual(m.elem_phase_correction_deg, zeros(1, tran.transducer.annular.elem_n), 'AbsTol', 1e-12);
            tc.verifyGreaterThan(min(m.metrics.r2), 0.999, 'geometric model scored against its own profiles');
        end

        function test_loader_applies_static_model_fields(tc)
            tran = tc.equip_param.trans.CTX_250_001;
            n    = tran.transducer.annular.elem_n;
            lib_dir = tempname(); mkdir(lib_dir);
            data.meta.tran_serial = 'CTX_250_001';
            data.global_model.elem_phase_correction_deg = zeros(1, n);
            data.global_model.elem_amp_weight   = [0.8 1.0 1.1 1.1];
            data.global_model.steering_sound_speed = 1480;
            data.global_model.amplitude_scaling.i30 = 100000;
            data.global_model.amplitude_by_depth = struct('depths_ep_mm', [20 40], ...
                'elem_p_amp', [80000 120000], 'isppa_wcm2', 30);
            yaml.dumpFile(fullfile(lib_dir, 'CTX_250_001.yaml'), data);
            out = load_transducer_from_library('CTX_250_001', 30, 30, tc.equip_param, lib_dir);
            tran.n_elem = n;
            tc.verifyEqual(mod(out.transducer.annular.elem_phase_deg - compute_phases(1480, tran, 30) + 180, 360) - 180, ...
                zeros(1, n), 'AbsTol', 1e-9, 'steering at steering_sound_speed');
            tc.verifyEqual(out.transducer.annular.elem_p_amp, round(100000 * [0.8 1.0 1.1 1.1]), ...
                'amplitude interpolated in depth and weighted per element');
            out120 = load_transducer_from_library('CTX_250_001', 30, 120, tc.equip_param, lib_dir);
            tc.verifyEqual(out120.transducer.annular.elem_p_amp, round(200000 * [0.8 1.0 1.1 1.1]), ...
                'AbsTol', 1, 'sqrt(I) scaling');
            rmdir(lib_dir, 's');
        end
        function test_nearfield_excess_detected(tc)
            % measured near field suppressed → the (unchanged) model is hotter there
            tran = tc.equip_param.trans.CTX_250_001;
            prof = tc.synthetic_profiles(tran, struct(), [20 35]);
            nf = prof.z_ep_mm >= 2 & prof.z_ep_mm < 10;
            prof.I(nf, :) = 0.01 * prof.I(nf, :);
            m = fit_transducer_model(tran, prof, struct('EvaluateOnly', true));
            tc.verifyTrue(all(isfinite(m.metrics.nearfield_excess)));
            if max(prof.I(nf, :) / 0.01, [], 'all') > 0.2 * max(prof.I(~nf, :), [], 'all')
                tc.verifyGreaterThan(max(m.metrics.nearfield_excess), 0.1);
                tc.verifyTrue(any(strcmp(m.flags, 'nearfield_excess')));
            end
        end
    end

    methods (Test, TestTags = {'fit'})

        function test_recovers_synthetic_static_model(tc)
            tran = tc.equip_param.trans.CTX_250_001;
            truth.elem_phase_correction_deg = [0 25 -30 40];
            truth.elem_amp_weight           = [0.85 1.05 1.15 0.95];
            truth.z_shift_mm                = 1.2;
            depths = [16 20 25 31 36 40 45 50 55 60];
            prof = tc.synthetic_profiles(tran, truth, depths);
            m = fit_transducer_model(tran, prof, struct('FreeParams', {{'phase', 'amp', 'z_shift'}}, 'NStarts', 8));
            tc.verifyGreaterThan(min(m.metrics.r2), 0.995);
            tc.verifyEqual(m.z_shift_mm, truth.z_shift_mm, 'AbsTol', 0.2);
            d = angle(exp(1i * (m.elem_phase_correction_deg - truth.elem_phase_correction_deg) * pi/180)) * 180/pi;
            tc.verifyLessThan(max(abs(d - mean(d))), 5, 'phases recovered up to a common offset');
        end
    end

    methods (Access = private)
        function prof = synthetic_profiles(~, tran, truth, depths)
            % Measured-like profiles from the model itself (exit-plane grid, 1 mm)
            ann = tran.transducer.annular;
            n   = ann.elem_n;
            ph  = zeros(1, n); g = ones(1, n); dz = 0;
            if isfield(truth, 'elem_phase_correction_deg'); ph = truth.elem_phase_correction_deg * pi/180; end
            if isfield(truth, 'elem_amp_weight');           g  = truth.elem_amp_weight; end
            if isfield(truth, 'z_shift_mm');                dz = truth.z_shift_mm; end
            z_ep = (0:1:110)';
            H    = axial_transfer_matrix(tran, z_ep + ann.curv_radius_mm - ann.dist_geom_ep_mm + dz, 1500, 994);
            dist = steering_path_lengths(tran, depths);
            k    = 2 * pi * tran.transducer.freq_hz / 1500;
            I    = zeros(numel(z_ep), numel(depths));
            for j = 1:numel(depths)
                p = H * (0.05 * g(:) .* exp(1i * (k * dist(j, :)' + ph(:))));
                I(:, j) = abs(p).^2 / (2 * 994 * 1500) * 1e-4;
            end
            prof = struct('depths_ep_mm', depths, 'z_ep_mm', z_ep, 'I', I);
        end
    end
end
