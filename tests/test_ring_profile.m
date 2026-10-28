classdef test_ring_profile < matlab.unittest.TestCase
% TEST_RING_PROFILE  Within-ring velocity profile: geometry, fit, library round trip.
%
% Analytical only (no k-Wave); uses the real equipment YAMLs.
%
%   results = runtests('tests/test_ring_profile.m');

    properties
        Equip
        Tran        % CTX_500_026 (4 rings)
        Tmp
    end

    methods (TestClassSetup)
        function setup(tc)
            root = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(root, 'functions'))); addpath(genpath(fullfile(root, 'external')));
            tc.Equip = load_equipment_config();
            tc.Tran  = tc.Equip.trans.CTX_500_026;
            tc.Tmp   = tempname; mkdir(tc.Tmp);
        end
    end

    methods (TestClassTeardown)
        function cleanup(tc)
            if isfolder(tc.Tmp); rmdir(tc.Tmp, 's'); end
        end
    end

    methods (Test)

        function test_uniform_profile_reproduces_rings(tc)
            tr = tc.Tran; a = tr.transducer.annular;
            trp = ring_profile_transducer(tr, struct('n_sub', 4));
            tc.verifyEqual(trp.transducer.annular.elem_n, 4 * a.elem_n);
            tc.verifyEqual(trp.n_elem, tr.n_elem, 'channel count kept');
            tc.verifyEqual(trp.transducer.annular.elem_od_mm(4:4:end), a.elem_od_mm, 'AbsTol', 1e-9);
            z = (5:0.5:120)' + a.curv_radius_mm - a.dist_geom_ep_mm;
            H0 = axial_transfer_matrix(tr, z, 1500, 994, 'rayleigh');
            Hs = axial_transfer_matrix(trp, z, 1500, 994, 'rayleigh');
            Hr = Hs * kron(eye(a.elem_n), ones(4, 1));          % sum the sub-rings of each ring
            tc.verifyLessThan(norm(Hr - H0) / norm(H0), 1e-3);
            tc.verifyEqual(ring_profile_transducer(trp, struct('n_sub', 4)), trp, 'idempotent');
        end

        function test_sub_rings_steer_with_their_ring(tc)
            trp = ring_profile_transducer(tc.Tran, struct('n_sub', 3));
            d = steering_path_lengths(trp, 50);
            tc.verifyEqual(d(1:3), repmat(d(1), 1, 3), 'AbsTol', 1e-12);
            ph = compute_phases(1500, trp, 50);
            tc.verifyEqual(ph(4:6), repmat(ph(4), 1, 3), 'AbsTol', 1e-9);
        end

        function test_fit_recovers_a_profile_model(tc)
            tr = tc.Tran; a = tr.transducer.annular; M = 4;
            rp = struct('n_sub', M, 'weight_abs', repmat([1.3 1.1 0.8 0.6], 1, a.elem_n), ...
                'weight_phase_deg', repmat([0 10 25 45], 1, a.elem_n));
            trp = ring_profile_transducer(tr, rp);
            z = (0:0.5:110)'; d = [40 50 60 70];
            H = axial_transfer_matrix(trp, z + a.curv_radius_mm - a.dist_geom_ep_mm, 1500, 994, 'rayleigh');
            rng(3); I = zeros(numel(z), numel(d));
            for k = 1:numel(d)
                c = exp(0.2 * randn(a.elem_n, 1)) .* exp(1i * 2 * pi * rand(a.elem_n, 1));
                I(:, k) = abs(H * c(kron(1:a.elem_n, ones(1, M)))).^2;
            end
            [fit, info] = fit_ring_profile(tr, struct('depths_ep_mm', d, 'z_ep_mm', z, 'I', I), ...
                struct('NSub', M, 'Smooth', 0.3, 'NStarts', 6, 'Iterations', 5));
            tc.verifyGreaterThan(median(info.r2_profile), 0.995);
            tc.verifyGreaterThan(median(info.r2_profile), median(info.r2_uniform));
            tc.verifyEqual(numel(fit.weight_abs), a.elem_n * M);
            area = trp.transducer.annular.elem_od_mm.^2 - trp.transducer.annular.elem_id_mm.^2;
            u = fit.weight_abs .* exp(1i * fit.weight_phase_deg * pi / 180);
            for i = 1:a.elem_n
                e = (i - 1) * M + (1:M);
                tc.verifyEqual(sum(area(e) .* u(e)) / sum(area(e)), 1, 'AbsTol', 1e-3, 'ring mean weight = 1');
            end
        end

        function test_field_matrix_matches_axis_and_converges(tc)
            tr = tc.Tran; a = tr.transducer.annular; off = a.curv_radius_mm - a.dist_geom_ep_mm;
            z = (5:2:90)';
            H0 = axial_transfer_matrix(tr, z + off, 1500, 994, 'rayleigh');
            H1 = field_transfer_matrix(tr, zeros(size(z)), z + off, 1500, 994);
            tc.verifyLessThan(norm(H1 - H0) / norm(H0), 1e-10, 'on axis = AXIAL_TRANSFER_MATRIX');
            r = [2 8 20 30]'; zz = [40 40 15 5]' + off;
            Ha = field_transfer_matrix(tr, r, zz, 1500, 994);
            Hb = field_transfer_matrix(tr, r, zz, 1500, 994, 4096);
            tc.verifyLessThan(norm(Ha - Hb) / norm(Hb), 1e-3, 'default azimuth sampling has converged');
        end

        function test_offaxis_data_constrain_the_profile(tc)
            % A transverse plane added to the axial profiles brings the weights closer to the truth
            tr = tc.Tran; a = tr.transducer.annular; M = 3; off = a.curv_radius_mm - a.dist_geom_ep_mm;
            rp = struct('n_sub', M, 'weight_abs', repmat([1.25 1 0.75], 1, a.elem_n), ...
                'weight_phase_deg', repmat([0 15 35], 1, a.elem_n));
            trp = ring_profile_transducer(tr, rp);
            z = (0:0.5:100)'; d = [40 55 70]; rr = (0:0.5:20)';
            H  = axial_transfer_matrix(trp, z + off, 1500, 994, 'rayleigh');
            Hr = field_transfer_matrix(trp, rr, 30 * ones(size(rr)) + off, 1500, 994);
            rng(4); I = zeros(numel(z), numel(d)); C = cell(1, numel(d));
            for k = 1:numel(d)
                C{k} = exp(0.2 * randn(a.elem_n, 1)) .* exp(1i * 2 * pi * rand(a.elem_n, 1));
                I(:, k) = abs(H * C{k}(kron(1:a.elem_n, ones(1, M)))).^2;
            end
            Ir = abs(Hr * C{1}(kron(1:a.elem_n, ones(1, M)))).^2;
            prof = struct('depths_ep_mm', d, 'z_ep_mm', z, 'I', I);
            o = struct('NSub', M, 'Smooth', 0.1, 'NStarts', 6, 'Iterations', 5);
            f_ax = fit_ring_profile(tr, prof, o);
            o.OffAxis = struct('depth_ep_mm', d(1), 'r_mm', rr, 'z_ep_mm', 30 * ones(size(rr)), 'I', Ir);
            [f_oa, info] = fit_ring_profile(tr, prof, o);
            area = trp.transducer.annular.elem_od_mm.^2 - trp.transducer.annular.elem_id_mm.^2;
            ut = rp.weight_abs .* exp(1i * rp.weight_phase_deg * pi / 180);
            for i = 1:a.elem_n; e = (i - 1) * M + (1:M); ut(e) = ut(e) / (sum(area(e) .* ut(e)) / sum(area(e))); end
            err = @(f) median(abs(f.weight_abs .* exp(1i * f.weight_phase_deg * pi / 180) - ut));
            tc.verifyLessThan(err(f_oa), err(f_ax), 'off-axis data reduce the weight error');
            tc.verifyGreaterThan(info.offaxis(1).r2_profile, 0.99);
        end

        function test_plane_scans_read_and_rescale(tc)
            % READ_OFFAXIS_SCANS: CSV and FUS Metrology export, radial reduction, rescaling
            Ir = @(r) 20 * exp(-(r / 3).^2);                     % radially symmetric plane [W/cm²]
            [x, y] = meshgrid(-10:0.5:10); x = x(:); y = y(:);
            [xs, zs] = meshgrid(-5:0.5:10, 20:0.5:50); xs = xs(:); zs = zs(:);
            T = table([40 * ones(numel(x), 1); 40 * ones(numel(xs), 1)], [x; xs], [y; 0 * xs], ...
                [35 * ones(numel(x), 1); zs], [Ir(hypot(x, y)); Ir(abs(xs)) .* exp(-((zs - 40) / 8).^2)], ...
                'VariableNames', {'focus_ep_mm', 'x_mm', 'y_mm', 'z_ep_mm', 'intensity_wcm2'});
            T.label = [repmat("XY 35", numel(x), 1); repmat("XZ", numel(xs), 1)];
            csv = fullfile(tc.Tmp, 'planes.csv'); writetable(T, csv);
            oa = read_offaxis_scans(csv, struct('RMaxMm', 8));
            tc.verifyEqual(sort({oa.kind}), {'axial', 'transverse'});
            t = oa(strcmp({oa.kind}, 'transverse'));
            tc.verifyEqual(t.z_ep_mm(1), 35); tc.verifyEqual(t.depth_ep_mm, 40);
            tc.verifyEqual(t.I(1), Ir(0), 'RelTol', 1e-9, 'axis bin = the centre sample');
            tc.verifyEqual(t.I, Ir(t.r_mm), 'AbsTol', 0.5, 'azimuthal average of a symmetric plane (0.5 mm bins)');
            a = oa(strcmp({oa.kind}, 'axial'));
            tc.verifyLessThanOrEqual(max(a.r_mm), 8); tc.verifyEqual(min(a.z_ep_mm), 20);
            tc.verifyEqual(read_offaxis_scans(csv, struct('Kind', 'transverse')).label, 'XY 35');
            % 3D scan, axis off the origin: (r, z) map around the centroid
            [xv, yv, zv] = ndgrid(-8:0.5:8, -8:0.5:8, 30:1:50); xv = xv(:) + 1; yv = yv(:) - 0.5; zv = zv(:);
            Iv = Ir(hypot(xv - 1, yv + 0.5)) .* exp(-((zv - 40) / 8).^2);
            V = table(40 * ones(size(xv)), xv, yv, zv, Iv, ...
                'VariableNames', {'focus_ep_mm', 'x_mm', 'y_mm', 'z_ep_mm', 'intensity_wcm2'});
            vcsv = fullfile(tc.Tmp, 'volume.csv'); writetable(V, vcsv);
            ov = read_offaxis_scans(vcsv, struct('Centre', 'centroid', 'RMaxMm', 6));
            tc.verifyEqual(ov.kind, 'volume');
            tc.verifyEqual(numel(unique(ov.z_ep_mm)), 21);
            at40 = ov.z_ep_mm == 40;
            tc.verifyEqual(ov.I(at40), Ir(ov.r_mm(at40)), 'AbsTol', 0.6, 'radial profile at each z');
            tc.verifyError(@() read_offaxis_scans(vcsv, struct('Kind', 'transverse')), 'read_offaxis_scans:none');
            % FUS Metrology export with its own axial profiles: planes rescaled to the fitted axial peak
            z = (0:0.5:80)'; ax = 30 * exp(-((z - 40) / 6).^2);
            pl = struct('type', 'Transverse XY - Holography', 'focus', 40, 'x', x, 'y', y, 'z', 35 + 0 * x, ...
                        'ISPPA', struct('data', single(Ir(hypot(x, y))'), 'unit', 'W/cm2'));
            axs = struct('type', 'Axial', 'focus', 40, 'x', 0 * z, 'y', 0 * z, 'z', z, ...
                         'ISPPA', struct('data', single(ax'), 'unit', 'W/cm2'));
            data = {axs, pl}; mat = fullfile(tc.Tmp, 'FUS_Metrology_Export_BableBrain.mat'); save(mat, 'data');
            meas = struct('depths_ep_mm', 39.6, 'z_ep_mm', z, 'I', 0.5 * ax);   % half the export's level
            o1 = read_offaxis_scans(mat);
            o2 = read_offaxis_scans(tc.Tmp, struct('Axial', meas, 'Kind', 'transverse'));
            tc.verifyEqual(o1.scale, 1);
            tc.verifyEqual(o2(1).scale, 0.5, 'RelTol', 1e-6);
            tc.verifyEqual(o2(1).I, 0.5 * o1.I, 'RelTol', 1e-6);
        end

        function test_plane_file_path_equals_struct(tc)
            % fit_ring_profile reads a plane-scan path itself; same result as the struct
            tr = tc.Tran; a = tr.transducer.annular; off = a.curv_radius_mm - a.dist_geom_ep_mm;
            z = (0:0.5:100)'; d = [40 55];
            H = axial_transfer_matrix(tr, z + off, 1500, 994, 'rayleigh');
            I = [abs(H * ones(a.elem_n, 1)).^2, abs(H * exp(1i * (0:a.elem_n - 1)' * 0.4)).^2];
            prof = struct('depths_ep_mm', d, 'z_ep_mm', z, 'I', I);
            [x, y] = meshgrid(-6:0.5:6); x = x(:); y = y(:);
            Hr = field_transfer_matrix(tr, hypot(x, y), 35 * ones(size(x)) + off, 1500, 994);
            T = table(40 * ones(size(x)), x, y, 35 * ones(size(x)), abs(Hr * ones(a.elem_n, 1)).^2, ...
                'VariableNames', {'focus_ep_mm', 'x_mm', 'y_mm', 'z_ep_mm', 'intensity_wcm2'});
            csv = fullfile(tc.Tmp, 'plane_fit.csv'); writetable(T, csv);
            o = struct('NSub', 2, 'NStarts', 2, 'Iterations', 2, 'OffAxisOptions', struct('RMaxMm', 6));
            [rp1, i1] = fit_ring_profile(tr, prof, setfield(o, 'OffAxis', csv)); %#ok<SFLD>
            oa = read_offaxis_scans(csv, struct('RMaxMm', 6, 'Axial', prof, 'NearFieldSkipMm', 4));
            [rp2, i2] = fit_ring_profile(tr, prof, setfield(rmfield(o, 'OffAxisOptions'), 'OffAxis', oa)); %#ok<SFLD>
            tc.verifyEqual(rp1.weight_abs, rp2.weight_abs, 'AbsTol', 1e-9);
            tc.verifyEqual(i1.offaxis.r2_profile, i2.offaxis.r2_profile, 'AbsTol', 1e-9);
            tc.verifyEqual(rp1.offaxis_source, csv);
            tc.verifyEqual(rp2.offaxis_source, '');
            tc.verifyGreaterThan(rp1.offaxis_r2, 0.9);
        end

        function test_drive_applies_weights_once(tc)
            trp = ring_profile_transducer(tc.Tran, struct('n_sub', 2, 'weight_abs', 0.5 * ones(1, 8), ...
                'weight_phase_deg', 30 * ones(1, 8)));
            ann = ring_profile_drive(trp.transducer.annular, zeros(1, 8), 1000 * ones(1, 8));
            tc.verifyEqual(ann.elem_phase_deg, 30 * ones(1, 8), 'AbsTol', 1e-9);
            tc.verifyEqual(ann.elem_p_amp, 500 * ones(1, 8));
            tc.verifyFalse(isfield(ann, 'elem_weight'));
            ann0 = ring_profile_drive(tc.Tran.transducer.annular, [1 2 3 4], [5 6 7 8]);
            tc.verifyEqual(ann0.elem_phase_deg, [1 2 3 4]); tc.verifyEqual(ann0.elem_p_amp, [5 6 7 8]);
        end

        function test_library_round_trip(tc)
            % calibrate_library_entry with a given profile, then load as a study would
            tr = tc.Tran; a = tr.transducer.annular; M = 2;
            rp = struct('n_sub', M, 'weight_abs', repmat([1.1 0.9], 1, a.elem_n), ...
                'weight_phase_deg', repmat([0 20], 1, a.elem_n), 'smooth', 0.3, 'near_field_skip_mm', 4);
            trp = ring_profile_transducer(tr, rp);
            z = (0:0.5:110)'; d = [45 60];
            H = axial_transfer_matrix(trp, z + a.curv_radius_mm - a.dist_geom_ep_mm, 1500, 994, 'rayleigh');
            I = zeros(numel(z), numel(d));
            for k = 1:numel(d)
                ph = compute_phases(1500, trp, d(k)) * pi / 180;
                I(:, k) = 30 * abs(H * exp(1i * ph(:))).^2 / max(abs(H * exp(1i * ph(:))).^2);
            end
            meas = struct('depths_ep_mm', d, 'z_ep_mm', z, 'I', I);
            lib = struct('meta', struct('tran_serial', 'CTX_500_026', 'driving_system', 'SC_105_010'), ...
                'global_model', struct('elem_phase_correction_deg', zeros(1, a.elem_n)));
            lib = calibrate_library_entry(lib, tr, meas, struct('RingProfile', rp, ...
                'FitOptions', struct('NStarts', 2), 'Intensity', 30));
            tc.verifyEqual(lib.global_model.ring_profile, rp);
            yaml.dumpFile(fullfile(tc.Tmp, 'CTX_500_026_SC_105_010.yaml'), lib);
            p = load_transducer_from_library('CTX_500_026_SC_105_010', 50, 30, tc.Equip, tc.Tmp);
            ann = p.transducer.annular;
            tc.verifyEqual(ann.elem_n, a.elem_n * M);
            tc.verifyEqual(numel(ann.elem_phase_deg), a.elem_n * M);
            tc.verifyEqual(numel(ann.elem_p_amp), a.elem_n * M);
            tc.verifyFalse(isfield(ann, 'elem_weight'), 'weights folded into the drive');
            tc.verifyEqual(p.transducer.tx_calibration.ring_profile.n_sub, M);
            % drive ratio within a ring follows the profile
            tc.verifyEqual(ann.elem_p_amp(2) / ann.elem_p_amp(1), 0.9 / 1.1, 'RelTol', 1e-2);
            r = validate_transducer_library(tc.Tmp, fullfile(fileparts(fileparts(mfilename('fullpath'))), 'config', 'equipment'));
            tc.verifyFalse(any(contains(string([r.issues]), 'elements')), 'validator accepts sub-ring arrays');
        end

    end
end
