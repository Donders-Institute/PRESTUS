classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_calibration_session < matlab.unittest.TestCase
% TEST_CALIBRATION_SESSION  The interactive calibration steps shared by CALIBRATION_GUI and the
% PRESTUS Nexus: CALIBRATION_DESIRED_PROFILES, CALIBRATION_FIT, CALIBRATION_EVALUATE,
% CALIBRATION_SAVE_ENTRY (meta.fit). Synthetic "hydrophone" profiles are the geometric model
% of a shipped annular transducer at three depths, so no measurement data is needed.

    properties
        Tran
        Serial
        Meas
        Dir
    end

    methods (TestClassSetup)
        function setup(tc)
            [tc.Serial, tc.Tran] = steerable_annular_serial();
            d0 = tc.Tran.transducer.annular.dist_geom_ep_mm;
            z = (0:0.25:1.4 * d0)';
            dm = round(d0 + [-10 0 10], 2);
            I = zeros(numel(z), 3);
            for j = 1:3
                m = fit_transducer_model(tc.Tran, struct('depths_ep_mm', dm(j), 'z_ep_mm', z, ...
                    'I', ones(numel(z), 1)), struct('EvaluateOnly', true, 'NearFieldSkipMm', 0));
                I(:, j) = 30 * interp1(m.z_ep_mm, m.I_model, z, 'linear', 0);
            end
            tc.Meas = struct('depths_ep_mm', dm, 'z_ep_mm', z, 'I', I);
        end
    end

    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tempname; mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
        end
    end

    methods (Test, TestTags = {'unit', 'calibration'})

        function test_desired_profiles_without_and_with_measurements(tc)
            d0 = tc.Meas.depths_ep_mm(2);
            p = calibration_desired_profiles(tc.Tran, [], [d0 + 5, d0 - 5, d0 - 5.001]);
            tc.verifyEqual(p.depths_ep_mm, round([d0 - 5, d0 + 5], 2));
            tc.verifyEqual(p.profile_type, {'unavailable', 'unavailable'});
            tc.verifyTrue(all(isnan(p.I), 'all'));
            p = calibration_desired_profiles(tc.Tran, tc.Meas, [d0, d0 + 5, d0 + 30], ...
                struct('MaxExtrapolationMm', 0));
            tc.verifyEqual(p.profile_type, {'measured', 'interpolated', 'unavailable'});
            tc.verifyEqual(p.n_z_scan, numel(tc.Meas.z_ep_mm));
            tc.verifyGreaterThan(p.z_ep_mm(end), tc.Meas.z_ep_mm(end));   % model range beyond the scan
            tc.verifyTrue(all(isnan(p.I(p.n_z_scan + 1:end, :)), 'all'));
        end

        function test_fit_skips_unavailable_and_fits_only_the_scan(tc)
            d0 = tc.Meas.depths_ep_mm(2);
            sel = calibration_desired_profiles(tc.Tran, tc.Meas, [d0 - 5, d0 + 40]);
            [fit, msg] = calibration_fit(tc.Tran, sel, struct('FitOptions', struct('NStarts', 2), ...
                'TranSerial', tc.Serial, 'CalibrateAt', 'desired', 'Tool', 'test'));
            tc.verifyEqual(fit.prof.depths_ep_mm, round(d0 - 5, 2));
            tc.verifyEqual(numel(fit.prof.z_ep_mm), numel(tc.Meas.z_ep_mm));
            tc.verifyFalse(isfield(fit.prof, 'n_z_scan'));
            tc.verifySubstring(msg, 'No profile at');
            tc.verifyGreaterThan(fit.r2, 0.9);
            tc.verifyEqual(fit.record.calibrate_at, 'desired');
            tc.verifyEqual(fit.record.tool, 'test');
            tc.verifyFalse(isfield(fit.record.fit_options, 'Prior'));
            tc.verifyTrue(isfield(fit.entry, 'refinement'));
            tc.verifyEqual(fit.entry.meta.tran_serial, tc.Serial);
        end

        function test_fit_errors(tc)
            p = tc.Meas; p.profile_type = repmat({'unavailable'}, 1, 3);
            tc.verifyError(@() calibration_fit(tc.Tran, p), 'calibration_fit:noProfiles');
            tc.verifyError(@() calibration_fit(tc.Tran, tc.Meas, struct('Mode', 'static')), 'calibration_fit:mode');   % removed
            tc.verifyError(@() calibration_fit(tc.Tran, tc.Meas, struct('Mode', 'joint')), 'calibration_fit:mode');
        end

        function test_evaluate_measured_and_predicted(tc)
            fo = struct('NStarts', 2);
            [res, msg] = calibration_evaluate(tc.Tran, tc.Meas, struct('FitOptions', fo));
            tc.verifyEqual(fieldnames(res), {'geometric'});
            tc.verifyGreaterThan(res.geometric.r2, 0.99);              % the "measurements" are geometric
            tc.verifySubstring(msg{1}, 'Geometric median R²');
            fit = calibration_fit(tc.Tran, tc.Meas, struct('FitOptions', fo, 'TranSerial', tc.Serial));
            d0 = tc.Meas.depths_ep_mm(2);
            p = calibration_desired_profiles(tc.Tran, tc.Meas, [d0 - 5, d0, d0 + 30]);
            [res, msg] = calibration_evaluate(tc.Tran, p, struct('Predicted', true, 'Fit', fit, ...
                'FitOptions', fo, 'MeasuredZmaxMm', max(tc.Meas.z_ep_mm)));
            tc.verifyEqual(fieldnames(res), {'geometric'; 'fit'});
            tc.verifyEqual(numel(res.fit.focus), 3);
            tc.verifyTrue(isnan(res.fit.r2(3)));                       % no profile there
            tc.verifyEqual({res.fit.settings.phase_model}, {'refined', 'refined', 'geometric'});   % no static model without a prior
            tc.verifySubstring(msg{1}, 'Predicted 3 depths');
            tc.verifySubstring(msg{2}, 'Fit (unsaved) at the desired depths');
        end

        function test_save_entry_records_fit_and_protects_shipped_library(tc)
            d0 = tc.Meas.depths_ep_mm(2);
            sel = calibration_desired_profiles(tc.Tran, tc.Meas, d0 - 5);
            fit = calibration_fit(tc.Tran, sel, struct('FitOptions', struct('NStarts', 2), ...
                'TranSerial', tc.Serial, 'DrivingSystem', 'DS1', 'IsppaWcm2', 20, 'CalibrateAt', 'desired', ...
                'ProfileOptions', struct('Align', 'peak')));
            [file, pcsv] = calibration_save_entry(fit.entry, fit, tc.Dir, ...
                struct('Kwave', struct('grid', '3D')));
            tc.verifyEqual(file, fullfile(tc.Dir, [tc.Serial '_DS1.yaml']));
            tc.verifyTrue(isfile(pcsv));                               % interpolated profile next to it
            lib = yaml.loadFile(file, 'ConvertToArray', true);
            tc.verifyEqual(char(lib.meta.calibration_type), 'user_empirical');
            tc.verifyEqual(char(lib.meta.fit.calibrate_at), 'desired');
            tc.verifyEqual(double(lib.meta.fit.isppa_wcm2), 20);
            tc.verifyEqual(char(lib.meta.fit.profile_options.Align), 'peak');
            tc.verifyEqual(char(lib.meta.fit.kwave.grid), '3D');
            tc.verifyEqual(cellstr(lib.meta.fit.fit_options.FreeParams), {'phase', 'amp'});
            tc.verifyError(@() calibration_save_entry(fit.entry, fit, tc.Dir, ...
                struct('ShippedLibraryPath', tc.Dir)), 'calibration_save_entry:shipped');
            calibration_save_entry(fit.entry, fit, tc.Dir, struct('ShippedLibraryPath', tc.Dir, 'AllowShipped', true));
        end

        function test_writers_stamp_the_same_provenance(tc)
            % calibration_save_entry and save_library_entry store the version
            % string in prestus_version and the commit hash in prestus_commit.
            [hash, ver] = prestus_version();
            lib0 = struct('meta', struct('tran_serial', 'PROV_TEST'));
            fit = struct('mode', 'per_depth', 'record', struct('mode', 'per_depth'), ...
                'prof', struct('profile_type', {{'measured'}}));
            file = calibration_save_entry(lib0, fit, tc.Dir);
            a = yaml.loadFile(file, 'ConvertToArray', true);
            report = struct('depth_ep_mm', 50, 'r2_library', 0.9, 'r2_refined', 0.95);
            file_b = fullfile(tc.Dir, 'b', 'PROV_TEST.yaml');
            save_library_entry(lib0, file_b, 'per_depth', report);
            b = yaml.loadFile(file_b, 'ConvertToArray', true);
            for m = {a.meta, b.meta}
                tc.verifyEqual(string(m{1}.prestus_version), string(ver));
                tc.verifyEqual(string(m{1}.prestus_commit), string(hash));
            end
        end

        function test_focus_beyond_ignores_the_near_field(tc)
            z = (0:0.5:60)';
            I = exp(-((z - 40) / 4).^2) + 2 * exp(-((z - 3) / 1).^2);  % near-field peak at 3 mm
            tc.verifyEqual(calibration_focus_beyond(I, z, 10), 40, 'AbsTol', 0.5);
            tc.verifyLessThan(calibration_focus_beyond(I, z, 0), 10);
        end
    end
end
