classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_calibration_fit_defaults < matlab.unittest.TestCase
% TEST_CALIBRATION_FIT_DEFAULTS  CALIBRATION_FIT fills fit options that are not
% given from the 'rayleigh_phase_gain' preset (CALIBRATION_FIT_OPTIONS), not
% from FIT_TRANSDUCER_MODEL's own defaults; options that are given win.
% Synthetic geometric profile of a shipped annular transducer; no data.
%
%   results = runtests('tests/test_calibration_fit_defaults.m');

    properties
        Tran
        Meas
    end

    methods (TestClassSetup)
        function setup(tc)
            [~, tc.Tran] = steerable_annular_serial();
            d0 = tc.Tran.transducer.annular.dist_geom_ep_mm;
            z  = (0:0.25:1.4 * d0)';
            m  = fit_transducer_model(tc.Tran, struct('depths_ep_mm', d0, 'z_ep_mm', z, ...
                'I', ones(numel(z), 1)), struct('EvaluateOnly', true, 'NearFieldSkipMm', 0));
            tc.Meas = struct('depths_ep_mm', d0, 'z_ep_mm', z, ...
                'I', 30 * interp1(m.z_ep_mm, m.I_model, z, 'linear', 0));
        end
    end

    methods (Test)
        function test_missing_options_come_from_the_preset(tc)
            fit = calibration_fit(tc.Tran, tc.Meas, struct('FitOptions', struct('NStarts', 2)));
            preset = calibration_fit_options('rayleigh_phase_gain');
            rec = fit.record.fit_options;
            for f = {'ForwardModel', 'Domain', 'HuberDelta', 'Lambda', 'NearFieldSkipMm'}
                tc.verifyTrue(isfield(rec, f{1}), sprintf('%s not set from the preset', f{1}));
                if isfield(rec, f{1})
                    tc.verifyEqual(rec.(f{1}), preset.(f{1}), f{1});
                end
            end
            tc.verifyEqual(rec.NStarts, 2);
        end

        function test_given_options_win(tc)
            % Negative control: explicit options are kept
            fo  = struct('NStarts', 2, 'Domain', 'amplitude', 'HuberDelta', 0.05);
            fit = calibration_fit(tc.Tran, tc.Meas, struct('FitOptions', fo));
            tc.verifyEqual(fit.record.fit_options.Domain, 'amplitude');
            tc.verifyEqual(fit.record.fit_options.HuberDelta, 0.05);
        end
    end
end
