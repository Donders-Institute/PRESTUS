classdef test_medium < matlab.unittest.TestCase
% TEST_MEDIUM  Tissue label -> medium property assignment (functions/medium).
%
% MEDIUM_SETUP on a tiny synthetic label volume (one slab per tissue of
% config_default.yaml, plus unlabelled voxels): each tissue gets the acoustic
% and thermal values of config_default.yaml (read here directly with
% yaml.loadFile, independent of LOAD_PARAMETERS), unlabelled voxels become
% water, perfusion is converted to a coefficient, attenuation is remapped to
% y = 2 at the source frequency (FITPOWERLAWPARAMSMULTI), and the
% thermal-only fields are split into medium_plus. Also GET_ALPHA_COEFF.
% No pseudo-CT (pct.enabled = 0); pCT mappings are covered by
% test_pct_density_calibration.
%
%   results = runtests('tests/test_medium.m');

    properties
        Repo
        Cfg     % config_default.yaml as parsed by yaml.loadFile
    end

    methods (TestClassSetup)
        function add_paths(tc)
            tc.Repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(tc.Repo, 'functions')));
            addpath(genpath(fullfile(tc.Repo, 'external')));
            addpath(fullfile(tc.Repo, 'config'));
            tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                {'prestus:toolboxDefault', 'prestus:noSimNIBS', 'prestus:missingCharm'}));
            tc.Cfg = yaml.loadFile(fullfile(tc.Repo, 'config', 'config_default.yaml'), ...
                'ConvertToArray', true);
        end
    end

    methods (Test, TestTags = {'medium_setup'})

        function test_properties_match_config_per_tissue(tc)
            p = medium_params();
            [masks, tissues] = slab_masks(p);
            [km, plus] = quiet_medium_setup(p, masks);
            for t = 1:numel(tissues)
                ref = tc.Cfg.medium_properties.(tissues{t});
                v = masks == t;
                tc.verifyEqual(unique(km.sound_speed(v)), ref.sound_speed, tissues{t});
                tc.verifyEqual(unique(km.density(v)), ref.density, tissues{t});
                tc.verifyEqual(unique(km.thermal_conductivity(v)), ref.thermal_conductivity, tissues{t});
                tc.verifyEqual(unique(km.specific_heat(v)), ref.specific_heat_capacity, tissues{t});
                tc.verifyEqual(unique(km.perfusion_coeff(v)), ...
                    ref.perfusion / 60 * ref.density * 1e-6, 'RelTol', 1e-12, tissues{t});
                tc.verifyEqual(unique(plus.absorption_fraction(v)), ref.absorption_fraction, tissues{t});
                tc.verifyEqual(unique(plus.temp_0(v)), tc.Cfg.thermal.temp_0.(tissues{t}), tissues{t});
            end
        end

        function test_attenuation_remapped_to_y2_at_source_frequency(tc)
            % alpha_fit = a0 * f_MHz^(y - 2): same absorption at f with y = 2
            p = medium_params();
            [masks, tissues] = slab_masks(p);
            km = quiet_medium_setup(p, masks);
            tc.verifyEqual(km.alpha_power, 2);
            f_mhz = p.transducer.freq_hz / 1e6;
            for t = 1:numel(tissues)
                ref = tc.Cfg.medium_properties.(tissues{t});
                tc.verifyEqual(unique(km.alpha_coeff(masks == t)), ...
                    ref.alpha_coeff * f_mhz^(ref.alpha_power - 2), 'RelTol', 1e-9, tissues{t});
            end
        end

        function test_unlabelled_voxels_become_water(tc)
            p = medium_params();
            masks = zeros(p.grid.dims);
            [km, plus] = quiet_medium_setup(p, masks);
            w = tc.Cfg.medium_properties.water;
            tc.verifyEqual(unique(km.sound_speed), w.sound_speed);
            tc.verifyEqual(unique(km.density), w.density);
            tc.verifyEqual(unique(plus.temp_0), tc.Cfg.thermal.temp_0.water);
            tc.verifyFalse(any(isnan(km.sound_speed(:))));
        end

        function test_thermal_only_fields_split_off(tc)
            % k-Wave rejects temp_0 / absorption_fraction in the acoustic medium
            p = medium_params();
            [km, plus] = quiet_medium_setup(p, zeros(p.grid.dims));
            tc.verifyFalse(isfield(km, 'temp_0'));
            tc.verifyFalse(isfield(km, 'absorption_fraction'));
            tc.verifyEqual(sort(fieldnames(plus)), sort({'temp_0'; 'absorption_fraction'}));
            tc.verifySize(km.sound_speed, p.grid.dims);
        end

        function test_no_remap_requires_single_power_law(tc)
            p = medium_params();
            p.fit_alpha_power = 0;
            masks = slab_masks(p);
            tc.verifyError(@() quiet_medium_setup(p, masks), 'medium_setup:alphaPowerNotUniform');
        end

        function test_no_remap_keeps_alpha_for_uniform_power(tc)
            p = medium_params();
            p.fit_alpha_power = 0;
            km = quiet_medium_setup(p, zeros(p.grid.dims));         % water only, y = 2
            tc.verifyEqual(unique(km.alpha_coeff), tc.Cfg.medium_properties.water.alpha_coeff);
            tc.verifyEqual(km.alpha_power, tc.Cfg.medium_properties.water.alpha_power);
        end

    end

    methods (Test, TestTags = {'absorption'})

        function test_absorption_coeff_converted_to_np_per_m(tc)
            % absorption_coeff [dB/(cm MHz^y)] with y = alpha_power, in Np/m at f
            p = medium_params();
            p.thermal.absorption_model = 'coefficient';
            [masks, tissues] = slab_masks(p);
            [~, plus] = quiet_medium_setup(p, masks);
            tc.verifyEqual(sort(fieldnames(plus)), sort({'temp_0'; 'absorption_np'}));
            f_mhz = p.transducer.freq_hz / 1e6;
            for t = 1:numel(tissues)
                ref = tc.Cfg.medium_properties.(tissues{t});
                tc.verifyEqual(unique(plus.absorption_np(masks == t)), ...
                    ref.absorption_coeff * f_mhz^ref.alpha_power * 100 / (20*log10(exp(1))), ...
                    'RelTol', 1e-12, tissues{t});
            end
        end

        function test_absorption_coeff_matches_fraction_where_equivalent(tc)
            % fraction 1 and absorption_coeff = alpha_coeff give the same heat source
            p = medium_params();
            [masks, tissues] = slab_masks(p);
            [km_f, plus_f] = quiet_medium_setup(p, masks);
            p.thermal.absorption_model = 'coefficient';
            [~, plus_c] = quiet_medium_setup(p, masks);
            w = 2*pi*p.transducer.freq_hz;
            alpha_np_f = db2neper(km_f.alpha_coeff .* plus_f.absorption_fraction, km_f.alpha_power) .* w.^km_f.alpha_power;
            for t = 1:numel(tissues)
                ref = tc.Cfg.medium_properties.(tissues{t});
                if ref.absorption_fraction == 1 && ref.absorption_coeff == ref.alpha_coeff
                    v = masks == t;
                    tc.verifyEqual(plus_c.absorption_np(v), alpha_np_f(v), 'RelTol', 1e-9, tissues{t});
                end
            end
        end

        function test_absorption_coeff_required_in_coefficient_mode(tc)
            p = medium_params();
            p.thermal.absorption_model = 'coefficient';
            p.medium_properties.brain = rmfield(p.medium_properties.brain, 'absorption_coeff');
            masks = slab_masks(p);
            tc.verifyError(@() quiet_medium_setup(p, masks), 'medium_setup:missingAbsorptionCoeff');
        end

        function test_absorption_model_resolution(tc)
            p = medium_params();
            tc.verifyEqual(thermal_absorption_model(p), 'fraction');      % shipped default
            p.thermal.absorption_model = 'coefficient';
            [m, forced] = thermal_absorption_model(p);
            tc.verifyEqual(m, 'coefficient'); tc.verifyFalse(forced);
            p.pct.enabled = 1; p.pct.mapping_attenuation = 'none';       % fixed skull values
            tc.verifyEqual(thermal_absorption_model(p), 'coefficient');
            p.pct.mapping_attenuation = 'k-plan';                         % mapped skull attenuation
            [m, forced] = thermal_absorption_model(p);
            tc.verifyEqual(m, 'fraction'); tc.verifyTrue(forced);
            p.thermal.absorption_model = 'other';
            tc.verifyError(@() thermal_absorption_model(p), 'thermal_absorption_model:unknownModel');
        end

        function test_fraction_mode_cache_signature_ignores_absorption_coeff(tc)
            % caches written before absorption_coeff existed stay valid in fraction mode
            p = medium_params();
            sig = cache_signature(p, 'medium');
            tc.verifyFalse(isfield(sig, 'absorption_model'));
            tc.verifyFalse(isfield(sig.medium_properties.brain, 'absorption_coeff'));
            p.thermal.absorption_model = 'coefficient';
            sig = cache_signature(p, 'medium');
            tc.verifyEqual(sig.absorption_model, 'coefficient');
            tc.verifyTrue(isfield(sig.medium_properties.brain, 'absorption_coeff'));
            sig = cache_signature(p, 'acoustic');
            tc.verifyFalse(isfield(sig.medium_properties.brain, 'absorption_coeff'));
        end

    end

    methods (Test, TestTags = {'attenuation'})

        function test_power_law_fit_identity_for_y2(tc)
            a = fitPowerLawParamsMulti(0.5, 2, 1500, 500e3, 2);
            tc.verifyEqual(a, 0.5, 'RelTol', 1e-9);
        end

        function test_power_law_fit_rejects_y_ref_1(tc)
            tc.verifyError(@() fitPowerLawParamsMulti(0.5, 1, 1500, 500e3, 1), ?MException);
        end

        function test_get_alpha_coeff(tc)
            % KNOWN BUG (2026-10): get_alpha_coeff.m:28 reads
            % load_parameters().medium.(name); the tissue table is
            % parameters.medium_properties and there is no top-level 'medium'
            % key, so every call errors. No caller in functions/ uses it.
            ref = tc.Cfg.medium_properties.brain;
            try
                a = get_alpha_coeff('brain', 500e3);
            catch ME
                tc.verifyFail(sprintf(['Known bug: get_alpha_coeff reads load_parameters().medium ' ...
                    'instead of .medium_properties (%s)'], ME.message));
                return
            end
            tc.verifyEqual(a, ref.alpha_coeff * 0.5^ref.alpha_power, 'RelTol', 1e-12);
        end

    end

end

% =========================================================================

function p = medium_params()
    p = load_parameters(struct('simulation', struct('debug', 0, 'interactive', 0)));
    p.subject_id = 1;
    p.simulation.medium = 'layered';
    p.grid.dims = [6 4 5];
    p.grid.resolution_mm = 1;
    p.transducer.freq_hz = 500e3;
    p.pct.enabled = 0;
    p.headmodel.smooth_properties = false;
end

function [masks, tissues] = slab_masks(p)
% One x-slab per medium index (medium_masks hold indices into medium_properties).
    tissues = fieldnames(p.medium_properties);
    masks = zeros(p.grid.dims);
    for t = 1:min(numel(tissues), p.grid.dims(1))
        masks(t, :, :) = t;
    end
    assert(numel(tissues) <= p.grid.dims(1), 'grid too small for all tissues');
end

function [km, plus] = quiet_medium_setup(p, masks)
    [~, km, plus] = evalc_medium(p, masks);
end

function [txt, km, plus] = evalc_medium(p, masks)
    txt = evalc('[km, plus] = medium_setup(p, masks, struct());');
end
