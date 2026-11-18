classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_medium < matlab.unittest.TestCase
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
            tc.verifyEqual(sort(fieldnames(plus)), ...
                sort({'temp_0'; 'absorption_np'; 'absorption_fraction'; 'absorption_from_coeff'}));
            tc.verifyFalse(plus.absorption_from_coeff);
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
            p = use_coefficients(medium_params());
            [masks, tissues] = slab_masks(p);
            [~, plus] = quiet_medium_setup(p, masks);
            tc.verifyTrue(plus.absorption_from_coeff);
            f_mhz = p.transducer.freq_hz / 1e6;
            for t = 1:numel(tissues)
                ref = p.medium_properties.(tissues{t});
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
            tc.verifyFalse(plus_f.absorption_from_coeff);
            p = use_coefficients(p);
            [~, plus_c] = quiet_medium_setup(p, masks);
            w = 2*pi*p.transducer.freq_hz;
            alpha_np_f = db2neper(km_f.alpha_coeff .* plus_f.absorption_fraction, km_f.alpha_power) .* w.^km_f.alpha_power;
            for t = 1:numel(tissues)
                ref = tc.Cfg.medium_properties.(tissues{t});
                if ref.absorption_fraction == 1
                    v = masks == t;
                    tc.verifyEqual(plus_c.absorption_np(v), alpha_np_f(v), 'RelTol', 1e-9, tissues{t});
                end
            end
        end

        function test_mixed_fields_keep_coefficient_and_fraction(tc)
            % brain as fraction, other tissues as coefficient: both maps are
            % effective for every tissue (coefficient = fraction x attenuation)
            p = use_coefficients(medium_params());
            p.medium_properties.brain.absorption_fraction = 0.77;
            p.medium_properties.brain.absorption_coeff = [];
            [masks, tissues] = slab_masks(p);
            [~, plus] = quiet_medium_setup(p, masks);
            tc.verifyTrue(plus.absorption_from_coeff);
            db2np = @(a, y) a * (p.transducer.freq_hz / 1e6)^y * 100 / (20*log10(exp(1)));
            for t = 1:numel(tissues)
                ref = p.medium_properties.(tissues{t});
                v = masks == t;
                if strcmp(tissues{t}, 'brain')
                    np = db2np(0.77 * ref.alpha_coeff, ref.alpha_power);
                    frac = 0.77;
                else
                    np = db2np(ref.absorption_coeff, ref.alpha_power);
                    frac = ref.absorption_coeff / ref.alpha_coeff;
                end
                tc.verifyEqual(unique(plus.absorption_np(v)), np, 'RelTol', 1e-12, tissues{t});
                tc.verifyEqual(unique(plus.absorption_fraction(v)), frac, 'RelTol', 1e-12, tissues{t});
            end
        end

        function test_absorption_properties_report_both(tc)
            p = use_coefficients(medium_params(), {'skull'});
            b = tc.Cfg.medium_properties.brain;
            s = p.medium_properties.skull;
            a = thermal_absorption_properties(p, 'brain');
            tc.verifyEqual(a.model, 'fraction');
            tc.verifyEqual(a.absorption_coeff, b.absorption_fraction * b.alpha_coeff, 'RelTol', 1e-12);
            a = thermal_absorption_properties(p, 'skull');
            tc.verifyEqual(a.model, 'coefficient');
            tc.verifyEqual(a.absorption_fraction, s.absorption_coeff / s.alpha_coeff, 'RelTol', 1e-12);
            tc.verifyFalse(a.attenuation_mapped);
        end

        function test_absorption_model_follows_set_field(tc)
            p = medium_params();
            tc.verifyEqual(thermal_absorption_model(p, 'brain'), 'fraction');   % shipped default
            p.medium_properties.brain.absorption_fraction = NaN;                % NaN is unset: fraction 1
            tc.verifyEqual(thermal_absorption_property(p, 'brain', 'absorption_fraction'), 1);
            p = use_coefficients(p);
            [m, forced] = thermal_absorption_model(p, 'brain');
            tc.verifyEqual(m, 'coefficient'); tc.verifyFalse(forced);
            p.pct.enabled = 1; p.pct.mapping_attenuation = 'none';             % fixed skull values
            tc.verifyEqual(thermal_absorption_model(p, 'skull'), 'coefficient');
            p.pct.mapping_attenuation = 'k-plan';                               % mapped skull attenuation
            tc.verifyEqual(thermal_absorption_model(p, 'brain'), 'coefficient'); % only the skull is mapped
            [m, forced] = thermal_absorption_model(p, 'skull');
            tc.verifyEqual(m, 'fraction'); tc.verifyTrue(forced);
            a = thermal_absorption_properties(p, 'skull');                      % coefficient applied as fraction
            tc.verifyEqual(a.absorption_fraction, 3.78 / p.medium_properties.skull.alpha_coeff, 'RelTol', 1e-12);
            p.medium_properties.brain.absorption_fraction = 0.5;                % both set after loading
            tc.verifyError(@() thermal_absorption_model(p, 'brain'), 'thermal_absorption_model:bothSet');
        end

        function test_config_layer_clears_other_absorption_field(tc)
            dir = tc.createTemporaryFolder();
            write_yaml(fullfile(dir, 'coeff.yaml'), 'medium_properties:\n  skull:\n    absorption_coeff: 3.78\n');
            write_yaml(fullfile(dir, 'frac.yaml'), 'medium_properties:\n  skull:\n    absorption_fraction: 0.3\n');
            write_yaml(fullfile(dir, 'both.yaml'), ...
                'medium_properties:\n  skull:\n    absorption_coeff: 3.78\n    absorption_fraction: 0.3\n');
            p = load_parameters('coeff.yaml', dir);                             % clears the default fraction
            tc.verifyEqual(p.medium_properties.skull.absorption_coeff, 3.78);
            tc.verifyEmpty(p.medium_properties.skull.absorption_fraction);
            tc.verifyEqual(p.medium_properties.brain.absorption_fraction, tc.Cfg.medium_properties.brain.absorption_fraction);
            tc.verifyEqual(thermal_absorption_model(p, 'skull'), 'coefficient');
            p = load_parameters('coeff.yaml', dir, 'frac.yaml', dir);           % a later layer switches back
            tc.verifyEqual(p.medium_properties.skull.absorption_fraction, 0.3);
            tc.verifyEmpty(p.medium_properties.skull.absorption_coeff);
            p = load_parameters('coeff.yaml', dir, {'medium_properties.skull.absorption_fraction', 0.3}); % overrides
            tc.verifyEqual(thermal_absorption_model(p, 'skull'), 'fraction');
            tc.verifyEmpty(p.medium_properties.skull.absorption_coeff);
            tc.verifyError(@() load_parameters('both.yaml', dir), 'medium_absorption_layer:bothSet');
        end

        function test_absorption_coeff_must_be_a_number(tc)
            p = use_coefficients(medium_params());
            p.medium_properties.brain.absorption_coeff = -1;
            masks = slab_masks(p);
            tc.verifyError(@() quiet_medium_setup(p, masks), 'medium_setup:badAbsorptionCoeff');
        end

        function test_coefficient_without_attenuation_stays_finite(tc)
            % alpha_coeff 0: the fraction map is 0 without absorption and 1
            % with it (which warns), never Inf or NaN
            p = use_coefficients(medium_params(), {'brain', 'skin'});
            p.medium_properties.brain.alpha_coeff = 0;
            p.medium_properties.brain.absorption_coeff = 0;
            p.medium_properties.skin.alpha_coeff = 0;
            [masks, tissues] = slab_masks(p);
            [km, plus] = quiet_medium_setup(p, masks); %#ok<ASGLU>
            tc.verifyTrue(all(isfinite(plus.absorption_fraction(:))));
            tc.verifyEqual(unique(plus.absorption_fraction(masks == find(strcmp(tissues, 'brain')))), 0);
            tc.verifyEqual(unique(plus.absorption_fraction(masks == find(strcmp(tissues, 'skin')))), 1);
            a = thermal_absorption_properties(p, 'brain');
            tc.verifyTrue(isnan(a.absorption_fraction), 'undefined fraction is reported as NaN');
        end

        function test_unit_conversion_db_to_np(tc)
            % 1 dB/(cm MHz) at 1 MHz = 100 / 8.686 Np/m; y scales with frequency
            tc.verifyEqual(db_cm_mhz_to_np_m(1, 1e6, 1), 100 / (20*log10(exp(1))), 'RelTol', 1e-12);
            tc.verifyEqual(db_cm_mhz_to_np_m(0.59, 500e3, 1.2), ...
                0.59 * 0.5^1.2 * 100 / (20*log10(exp(1))), 'RelTol', 1e-12);
            tc.verifyEqual(db_cm_mhz_to_np_m([1 2], 2e6, [1 2]), [2 8] * 100 / (20*log10(exp(1))), 'RelTol', 1e-12);
        end

        function test_cache_signature_keeps_the_set_absorption_field(tc)
            % caches written before absorption_coeff existed stay valid for fractions
            p = medium_params();
            sig = cache_signature(p, 'medium');
            tc.verifyFalse(isfield(sig.medium_properties.brain, 'absorption_coeff'));
            tc.verifyTrue(isfield(sig.medium_properties.brain, 'absorption_fraction'));
            p = use_coefficients(p, {'skull'});
            sig = cache_signature(p, 'medium');
            tc.verifyFalse(isfield(sig.medium_properties.brain, 'absorption_coeff'));
            tc.verifyTrue(isfield(sig.medium_properties.skull, 'absorption_coeff'));
            tc.verifyFalse(isfield(sig.medium_properties.skull, 'absorption_fraction'));
            sig = cache_signature(p, 'acoustic');
            tc.verifyFalse(isfield(sig.medium_properties.skull, 'absorption_coeff'));
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
            ref = tc.Cfg.medium_properties.brain;
            a = get_alpha_coeff('brain', 500e3);
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

function p = use_coefficients(p, tissues)
% absorption_coeff instead of absorption_fraction: 3.78 for the skull layers,
% alpha_coeff (all attenuation absorbed) for the others
    if nargin < 2
        tissues = fieldnames(p.medium_properties);
    end
    for t = 1:numel(tissues)
        m = p.medium_properties.(tissues{t});
        if startsWith(tissues{t}, 'skull')
            m.absorption_coeff = 3.78;
        else
            m.absorption_coeff = m.alpha_coeff;
        end
        m.absorption_fraction = [];
        p.medium_properties.(tissues{t}) = m;
    end
end

function v = thermal_absorption_property(p, tissue, field)
    a = thermal_absorption_properties(p, tissue);
    v = a.(field);
end

function write_yaml(file, text)
    fid = fopen(file, 'w');
    fprintf(fid, text);
    fclose(fid);
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
