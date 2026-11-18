classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_pct_mapping_curves < matlab.unittest.TestCase
% TEST_PCT_MAPPING_CURVES  HU -> property curves of the configured mappings (PCT_MAPPING_CURVES).
%
% The curves use PCT_MAP_PROPERTIES, the same calls MEDIUM_SETUP makes, so a
% curve point equals what a skull voxel of that HU is given in a run.
%
% Run with:   results = runtests('tests/test_pct_mapping_curves.m');

    methods (Test)
        function test_kplan_curves(tc)
            p = load_parameters();
            p.pct.mapping_density = 'k-plan';
            p.pct.mapping_soundspeed = 'k-plan';
            p.pct.mapping_attenuation = 'k-plan';
            c = pct_mapping_curves(p, [0 500 1000 1500]);
            tc.verifySize(c.density, [1 4]);
            tc.verifyTrue(all(diff(c.density) > 0), 'density rises with HU');
            tc.verifyGreaterThanOrEqual(min(c.density), p.medium_properties.water.density - 1e-6);
            tc.verifyEqual(c.mapping.density, 'k-plan');
            tc.verifyTrue(isfield(c.mapping, 'density_calibration'));
            tc.verifyEqual(c.reference.skull.density, p.medium_properties.skull.density);
        end

        function test_none_keeps_uniform_skull(tc)
            p = load_parameters();
            p.pct.mapping_density = 'none';
            p.pct.mapping_soundspeed = 'none';
            p.pct.mapping_attenuation = 'none';
            c = pct_mapping_curves(p, [0 1000]);
            tc.verifyEqual(c.sound_speed, repmat(p.medium_properties.skull.sound_speed, 1, 2), 'AbsTol', 1e-9);
        end

        function test_density_mappings_rise_with_hu(tc)
            % Every density mapping gives a density that does not fall with HU
            % within the bone range, and stays at or above water density.
            p = load_parameters();
            hu = [400 800 1200 1600];
            maps = {'k-plan', 'k-wave', 'marsac', 'marsac-babelbrain', 'aubry', 'air-tissue'};
            for i = 1:numel(maps)
                p.pct.mapping_density = maps{i};
                c = pct_mapping_curves(p, hu);
                tc.verifyTrue(all(diff(c.density) >= 0), [maps{i} ': density falls with HU']);
                tc.verifyGreaterThanOrEqual(min(c.density), p.medium_properties.water.density - 10, maps{i});
            end
        end

        function test_marsac_babelbrain_density(tc)
            % BabelBrain's HUtoDensityMarsac: 1000 + 1700 HU / max HU, HU raised to 300.
            p = load_parameters();
            p.pct.mapping_density = 'marsac-babelbrain';
            c = pct_mapping_curves(p, [100 300 1000 2000]);
            hu = [300 300 1000 2000];
            tc.verifyEqual(c.density, 1000 + 1700 * hu / 2000, 'AbsTol', 1e-9);
        end

        function test_matches_map_properties(tc)
            p = load_parameters();
            hu = [200 900 1800];
            c = pct_mapping_curves(p, hu);
            props = pct_map_properties(p, hu, 1:3);
            tc.verifyEqual(c.alpha_coeff, double(props.alpha_coeff), 'AbsTol', 1e-9);
        end

        function test_webb_coefficients(tc)
            % pct.webb_*: unset and [] give the built-in coefficients; set
            % values enter c = slope * HU + 1000 * intercept and
            % alpha = alpha0 * exp(c_att * HU) [Np/cm] -> dB/cm, alpha_power = beta.
            p = load_parameters();
            p.pct.mapping_soundspeed = 'webb';
            p.pct.mapping_attenuation = 'webb';
            hu = [400 900 1600];
            keys = {'webb_sos_slope', 'webb_sos_intercept', 'webb_att_alpha0', 'webb_att_beta', 'webb_att_c'};
            base = pct_mapping_curves(p, hu);
            tc.verifyEqual(base.sound_speed, 0.75 * hu + 1320, 'AbsTol', 1e-6);

            q = p;
            for k = keys
                q.pct.(k{1}) = [];
            end
            empty = pct_mapping_curves(q, hu);
            tc.verifyEqual(empty.sound_speed, base.sound_speed);
            tc.verifyEqual(empty.alpha_coeff, base.alpha_coeff);
            tc.verifyEqual(empty.alpha_power, base.alpha_power);

            vals = {1.0, 1.2, 20, 1.1, -0.001};
            for i = 1:numel(keys)
                q.pct.(keys{i}) = vals{i};
            end
            custom = pct_mapping_curves(q, hu);
            tc.verifyEqual(custom.sound_speed, 1.0 * hu + 1200, 'AbsTol', 1e-6);
            tc.verifyEqual(custom.alpha_coeff, 20 * exp(-0.001 * hu) * 8.686, 'RelTol', 1e-6);
            tc.verifyEqual(custom.alpha_power, repmat(1.1, 1, numel(hu)), 'AbsTol', 1e-9);
            tc.verifyNotEqual(custom.alpha_coeff, base.alpha_coeff);
        end
    end
end
