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

        function test_matches_map_properties(tc)
            p = load_parameters();
            hu = [200 900 1800];
            c = pct_mapping_curves(p, hu);
            props = pct_map_properties(p, hu, 1:3);
            tc.verifyEqual(c.alpha_coeff, double(props.alpha_coeff), 'AbsTol', 1e-9);
        end
    end
end
