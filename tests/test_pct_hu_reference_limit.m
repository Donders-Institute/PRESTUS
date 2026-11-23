classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_pct_hu_reference_limit < matlab.unittest.TestCase
% TEST_PCT_HU_REFERENCE_LIMIT  Max-normalised pCT mappings stay within the
% water-to-bone range when skull HU lie above the reference (percentile < 100)
% or below 0, for every mapping that scales by the reference.
%
% Run with:   results = runtests('tests/test_pct_hu_reference_limit.m');

    properties (TestParameter)
        quantity = {'density', 'soundspeed', 'attenuation'}
    end

    methods (Test)
        function test_porosity_limited(tc)
            phi = medium_pct_porosity(params(1000), [-200 0 500 1000 1500]);
            tc.verifyEqual(phi, [1 1 0.5 0 0], 'AbsTol', 1e-12);
        end

        function test_aubry_above_reference_is_densest_bone(tc, quantity)
            % HU 1500 and 3000 above the reference 1000 map like 1000
            v = map_aubry(quantity, [1000 1500 3000]);
            tc.verifyEqual(v(2:3), [v(1) v(1)], 'AbsTol', 1e-9);
        end

        function test_aubry_below_zero_is_water(tc, quantity)
            v = map_aubry(quantity, [0 -300]);
            tc.verifyEqual(v(2), v(1), 'AbsTol', 1e-9);
        end

        function test_aubry_within_range_unchanged(tc)
            % negative control: inside [0, reference] the mixture model is linear in HU
            v = map_aubry('density', [0 500 1000]);
            tc.verifyEqual(v, [1000 1450 1900], 'AbsTol', 1e-9);
        end

        function test_marsac_babelbrain_limited(tc)
            p = params(1000);
            hu = [300 1000 2500];
            d = medium_pct_density(p, zeros(size(hu)), hu, 1:numel(hu), 'marsac-babelbrain');
            tc.verifyEqual(d, [1000 + 1700 * 0.3, 2700, 2700], 'AbsTol', 1e-9);
        end

        function test_percentile_100_limits_nothing(tc)
            % without a stored reference the largest HU is the reference
            p = params([]);
            hu = [0 500 2000];
            phi = medium_pct_porosity(p, hu);
            tc.verifyEqual(phi, 1 - hu / 2000, 'AbsTol', 1e-12);
        end
    end
end

function p = params(hu_ref)
p = struct();
p.simulation.debug = 0;
p.pct.hu_reference_value = hu_ref;
p.medium_properties.water = struct('density', 1000, 'sound_speed', 1500);
p.medium_properties.skull = struct('density', 1900, 'sound_speed', 3100, 'alpha_power', 1);
end

function v = map_aubry(quantity, hu)
p = params(1000);
idx = 1:numel(hu);
z = zeros(size(hu));
switch quantity
    case 'density'
        v = medium_pct_density(p, z, hu, idx, 'aubry');
    case 'soundspeed'
        v = medium_pct_soundspeed(p, z, z + 1500, hu, idx, 'aubry');
    case 'attenuation'
        v = medium_pct_attenuation(p, z, z, hu, idx, 'aubry', z + 1500);
end
end
