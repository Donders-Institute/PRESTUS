classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_babelbrain_insitu_isppa < matlab.unittest.TestCase
% TEST_BABELBRAIN_INSITU_ISPPA  Unit tests for the Isppa handed to BabelBrain thermal.
%
% BabelBrain's BHTE normalises the field to max p^2/(2*rho*c) inside the brain,
% so the value PRESTUS passes must be the in-brain peak of the scaled field —
% not the free-water calibration ISPPA (which ignores transcranial loss).
%
%   Run with:   results = runtests('tests/test_babelbrain_insitu_isppa.m');

    methods (Test, TestTags = {'unit', 'quick', 'nogpu'})

        function test_uses_brain_peak_not_global_peak(tc)
            [p, sd, masks] = tc.setup_field();
            % Global peak sits in the non-brain half; brain peak is lower.
            isppa = babelbrain_insitu_isppa(p, sd, struct(), masks);
            rho_c = p.medium_properties.brain.sound_speed * p.medium_properties.brain.density;
            tc.verifyEqual(isppa, (2e5)^2 / (2 * rho_c) * 1e-4, 'RelTol', 1e-12);
        end

        function test_prefers_babelbrain_material_map(tc)
            [p, sd, masks] = tc.setup_field();
            % BabelBrain map labels the high-pressure half as brain (label 4).
            sd.babelbrain_material_map = zeros(size(masks));
            sd.babelbrain_material_map(:, :, 1:5) = 4;
            isppa = babelbrain_insitu_isppa(p, sd, struct(), masks);
            rho_c = p.medium_properties.brain.sound_speed * p.medium_properties.brain.density;
            tc.verifyEqual(isppa, (5e5)^2 / (2 * rho_c) * 1e-4, 'RelTol', 1e-12);
        end

        function test_uses_voxelwise_rho_c(tc)
            [p, sd, masks] = tc.setup_field();
            km.sound_speed = 1500 * ones(size(masks));
            km.density     = 1000 * ones(size(masks));
            isppa = babelbrain_insitu_isppa(p, sd, km, masks);
            tc.verifyEqual(isppa, (2e5)^2 / (2 * 1500 * 1000) * 1e-4, 'RelTol', 1e-12);
        end

        function test_no_brain_errors(tc)
            [p, sd, masks] = tc.setup_field();
            masks(:) = 1;   % water only
            tc.verifyError(@() babelbrain_insitu_isppa(p, sd, struct(), masks), ...
                'babelbrain_insitu_isppa:noBrainMask');
        end

    end

    methods (Access = private)

        function [p, sd, masks] = setup_field(~)
            p.medium_properties.water = struct('sound_speed', 1500, 'density', 994);
            p.medium_properties.brain = struct('sound_speed', 1546, 'density', 1046);
            masks = ones(6, 6, 10);           % 1 = water
            masks(:, :, 6:10) = 2;            % 2 = brain
            sd.p_max_all = 5e5 * ones(size(masks));
            sd.p_max_all(:, :, 6:10) = 1e5;
            sd.p_max_all(3, 3, 8) = 2e5;      % brain peak
        end

    end

end
