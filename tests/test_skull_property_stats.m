classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_skull_property_stats < matlab.unittest.TestCase
% TEST_SKULL_PROPERTY_STATS  Modelled skull properties in all skull voxels and on the beam path.
%
% A 40 x 40 x 60 grid (1 mm) with a 6 mm skull slab (layers skull_cortical
% and skull_trabecular) crossed at a right angle by the beam: the path
% through skull is 6 mm, the medians are the known slab values, a second
% transducer beside the slab adds a zero-length path, and a model without
% skull layers or without positions gives an empty struct or no path fields.
% The second output holds the full distributions written for the Nexus.
%
%   results = runtests('tests/test_skull_property_stats.m');

    methods (Test)

        function test_all_voxels_and_beam_path(tc)
            [p, masks, c, rho, att] = slab_model();
            s = skull_property_stats(p, masks, c, rho, att, 1);
            % all skull voxels: 4 cortical slices (c 2800) and 2 trabecular (c 2300)
            tc.verifyEqual(s.skull_c_median, 2800);
            tc.verifyEqual(s.skull_c_p10, 2300);
            tc.verifyEqual(s.skull_rho_median, 1900);
            tc.verifyEqual(s.skull_atten_median, 8);
            % beam crosses the slab at a right angle: 6 mm (half-voxel samples, +-1 sample)
            tc.verifyEqual(s.skull_path_mm, 6, 'AbsTol', 0.5);
            tc.verifyEqual(s.skull_path_c_median, 2800, 'AbsTol', 1e-9);
            tc.verifyTrue(all(isfield(s, {'skull_path_rho_median', 'skull_path_atten_median'})));
        end

        function test_full_distributions_for_the_nexus(tc)
            [p, masks, c, rho, att] = slab_model();
            [~, d] = skull_property_stats(p, masks, c, rho, att, 1);
            tc.verifyEqual(d.layers, {'skull_cortical', 'skull_trabecular'});
            tc.verifyEqual(d.n_voxels, 6 * 40 * 40);
            tc.verifyEqual(d.properties.sound_speed.median, 2800);
            tc.verifyEqual(d.properties.alpha_coeff.n, 6 * 40 * 40);
            tc.verifyTrue(isfield(d.properties.density, 'histogram'));
            tc.verifyEqual(d.beam_path.path_mm, 6, 'AbsTol', 0.5);
            tc.verifyEqual(d.beam_path.properties.sound_speed.n, d.beam_path.n_samples);
            % attenuation at the source frequency: alpha_coeff * f_MHz^alpha_power (2 MHz, y = 1)
            p.transducer.freq_hz = 2e6;
            [s2, d2] = skull_property_stats(p, masks, c, rho, att, 1);
            tc.verifyEqual(s2.skull_atten_median, 16);
            tc.verifyEqual(d2.properties.alpha_coeff.median, 8);
        end

        function test_second_transducer_beside_the_slab(tc)
            [p, masks, c, rho, att] = slab_model();
            p.transducer(2) = p.transducer(1);
            p.transducer(2).trans_pos = [20 20 2];      % beam along z, never reaches the slab
            p.transducer(2).focus_pos = [20 20 15];
            s = skull_property_stats(p, masks, c, rho, att, 1);
            tc.verifyEqual(numel(s.skull_path_mm), 2);
            tc.verifyEqual(s.skull_path_mm(2), 0);
        end

        function test_no_skull_or_no_positions(tc)
            [p, masks, c, rho, att] = slab_model();
            p2 = p;
            p2.layers = rmfield(p2.layers, {'skull_cortical', 'skull_trabecular'});
            tc.verifyEmpty(fieldnames(skull_property_stats(p2, masks, c, rho, att, 1)));
            p3 = p;
            p3.transducer = rmfield(p3.transducer, 'trans_pos');
            s = skull_property_stats(p3, masks, c, rho, att, 1);
            tc.verifyTrue(isfield(s, 'skull_c_median'));
            tc.verifyFalse(isfield(s, 'skull_path_mm'));
        end

    end
end

function [p, masks, c, rho, att] = slab_model()
    p = struct();
    p.grid = struct('resolution_mm', 1);
    p.medium_properties = struct('water', struct(), 'brain', struct(), ...
        'skull_cortical', struct(), 'skull_trabecular', struct());
    p.layers = struct('water', 0, 'brain', 1, 'skull_cortical', 2, 'skull_trabecular', 3);
    p.transducer = struct('freq_hz', 1e6, 'trans_pos', [20 20 5], 'focus_pos', [20 20 50]);
    masks = ones(40, 40, 60) * 2;                     % brain (index 2 in medium_properties)
    masks(:, :, 1:19) = 1;                            % water
    masks(:, :, [20 21 24 25]) = 3;                   % cortical
    masks(:, :, [22 23]) = 4;                         % trabecular
    c = 1500 * ones(size(masks));   c(masks == 2) = 1560; c(masks == 3) = 2800; c(masks == 4) = 2300;
    rho = 1000 * ones(size(masks)); rho(masks == 3) = 1900; rho(masks == 4) = 1700;
    att = zeros(size(masks));       att(masks == 3) = 8;    att(masks == 4) = 6;
end
