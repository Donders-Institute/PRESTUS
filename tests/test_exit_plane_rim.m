classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_exit_plane_rim < matlab.unittest.TestCase
% TEST_EXIT_PLANE_RIM  A missing dist_geom_ep_mm puts the exit plane at the bowl rim.
%
% VALIDATE_ANNULAR_TRANSDUCER fills dist_geom_ep_mm = sqrt(R^2 - (D/2)^2) at
% load (D = largest element outer diameter), and EXIT_PLANE_DEPTH_MM's fallback
% for structs that bypass the loader gives the same apex-to-exit-plane distance
% (outer_diameter_mm for matrix arrays). Given values are kept; flat elements
% have the exit plane at the face.
%
%   results = runtests('tests/test_exit_plane_rim.m');

    properties (Constant)
        Annular = struct('elem_n', 4, 'elem_id_mm', [0 32.9 46.0 56.0], ...
            'elem_od_mm', [32.3 45.4 55.4 64], 'curv_radius_mm', 63.2, ...
            'elem_phase_deg', [0 0 0 0], 'elem_p_amp', 1e5)
    end

    methods (Test, TestTags = {'placement'})

        function test_annular_fill_matches_rim_fallback(tc)
            g = tc.Annular;
            R = g.curv_radius_mm;
            rim_depth = tc.fallback_depth(g);                   % today's placement fallback
            tr = validate_annular_transducer(struct('type', 'annular', 'annular', g), 1);
            tc.verifyEqual(tr.annular.dist_geom_ep_mm, sqrt(R^2 - 32^2), 'AbsTol', 1e-12);
            tc.verifyEqual(annular_bowl_depth_mm(tr.annular), rim_depth, 'AbsTol', 1e-12);
            tc.verifyEqual(exit_plane_depth_mm(tr.annular), rim_depth, 'AbsTol', 1e-12);
            tc.verifyWarningFree(@() exit_plane_depth_mm(tr.annular));
        end

        function test_annular_fill_prints_notice(tc)
            out = evalc('validate_annular_transducer(struct(''type'', ''annular'', ''annular'', tc.Annular), 1);');
            tc.verifySubstring(out, 'derived from the bowl rim');
        end

        function test_given_value_is_kept(tc)
            g = tc.Annular;
            g.dist_geom_ep_mm = 52.38;
            tr = validate_annular_transducer(struct('type', 'annular', 'annular', g), 1);
            tc.verifyEqual(tr.annular.dist_geom_ep_mm, 52.38);
            tc.verifyEqual(exit_plane_depth_mm(tr.annular), 63.2 - 52.38, 'AbsTol', 1e-12);
        end

        function test_flat_annular_exit_plane_is_the_face(tc)
            g = struct('elem_n', 1, 'elem_id_mm', 0, 'elem_od_mm', 20, 'curv_radius_mm', Inf, ...
                'elem_phase_deg', 0, 'elem_p_amp', 1e5);
            tr = validate_annular_transducer(struct('type', 'annular', 'annular', g), 1);
            tc.verifyEqual(tr.annular.dist_geom_ep_mm, Inf);
            tc.verifyEqual(exit_plane_depth_mm(tr.annular), 0);
        end

        function test_matrix_fallback_uses_outer_diameter(tc)
            % as validate_matrix_transducer: D = outer_diameter_mm, not the element size
            g = struct('curv_radius_mm', 72, 'outer_diameter_mm', 60, 'elem_od_mm', 4);
            d = tc.fallback_depth(g);
            tc.verifyEqual(d, 72 - sqrt(72^2 - 30^2), 'AbsTol', 1e-12);
            g.dist_geom_ep_mm = sqrt(72^2 - 30^2);              % the value the loader fills
            tc.verifyEqual(exit_plane_depth_mm(g), d, 'AbsTol', 1e-12);
        end
    end

    methods (Access = private)
        function d = fallback_depth(tc, g)
            d = tc.verifyWarning(@() exit_plane_depth_mm(g), 'prestus:placement:rimExitPlane');
        end
    end
end
