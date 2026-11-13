classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_matrix_transducer < matlab.unittest.TestCase
% TEST_MATRIX_TRANSDUCER  Unit tests for matrix-array geometry options (issue #114)
% and transducer markers outside the plotted image (issue #130).
%
%   results = runtests('tests/test_matrix_transducer.m');


    methods (Test)

        %% Clover layout
        function test_clover_default_azimuth_follows_n_leaves(tc)
            for n = [2 3 4]
                [pos, ~, P] = tc.make_clover(n, struct());
                az = sort(mod(rad2deg(atan2(P.cy(:,2), P.cy(:,1))), 360));
                tc.verifyEqual(diff(az), repmat(360/n, n-1, 1), 'AbsTol', 1e-6, ...
                    sprintf('n_leaves = %d', n));
                tc.verifySize(pos, [3, n * P.n_base]);
            end
        end

        function test_clover_leaves_converge_on_parent_centre(tc)
            [pos, aim, P] = tc.make_clover(3, struct());
            % Leaf ROC == ROC_parent: every element on the parent sphere,
            % aiming at its centre
            d = vecnorm(pos - P.pc, 2, 1);
            tc.verifyEqual(d, repmat(P.roc * 1e-3, size(d)), 'AbsTol', 1e-9);
            tc.verifyEqual(aim, repmat(P.pc, 1, size(aim, 2)), 'AbsTol', 1e-12);
        end

        function test_clover_default_elevation_matches_aperture(tc)
            [~, ~, P] = tc.make_clover(3, struct());
            % Adjacent leaf apices one aperture diameter apart on the parent bowl
            chord = norm(P.apex(1,:) - P.apex(2,:));
            tc.verifyEqual(chord, P.aperture, 'RelTol', 1e-9);
        end

        function test_clover_overrides(tc)
            [~, ~, P] = tc.make_clover(2, struct('az_sep_deg', 90, 'elevation_deg', 20));
            az = mod(rad2deg(atan2(P.cy(:,2), P.cy(:,1))), 360);
            tc.verifyEqual(abs(diff(az)), 90, 'AbsTol', 1e-6);
            tilt = acosd(dot(P.pc' - P.apex(1,:)*1e-3, [0 0 1]) / norm(P.pc' - P.apex(1,:)*1e-3));
            tc.verifyEqual(tilt, 20, 'AbsTol', 1e-6);
        end

        function test_single_leaf_stays_on_axis(tc)
            [pos, ~, P] = tc.make_clover(1, struct());
            tc.verifyEqual(pos, P.base, 'AbsTol', 1e-12);
        end

        function test_clover_aims_at_own_leaf_focus_when_radii_differ(tc)
            [pos, aim, P] = tc.make_clover(3, struct(), 50);   % leaf ROC 50, parent 64
            leaf1 = 1:P.n_base;
            % Rigid leaf: its elements stay equidistant (ROC_leaf) from the aim point
            d = vecnorm(pos(:, leaf1) - aim(:, leaf1), 2, 1);
            tc.verifyEqual(d, repmat(50e-3, size(d)), 'AbsTol', 1e-9);
            tc.verifyGreaterThan(norm(aim(:,1) - P.pc), 1e-3);
        end

        function test_clover_that_does_not_fit_errors(tc)
            tc.verifyError(@() tc.make_clover(40, struct()), ?MException);
        end

        %% Element aperture
        function test_elem_diameter_equal_area_and_override(tc)
            m = struct('elem_shape', 'disc', 'elem_width_mm', 3, 'elem_height_mm', 4, ...
                       'curv_radius_mm', 60);
            tc.verifyEqual(pi * (matrix_elem_diameter(m)/2)^2, 12e-6, 'RelTol', 1e-12);

            m.elem_shape = 'bowl';
            a = matrix_elem_diameter(m) / 2;  r = 60e-3;
            tc.verifyEqual(2*pi*r*(r - sqrt(r^2 - a^2)), 12e-6, 'RelTol', 1e-9);

            m.elem_diameter_mm = 5;
            tc.verifyEqual(matrix_elem_diameter(m), 5e-3);
        end

        %% Phase override
        function test_phase_override_replaces_geometric_phases(tc)
            [kgrid, tr, elem_pos_m, tp, fp] = tc.make_karray_inputs();
            [~, tr_geo] = create_matrix_karray(kgrid, kWaveArray, tc.params(), tr, elem_pos_m, tp, fp);
            tc.verifyNotEqual(tr_geo.matrix.elem_phase_deg, zeros(1, 4));

            tr.matrix.elem_phase_override_deg = [0 90 180 450];
            [~, tr_ovr] = create_matrix_karray(kgrid, kWaveArray, tc.params(), tr, elem_pos_m, tp, fp);
            tc.verifyEqual(tr_ovr.matrix.elem_phase_deg, [0 90 180 90], 'AbsTol', 1e-10);
        end

        function test_phase_override_wrong_length_errors(tc)
            [kgrid, tr, elem_pos_m, tp, fp] = tc.make_karray_inputs();
            tr.matrix.elem_phase_override_deg = [0 90];
            tc.verifyError(@() create_matrix_karray(kgrid, kWaveArray, tc.params(), tr, elem_pos_m, tp, fp), ...
                'create_matrix_karray:phaseOverrideSize');
        end

        %% Issue #130: transducer outside the plotted image
        function test_plot_overlay_allows_transducer_outside_image(tc)
            img = rand(20, 20, 20);
            p.grid.resolution_mm = 1;
            p.modules.run_grid_setup = 0;
            p.transducer = struct('type', 'annular', 'annular', struct('curv_radius_mm', 10));
            [~,~,~,~,~,~,~,h] = plot_overlay(img, img, zeros(size(img)), p, {'y', 10}, ...
                [10 10 -15], [10 10 10], [10 10 10], 'show_colorbar', 0);
            close(h);
            tc.verifyError(@() plot_overlay(img, img, zeros(size(img)), p, {'y', 10}, ...
                [10 10 5], [10 10 0], [10 10 10]), ?MException);
        end
    end

    methods (Static, Access = private)

        function [pos, aim, P] = make_clover(n_leaves, clover_extra, roc_leaf)
            if nargin < 3; roc_leaf = 64; end
            P.roc = 64;  P.aperture = 30;
            % Small leaf: 5 elements on the leaf bowl, apex at trans_pos, beam +Z
            tp = [0.01; 0.02; 0.005];
            r  = [0 5 5 5 5];  th = [0 0 pi/2 pi 3*pi/2];
            x = r .* cos(th);  y = r .* sin(th);
            z = roc_leaf - sqrt(roc_leaf^2 - r.^2);
            P.base = tp + [x; y; z] * 1e-3;
            P.n_base = numel(r);
            P.pc = tp + [0; 0; P.roc * 1e-3];

            m = struct('outer_diameter_mm', P.aperture, 'curv_radius_mm', roc_leaf, ...
                       'elem_shape', 'rect');
            m.clover = struct('n_leaves', n_leaves, 'ROC_parent', P.roc);
            f = fieldnames(clover_extra);
            for k = 1:numel(f); m.clover.(f{k}) = clover_extra.(f{k}); end
            params.simulation.debug = 0;

            [pos, aim] = create_clover_array(params, m, P.base, tp, tp + [0; 0; 0.05]);

            % Leaf apices (first element of each leaf is the leaf apex)
            P.apex = (pos(:, 1:P.n_base:end))' * 1e3;            % [n x 3] mm
            P.cy   = P.apex - P.pc' * 1e3;                        % relative to parent centre
        end

        function p = params()
            p.medium_properties.water.sound_speed = 1500;
            p.simulation.debug = 0;
        end

        function [kgrid, tr, elem_pos_m, tp, fp] = make_karray_inputs()
            N = 32;  dx = 1e-3;
            kgrid = kWaveGrid(N, dx, N, dx, N, dx);
            tp = [16 16 4];  fp = [18 15 24];
            tp_m = [kgrid.x_vec(tp(1)); kgrid.y_vec(tp(2)); kgrid.z_vec(tp(3))];
            elem_pos_m = tp_m + [-3 3 -3 3; -3 -3 3 3; 0.2 0.2 0.2 0.2] * 1e-3;
            tr.freq_hz = 500e3;
            tr.matrix = struct('elem_shape', 'rect', 'elem_width_mm', 2, 'elem_height_mm', 2, ...
                'curv_radius_mm', 40, 'elem_p_amp', 1e5, 'is_clover_setup', false);
        end
    end
end
