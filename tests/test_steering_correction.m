classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'unit'}) test_steering_correction < matlab.unittest.TestCase
% TEST_STEERING_CORRECTION  Steering correction of matrix arrays, flat arrays and
% the matrix Rayleigh model, on synthetic data.
%
%   results = runtests('tests/test_steering_correction.m');

    properties
        repo_root
        tmp
    end

    methods (TestClassSetup)
        function setup_paths(tc)
            tc.repo_root = fileparts(fileparts(mfilename('fullpath')));
            tc.tmp = tempname; mkdir(tc.tmp);
        end
    end

    methods (TestClassTeardown)
        function cleanup_tmp(tc)
            if isfolder(tc.tmp); rmdir(tc.tmp, 's'); end
        end
    end

    methods (Test)

        %% Fit
        function test_quadratic_recovers_smooth_map(tc)
            [X, sim] = tc.targets(60);
            meas = tc.k_true(X) .* sim;
            sc = fit_steering_correction(X, sim, meas);
            tc.verifyEqual(sc.model, 'quadratic');
            Xt = (X(1:30, :) + X(31:60, :)) / 2;   % unmeasured targets inside the measured range
            tc.verifyEqual(steering_correction_factor(sc, Xt), tc.k_true(Xt), 'AbsTol', 1e-9);
            tc.verifyLessThan(sc.cv.quadratic_med_pct, 1e-6);
        end

        function test_auto_keeps_scale_without_structure(tc)
            [X, sim] = tc.targets(60);
            rng(3); meas = 0.7 * sim .* (1 + 0.03 * randn(size(sim)));
            sc = fit_steering_correction(X, sim, meas);
            tc.verifyEqual(sc.model, 'scale');
            tc.verifyEqual(steering_correction_factor(sc, X(1:3, :)), repmat(sc.scale, 3, 1));
            tc.verifyEqual(sc.scale, 0.7, 'RelTol', 0.02);
        end

        function test_auto_needs_min_targets(tc)
            [X, sim] = tc.targets(12);
            sc = fit_steering_correction(X, sim, tc.k_true(X) .* sim);
            tc.verifyEqual(sc.model, 'scale');   % 12 < min_targets (15)
            sc = fit_steering_correction(X, sim, tc.k_true(X) .* sim, struct('model', 'quadratic'));
            tc.verifyEqual(sc.model, 'quadratic');
            tc.verifyError(@() fit_steering_correction(X(1:8, :), sim(1:8), sim(1:8), ...
                struct('model', 'quadratic')), 'fit_steering_correction:tooFew');
        end

        function test_factor_does_not_extrapolate(tc)
            [X, sim] = tc.targets(60);
            sc = fit_steering_correction(X, sim, tc.k_true(X) .* sim);
            far = sc.bounds_mm(2, :) + [0 0 30];
            edge = sc.bounds_mm(2, :);
            tc.verifyWarning(@() steering_correction_factor(sc, far), 'steering_correction_factor:outside');
            warning('off', 'steering_correction_factor:outside');
            c = onCleanup(@() warning('on', 'steering_correction_factor:outside'));
            tc.verifyEqual(steering_correction_factor(sc, far), steering_correction_factor(sc, edge), 'AbsTol', 1e-12);
        end

        function test_build_and_library_round_trip(tc)
            [X, sim] = tc.targets(40);
            T = table(X(:, 1), X(:, 2), X(:, 3), tc.k_true(X) .* sim / 1e6 * 2, sim / 1e6, repmat(20, 40, 1), ...
                'VariableNames', {'target_x_mm', 'target_y_mm', 'target_z_mm', 'peak_pressure_mpa', ...
                'model_pressure_mpa', 'drive_v'});
            T.drive_v(1:2:end) = 40;                       % half the rows at twice the voltage
            T.peak_pressure_mpa(1:2:end) = 2 * T.peak_pressure_mpa(1:2:end);
            csv = fullfile(tc.tmp, 'meas.csv'); writetable(T, csv);
            out = fullfile(tc.tmp, 'EXAMPLE_MATRIX.yaml');
            tr = tc.flat_tr();
            sc = build_steering_correction(struct(), tr, csv, struct('out_file', out, 'serial', 'EXAMPLE_MATRIX'));
            tc.verifyEqual(sc.model, 'quadratic');
            tc.verifyEqual(sc.model_source, 'file');
            tc.verifyEqual(sc.drive_ref_v, 30);
            lib = yaml.loadFile(out, 'ConvertToArray', true);
            % median voltage 30 V: rows scaled by 30/20 and 30/40 give measured / modelled = 3 k
            tc.verifyEqual(steering_correction_factor(lib.steering_correction, X(1:5, :)), ...
                3 * tc.k_true(X(1:5, :)), 'RelTol', 1e-6);
        end

        %% Flat arrays and element rotation
        function test_flat_layout_faces_beam_axis(tc)
            L = matrix_element_layout(struct(), tc.flat_tr());
            tc.verifyEqual(L.ez, repmat([0; 0; 1], 1, L.n), 'AbsTol', 1e-12);
            tc.verifyEqual(L.pos_m(3, :), zeros(1, L.n), 'AbsTol', 1e-12);
        end

        function test_element_rotation_turns_rectangle(tc)
            tr = tc.flat_tr();
            tr.matrix.elem_rot_deg = 90 * ones(1, 16);
            tr.matrix.matrix_shape.type = 'define_here';
            L0 = matrix_element_layout(struct(), tc.flat_tr());
            L = matrix_element_layout(struct(), tr);
            tc.verifyEqual(L.ex, L0.ey, 'AbsTol', 1e-12);
            tc.verifyEqual(L.ey, -L0.ex, 'AbsTol', 1e-12);
        end

        function test_flat_file_layout_and_rotation_column(tc)
            f = fullfile(tc.tmp, 'flat_coords.csv');
            writetable(table([-3; 3], [0; 0], [0; 0], [0; 90], 'VariableNames', {'X', 'Y', 'Z', 'rot'}), f);
            tr = tc.flat_tr();
            tr.matrix.matrix_shape = struct('type', 'extract_from_file', 'extract_from_file', struct( ...
                'file_path', f, 'start_row', 1, 'start_col', 1, 'rot_col', 4, 'elem_n', 2, ...
                'select_random_subset', false, 'project_on_new_ROC', false));
            L = matrix_element_layout(struct(), tr);
            tc.verifyEqual(L.pos_m, [-3 3; 0 0; 0 0] * 1e-3, 'AbsTol', 1e-12);
            tc.verifyEqual(L.ex(:, 2), [0; 1; 0], 'AbsTol', 1e-12);
        end

        %% Rayleigh model
        function test_rayleigh_disc_on_axis_matches_piston(tc)
            f = 500e3; c = 1500; k = 2 * pi * f / c; a = 5e-3;
            L = struct('pos_m', [0; 0; 0], 'ex', [1; 0; 0], 'ey', [0; 1; 0], 'ez', [0; 0; 1], ...
                'shape', 'disc', 'width_m', [], 'height_m', [], 'diameter_m', 2 * a, 'n', 1);
            z = (20:10:80) * 1e-3;
            p = matrix_rayleigh_pressure(L, f, c, 1e5, 0, [zeros(2, numel(z)); z], struct('patch_per_lambda', 30));
            p_exact = 2 * 1e5 * abs(sin(k / 2 * (sqrt(z.^2 + a^2) - z)));
            tc.verifyEqual(abs(p(:))', p_exact, 'RelTol', 0.02);
        end

        function test_focal_peak_lands_at_steered_target(tc)
            tr = tc.flat_tr();
            T = matrix_focal_peak(struct('medium_properties', struct('water', struct('sound_speed', 1500))), ...
                tr, [0 0 30; 3 0 30]);
            tc.verifyLessThanOrEqual(abs(T.peak_x_mm - T.target_x_mm), 1);
            tc.verifyLessThanOrEqual(abs(T.peak_y_mm - T.target_y_mm), 1);
            tc.verifyGreaterThan(T.peak_pressure_pa(1), T.peak_pressure_pa(2));   % steering off-axis costs pressure
        end

        %% Applied to the k-Wave source
        function test_karray_drive_scaled_by_correction(tc)
            kgrid = kWaveGrid(48, 0.5e-3, 48, 0.5e-3, 64, 0.5e-3);
            params = struct('medium_properties', struct('water', struct('sound_speed', 1500)), ...
                'simulation', struct('debug', 0));
            tr = tc.flat_tr();
            pos = [-1.5 1.5 -1.5 1.5; -1.5 -1.5 1.5 1.5; 0 0 0 0] * 1e-3 + [kgrid.x_vec(24); kgrid.y_vec(24); kgrid.z_vec(5)];
            tr.matrix.steering_correction = struct('model', 'scale', 'scale', 0.8);
            [~, tr1] = create_matrix_karray(kgrid, kWaveArray(), params, tr, pos, [24 24 5], [24 24 50]);
            tc.verifyEqual(tr1.matrix.steering_k, 0.8);
            tc.verifyEqual(tr1.matrix.elem_p_amp, 0.8 * 1e5 * ones(1, 4), 'RelTol', 1e-12);
            tr.matrix.use_steering_correction = false;
            [~, tr2] = create_matrix_karray(kgrid, kWaveArray(), params, tr, pos, [24 24 5], [24 24 50]);
            tc.verifyEqual(tr2.matrix.elem_p_amp, 1e5 * ones(1, 4));
        end
    end

    methods (Static)
        function [X, sim] = targets(n, seed)
            % Steering targets in a cone in front of the array and a smooth modelled peak
            if nargin < 2; seed = 1; end
            rng(seed);
            z = 30 + 50 * rand(n, 1); r = (5 + 0.3 * (z - 30)) .* sqrt(rand(n, 1)); t = 2 * pi * rand(n, 1);
            X = [r .* cos(t), r .* sin(t), z];
            sim = 1e6 * (1.2 - 0.005 * (z - 30) - 0.004 * r);
        end

        function k = k_true(X)
            k = 0.7 + 0.002 * X(:, 1) - 0.001 * X(:, 2) + 1e-5 * (X(:, 3) - 55).^2 - 2e-5 * X(:, 1).^2;
        end

        function tr = flat_tr()
            % A generic flat 4 x 4 array of 3 x 3 mm elements at 500 kHz
            tr = struct('freq_hz', 500e3);
            tr.matrix = struct('elem_shape', 'rect', 'elem_width_mm', 3, 'elem_height_mm', 3, ...
                'elem_p_amp', 1e5, 'outer_diameter_mm', 30, 'is_curved', false, 'curv_radius_mm', inf, ...
                'is_clover_setup', false);
            tr.matrix.matrix_shape = struct('type', 'define_here', 'define_here', struct('grid_shape', ...
                struct('type', 'rect', 'rect', struct('elem_n_row', 4, 'elem_n_col', 4, ...
                'elem_spacing_height_mm', 0.5, 'elem_spacing_width_mm', 0.5, 'sparsity_factor', 1))));
        end
    end
end
