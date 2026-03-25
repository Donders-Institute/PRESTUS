classdef test_helper < matlab.unittest.TestCase
% TEST_HELPER  Unit tests for functions/helper/
%
%   Run with:   results = runtests('tests/test_helper.m');

    % ------------------------------------------------------------------ %
    %% round_if_integer
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'round_if_integer'})

        function test_round_if_integer_exact(tc)
            tc.verifyEqual(round_if_integer(3.0, 'err'), 3);
        end

        function test_round_if_integer_within_tolerance(tc)
            tc.verifyEqual(round_if_integer(3 + 1e-7, 'err'), 3);
        end

        function test_round_if_integer_array(tc)
            % round_if_integer asserts element-wise; call per element
            tc.verifyEqual(round_if_integer(1.0, 'err'), 1);
            tc.verifyEqual(round_if_integer(2.0, 'err'), 2);
        end

        function test_round_if_integer_non_integer_errors(tc)
            % Any error should be thrown for a non-integer input
            tc.verifyError(@() round_if_integer(3.5, 'not an integer'), ?MException);
        end

    end

    % ------------------------------------------------------------------ %
    %% find_min_factor
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'find_min_factor'})

        function test_find_min_factor_power_of_two(tc)
            % 128 = 2^7; should win over nearby primes in range [120,130]
            result = find_min_factor(120, 130);
            tc.verifyEqual(result, 128);
        end

        function test_find_min_factor_single_element_range(tc)
            result = find_min_factor(16, 16);
            tc.verifyEqual(result, 16);
        end

        function test_find_min_factor_returns_value_in_range(tc)
            lo = 50; hi = 60;
            result = find_min_factor(lo, hi);
            tc.verifyGreaterThanOrEqual(result, lo);
            tc.verifyLessThanOrEqual(result, hi);
        end

    end

    % ------------------------------------------------------------------ %
    %% get_crop_dims
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'get_crop_dims'})

        function test_get_crop_dims_single_voxel(tc)
            img = zeros(10, 10, 10);
            img(5, 5, 5) = 1;
            margin = 2;
            [mn, mx, sz] = get_crop_dims(img, margin);
            tc.verifyEqual(mn, [3, 3, 3]);
            tc.verifyEqual(mx, [7, 7, 7]);
            tc.verifyEqual(sz, [5, 5, 5]);
        end

        function test_get_crop_dims_cuboid_object(tc)
            img = zeros(20, 20, 20);
            img(3:7, 4:8, 5:9) = 1;
            margin = 1;
            [mn, mx, ~] = get_crop_dims(img, margin);
            tc.verifyEqual(mn, [2, 3, 4]);
            tc.verifyEqual(mx, [8, 9, 10]);
        end

        function test_get_crop_dims_grid_size_consistent(tc)
            img = zeros(15, 15, 15);
            img(6:10, 6:10, 6:10) = 1;
            [mn, mx, sz] = get_crop_dims(img, 0);
            tc.verifyEqual(sz, mx - mn + 1);
        end

    end

    % ------------------------------------------------------------------ %
    %% masked_max_3d
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'masked_max_3d'})

        function test_masked_max_3d_known_location(tc)
            vol = zeros(5, 5, 5);
            vol(2, 3, 4) = 99;
            mask = ones(5, 5, 5);
            [val, Ix, Iy, Iz] = masked_max_3d(vol, mask);
            tc.verifyEqual(val, 99);
            tc.verifyEqual([Ix, Iy, Iz], [2, 3, 4]);
        end

        function test_masked_max_3d_mask_excludes_true_max(tc)
            vol = zeros(5, 5, 5);
            vol(1, 1, 1) = 100;  % true max, but outside mask
            vol(3, 3, 3) = 50;   % max within mask
            mask = zeros(5, 5, 5);
            mask(3, 3, 3) = 1;
            [val, ~, ~, ~] = masked_max_3d(vol, mask);
            tc.verifyEqual(val, 50);
        end

    end

    % ------------------------------------------------------------------ %
    %% charm_seg_labels
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'charm_seg_labels'})

        function test_charm_seg_labels_returns_struct(tc)
            labels = charm_seg_labels();
            tc.verifyClass(labels, 'struct');
        end

        function test_charm_seg_labels_expected_fields(tc)
            labels = charm_seg_labels();
            expected = {'wm','gm','csf','skull','skin','skull_cortical','skull_trabecular'};
            for i = 1:numel(expected)
                tc.verifyTrue(isfield(labels, expected{i}), ...
                    sprintf('Missing field: %s', expected{i}));
            end
        end

        function test_charm_seg_labels_values(tc)
            labels = charm_seg_labels();
            tc.verifyEqual(labels.wm,               1);
            tc.verifyEqual(labels.gm,               2);
            tc.verifyEqual(labels.csf,              3);
            tc.verifyEqual(labels.skull,            4);
            tc.verifyEqual(labels.skin,             5);
            tc.verifyEqual(labels.skull_cortical,   7);
            tc.verifyEqual(labels.skull_trabecular, 8);
        end

    end

    % ------------------------------------------------------------------ %
    %% get_flhm_center_position
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'get_flhm_center_position'})

        function test_flhm_symmetric_gaussian(tc)
            x = -10:0.1:10;
            y = exp(-x.^2);  % Gaussian centred at 0
            [centre, ~] = get_flhm_center_position(x, y);
            tc.verifyEqual(centre, 0, 'AbsTol', 0.1);
        end

        function test_flhm_shifted_gaussian(tc)
            x = 0:0.1:20;
            y = exp(-(x - 10).^2);  % Gaussian centred at 10
            [centre, ~] = get_flhm_center_position(x, y);
            tc.verifyEqual(centre, 10, 'AbsTol', 0.1);
        end

    end

    % ------------------------------------------------------------------ %
    %% cast_struct
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'cast_struct'})

        function test_cast_struct_single_level(tc)
            s.a = double(1.0);
            s.b = double(2.0);
            s = cast_struct(s, 'single');
            tc.verifyClass(s.a, 'single');
            tc.verifyClass(s.b, 'single');
        end

        function test_cast_struct_nested(tc)
            s.x = double(1.0);
            s.sub.y = double(2.0);
            s = cast_struct(s, 'single');
            tc.verifyClass(s.x, 'single');
            tc.verifyClass(s.sub.y, 'single');
        end

        function test_cast_struct_preserves_non_numeric(tc)
            s.label = 'hello';
            s.val = double(3.0);
            s = cast_struct(s, 'single');
            tc.verifyClass(s.label, 'char');   % untouched
            tc.verifyClass(s.val,   'single'); % cast
        end

    end

    % ------------------------------------------------------------------ %
    %% get_xyz_mesh
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'get_xyz_mesh'})

        function test_get_xyz_mesh_size(tc)
            img = zeros(3, 4, 5);
            mesh = get_xyz_mesh(img);
            tc.verifyEqual(size(mesh), [3*4*5, 3]);
        end

        function test_get_xyz_mesh_range(tc)
            img = zeros(3, 4, 5);
            mesh = get_xyz_mesh(img);
            tc.verifyEqual(min(mesh(:,1)), 1);
            tc.verifyEqual(max(mesh(:,1)), 3);
            tc.verifyEqual(max(mesh(:,3)), 5);
        end

    end

    % ------------------------------------------------------------------ %
    %% zip_fields / subset_fields
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'struct_utils'})

        function test_zip_fields_order(tc)
            s.alpha = 1;
            s.beta  = 2;
            c = zip_fields(s);
            tc.verifyEqual(c{1}, 'alpha');
            tc.verifyEqual(c{2}, 1);
            tc.verifyEqual(c{3}, 'beta');
            tc.verifyEqual(c{4}, 2);
        end

        function test_subset_fields_copies_requested(tc)
            s.a = 1; s.b = 2; s.c = 3;
            r = subset_fields(s, {'a', 'c'});
            tc.verifyTrue(isfield(r, 'a'));
            tc.verifyTrue(isfield(r, 'c'));
            tc.verifyFalse(isfield(r, 'b'));
        end

    end

end
