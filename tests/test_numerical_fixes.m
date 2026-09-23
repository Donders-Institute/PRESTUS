classdef test_numerical_fixes < matlab.unittest.TestCase
% TEST_NUMERICAL_FIXES  Unit tests for smoothing edge padding (SMOOTH_IMG),
% phase-optimisation vector orientation (PHASE_OPTIMIZATION_ANNULUS_FULL_CURVE),
% jwave PML sizing (RESOLVE_JWAVE_PML_SIZE) and axial-profile beam direction
% (PLOT_AXIAL_PROFILE).
%
%   Run with:   results = runtests('tests/test_numerical_fixes.m');

    properties
        TmpDir
    end

    methods (TestMethodSetup)
        function make_tmp(tc)
            tc.TmpDir = tempname;
            mkdir(tc.TmpDir);
            tc.addTeardown(@() rmdir(tc.TmpDir, 's'));
        end
    end

    % ------------------------------------------------------------------ %
    %% smooth_img: replicate padding
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'smooth_img'})

        function test_gaussian_constant_volume_edges_preserved(tc)
            img = ones(16, 16, 16);
            out = smooth_img(img, 3, 1, 0, "gaussian");
            tc.verifyEqual(size(out), size(img));
            tc.verifyEqual(out, img, 'AbsTol', 1e-10);
        end

        function test_box_constant_volume_edges_preserved(tc)
            img = ones(16, 16, 16);
            out = smooth_img(img, 3, 1, 0, "box");
            tc.verifyEqual(size(out), size(img));
            tc.verifyEqual(out, img, 'AbsTol', 1e-10);
        end

        function test_binary_mask_touching_border_survives_threshold(tc)
            % A slab filling the lower half of the volume touches the border;
            % zero padding used to erode it at the faces after thresholding.
            img = zeros(16, 16, 16);
            img(:, :, 1:8) = 1;
            out = smooth_img(img, 3, 1, 0.5, "gaussian");
            tc.verifyTrue(all(out(:, :, 1:6) == 1, 'all'));
        end

    end

    % ------------------------------------------------------------------ %
    %% phase_optimization_annulus_full_curve: vector orientation
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'calibration'})

        function test_row_and_column_targets_give_same_error(tc)
            p = annular_params();
            phase = zeros(1, p.transducer.annular.elem_n);
            velocity = 0.02;
            axial_col = (1:0.5:100)';
            target = 0.5 + 0.01 * axial_col;          % arbitrary non-zero target

            err_col = phase_optimization_annulus_full_curve(phase, p, velocity, ...
                axial_col, target, 0, [10 90], 0);
            err_row = phase_optimization_annulus_full_curve(phase, p, velocity, ...
                axial_col, target', 0, [10 90], 0);
            err_all_row = phase_optimization_annulus_full_curve(phase, p, velocity, ...
                axial_col', target', 0, [10 90], 0);

            tc.verifyTrue(isscalar(err_col) && isfinite(err_col));
            tc.verifyEqual(err_row, err_col, 'RelTol', 1e-12);
            tc.verifyEqual(err_all_row, err_col, 'RelTol', 1e-12);
        end

        function test_unknown_forward_model_raises_clear_error(tc)
            % The output used to be named `error`, shadowing error() so this
            % branch crashed with "used before defined" instead of its message.
            p = annular_params();
            p.calibration.forward_model = 'bogus';
            tc.verifyError(@() phase_optimization_annulus_full_curve( ...
                zeros(1, 4), p, 0.02, (1:100)', ones(100, 1)), ?MException);
            try
                phase_optimization_annulus_full_curve(zeros(1, 4), p, 0.02, (1:100)', ones(100, 1));
            catch e
                tc.verifySubstring(e.message, 'unknown forward_model');
            end
        end

    end

    % ------------------------------------------------------------------ %
    %% resolve_jwave_pml_size: sized from c_max
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'jwave_pml'})

        function test_auto_uses_longest_wavelength(tc)
            pml_water = resolve_jwave_pml_size('auto', 3e5, 1500, 0.75e-3);
            pml_skull = resolve_jwave_pml_size('auto', 3e5, 2450, 0.75e-3);
            tc.verifyEqual(pml_water, ceil(2.5 * (1500/3e5) / 0.75e-3));
            tc.verifyEqual(pml_skull, ceil(2.5 * (2450/3e5) / 0.75e-3));
            tc.verifyGreaterThan(pml_skull, pml_water);
        end

        function test_auto_clamped_to_range(tc)
            tc.verifyEqual(resolve_jwave_pml_size('auto', 3e5, 2450, 0.1e-3), 40);
        end

        function test_numeric_passthrough(tc)
            tc.verifyEqual(resolve_jwave_pml_size(17, 3e5, 2450, 0.75e-3), 17);
        end

    end

    % ------------------------------------------------------------------ %
    %% plot_axial_profile: beam direction
    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'axial_profile'})

        function test_minus_z_beam_gives_positive_depths(tc)
            % Transducer on the high-z face, beam travelling towards low z.
            [p, Ipa, P, masks] = axial_profile_inputs(tc.TmpDir);
            exitplane_pos = [5 5 36];
            focus_pos     = [5 5 16];
            plot_axial_profile(Ipa, P, masks, p, exitplane_pos, focus_pos, focus_pos);

            T = readtable(axial_csv(p));
            d = T.depth_from_exit_plane_mm;
            tc.verifyEqual(d(16), 20 * p.grid.resolution_mm, 'AbsTol', 1e-9);   % focus 20 vox deep
            tc.verifyEqual(d(36), 0, 'AbsTol', 1e-9);                          % exit plane
            tc.verifyLessThan(d(40), 0);                                       % behind transducer
        end

        function test_plus_z_beam_unchanged(tc)
            [p, Ipa, P, masks] = axial_profile_inputs(tc.TmpDir);
            exitplane_pos = [5 5 5];
            focus_pos     = [5 5 25];
            plot_axial_profile(Ipa, P, masks, p, exitplane_pos, focus_pos, focus_pos);

            T = readtable(axial_csv(p));
            d = T.depth_from_exit_plane_mm;
            tc.verifyEqual(d, ((1:40)' - 5) * p.grid.resolution_mm, 'AbsTol', 1e-9);
        end

    end

end

% ====================================================================== %
%% Local helpers
% ====================================================================== %

function p = annular_params()
    % CTX_250_026 geometry (config/equipment/CTX_250_026.yaml)
    p = struct();
    p.transducer.freq_hz = 250e3;
    p.transducer.annular.elem_n         = 4;
    p.transducer.annular.elem_id_mm     = [0.0, 30.1788, 42.1388, 51.1088];
    p.transducer.annular.elem_od_mm     = [29.62, 41.58, 50.55, 57.94];
    p.transducer.annular.curv_radius_mm = 62.94;
    p.transducer.annular.dist_geom_ep_mm = 52.38;
    p.medium_properties.water.sound_speed = 1500;
    p.medium_properties.water.density     = 994;
    p.calibration.forward_model = 'oneil';
end

function [p, Ipa, P, masks] = axial_profile_inputs(out_dir)
    p = struct();
    p.grid.resolution_mm = 0.5;
    p.medium_properties = struct('water', struct(), 'brain', struct());
    p.io.dir_img       = out_dir;
    p.io.output_affix  = '_test';
    p.subject_id       = 1;
    p.simulation.medium = 'water';
    Ipa   = rand(9, 9, 40);
    P     = rand(9, 9, 40) * 1e5;
    masks = ones(9, 9, 40);
end

function f = axial_csv(p)
    f = fullfile(p.io.dir_img, sprintf('sub-%03d_%s_axial_profile%s.csv', ...
        p.subject_id, p.simulation.medium, p.io.output_affix));
end
