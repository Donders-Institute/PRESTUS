classdef test_axisymmetric_phantom < matlab.unittest.TestCase
% TEST_AXISYMMETRIC_PHANTOM  Axisymmetric and phantom grid bookkeeping on small synthetic arrays.
%
% Covers AXISYMMETRIC_OUTPUT_MODE (2D/3D form of the expansion),
% CONVERT_AXISYMMETRIC_TO_2D with the beam axis off centre,
% GRID_TRANSDUCER_LOCATION + GRID_AXISYMMETRY position handling,
% PHANTOM_VOXEL_ORDER (output NIfTIs in the phantom file's voxel order and
% world space) and GRID_TISSUE_SETUP for phantoms (isotropy, plane cut from
% a 3D phantom). No k-Wave run, no subject data.
%
%   results = runtests('tests/test_axisymmetric_phantom.m');

    properties
        Dir
    end

    methods (TestClassSetup)
        function add_paths(~)
            repo = fileparts(fileparts(mfilename('fullpath')));
            addpath(genpath(fullfile(repo, 'functions')));
            addpath(fullfile(repo, 'config'));
            addpath(fullfile(repo, 'tests', 'fixtures'));
        end
    end

    methods (TestMethodSetup)
        function make_dir(tc)
            tc.Dir = tempname;
            mkdir(tc.Dir);
            tc.addTeardown(@() rmdir(tc.Dir, 's'));
        end
    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'axisymmetric'})

        function test_output_mode(tc)
            p = axisym_params([12 4], [12 8], [1 5], [9 5]);
            p.simulation.medium = 'phantom';
            tc.verifyEqual(axisymmetric_output_mode(p), '3d');
            p.simulation.medium = 'water';
            tc.verifyEqual(axisymmetric_output_mode(p), '2d');
            p.modules.run_heating_sims = 1;
            tc.verifyEqual(axisymmetric_output_mode(p), '3d');
            p.grid.axisymmetric_output = '2d';
            tc.verifyError(@() axisymmetric_output_mode(p), 'acoustic_wrapper:axisymOutput2dHeating');
            p.grid.axisymmetric_output = 'x';
            tc.verifyError(@() axisymmetric_output_mode(p), 'acoustic_wrapper:axisymOutput');
            p.grid.dims = [8 8 12];                         % already expanded
            tc.verifyEmpty(axisymmetric_output_mode(p));
            p.grid.dims = [12 4]; p.grid.axisymmetric = 0;
            tc.verifyEmpty(axisymmetric_output_mode(p));
        end

        function test_convert_2d_axis_right_of_centre(tc)
            % Nlateral 10, axis at 8: half-plane has Nr = 3 < 8 columns.
            % The left side beyond the simulated radius repeats the outer
            % column, as radialExpand2DTo3D fills the 3D expansion.
            half = reshape(1:12, 4, 3);
            p = axisym_params([4 3], [4 10], [1 8], [3 8]);
            [sd, p2, seg] = convert_2d(p, half);
            tc.verifyEqual(size(sd.p_max_all), [10 4]);
            tc.verifyEqual(p2.grid.dims, [10 4]);
            tc.verifyEqual(p2.transducer(1).trans_pos, [8 1]);
            vol = radialExpand2DTo3D(half, 'linear', 10, 8);
            tc.verifyEqual(sd.p_max_all, squeeze(vol(8, :, :)), 'AbsTol', 1e-12);
            tc.verifyEqual(seg, uint8(sd.p_max_all));
        end

        function test_convert_2d_axis_left_of_centre(tc)
            half = reshape(1:20, 4, 5);                       % Nlateral 8, axis 4
            p = axisym_params([4 5], [4 8], [1 4], [3 4]);
            sd = convert_2d(p, half);
            vol = radialExpand2DTo3D(half, 'linear', 8, 4);
            tc.verifyEqual(sd.p_max_all, squeeze(vol(4, :, :)), 'AbsTol', 1e-12);
        end

        function test_focus_pos_stays_row_and_beam_axis_inferred(tc)
            % focus_pos(1) > focus_pos(2) used to turn focus_pos into a column,
            % which GRID_AXISYMMETRY then did not flip with trans_pos
            p = make_minimal_parameters(tc.Dir);
            p.simulation.medium = 'phantom';
            p.grid.axisymmetric = 1;
            p.grid.dims = [140 240];
            p.transducer = p.transducer(1);
            p.transducer.trans_pos = [70 1];
            p.transducer.focus_pos = [70 50];
            [~, p] = evalc('grid_transducer_location(p, struct())');
            tc.verifyEqual(p.transducer.focus_pos, [70 50]);
            z = zeros(140, 240);
            p = grid_axisymmetry(p, z, z, [], z);
            tc.verifyTrue(p.grid.axisym_transposed);
            tc.verifyEqual(p.transducer.trans_pos_bilateral, [1 70]);
            tc.verifyEqual(p.transducer.focus_pos_bilateral, [50 70]);
            tc.verifyEqual(p.transducer.focus_pos, [50 1]);
            tc.verifyEqual(p.grid.dims, [240 71]);
        end

    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'phantom_nifti'})

        function test_voxel_order_slices_at_beam_axis(tc)
            % Axis at 5 of 8: the centre (4) is off the beam axis
            half = rand(12, 4);
            p = axisym_params([8 12], [12 8], [1 5], [9 5]);
            p.grid.axisym_transposed = true;
            vol = radialExpand2DTo3D(half, 'linear', 8, 5);
            [out, tf] = phantom_voxel_order(vol, p, struct());
            tc.verifyEqual(out, squeeze(vol(5, :, :)));
            tc.verifyEqual(out(5, :), half(:, 1)', 'AbsTol', 1e-12);
            tc.verifyEmpty(tf);
        end

        function test_voxel_order_itrusst_layout_unchanged(tc)
            % 2D phantom [140 x 240], dim 1 lateral, axis at 70: as before
            p = axisym_params([140 140 240], [240 140], [1 70], [100 70]);
            p.grid.axisym_transposed = true;
            vol = rand(140, 140, 240);
            planimg = phantom_planimg([140 240]);
            [out, tf] = phantom_voxel_order(vol, p, planimg);
            tc.verifyEqual(out, squeeze(vol(round(140/2), :, :)));
            tc.verifyEqual(tf.T, planimg.phantom_header.Transform.T);
        end

        function test_voxel_order_axial_first_phantom(tc)
            % [Naxial x Nlateral] phantom: GRID_AXISYMMETRY kept its order
            p = axisym_params([8 8 12], [12 8], [1 5], [9 5]);
            p.grid.axisym_transposed = false;
            vol = rand(8, 8, 12);
            out = phantom_voxel_order(vol, p, phantom_planimg([12 8]));
            tc.verifyEqual(out, squeeze(vol(5, :, :)).');
            % the mirrored 2D plane likewise
            plane = rand(8, 12);
            tc.verifyEqual(phantom_voxel_order(plane, p, phantom_planimg([12 8])), plane.');
        end

        function test_voxel_order_transposed_phantom(tc)
            % Phantom file [12 x 8] transposed to default_dims [8 12]
            p = axisym_params([8 8 12], [12 8], [1 5], [9 5]);
            p.grid.axisym_transposed = true;
            planimg = phantom_planimg([12 8]);
            planimg.phantom_transposed = true;
            vol = rand(8, 8, 12);
            out = phantom_voxel_order(vol, p, planimg);
            tc.verifyEqual(out, squeeze(vol(5, :, :)).');
            % non-axisymmetric 2D run on the transposed phantom
            q = p; q.grid.axisymmetric = 0;
            g = rand(8, 12);
            tc.verifyEqual(phantom_voxel_order(g, q, planimg), g.');
        end

        function test_voxel_order_plane_of_3d_phantom(tc)
            % Plane y = 3 of a [6 x 5 x 7] phantom: voxel (i, j) of the output
            % lands where voxel (i, 3, j) of the phantom does
            p = axisym_params([6 6 7], [7 6], [1 3], [5 3]);
            p.grid.axisym_transposed = true;
            T = [0.5 0 0 0; 0 0.6 0.1 0; 0 0 0.7 0; 10 20 30 1];
            planimg = phantom_planimg([6 5 7], T);
            planimg.phantom_slice = struct('dim', 2, 'index', 3);
            [out, tf] = phantom_voxel_order(rand(6, 6, 7), p, planimg);
            tc.verifyEqual(size(out), [6 7]);
            for ij = [1 1; 4 2; 6 7]'
                world_out = [ij(1)-1, ij(2)-1, 0, 1] * tf.T;
                world_ph  = [ij(1)-1, 3-1, ij(2)-1, 1] * T;
                tc.verifyEqual(world_out, world_ph, 'AbsTol', 1e-12);
            end
        end

    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'phantom_grid'})

        function test_3d_phantom_plane_through_transducer(tc)
            img = zeros(6, 5, 7); img(:, 4, :) = 1;           % brain only in y = 4
            p = phantom_params(tc.Dir, img, [1 1 1]);
            p.transducer.trans_pos = [3 4 1];
            p.transducer.focus_pos = [3 4 6];
            [~, p2, ~, seg, ~, ~, planimg] = evalc('grid_tissue_setup(p)');
            tc.verifyEqual(planimg.phantom_slice, struct('dim', 2, 'index', 4));
            tc.verifyEqual(p2.grid.dims, [6 7]);
            tc.verifyTrue(all(seg(:) == 1));
            tc.verifyEqual(p2.transducer.trans_pos, [3 1]);
            tc.verifyEqual(p2.transducer.focus_pos, [3 6]);
        end

        function test_3d_phantom_even_centre_warns(tc)
            p = phantom_params(tc.Dir, zeros(6, 4, 7), [1 1 1]);
            p.transducer.trans_pos = [3 1];
            p.transducer.focus_pos = [3 6];
            tc.verifyWarning(@() run_quiet(@() grid_tissue_setup(p)), ...
                'grid_tissue_setup:ambiguousPhantomCentre');
        end

        function test_anisotropic_phantom_errors(tc)
            p = phantom_params(tc.Dir, zeros(6, 5, 7), [1 1 2]);
            tc.verifyError(@() run_quiet(@() grid_tissue_setup(p)), 'grid_tissue_setup:anisotropicPhantom');
        end

        function test_2d_phantom_transpose_recorded(tc)
            p = phantom_params(tc.Dir, zeros(12, 8), [1 1]);
            p.grid.default_dims = [8 12];
            [~, p2, ~, ~, ~, ~, planimg] = evalc('grid_tissue_setup(p)');
            tc.verifyEqual(p2.grid.dims, [8 12]);
            tc.verifyTrue(planimg.phantom_transposed);
            tc.verifyEmpty(planimg.phantom_slice);
        end

    end

end

% =========================================================================
%% Local helpers
% =========================================================================

function p = axisym_params(dims, bilateral_dims, tp_bilateral, fp_bilateral)
% Parameters of an axisymmetric run after GRID_AXISYMMETRY.
    p = struct();
    p.simulation = struct('medium', 'phantom');
    p.modules = struct('run_heating_sims', 0);
    p.grid = struct('dims', dims, 'axisymmetric', 1, 'resolution_mm', 1, ...
        'axisym_bilateral_dims', bilateral_dims);
    p.transducer = struct('trans_pos', [tp_bilateral(1) 1], 'focus_pos', [fp_bilateral(1) 1], ...
        'trans_pos_bilateral', tp_bilateral, 'focus_pos_bilateral', fp_bilateral);
end

function [sd, p, seg] = convert_2d(p, half)
    sd = struct('p_final', half, 'p_max_all', half);
    med = struct('sound_speed', half, 'density', 1000);
    src = struct('p_mask', half > 6);
    [sd, p, seg] = convert_axisymmetric_to_2d(sd, p, uint8(half), uint8(half), med, src, half);
end

function planimg = phantom_planimg(sz, T)
% planimg with the header of a phantom NIfTI of size sz.
    if nargin < 2
        T = diag([1 1 1 1]); T(4, 1:3) = [-5 -6 -7];
    end
    hdr = struct('ImageSize', sz, 'PixelDimensions', ones(1, numel(sz)), 'Transform', affine3d(T));
    planimg = struct('phantom_header', hdr, 'phantom_transposed', false, 'phantom_slice', []);
end

function p = phantom_params(root, img, pix)
% Defaults for an axisymmetric phantom run on img, written as the subject's
% final_tissues.nii.gz with voxel size pix [mm].
    p = make_minimal_parameters(root);
    p.simulation.medium = 'phantom';
    p.grid.axisymmetric = 1;
    p.grid.resolution_mm = pix(1);
    p.grid.default_dims = [size(img, 1), size(img, ndims(img))];
    p.transducer = p.transducer(1);
    m2m = fullfile(root, sprintf('m2m_sub-%03d', p.subject_id));
    mkdir(m2m);
    nii = fullfile(m2m, 'final_tissues.nii');
    niftiwrite(uint8(img), nii);
    info = niftiinfo(nii);
    info.PixelDimensions = pix;
    niftiwrite(uint8(img), nii, info);
    [status, msg] = system(sprintf('gzip -f "%s"', nii));
    assert(status == 0, msg);
end

function run_quiet(fn) %#ok<INUSD>
% Run fn with its console output captured.
    evalc('fn();');
end
