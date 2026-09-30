classdef test_transducer_placement < matlab.unittest.TestCase
% TEST_TRANSDUCER_PLACEMENT  Tests for preproc_transducer_placement dispatch.
%
%   Manual-mode tests run without any data files. Localite-mode tests
%   require the example XML fixtures and are skipped when absent.
%
%   Run with:
%     results = runtests('tests/test_transducer_placement.m');

    properties (Constant)
        % Path to a Localite InstrumentMarker XML with known 4×4 matrix entries.
        % Update this path if you move the example data.
        LocaliteXML = '/Volumes/2425122.01/v05/data/localite_examples/amydacc/old_format/InstrumentMarker20250401163545816.xml';
    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'placement', 'manual'})

        function test_manual_mode_is_passthrough(tc)
            % trans_pos and focus_pos must be unchanged after manual placement
            p = make_placement_parameters();
            original_trans = p.transducer(1).trans_pos;
            original_focus = p.transducer(1).focus_pos;

            p_out = preproc_transducer_placement(p);

            tc.verifyEqual(p_out.transducer(1).trans_pos, original_trans, ...
                'Manual mode must not modify trans_pos');
            tc.verifyEqual(p_out.transducer(1).focus_pos, original_focus, ...
                'Manual mode must not modify focus_pos');
        end

        function test_default_mode_is_manual(tc)
            % Omitting placement field should behave identically to manual mode
            p = make_placement_parameters();
            p = rmfield(p, 'placement');
            original_trans = p.transducer(1).trans_pos;

            p_out = preproc_transducer_placement(p);

            tc.verifyEqual(p_out.transducer(1).trans_pos, original_trans, ...
                'Default (no placement field) should be a passthrough');
        end

        function test_unknown_mode_errors(tc)
            % An unrecognised placement mode must raise an error
            p = make_placement_parameters();
            p.placement.mode = 'nonsense';
            tc.verifyError(@() preproc_transducer_placement(p), '');
        end

        function test_manual_explicit_voxel_is_passthrough(tc)
            p = make_placement_parameters();
            p.placement.manual.coordinate_space = 'voxel';
            p_out = preproc_transducer_placement(p);
            tc.verifyEqual(p_out.transducer(1).trans_pos, p.transducer(1).trans_pos);
            tc.verifyEqual(p_out.transducer(1).focus_pos, p.transducer(1).focus_pos);
        end

        function test_manual_unknown_coordinate_space_errors(tc)
            p = make_placement_parameters();
            p.placement.manual.coordinate_space = 'lps';
            tc.verifyError(@() preproc_transducer_placement(p), ...
                'PRESTUS:placement:badCoordinateSpace');
        end

        function test_manual_ras_converts_per_transducer(tc)
            % RAS+ mm positions are converted via the T1 affine, per transducer
            [p, A] = make_ras_placement_parameters(tc);
            ras_trans = [p.transducer(1).trans_pos; p.transducer(2).trans_pos];
            ras_focus = [p.transducer(1).focus_pos; p.transducer(2).focus_pos];

            p_out = preproc_transducer_placement(p);

            for ti = 1:2
                exp_trans = round(A \ [ras_trans(ti,:) 1]') + 1;   % 1-based voxel
                exp_focus = round(A \ [ras_focus(ti,:) 1]') + 1;
                tc.verifyEqual(p_out.transducer(ti).trans_pos, exp_trans(1:3)', ...
                    sprintf('Transducer %d trans_pos', ti));
                tc.verifyEqual(p_out.transducer(ti).focus_pos, exp_focus(1:3)', ...
                    sprintf('Transducer %d focus_pos', ti));
            end
            tc.verifyEqual(p_out.placement.manual.coordinate_space, 'voxel', ...
                'Positions are voxels after conversion');
            tc.verifyEqual(p_out.placement.manual.input_trans_pos, ras_trans);
            tc.verifyEqual(p_out.placement.manual.input_focus_pos, ras_focus);
        end

        function test_manual_ras_defers_focal_distance_fallback(tc)
            % The distance fallback must not run on mm positions at load time,
            % but must be filled once positions are voxels
            [p, ~] = make_ras_placement_parameters(tc);
            % As after load_parameters: fields present but empty
            [p.transducer.focal_distance_bowl]   = deal([]);
            [p.transducer.focal_distance_ep]     = deal([]);
            [p.transducer.focal_distance_offset] = deal([]);

            p_loaded = focal_distance_calculation(p);
            tc.verifyEmpty(p_loaded.transducer(1).focal_distance_bowl, ...
                'Fallback must be skipped while positions are RAS mm');

            p_out = preproc_transducer_placement(p_loaded);
            tc.verifyNotEmpty(p_out.transducer(1).focal_distance_bowl, ...
                'Fallback must be filled after conversion to voxels');
        end

        function test_localite_missing_file_errors(tc)
            % Localite mode with no file and no path.localite must raise an error
            p = make_placement_parameters();
            p.placement.mode = 'localite';
            p.placement.localite = struct();   % empty — no file, no path.localite
            tc.verifyError(@() preproc_transducer_placement(p), 'PRESTUS:localite:noSource');
        end

    end

    % ------------------------------------------------------------------ %
    methods (Test, TestTags = {'placement', 'localite'})

        function test_localite_matrix_to_positions_identity(tc)
            % localite_matrix_to_positions with an identity translation should
            % round-trip: a point at a known RAS location maps to a predictable voxel.
            %
            % Uses a synthetic t1_header with a 1mm isotropic voxel grid so that
            % RAS mm ≈ voxel index (modulo offset from image origin).

            p = make_placement_parameters();
            % Matrix origin = bowl rear centre. tracker_to_bowl_mm is the
            % highest-priority offset; reference_distance_mm (legacy) is only
            % used when the transducer geometry cannot give one.
            p.placement.localite.tracker_to_bowl_mm = 0;

            % Build a minimal synthetic NIfTI header: 1 mm iso, no rotation
            t1_header          = struct();
            t1_header.ImageSize = [200, 200, 200];
            t1_header.Transform = struct();
            % affine: RAS = vox0 - [100 100 100] (0-based) → 1-based vox = RAS + 101
            T = eye(4);
            T(1:3, 4) = -[100; 100; 100];
            t1_header.Transform.T = T';   % niftiinfo stores T transposed

            % Identity-like 4×4 Localite matrix: translation = [10; 20; 30; 1]
            % direction column = [0 0 1 0]' (pointing along z)
            coord_matrix = eye(4);
            coord_matrix(1:3, 4) = [10; 20; 30];   % RAS position (mm)
            coord_matrix(1:3, 1) = [0; 0; 1];       % direction

            [trans_pos, focus_pos, trans_pos_ras, ~] = ...
                localite_matrix_to_positions(coord_matrix, t1_header, p);

            % With tracker_to_bowl_mm = 0: trans_pos_ras = [10 20 30]
            tc.verifyEqual(trans_pos_ras, [10, 20, 30], 'AbsTol', 1e-6, ...
                'RAS position should match matrix translation column');

            % 1-based voxel = RAS + [101 101 101] = [111 121 131]
            tc.verifyEqual(trans_pos, [111, 121, 131], 'AbsTol', 0.5, ...
                'Voxel position should be RAS + image origin offset (1-based)');

            % Focus is offset along direction by focal_distance_bowl
            focal_dist = p.transducer(1).focal_distance_bowl;
            expected_focus_ras = [10, 20, 30 + focal_dist];
            tc.verifyEqual(double(focus_pos(3)), expected_focus_ras(3) + 101, ...
                'AbsTol', 1.5, ...
                'Focus z-voxel should equal (focus z-RAS + 101) to within 1.5 vox');
        end

        function test_localite_xml_parse_returns_valid_positions(tc)
            % position_transducer_localite should return finite, non-zero positions
            % when given a real InstrumentMarker XML with actual matrix data.
            tc.assumeTrue(isfile(tc.LocaliteXML), ...
                ['Localite fixture not found — skipping. ', ...
                 'Set LocaliteXML to a valid InstrumentMarker XML path.']);

            p = make_placement_parameters();
            p.placement.localite.reference_distance_mm = 0;

            % Synthetic 1mm iso header centred at [100 100 100]
            t1_header           = struct();
            t1_header.ImageSize = [200, 200, 200];
            t1_header.Transform = struct();
            T = eye(4);
            T(1:3, 4) = -[100; 100; 100];
            t1_header.Transform.T = T';

            [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
                position_transducer_localite(tc.LocaliteXML, t1_header, p);

            tc.verifySize(trans_pos, [1, 3], 'trans_pos must be [1×3]');
            tc.verifySize(focus_pos, [1, 3], 'focus_pos must be [1×3]');
            tc.verifyTrue(all(isfinite(trans_pos)),  'trans_pos must be finite');
            tc.verifyTrue(all(isfinite(focus_pos)),  'focus_pos must be finite');
            tc.verifyTrue(all(isfinite(trans_pos_ras)), 'trans_pos_ras must be finite');
            tc.verifyTrue(all(isfinite(focus_pos_ras)), 'focus_pos_ras must be finite');
            % RAS position should be non-trivial (the XML has real matrix values)
            tc.verifyGreaterThan(norm(trans_pos_ras), 1, ...
                'trans_pos_ras should be a non-trivial RAS coordinate');
        end

    end

end

% =========================================================================
%% Fixture helper
% =========================================================================

function p = make_placement_parameters()
% Minimal parameters struct for placement tests (no simulation, no data files).
    p = struct();
    p.subject_id = 1;
    p.platform   = 'matlab';
    p.placement.mode = 'manual';

    % Minimal transducer with geometry required by localite_matrix_to_positions
    p.transducer(1).type            = 'annular';
    p.transducer(1).trans_pos       = [36, 36, 5];
    p.transducer(1).focus_pos       = [36, 36, 64];
    p.transducer(1).annular.curv_radius_mm   = 63.2;
    p.transducer(1).annular.dist_geom_ep_mm  = 52.38;
    p.transducer(1).freq_hz         = 250e3;

    % focal_distance_bowl is expected by localite_matrix_to_positions;
    % pre-compute it here so tests that call that function directly work
    % without running the full focal_distance_calculation step.
    p.transducer(1).focal_distance_bowl = ...
        p.transducer(1).annular.curv_radius_mm;  % simplified: bowl depth ≈ radius
end

function [p, A] = make_ras_placement_parameters(tc)
% Two transducers in RAS+ mm, plus a synthetic T1 with a non-trivial affine
% (2 mm voxels, shifted origin) written to a temporary folder.
    tmp = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
    A = [2 0 0 -40; 0 2 0 -50; 0 0 2 -30; 0 0 0 1];   % voxel -> RAS+ (column form)
    t1_file = fullfile(tmp, 'sub-001_T1w.nii');
    niftiwrite(zeros(40, 50, 40, 'single'), t1_file);
    info = niftiinfo(t1_file);
    info.Transform = affine3d(A');
    info.TransformName = 'Sform';
    niftiwrite(zeros(40, 50, 40, 'single'), t1_file, info);

    p = make_placement_parameters();
    p.path.anat = tmp;
    p.path.t1_pattern = 'sub-%1$03d_T1w.nii*';
    p.grid.resolution_mm = 1;
    p.placement.manual.coordinate_space = 'ras';
    p.transducer(1).focal_distance_ep     = [];   % as after load_parameters
    p.transducer(1).focal_distance_offset = [];
    p.transducer(2) = p.transducer(1);
    p.transducer(1).trans_pos = [-20.3, 10.1, 30.0];
    p.transducer(1).focus_pos = [-5.0, 12.2, 8.9];
    p.transducer(2).trans_pos = [15.0, -20.4, 25.2];
    p.transducer(2).focus_pos = [4.1, -3.0, 6.0];
end
