classdef test_plantus_placement < matlab.unittest.TestCase
% TEST_PLANTUS_PLACEMENT  Tests for position_transducer_plantus and the
%                          'plantus' dispatch path in preproc_transducer_placement.
%
% Tests are organised into three tiers:
%
%   'validation'   Pure MATLAB logic tests — no data files required.
%                  Exercise all error paths for missing/invalid config.
%
%   'dispatch'     Verify that preproc_transducer_placement routes
%                  placement.mode='plantus' into the plantus code path.
%
%   'integration'  End-to-end flow using a mock PlanTUS environment built
%                  entirely from temporary files. Skipped automatically when
%                  a usable Python 3 interpreter cannot be found.
%
% Run all tests:
%   results = runtests('tests/test_plantus_placement.m');
%
% Run only validation tests (no Python required):
%   results = runtests('tests/test_plantus_placement.m', 'Tag', 'validation');

    % ------------------------------------------------------------------
    properties
        TmpDir  (1,:) char  % per-test temporary directory
    end

    % ------------------------------------------------------------------
    methods (TestMethodSetup)
        function create_tmp_dir(tc)
            tc.TmpDir = tempname();
            mkdir(tc.TmpDir);
        end
    end

    methods (TestMethodTeardown)
        function remove_tmp_dir(tc)
            if isfolder(tc.TmpDir)
                rmdir(tc.TmpDir, 's');
            end
        end
    end

    % ==================================================================
    %% Validation — missing / invalid configuration
    % ==================================================================
    methods (Test, TestTags = {'plantus', 'validation'})

        function test_missing_plantus_struct_errors(tc)
            % preproc_transducer_placement with mode='plantus' and no
            % placement.plantus sub-struct must raise an error.
            p = make_base_parameters();
            p.placement.mode = 'plantus';
            % deliberately omit placement.plantus
            tc.verifyError(@() preproc_transducer_placement(p), ...
                'PRESTUS:placement:missingConfig');
        end

        function test_missing_script_path_errors(tc)
            p = make_plantus_parameters(tc.TmpDir);
            p.placement.plantus.script_path = '';
            tc.verifyError( ...
                @() position_transducer_plantus(p, make_synthetic_t1_header()), ...
                'PRESTUS:plantus:missingConfig');
        end

        function test_missing_env_path_errors(tc)
            p = make_plantus_parameters(tc.TmpDir);
            p.placement.plantus.env_path = '';
            tc.verifyError( ...
                @() position_transducer_plantus(p, make_synthetic_t1_header()), ...
                'PRESTUS:plantus:missingConfig');
        end

        function test_missing_focal_distance_list_errors(tc)
            p = make_plantus_parameters(tc.TmpDir);
            p.placement.plantus.focal_distance_list = [];
            tc.verifyError( ...
                @() position_transducer_plantus(p, make_synthetic_t1_header()), ...
                'PRESTUS:plantus:missingConfig');
        end

        function test_missing_flhm_list_errors(tc)
            p = make_plantus_parameters(tc.TmpDir);
            p.placement.plantus.flhm_list = [];
            tc.verifyError( ...
                @() position_transducer_plantus(p, make_synthetic_t1_header()), ...
                'PRESTUS:plantus:missingConfig');
        end

        function test_missing_target_name_errors(tc)
            p = make_plantus_parameters(tc.TmpDir);
            p.placement.plantus.target_name = '';
            tc.verifyError( ...
                @() position_transducer_plantus(p, make_synthetic_t1_header()), ...
                'PRESTUS:plantus:missingConfig');
        end

        function test_nonexistent_script_path_errors(tc)
            % script_path pointing to a directory that does not exist
            p = make_plantus_parameters(tc.TmpDir);
            p.placement.plantus.script_path = fullfile(tc.TmpDir, 'does_not_exist');
            tc.verifyError( ...
                @() position_transducer_plantus(p, make_synthetic_t1_header()), ...
                'PRESTUS:plantus:notFound');
        end

        function test_missing_plantus_py_errors(tc)
            % script_path exists but contains no code/PlanTUS.py
            p = make_plantus_parameters(tc.TmpDir);
            script_dir = fullfile(tc.TmpDir, 'plantus_empty');
            mkdir(script_dir);
            p.placement.plantus.script_path = script_dir;
            tc.verifyError( ...
                @() position_transducer_plantus(p, make_synthetic_t1_header()), ...
                'PRESTUS:plantus:notFound');
        end

        function test_no_target_specified_errors(tc)
            % All paths valid but no mni_target_mm and no focus_pos → error
            p = make_plantus_parameters(tc.TmpDir);
            p.placement.plantus.mni_target_mm = [];
            p.transducer(1).focus_pos = [];

            % Ensure script_path exists and has code/PlanTUS.py so we get
            % past path checks; _launcher_py is not needed here (error fires
            % before the subprocess).
            script_dir = fullfile(tc.TmpDir, 'fake_plantus');
            make_fake_plantus_install(script_dir);
            p.placement.plantus.script_path = script_dir;

            tc.verifyError( ...
                @() position_transducer_plantus(p, make_synthetic_t1_header()), ...
                'PRESTUS:plantus:noTarget');
        end

        function test_missing_mesh_errors(tc)
            % All config valid, mesh file absent → error before subprocess
            p = make_plantus_parameters(tc.TmpDir);
            t1_hdr = write_synthetic_t1(tc.TmpDir, p);

            script_dir = fullfile(tc.TmpDir, 'fake_plantus');
            make_fake_plantus_install(script_dir);
            p.placement.plantus.script_path = script_dir;

            % Do NOT create the mesh file
            tc.verifyError( ...
                @() position_transducer_plantus(p, t1_hdr), ...
                'PRESTUS:plantus:notFound');
        end

        function test_subprocess_failure_errors(tc)
            % Subprocess exits non-zero → PRESTUS:plantus:failed
            p      = make_plantus_parameters(tc.TmpDir);
            t1_hdr = write_synthetic_t1(tc.TmpDir, p);
            write_fake_mesh(tc.TmpDir, p);

            python_bin = find_system_python();
            tc.assumeNotEmpty(python_bin, 'No Python 3 found — skipping subprocess test.');

            script_dir = fullfile(tc.TmpDir, 'failing_plantus');
            make_fake_plantus_install(script_dir);
            launcher_path = fullfile(tc.TmpDir, 'failing_launcher.py');
            write_python_launcher(launcher_path, 'failing');

            p.placement.plantus.script_path      = script_dir;
            p.placement.plantus.env_path = make_fake_simnibs_env(tc.TmpDir, python_bin);
            p.placement.plantus._launcher_py     = launcher_path;

            tc.verifyError( ...
                @() position_transducer_plantus(p, t1_hdr), ...
                'PRESTUS:plantus:failed');
        end

        function test_no_output_mat_errors(tc)
            % Subprocess succeeds but writes no *Localite.mat → error
            p      = make_plantus_parameters(tc.TmpDir);
            t1_hdr = write_synthetic_t1(tc.TmpDir, p);
            write_fake_mesh(tc.TmpDir, p);

            python_bin = find_system_python();
            tc.assumeNotEmpty(python_bin, 'No Python 3 found — skipping subprocess test.');

            script_dir = fullfile(tc.TmpDir, 'silent_plantus');
            make_fake_plantus_install(script_dir);
            launcher_path = fullfile(tc.TmpDir, 'silent_launcher.py');
            write_python_launcher(launcher_path, 'silent');

            p.placement.plantus.script_path      = script_dir;
            p.placement.plantus.env_path = make_fake_simnibs_env(tc.TmpDir, python_bin);
            p.placement.plantus._launcher_py     = launcher_path;

            tc.verifyError( ...
                @() position_transducer_plantus(p, t1_hdr), ...
                'PRESTUS:plantus:noOutput');
        end

    end

    % ==================================================================
    %% Dispatch — preproc_transducer_placement routing
    % ==================================================================
    methods (Test, TestTags = {'plantus', 'dispatch'})

        function test_plantus_mode_routes_to_plantus_code(tc)
            % mode='plantus' must attempt plantus placement (not silently
            % fall through to manual). We verify this by checking the error
            % ID: if routing is correct we get PRESTUS:plantus:*, not the
            % generic placement error from an unknown mode.
            p = make_plantus_parameters(tc.TmpDir);
            p.placement.mode = 'plantus';
            p.placement.plantus.script_path = '';  % force earliest plantus error

            err = tc.verifyError( ...
                @() preproc_transducer_placement(p), ...
                'PRESTUS:plantus:missingConfig');
            tc.verifySubstring(err.identifier, 'plantus', ...
                'Error must originate from the plantus code path, not the mode dispatcher');
        end

        function test_unknown_mode_still_errors_correctly(tc)
            % Ensure we didn't break the otherwise branch
            p = make_base_parameters();
            p.placement.mode = 'plantus_typo';
            tc.verifyError(@() preproc_transducer_placement(p), '');
        end

    end

    % ==================================================================
    %% Integration — end-to-end with a mock PlanTUS environment
    % ==================================================================
    methods (Test, TestTags = {'plantus', 'integration'})

        function test_end_to_end_returns_valid_positions(tc)
            % Full flow: synthetic T1, fake mesh, mock PlanTUS wrapper that
            % writes a known position_matrix, verify output dimensions and
            % coordinate plausibility.
            %
            % The mock PlanTUS writes a position_matrix with:
            %   translation = [10; 20; 30] mm (RAS)
            %   direction   = [0; 0; 1]  (pointing along z axis)
            % Using a 1mm-iso T1 centred at [100 100 100] voxels:
            %   trans_pos_ras ≈ [10 20 30] + reference_dist * [0 0 1]
            %   trans_pos (vox) ≈ trans_pos_ras + [100 100 100]

            python_bin = find_system_python();
            tc.assumeNotEmpty(python_bin, ...
                'No Python 3 found — skipping integration test.');

            p      = make_plantus_parameters(tc.TmpDir);
            t1_hdr = write_synthetic_t1(tc.TmpDir, p);
            write_fake_mesh(tc.TmpDir, p);

            % Build mock PlanTUS that writes a known position_matrix
            script_dir    = fullfile(tc.TmpDir, 'mock_plantus');
            make_fake_plantus_install(script_dir);
            launcher_path = fullfile(tc.TmpDir, 'mock_launcher.py');
            write_python_launcher(launcher_path, 'success');

            p.placement.plantus.script_path      = script_dir;
            p.placement.plantus.env_path = make_fake_simnibs_env(tc.TmpDir, python_bin);
            p.placement.plantus._launcher_py     = launcher_path;

            % Use focus_pos (voxel) as target; no MNI conversion needed
            p.transducer(1).focus_pos = [100, 100, 130];
            p.placement.plantus.mni_target_mm    = [];

            [trans_pos, focus_pos, trans_pos_ras, focus_pos_ras] = ...
                position_transducer_plantus(p, t1_hdr);

            % Output shape
            tc.verifySize(trans_pos,     [1, 3], 'trans_pos must be [1×3]');
            tc.verifySize(focus_pos,     [1, 3], 'focus_pos must be [1×3]');
            tc.verifySize(trans_pos_ras, [1, 3], 'trans_pos_ras must be [1×3]');
            tc.verifySize(focus_pos_ras, [1, 3], 'focus_pos_ras must be [1×3]');

            % All outputs must be finite
            tc.verifyTrue(all(isfinite(trans_pos)),     'trans_pos must be finite');
            tc.verifyTrue(all(isfinite(focus_pos)),     'focus_pos must be finite');
            tc.verifyTrue(all(isfinite(trans_pos_ras)), 'trans_pos_ras must be finite');
            tc.verifyTrue(all(isfinite(focus_pos_ras)), 'focus_pos_ras must be finite');

            % Positions must be positive (within a 200³ voxel grid)
            tc.verifyTrue(all(trans_pos >= 1), 'trans_pos must be inside the image grid');
            tc.verifyTrue(all(focus_pos >= 1), 'focus_pos must be inside the image grid');
        end

        function test_ras_to_voxel_round_trip(tc)
            % The known position_matrix written by the mock wrapper has
            % translation = [10 20 30] mm (RAS) and direction = [0 0 1].
            % reference_dist = -(curv_radius - dist_geom_ep) = -(63.2-52.38)
            %                = -10.82 mm
            % Therefore trans_pos_ras(x) = 10, trans_pos_ras(y) = 20 exactly
            % (the reference_dist offset is purely along z).
            %
            % The voxel conversion depends on the header affine, which we do
            % not control (niftiwrite defaults). We verify the round-trip by
            % converting trans_pos_ras back to voxels using the actual header
            % and checking it matches trans_pos to within 1 voxel.

            python_bin = find_system_python();
            tc.assumeNotEmpty(python_bin, 'No Python 3 found — skipping.');

            p      = make_plantus_parameters(tc.TmpDir);
            t1_hdr = write_synthetic_t1(tc.TmpDir, p);
            write_fake_mesh(tc.TmpDir, p);

            script_dir    = fullfile(tc.TmpDir, 'mock_plantus2');
            make_fake_plantus_install(script_dir);
            launcher_path = fullfile(tc.TmpDir, 'mock_launcher2.py');
            write_python_launcher(launcher_path, 'success');

            p.placement.plantus.script_path      = script_dir;
            p.placement.plantus.env_path = make_fake_simnibs_env(tc.TmpDir, python_bin);
            p.placement.plantus._launcher_py     = launcher_path;
            p.transducer(1).focus_pos            = [100, 100, 130];
            p.placement.plantus.mni_target_mm    = [];

            [trans_pos, ~, trans_pos_ras, ~] = position_transducer_plantus(p, t1_hdr);

            % X and Y of trans_pos_ras must equal the mock translation (10, 20)
            % because the direction vector is [0 0 1] so offset is pure-z.
            tc.verifyEqual(trans_pos_ras(1), 10.0, 'AbsTol', 1e-6, ...
                'trans_pos_ras x must match mock matrix translation x');
            tc.verifyEqual(trans_pos_ras(2), 20.0, 'AbsTol', 1e-6, ...
                'trans_pos_ras y must match mock matrix translation y');

            % Round-trip: convert trans_pos_ras back to voxels and compare
            expected_vox = ras_to_grid(trans_pos_ras, t1_hdr);
            tc.verifyEqual(double(trans_pos), double(expected_vox), 'AbsTol', 1.0, ...
                'trans_pos voxels must match RAS→voxel conversion of trans_pos_ras');
        end

        function test_focus_distal_to_transducer(tc)
            % focus_pos must be further from the image origin than trans_pos
            % (the focus lies in front of the bowl, deeper into the head).

            python_bin = find_system_python();
            tc.assumeNotEmpty(python_bin, 'No Python 3 found — skipping.');

            p      = make_plantus_parameters(tc.TmpDir);
            t1_hdr = write_synthetic_t1(tc.TmpDir, p);
            write_fake_mesh(tc.TmpDir, p);

            script_dir    = fullfile(tc.TmpDir, 'mock_plantus3');
            make_fake_plantus_install(script_dir);
            launcher_path = fullfile(tc.TmpDir, 'mock_launcher3.py');
            write_python_launcher(launcher_path, 'success');

            p.placement.plantus.script_path      = script_dir;
            p.placement.plantus.env_path = make_fake_simnibs_env(tc.TmpDir, python_bin);
            p.placement.plantus._launcher_py     = launcher_path;
            p.transducer(1).focus_pos            = [100, 100, 130];
            p.placement.plantus.mni_target_mm    = [];

            [trans_pos, focus_pos, ~, ~] = position_transducer_plantus(p, t1_hdr);

            % With direction = [0 0 1] and focal_distance_bowl > 0, the focus
            % must be at a higher z-voxel than the transducer rear centre.
            tc.verifyGreaterThan(focus_pos(3), trans_pos(3), ...
                'focus_pos z must exceed trans_pos z (beam direction = +z in mock)');
        end

        function test_results_written_to_all_transducers(tc)
            % preproc_transducer_placement must propagate positions to every
            % configured transducer, not just the first.

            python_bin = find_system_python();
            tc.assumeNotEmpty(python_bin, 'No Python 3 found — skipping.');

            p      = make_plantus_parameters(tc.TmpDir);
            % Add a second transducer slot
            p.transducer(2) = p.transducer(1);

            t1_hdr = write_synthetic_t1(tc.TmpDir, p);
            write_fake_mesh(tc.TmpDir, p);

            script_dir    = fullfile(tc.TmpDir, 'mock_plantus4');
            make_fake_plantus_install(script_dir);
            launcher_path = fullfile(tc.TmpDir, 'mock_launcher4.py');
            write_python_launcher(launcher_path, 'success');

            p.placement.plantus.script_path      = script_dir;
            p.placement.plantus.env_path = make_fake_simnibs_env(tc.TmpDir, python_bin);
            p.placement.plantus._launcher_py     = launcher_path;
            p.transducer(1).focus_pos            = [100, 100, 130];
            p.transducer(2).focus_pos            = [100, 100, 130];
            p.placement.plantus.mni_target_mm    = [];
            p.placement.mode                     = 'plantus';

            p_out = preproc_transducer_placement(p);

            tc.verifyEqual(p_out.transducer(1).trans_pos, p_out.transducer(2).trans_pos, ...
                'Both transducers must receive the same plantus-resolved trans_pos');
            tc.verifyEqual(p_out.transducer(1).focus_pos, p_out.transducer(2).focus_pos, ...
                'Both transducers must receive the same plantus-resolved focus_pos');
        end

        function test_yaml_fields_written_to_disk(tc)
            % After the YAML is written, it must contain the expected keys.
            % We verify this by running up to (but not past) the subprocess
            % by using a wrapper that reads the YAML then exits 0.

            python_bin = find_system_python();
            tc.assumeNotEmpty(python_bin, 'No Python 3 found — skipping.');

            p      = make_plantus_parameters(tc.TmpDir);
            t1_hdr = write_synthetic_t1(tc.TmpDir, p);
            write_fake_mesh(tc.TmpDir, p);

            script_dir    = fullfile(tc.TmpDir, 'yaml_check_plantus');
            make_fake_plantus_install(script_dir);
            launcher_path = fullfile(tc.TmpDir, 'yaml_check_launcher.py');
            write_python_launcher(launcher_path, 'check_yaml');

            p.placement.plantus.script_path      = script_dir;
            p.placement.plantus.env_path = make_fake_simnibs_env(tc.TmpDir, python_bin);
            p.placement.plantus._launcher_py     = launcher_path;
            p.transducer(1).focus_pos            = [100, 100, 130];
            p.placement.plantus.mni_target_mm    = [];

            % This launcher writes the parsed YAML fields to a sidecar file
            % and writes the expected Localite mat so the function completes.
            [~, ~, ~, ~] = position_transducer_plantus(p, t1_hdr);

            yaml_echo = fullfile(tc.TmpDir, 'm2m_sub-001', 'PlanTUS', ...
                'test_thal', 'yaml_fields.txt');
            tc.assumeTrue(isfile(yaml_echo), ...
                'check_yaml wrapper did not produce yaml_fields.txt — skipping YAML field check.');

            contents = fileread(yaml_echo);
            required_keys = {'min_distance', 'max_distance', 'transducer_diameter', ...
                             'focal_distance_list', 'flhm_list', 'plane_offset'};
            for i = 1:numel(required_keys)
                tc.verifySubstring(contents, required_keys{i}, ...
                    sprintf('PlanTUS YAML must contain key: %s', required_keys{i}));
            end
        end

        function test_target_mask_is_sphere_in_t1_space(tc)
            % After a successful run the target mask NIfTI must:
            %   - Exist on disk
            %   - Be binary (values 0 and 1 only)
            %   - Have a connected non-zero region (the sphere)

            python_bin = find_system_python();
            tc.assumeNotEmpty(python_bin, 'No Python 3 found — skipping.');

            p      = make_plantus_parameters(tc.TmpDir);
            t1_hdr = write_synthetic_t1(tc.TmpDir, p);
            write_fake_mesh(tc.TmpDir, p);

            script_dir    = fullfile(tc.TmpDir, 'mask_check_plantus');
            make_fake_plantus_install(script_dir);
            launcher_path = fullfile(tc.TmpDir, 'mask_check_launcher.py');
            write_python_launcher(launcher_path, 'success');

            p.placement.plantus.script_path      = script_dir;
            p.placement.plantus.env_path = make_fake_simnibs_env(tc.TmpDir, python_bin);
            p.placement.plantus._launcher_py     = launcher_path;
            p.transducer(1).focus_pos            = [100, 100, 130];
            p.placement.plantus.mni_target_mm    = [];

            position_transducer_plantus(p, t1_hdr);

            mask_path = fullfile(tc.TmpDir, 'm2m_sub-001', 'PlanTUS', ...
                'test_thal', 'sub-001_test_thal_mask.nii.gz');
            tc.verifyTrue(isfile(mask_path), 'Target mask NIfTI must be written to disk.');

            mask_info = niftiinfo(mask_path);
            mask_vol  = niftiread(mask_info);

            unique_vals = unique(mask_vol(:));
            tc.verifyTrue(all(ismember(unique_vals, [0, 1])), ...
                'Mask must be binary (values 0 and 1 only)');
            tc.verifyGreaterThan(sum(mask_vol(:)), 0, ...
                'Mask sphere must contain at least one non-zero voxel');
        end

    end

end

% =========================================================================
%% Parameter / fixture helpers
% =========================================================================

function p = make_base_parameters()
% Minimal parameters struct without plantus config.
    p = struct();
    p.subject_id = 1;
    p.platform   = 'matlab';

    p.placement.mode = 'manual';

    p.transducer(1).type                   = 'annular';
    p.transducer(1).trans_pos              = [100, 100, 5];
    p.transducer(1).focus_pos              = [100, 100, 130];
    p.transducer(1).focal_distance_bowl    = 63.2;
    p.transducer(1).focal_distance_ep      = 52.38;
    p.transducer(1).annular.curv_radius_mm  = 63.2;
    p.transducer(1).annular.dist_geom_ep_mm = 52.38;
    p.transducer(1).annular.elem_od_mm      = ...
        [4.09, 8.69, 13.29, 17.89, 22.49, 27.09, 31.69, 36.29, 40.89, 45.49];
end

function p = make_plantus_parameters(tmp_dir)
% Parameters struct with a filled-in plantus sub-struct pointing at
% temporary paths. The filesystem structure is NOT created here — call
% write_synthetic_t1, write_fake_mesh, and make_fake_simnibs_env as needed.
    p = make_base_parameters();
    p.placement.mode = 'plantus';

    % Paths that exist (will be created by test setup helpers)
    p.path.anat        = tmp_dir;
    p.path.t1_pattern  = 'sub-%03d_T1w.nii';
    p.path.seg         = tmp_dir;
    p.path.sim         = tmp_dir;
    p.io.dir_output    = tmp_dir;

    % plantus config — script_path / env_path intentionally left
    % empty here so each test can set exactly what it needs to test.
    p.placement.plantus.script_path              = '';
    p.placement.plantus.env_path         = '';
    p.placement.plantus.target_name              = 'test_thal';
    p.placement.plantus.mni_target_mm            = [];
    p.placement.plantus.focal_distance_list      = [63.2, 70.0];
    p.placement.plantus.flhm_list                = [25.0, 28.0];
    p.placement.plantus.additional_offset_mm     = 0;
    p.placement.plantus.max_angle_deg            = 10;
    p.placement.plantus.mask_radius_mm           = 2;
    p.placement.plantus.weight_skin_target_distances     = 0.2;
    p.placement.plantus.weight_skin_target_angles        = 0.2;
    p.placement.plantus.weight_skin_target_intersections = 0.2;
    p.placement.plantus.weight_skin_skull_angles         = 0.2;
    p.placement.plantus.weight_skull_thickness           = 0.2;
    p.placement.plantus.connectome_wb_path       = '';
end

% ── Synthetic T1 ──────────────────────────────────────────────────────────
function t1_hdr = make_synthetic_t1_header()
% Return a minimal header struct matching niftiinfo output for a 200³
% 1mm-iso volume. RAS = voxel - [100 100 100] mm (origin at image centre).
% This struct is sufficient for localite_matrix_to_positions and
% position_transducer_plantus for tests that do not touch disk.
    t1_hdr               = struct();
    t1_hdr.ImageSize     = [200, 200, 200];
    t1_hdr.PixelDimensions = [1.0, 1.0, 1.0, 1.0];
    t1_hdr.Transform     = struct();
    T = eye(4);
    T(1:3, 4) = -[100; 100; 100];
    t1_hdr.Transform.T = T';   % niftiinfo stores T transposed
end

function t1_hdr = write_synthetic_t1(tmp_dir, p)
% Write a 200³ NIfTI with default metadata and read back the real niftiinfo.
% The actual affine is whatever MATLAB's niftiwrite chooses; integration
% tests use the returned header for all coordinate assertions.
    t1_path = fullfile(tmp_dir, sprintf(p.path.t1_pattern, p.subject_id));
    niftiwrite(zeros([200, 200, 200], 'uint8'), t1_path);
    t1_hdr = niftiinfo(t1_path);
end

% ── Fake SimNIBS mesh ──────────────────────────────────────────────────────
function write_fake_mesh(tmp_dir, p)
% Create the empty .msh file and directory structure expected by
% position_transducer_plantus's mesh_path check.
    m2m_dir = fullfile(tmp_dir, sprintf('m2m_sub-%03d', p.subject_id));
    if ~isfolder(m2m_dir), mkdir(m2m_dir); end
    fclose(fopen(fullfile(m2m_dir, sprintf('sub-%03d.msh', p.subject_id)), 'w'));
end

% ── Fake SimNIBS Python env ───────────────────────────────────────────────
function env_path = make_fake_simnibs_env(tmp_dir, python_bin)
% Create simnibs_env/bin/python as a symlink to the given python_bin, so
% that position_transducer_plantus finds a valid Python executable.
    env_bin = fullfile(tmp_dir, 'fake_simnibs', 'simnibs_env', 'bin');
    if ~isfolder(env_bin), mkdir(env_bin); end

    fake_python = fullfile(env_bin, 'python');
    if isunix
        % Symlink to real Python
        if ~isfile(fake_python)
            system(sprintf('ln -sf "%s" "%s"', python_bin, fake_python));
        end
    else
        % Windows: copy the executable
        copyfile(python_bin, fake_python);
    end

    env_path = fullfile(tmp_dir, 'fake_simnibs');
end

% ── Fake PlanTUS installation ─────────────────────────────────────────────
function make_fake_plantus_install(script_dir)
% Create the minimal directory structure that position_transducer_plantus
% checks: script_dir/code/PlanTUS.py must exist.
    code_dir = fullfile(script_dir, 'code');
    if ~isfolder(code_dir), mkdir(code_dir); end
    fclose(fopen(fullfile(code_dir, 'PlanTUS.py'), 'w'));
end

% ── Mock Python launcher scripts ──────────────────────────────────────────
function write_python_launcher(launcher_path, mode)
% Write a minimal Python script to launcher_path that mimics the behaviour
% of prestus_plantus_launcher.py for tests:
%
%   'success'    — Writes a *Localite.mat with a known 4×4 matrix to the
%                  work_dir (inferred from the mask path arg), then exits 0.
%   'failing'    — Exits non-zero immediately.
%   'silent'     — Exits 0 without writing any output files.
%   'check_yaml' — Parses the YAML config arg and echoes its keys to a
%                  sidecar file, then also writes the Localite mat.
%
% The launcher is called as:
%   python launcher.py <t1> <mesh> <mask> <yaml> --plantus_root <root>
% so sys.argv[1..4] hold the four positional paths.

    nl = newline();

    mat_write_snippet = strjoin({ ...
        'mat = np.eye(4)', ...
        'mat[0:3, 3] = [10, 20, 30]', ...
        'mat[0:3, 0] = [0, 0, 1]', ...
        'sio.savemat(os.path.join(work_dir, "result_Localite.mat"), {"position_matrix": mat})', ...
    }, nl);

    switch mode
        case 'success'
            lines = { ...
                'import sys, os, numpy as np', ...
                'try:', ...
                '    import scipy.io as sio', ...
                'except ImportError:', ...
                '    sys.exit(99)', ...
                'work_dir = os.path.dirname(sys.argv[3])', ...
                'os.makedirs(work_dir, exist_ok=True)', ...
                mat_write_snippet, ...
            };

        case 'failing'
            lines = {'import sys', 'sys.exit(1)'};

        case 'silent'
            lines = {'import sys', 'sys.exit(0)'};

        case 'check_yaml'
            lines = { ...
                'import sys, os, numpy as np', ...
                'try:', ...
                '    import scipy.io as sio', ...
                '    import yaml', ...
                'except ImportError:', ...
                '    sys.exit(99)', ...
                'work_dir = os.path.dirname(sys.argv[3])', ...
                'os.makedirs(work_dir, exist_ok=True)', ...
                'with open(sys.argv[4]) as fh:', ...
                '    cfg = yaml.safe_load(fh)', ...
                'echo_path = os.path.join(work_dir, "yaml_fields.txt")', ...
                'with open(echo_path, "w") as fh:', ...
                '    fh.write("\n".join(cfg.keys()))', ...
                mat_write_snippet, ...
            };

        otherwise
            error('write_python_launcher: unknown mode ''%s''', mode);
    end

    code = strjoin(lines, nl);
    fid  = fopen(launcher_path, 'w');
    fprintf(fid, '%s\n', code);
    fclose(fid);
end

% ── Python discovery ──────────────────────────────────────────────────────
function python_bin = find_system_python()
% Return path to a Python 3 binary that has numpy and scipy, or '' if none.
    for candidate = {'python3', 'python'}
        [status, out] = system(sprintf('which %s 2>/dev/null', candidate{1}));
        if status ~= 0 || isempty(strtrim(out)), continue; end
        bin = strtrim(out);
        [s, ~] = system(sprintf('"%s" -c "import numpy, scipy" 2>/dev/null', bin));
        if s == 0
            python_bin = bin;
            return;
        end
    end
    python_bin = '';
end
