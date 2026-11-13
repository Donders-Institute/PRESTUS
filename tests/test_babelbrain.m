classdef (SharedTestFixtures = {prestus_paths}, TestTags = {'standalone'}) test_babelbrain < matlab.unittest.TestCase
% TEST_BABELBRAIN  Staged integration tests for the PRESTUS–BabelBrain backend.
%
% Covers five progressively deeper stages:
%
%   Stage 1 (paths)      — PRESTUS functions on path, Python accessible,
%                          BabelBrain submodule present.
%   Stage 2 (write)      — write_babelbrain_input produces a valid NIfTI with
%                          correct label counts and focal voxel.
%   Stage 3 (transducer) — babelbrain_transducer_map resolves a known preset
%                          and generates a custom CSV for a novel transducer.
%   Stage 4 (dry_run)    — Python bridge executes with dry_run=true; verifies
%                          BabelBrain imports and writes the input NIfTI.
%   Stage 5 (full)       — Full GPU simulation; checks p_max_all is present
%                          and physically plausible.
%
% Each stage skips gracefully when its prerequisites are absent.
%
% Environment variables
% ---------------------
%   BABELBRAIN_PYTHON   Path to the Python executable  OR  a conda/venv
%                       environment root directory.  Falls back to 'python3'.
%                       Examples:
%                         export BABELBRAIN_PYTHON=/opt/miniconda3/envs/babelbrain
%                         export BABELBRAIN_PYTHON=/usr/local/bin/python3
%   PRESTUS_DEMO_DATA   Root of the prestus_testdata directory (must contain
%                       sim_outputs/sub-009_tmp2/… and configs/).
%
% Run with:
%   results = runtests('tests/test_babelbrain.m');
%   table(results)
%
% Stages by tag:
%   'paths'      — Stage 1 only
%   'write'      — Stage 2
%   'transducer' — Stage 3
%   'dry_run'    — Stage 4  (requires BabelBrain importable)
%   'full'       — Stage 5  (requires GPU)
%   'quick'      — Stages 1–3
%   'nogpu'      — Stages 1–4

    % ------------------------------------------------------------------ %
    properties
        SimPath     % temporary output directory for this test run
        DemoRoot    % path to prestus_testdata/
        PyBin       % resolved Python executable path
        BbRoot      % path to BabelBrain submodule root
        MediumMasks % medium_masks volume loaded in Stage 2, reused later
        FocusPos    % focus position in grid space [x y z]
        TransPos    % transducer position in grid space [x y z]
        SubjectId   % subject number (integer)
    end

    % ================================================================== %
    %% CLASS SETUP / TEARDOWN
    % ================================================================== %
    methods (TestClassSetup)

        function classSetup(tc)
            % Require PRESTUS_DEMO_DATA — skip the whole class if absent.
            demo_root = getenv('PRESTUS_DEMO_DATA');
            tc.assumeNotEmpty(demo_root, ...
                'Set PRESTUS_DEMO_DATA to the prestus_testdata directory to run BabelBrain tests.');
            tc.assumeTrue(isfolder(demo_root), ...
                sprintf('PRESTUS_DEMO_DATA does not exist: %s', demo_root));
            tc.DemoRoot  = demo_root;
            tc.SubjectId = 9;
            tc.FocusPos  = [112 104 163];
            tc.TransPos  = [12  101 166];

            % Resolve Python binary — accepts executable path or env root.
            [py_bin, err] = test_babelbrain.resolve_python_bin(getenv('BABELBRAIN_PYTHON'));
            tc.assumeEmpty(err, err);
            tc.PyBin = py_bin;

            % Locate BabelBrain submodule.
            prestus_root = fileparts(fileparts(mfilename('fullpath')));
            bb_candidate = fullfile(prestus_root, 'external', 'BabelBrain');
            if isfolder(fullfile(bb_candidate, 'TranscranialModeling'))
                tc.BbRoot = bb_candidate;
            else
                tc.BbRoot = '';
            end

            % Temporary sim output directory.
            tc.SimPath = fullfile(tempdir, ...
                sprintf('prestus_bb_test_%s', datestr(now, 'yyyymmddHHMMSS'))); %#ok<TNOW1,DATST>
            mkdir(tc.SimPath);
        end

    end

    methods (TestClassTeardown)

        function classTeardown(tc)
            % SimPath is unset when classSetup was skipped (no demo data).
            if ~isempty(tc.SimPath) && isfolder(tc.SimPath)
                rmdir(tc.SimPath, 's');
            end
        end

    end

    % ================================================================== %
    %% STAGE 1 — PATH AND IMPORT CHECKS
    % ================================================================== %
    methods (Test, TestTags = {'paths', 'quick', 'nogpu'})

        function test_prestus_functions_on_path(tc)
            required = {'write_babelbrain_input', 'babelbrain_transducer_map', ...
                        'read_babelbrain_output', 'babelbrain_wrapper', ...
                        'load_parameters'};
            for k = 1:numel(required)
                tc.verifyNotEmpty(which(required{k}), ...
                    sprintf('PRESTUS function not found on path: %s', required{k}));
            end
        end

        function test_babelbrain_submodule_present(tc)
            tc.verifyTrue(~isempty(tc.BbRoot) && isfolder(tc.BbRoot), ...
                'BabelBrain submodule not found. Run: git submodule update --init external/BabelBrain');
        end

        function test_python_executable_found(tc)
            [st, ver] = system(sprintf('"%s" --version 2>&1', tc.PyBin));
            tc.verifyEqual(st, 0, ...
                sprintf('Python executable not usable (%s): %s', tc.PyBin, strtrim(ver)));
            fprintf('[test_babelbrain] Python: %s  (%s)\n', tc.PyBin, strtrim(ver));
        end

        function test_runner_script_present(tc)
            prestus_root = fileparts(fileparts(mfilename('fullpath')));
            runner = fullfile(prestus_root, 'functions', 'babelbrain', 'run_babelbrain_case.py');
            tc.verifyTrue(isfile(runner), ...
                sprintf('Python runner not found: %s', runner));
        end

    end

    % ================================================================== %
    %% STAGE 2 — write_babelbrain_input
    % ================================================================== %
    methods (Test, TestTags = {'write', 'quick', 'nogpu'})

        function test_write_babelbrain_input_creates_nifti(tc)
            [masks, p] = tc.load_masks_and_params();
            tc.assumeNotEmpty(masks, 'medium_masks not available — run demo pipeline first.');

            out_nii = write_babelbrain_input(masks, tc.FocusPos, p, ...
                tc.make_planimg(), tc.SimPath, 'M1_left');

            tc.verifyTrue(isfile(out_nii), ...
                sprintf('write_babelbrain_input did not create: %s', out_nii));
        end

        function test_write_babelbrain_input_label_mapping(tc)
            [masks, p] = tc.load_masks_and_params();
            tc.assumeNotEmpty(masks, 'medium_masks not available.');

            out_nii = write_babelbrain_input(masks, tc.FocusPos, p, ...
                tc.make_planimg(), tc.SimPath, 'label_check');

            vol = niftiread(out_nii);
            tc.verifyGreaterThanOrEqual(max(vol(:)), 4, ...
                'No brain voxels (label >= 4) found — label mapping likely failed.');
            tc.verifyTrue(any(vol(:) == 5), ...
                'Focal voxel (label 5) missing from output.');
        end

        function test_write_babelbrain_input_resolution(tc)
            [masks, p] = tc.load_masks_and_params();
            tc.assumeNotEmpty(masks, 'medium_masks not available.');

            out_nii = write_babelbrain_input(masks, tc.FocusPos, p, ...
                tc.make_planimg(), tc.SimPath, 'res_check');

            % The mask is resampled to BabelBrain's own native step
            % (SmallestSOS / freq / ppw), not written at grid.resolution_mm,
            % because BabelBrain sizes its domain from the voxel count alone.
            % Mirror write_babelbrain_input's auto-PPW (ppw_override = 0).
            smallest_sos = 1000.0 * 0.422 + 680.515;   % Pichardo shear floor
            mp_names = fieldnames(p.medium_properties);
            for mi = 1:numel(mp_names)
                mp = p.medium_properties.(mp_names{mi});
                if isfield(mp, 'sound_speed') && mp.sound_speed > 0
                    smallest_sos = min(smallest_sos, mp.sound_speed);
                end
            end
            freq_hz = p.transducer(1).freq_hz;
            ppw = round(smallest_sos / (freq_hz * p.grid.resolution_mm * 1e-3));
            expected_res_mm = (smallest_sos / freq_hz) / ppw * 1e3;

            info = niftiinfo(out_nii);
            tc.verifyEqual(info.PixelDimensions(1), expected_res_mm, 'RelTol', 1e-4, ...
                'NIfTI pixel dimensions do not match BabelBrain''s native resolution.');
        end

    end

    % ================================================================== %
    %% STAGE 3 — babelbrain_transducer_map
    % ================================================================== %
    methods (Test, TestTags = {'transducer', 'quick', 'nogpu'})

        function test_known_preset_resolved(tc)
            p = tc.make_r15473_params();
            [preset, csv_path] = babelbrain_transducer_map(p, struct(), tc.SimPath);
            tc.verifyEqual(preset, 'R15473', ...
                sprintf('Expected preset R15473, got ''%s''.', preset));
            tc.verifyEmpty(csv_path, ...
                'No custom CSV expected for a known preset.');
        end

        function test_novel_transducer_generates_csv(tc)
            p = tc.make_novel_transducer_params();
            [preset, csv_path] = babelbrain_transducer_map(p, struct(), tc.SimPath);
            tc.verifyEmpty(preset, ...
                'Novel transducer should not match any preset.');
            tc.verifyTrue(~isempty(csv_path) && isfile(csv_path), ...
                'Custom element CSV should have been generated for novel transducer.');
        end

    end

    % ================================================================== %
    %% STAGE 4 — Python bridge dry run
    % ================================================================== %
    methods (Test, TestTags = {'dry_run', 'nogpu'})

        function test_babelbrain_importable(tc)
            [st, ~] = system(sprintf('"%s" -c "import BabelViscoFDTD" 2>/dev/null', tc.PyBin));
            tc.assumeEqual(st, 0, ...
                'BabelViscoFDTD not importable — set BABELBRAIN_PYTHON to the env that has it installed.');
        end

        function test_dry_run_returns_struct(tc)
            tc.skip_if_babelbrain_missing();
            [masks, p] = tc.load_masks_and_params();
            tc.assumeNotEmpty(masks, 'medium_masks not available.');

            p = tc.inject_babelbrain_config(p, true);
            kw_stub = tc.make_kwave_stub(masks);

            sd = babelbrain_wrapper(p, masks, kw_stub, tc.FocusPos, tc.TransPos, tc.make_planimg());
            tc.verifyTrue(isstruct(sd), 'babelbrain_wrapper must return a struct.');
        end

        function test_dry_run_writes_input_nifti(tc)
            tc.skip_if_babelbrain_missing();
            [masks, p] = tc.load_masks_and_params();
            tc.assumeNotEmpty(masks, 'medium_masks not available.');

            p = tc.inject_babelbrain_config(p, true);
            kw_stub = tc.make_kwave_stub(masks);

            babelbrain_wrapper(p, masks, kw_stub, tc.FocusPos, tc.TransPos, tc.make_planimg());

            freq_kHz = round(p.transducer(1).freq_hz / 1e3);
            c_min    = 1500;
            ppw      = max(6, round(c_min / (p.transducer(1).freq_hz * p.grid.resolution_mm * 1e-3)));
            subj_str = sprintf('sub-%03d', tc.SubjectId);
            target   = strrep(p.io.output_affix, '_', '');
            expected = fullfile(tc.SimPath, subj_str, ...
                sprintf('%s_%dkHz_%dPPW_BabelViscoInput.nii.gz', target, freq_kHz, ppw));
            tc.verifyTrue(isfile(expected), ...
                sprintf('Expected BabelViscoInput NIfTI not found:\n  %s', expected));
        end

    end

    % ================================================================== %
    %% STAGE 5 — Full simulation (GPU required)
    % ================================================================== %
    methods (Test, TestTags = {'full'})

        function test_full_simulation_runs(tc)
            tc.skip_if_babelbrain_missing();
            [masks, p] = tc.load_masks_and_params();
            tc.assumeNotEmpty(masks, 'medium_masks not available.');

            p = tc.inject_babelbrain_config(p, false);
            kw_stub = tc.make_kwave_stub(masks);

            t0 = tic;
            sd = babelbrain_wrapper(p, masks, kw_stub, tc.FocusPos, tc.TransPos, tc.make_planimg());
            elapsed = toc(t0);

            tc.verifyTrue(isfield(sd, 'p_max_all'), 'sensor_data missing p_max_all.');
            tc.verifyNotEmpty(sd.p_max_all, 'p_max_all is empty.');
            p_max = max(sd.p_max_all(:));
            fprintf('[test_babelbrain] Full sim: %.1f s  peak=%.1f kPa\n', elapsed, p_max/1e3);
            tc.verifyGreaterThan(p_max, 1e3, ...
                sprintf('Peak pressure %.1f Pa is suspiciously low.', p_max));
        end

    end

    % ================================================================== %
    %% STATIC HELPERS
    % ================================================================== %
    methods (Static)

        function [py_bin, err] = resolve_python_bin(env_val)
        % RESOLVE_PYTHON_BIN  Resolve a Python path spec to an executable.
        %
        %   env_val can be:
        %     ''          — fall back to 'python3' on PATH
        %     '/path/to/python3'  — direct executable (file must exist)
        %     '/path/to/env'      — conda or venv root; tries bin/python3,
        %                           bin/python, then python.exe (Windows)
        %
        %   If the path cannot be resolved, the user is prompted interactively.
        %   Returns py_bin (char) and err (empty on success, message on failure).
            err = '';
            if isempty(env_val)
                py_bin = 'python3';
                return
            end
            [py_bin, err] = test_babelbrain.try_resolve(env_val, 'BABELBRAIN_PYTHON');
            if ~isempty(err)
                fprintf('[test_babelbrain] %s\n', err);
                py_bin = strtrim(input( ...
                    'Enter path to Python executable or conda/venv root: ', 's'));
                if isempty(py_bin)
                    err = 'No Python path provided — cannot continue.';
                    return
                end
                [py_bin, err] = test_babelbrain.try_resolve(py_bin, 'user input');
            end
        end

        function [py_bin, err] = try_resolve(env_val, var_name)
        % Try to resolve env_val to an executable path; return err on failure.
            err = '';
            if ~isfolder(env_val)
                if ~isfile(env_val)
                    py_bin = '';
                    err = sprintf('%s is not a file or directory: %s', var_name, env_val);
                    return
                end
                py_bin = env_val;
                return
            end
            candidates = { fullfile(env_val, 'bin', 'python3'), ...
                           fullfile(env_val, 'bin', 'python'), ...
                           fullfile(env_val, 'python.exe') };
            py_bin = '';
            for k = 1:numel(candidates)
                if isfile(candidates{k})
                    py_bin = candidates{k};
                    return
                end
            end
            err = sprintf(['%s points to a directory but no Python executable was found inside:\n' ...
                           '  %s\nSearched: bin/python3, bin/python, python.exe'], var_name, env_val);
        end

    end

    % ================================================================== %
    %% PRIVATE HELPERS
    % ================================================================== %
    methods (Access = private)

        function [masks, p] = load_masks_and_params(tc)
            masks_nii = fullfile(tc.DemoRoot, 'sim_outputs', 'sub-009_tmp2', 'nii', ...
                'sub-009_layered_T1w_bm_limFOV_d60_l120_medium_masks.nii.gz');
            if ~isfile(masks_nii)
                masks = [];
                p     = struct();
                return
            end
            masks = niftiread(masks_nii);

            p.subject_id                       = tc.SubjectId;
            p.grid.resolution_mm               = 2.0;
            p.transducer(1).freq_hz            = 300000;
            p.medium_properties.brain.sound_speed = 1540;
            p.io.output_affix                  = '_bm_bb';
            p.layers.water                     = [0 3 6 9 10];
            p.layers.brain                     = [1 2];
            p.layers.skin                      = [5];
            p.layers.skull                     = [4];
            p.layers.skull_cortical            = [7];
            p.layers.skull_trabecular          = [8];
        end

        function p = make_r15473_params(~)
            p.transducer(1).freq_hz                  = 300000;
            p.transducer(1).annular.elem_n            = 10;
            p.transducer(1).annular.curv_radius_mm    = 100.0;
            p.transducer(1).annular.elem_od_mm        = [21.1 28.8 35.0 40.4 45.3 49.7 53.9 57.8 61.5 65.0];
            p.transducer(1).annular.elem_id_mm        = [10.0 22.1 29.8 36.0 41.4 46.3 50.7 54.9 58.7 62.4];
        end

        function p = make_novel_transducer_params(~)
            p.transducer(1).freq_hz                  = 400000;
            p.transducer(1).annular.elem_n            = 3;
            p.transducer(1).annular.curv_radius_mm    = 55.0;
            p.transducer(1).annular.elem_od_mm        = [10 20 30];
            p.transducer(1).annular.elem_id_mm        = [0 11 21];
        end

        function planimg = make_planimg(~)
            planimg.t1_header    = [];
            planimg.transf       = [];
            planimg.inv_transf   = [];
            planimg.origin_ras_mm = [0 0 0];
        end

        function kw_stub = make_kwave_stub(~, masks)
            kw_stub.sound_speed = ones(size(masks)) * 1500;
        end

        function p = inject_babelbrain_config(tc, p, dry_run)
            prestus_root = fileparts(fileparts(mfilename('fullpath')));
            p.acoustic_backend                      = 'babelbrain';
            p.babelbrain.tx_system                  = 'R15473';
            p.babelbrain.computing_backend           = 1;   % CUDA
            p.babelbrain.ppw                        = 6;
            p.babelbrain.do_refocusing              = false;
            p.babelbrain.z_beyond_focus_mm          = 40;
            p.babelbrain.python_env                 = tc.PyBin;
            p.babelbrain.root                       = fullfile(prestus_root, 'external', 'BabelBrain');
            p.babelbrain.dry_run                    = dry_run;
            p.io.dir_cache                          = tc.SimPath;
            p.transducer(1).annular.elem_phase_rad  = ...
                [0.003 352.671 345.725 25.437 46.156 36.306 33.659 32.695 60.623 336.098] * pi/180;
            p.transducer(1).annular.elem_p_amp        = repmat(156587.0, 1, 10);
        end

        function skip_if_babelbrain_missing(tc)
            [st, ~] = system(sprintf('"%s" -c "import BabelViscoFDTD" 2>/dev/null', tc.PyBin));
            tc.assumeEqual(st, 0, ...
                'BabelViscoFDTD not importable — set BABELBRAIN_PYTHON to the env that has it installed.');
        end

    end

end
