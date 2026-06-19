classdef test_jwave < matlab.unittest.TestCase
% TEST_JWAVE  Compare jwave Helmholtz and k-Wave acoustic backends via the PRESTUS pipeline.
%
% Runs the same water simulation through both backends on a small grid and
% reports: wall-clock time, peak pressure at focus, focal FWHM, and spatial
% correlation of the pressure field. Fails tests only on hard correctness
% criteria (positive pressure, non-empty output, plausible correlation with
% k-Wave). Timing and quantitative differences are logged but not enforced.
%
% jwave uses the Helmholtz solver (frequency-domain GMRES with Fourier spectral
% collocation), matching the solver used in the ITRUSST intercomparison benchmarks.
%
% Environment variables
% ---------------------
%   JWAVE_PYTHON        Path to the Python executable  OR  a conda/venv
%                       environment root for CPU jwave runs.
%                       Falls back to system 'python3' if unset.
%                       Examples:
%                         export JWAVE_PYTHON=/opt/miniconda3/envs/jwave_cpu
%                         export JWAVE_PYTHON=/usr/local/bin/python3
%   JWAVE_PYTHON_GPU    Path to the Python executable  OR  a conda/venv
%                       environment root for GPU jwave runs (optional).
%                       If unset, GPU tests are skipped or fall back to CPU.
%                       Examples:
%                         export JWAVE_PYTHON_GPU=/opt/miniconda3/envs/jwave_gpu
%                         export JWAVE_PYTHON_GPU=/usr/local/cuda-env/bin/python3
%   JWAVE_USE_GPU       'true' | 'false' — override GPU auto-detection.
%                       Leave unset for automatic nvidia-smi probe.
%
% Run with:
%   results = runtests('tests/test_jwave.m');
%   table(results)              % summary
%
% To run only the comparison tests without the slow full-pipeline sweeps:
%   results = runtests('tests/test_jwave.m', 'Tag', 'quick');
%
% To run timing sweeps across grid sizes:
%   results = runtests('tests/test_jwave.m', 'Tag', 'timing');

    % ------------------------------------------------------------------ %
    properties
        SimPath     % root output directory for this test run
        RefData     % struct: kwave reference sensor_data, kgrid, parameters
        PyBinCpu    % resolved CPU Python executable path
        PyBinGpu    % resolved GPU Python executable path (may be empty)
    end

    properties (Constant)
        % Minimum acceptable Pearson r between jwave and k-Wave p_max_all.
        MIN_CORRELATION = 0.95;

        % Maximum acceptable relative error on peak pressure.
        REL_TOL = 0.05;

        % Grid sizes (voxels per side, isotropic) used in timing sweep.
        TIMING_GRID_SIZES = [32, 48, 64, 96];
    end

    % ------------------------------------------------------------------ %
    methods (TestClassSetup)
        function classSetup(tc)
            % Ensure fixtures directory is on the path when the test is run directly.
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(here, 'fixtures'));

            % Resolve CPU Python env — skip the whole class if jwave is not importable.
            [cpu_bin, err] = test_jwave.resolve_python_bin(getenv('JWAVE_PYTHON'), 'JWAVE_PYTHON');
            tc.assumeEmpty(err, err);
            [status, ~] = system(sprintf('"%s" -c "import jwave" 2>/dev/null', cpu_bin));
            tc.assumeEqual(status, 0, ...
                sprintf(['jwave not importable in: %s\n' ...
                         'Set JWAVE_PYTHON to the Python executable or conda/venv root ' ...
                         'that has jwave installed.'], cpu_bin));
            tc.PyBinCpu = cpu_bin;

            % Resolve GPU Python env — warn but do not skip if misconfigured.
            gpu_env_val = getenv('JWAVE_PYTHON_GPU');
            if ~isempty(gpu_env_val)
                [gpu_bin, gpu_err] = test_jwave.resolve_python_bin(gpu_env_val, 'JWAVE_PYTHON_GPU');
                if ~isempty(gpu_err)
                    warning('test_jwave:gpuEnvBadPath', '%s', gpu_err);
                    tc.PyBinGpu = '';
                else
                    [gpu_status, ~] = system(sprintf('"%s" -c "import jwave" 2>/dev/null', gpu_bin));
                    if gpu_status ~= 0
                        warning('test_jwave:gpuEnvMissing', ...
                            'JWAVE_PYTHON_GPU is set but jwave is not importable there: %s', gpu_bin);
                    end
                    tc.PyBinGpu = gpu_bin;
                end
            else
                tc.PyBinGpu = '';
            end

            tc.SimPath = fullfile(tempdir, ...
                sprintf('prestus_jwave_test_%s', datestr(now, 'yyyymmddHHMMSS'))); %#ok<TNOW1,DATST>
            mkdir(tc.SimPath);

            % RefData is populated lazily on first access via get_ref_data().
            % Tests that do not compare against k-Wave skip the reference sim.
            tc.RefData = [];
        end
    end

    methods (TestClassTeardown)
        function classTeardown(tc)
            if exist(tc.SimPath, 'dir')
                rmdir(tc.SimPath, 's');
            end
        end
    end

    % ================================================================== %
    %% CORRECTNESS TESTS
    % ================================================================== %
    methods (Test, TestTags = {'quick', 'correctness'})

        % -------------------------------------------------------------- %
        function test_jwave_runs_without_error(tc)
        % jwave Helmholtz backend should complete without throwing.
            p = tc.make_jwave_params('helm_basic', [48 48 80]);
            tc.verifyWarningFree(@() prestus_pipeline(p));
        end

        % -------------------------------------------------------------- %
        function test_jwave_returns_positive_pressure(tc)
        % Peak pressure field should be non-empty and contain positive values.
            p   = tc.make_jwave_params('helm_pos', [48 48 80]);
            out = prestus_pipeline(p);
            sd  = out.sensor_data;
            tc.verifyNotEmpty(sd.p_max_all,           'p_max_all is empty');
            tc.verifyGreaterThan(max(sd.p_max_all(:)), 0, 'p_max_all has no positive values');
        end

        % -------------------------------------------------------------- %
        function test_jwave_output_fields_present(tc)
        % jwave_wrapper must return the same field names as acoustic_simulation.
            p   = tc.make_jwave_params('helm_fields', [48 48 80]);
            out = prestus_pipeline(p);
            sd  = out.sensor_data;
            for f = {'p_max_all', 'p_min_all', 'p_final'}
                tc.verifyTrue(isfield(sd, f{1}), sprintf('Missing field: %s', f{1}));
            end
        end

        % -------------------------------------------------------------- %
        function test_jwave_p_min_negative(tc)
        % Peak negative pressure should be strictly negative (rarefaction).
            p   = tc.make_jwave_params('helm_pmin', [48 48 80]);
            out = prestus_pipeline(p);
            sd  = out.sensor_data;
            tc.verifyLessThan(min(sd.p_min_all(:)), 0, 'p_min_all has no negative values');
        end

        % -------------------------------------------------------------- %
        function test_axisymmetric_raises_error(tc)
        % jwave backend must refuse axisymmetric grids with a clear message.
            p = tc.make_jwave_params('helm_as', [1 48 80]);
            p.grid.axisymmetric = 1;
            tc.verifyError(@() prestus_pipeline(p), '');
        end

        % -------------------------------------------------------------- %
        function test_kWaveArray_overridden_with_warning(tc)
        % load_parameters should silently override use_kWaveArray=1 to 0 when
        % the jwave backend is selected, emitting a warning rather than erroring.
            p = tc.make_jwave_params('helm_kwa', [48 48 80]);
            p.grid.use_kWaveArray = 1;
            tc.verifyWarning( ...
                @() load_parameters(p), ...
                'prestus:jwaveKWaveArrayUnsupported');
            p_fixed = load_parameters(p);
            tc.verifyEqual(p_fixed.grid.use_kWaveArray, 0, ...
                'use_kWaveArray should be overridden to 0 for jwave backend');
        end

        % -------------------------------------------------------------- %
        function test_dry_run_skips_simulation(tc)
        % Dry-run flag should write HDF5 input but return empty pressure.
            p = tc.make_jwave_params('helm_dryrun', [48 48 80]);
            p.jwave.dry_run = true;
            out = prestus_pipeline(p);
            sd  = out.sensor_data;
            tc.verifyEmpty(sd.p_max_all, 'Dry-run should return empty p_max_all');
            tc.verifyTrue(isfile(sd.jwave_input_h5), 'Input HDF5 should exist after dry run');
        end

    end

    % ================================================================== %
    %% ENVIRONMENT SELECTION TESTS
    % ================================================================== %
    methods (Test, TestTags = {'quick', 'env'})

        % -------------------------------------------------------------- %
        function test_force_cpu_uses_cpu_env(tc)
        % use_gpu = false must select python_env regardless of GPU presence.
            p = tc.make_jwave_params('env_cpu', [48 48 80]);
            p.jwave.use_gpu = false;
            p.jwave.dry_run = true;
            tc.verifyWarningFree(@() prestus_pipeline(p));
        end

        % -------------------------------------------------------------- %
        function test_force_gpu_errors_without_gpu_env(tc)
        % use_gpu = true with python_env_gpu unset must throw.
            p = tc.make_jwave_params('env_gpu_noenv', [48 48 80]);
            p.jwave.use_gpu        = true;
            p.jwave.python_env_gpu = '';
            p.jwave.dry_run        = true;
            tc.verifyError(@() prestus_pipeline(p), '');
        end

        % -------------------------------------------------------------- %
        function test_force_gpu_uses_gpu_env(tc)
        % use_gpu = true with python_env_gpu set must select the GPU env.
        % Skipped when JWAVE_PYTHON_GPU is not configured.
            tc.assumeNotEmpty(tc.PyBinGpu, 'JWAVE_PYTHON_GPU not set — skipping GPU-env test.');
            p = tc.make_jwave_params('env_gpu_ok', [48 48 80]);
            p.jwave.use_gpu        = true;
            p.jwave.python_env_gpu = tc.PyBinGpu;
            p.jwave.dry_run        = true;
            tc.verifyWarningFree(@() prestus_pipeline(p));
        end

        % -------------------------------------------------------------- %
        function test_autodetect_warns_when_gpu_found_but_no_gpu_env(tc)
        % Auto-detect: if a GPU is present but python_env_gpu is empty,
        % jwave_wrapper should issue a warning and fall back to CPU env.
            [smi_status, ~] = system('nvidia-smi -L 2>/dev/null');
            tc.assumeEqual(smi_status, 0, 'No GPU detected — skipping this test.');
            p = tc.make_jwave_params('env_auto_nogpuenv', [48 48 80]);
            p.jwave.use_gpu        = [];   % auto
            p.jwave.python_env_gpu = '';
            p.jwave.dry_run        = true;
            tc.verifyWarning(@() prestus_pipeline(p), 'jwave_wrapper:noGpuEnv');
        end

    end

    % ================================================================== %
    %% QUANTITATIVE AGREEMENT WITH K-WAVE
    % ================================================================== %
    methods (Test, TestTags = {'quick', 'agreement'})

        % -------------------------------------------------------------- %
        function test_jwave_peak_pressure_plausible(tc)
        % jwave Helmholtz peak pressure should be within REL_TOL of k-Wave.
            ref_peak = max(tc.get_ref_data().sensor_data.p_max_all(:));

            p   = tc.make_jwave_params('helm_agree', [48 48 80]);
            out = prestus_pipeline(p);
            jw_peak = max(out.sensor_data.p_max_all(:));

            rel_err = abs(jw_peak - ref_peak) / ref_peak;
            fprintf('[test_jwave] Peak pressure: kwave=%.1f Pa  jwave_helm=%.1f Pa  rel_err=%.1f%%\n', ...
                ref_peak, jw_peak, rel_err * 100);

            tc.verifyLessThan(rel_err, tc.REL_TOL, ...
                sprintf('Helmholtz peak pressure relative error %.1f%% exceeds threshold', rel_err * 100));
        end

        % -------------------------------------------------------------- %
        function test_jwave_spatial_correlation(tc)
        % Pearson r between jwave and k-Wave p_max_all should exceed MIN_CORRELATION.
            ref = tc.get_ref_data().sensor_data.p_max_all(:);
            p   = tc.make_jwave_params('helm_corr', [48 48 80]);
            out = prestus_pipeline(p);
            jw  = out.sensor_data.p_max_all(:);

            r = corr(double(ref), double(jw));
            fprintf('[test_jwave] Spatial correlation (p_max_all): r = %.4f\n', r);

            tc.verifyGreaterThan(r, tc.MIN_CORRELATION, ...
                sprintf('Spatial correlation %.4f below threshold %.2f', r, tc.MIN_CORRELATION));
        end

        % -------------------------------------------------------------- %
        function test_jwave_focal_depth_within_1voxel(tc)
        % The focal depth (axial index of peak pressure) should agree with
        % k-Wave to within 1 voxel.
            ref_sd  = tc.get_ref_data().sensor_data;
            ref_axial = squeeze(max(max(ref_sd.p_max_all, [], 1), [], 2));
            [~, ref_z] = max(ref_axial);

            p   = tc.make_jwave_params('helm_focal', [48 48 80]);
            out = prestus_pipeline(p);
            jw_axial = squeeze(max(max(out.sensor_data.p_max_all, [], 1), [], 2));
            [~, jw_z] = max(jw_axial);

            fprintf('[test_jwave] Focal depth: kwave=%d vox  jwave_helm=%d vox\n', ref_z, jw_z);
            tc.verifyLessThanOrEqual(abs(ref_z - jw_z), 1, ...
                'Focal depth differs by more than 1 voxel from k-Wave');
        end

    end

    % ================================================================== %
    %% TIMING TESTS  (tagged 'timing' — not run in quick CI)
    % ================================================================== %
    methods (Test, TestTags = {'timing'})

        % -------------------------------------------------------------- %
        function test_timing_sweep(tc)
        % Run kwave_cpu and jwave_helm across TIMING_GRID_SIZES and print a table.
        % No assertions — informational only.
            sizes    = tc.TIMING_GRID_SIZES;
            backends = {'kwave_cpu', 'jwave_helm'};
            results  = nan(numel(sizes), numel(backends));

            for gi = 1:numel(sizes)
                N = sizes(gi);
                grid_dims = [N, N, round(N * 1.6)];

                for bi = 1:numel(backends)
                    be = backends{bi};
                    try
                        switch be
                            case 'kwave_cpu'
                                p = tc.make_params(sprintf('timing_kw_%d', N), grid_dims);
                            case 'jwave_helm'
                                p = tc.make_jwave_params(sprintf('timing_helm_%d', N), grid_dims);
                        end
                        t0 = tic;
                        prestus_pipeline(p);
                        results(gi, bi) = toc(t0);
                    catch ME
                        fprintf('[test_jwave] %s N=%d FAILED: %s\n', be, N, ME.message);
                    end
                end
            end

            fprintf('\n=== jwave vs k-Wave timing (wall-clock seconds) ===\n');
            fprintf('%-12s', 'Grid (N)');
            for bi = 1:numel(backends)
                fprintf('  %-14s', backends{bi});
            end
            fprintf('\n');
            for gi = 1:numel(sizes)
                fprintf('%-12d', sizes(gi));
                for bi = 1:numel(backends)
                    if isnan(results(gi, bi))
                        fprintf('  %-14s', 'FAILED');
                    else
                        fprintf('  %-14.1f', results(gi, bi));
                    end
                end
                fprintf('\n');
            end
            fprintf('\n');
        end

    end

    % ================================================================== %
    %% STATIC HELPERS
    % ================================================================== %
    methods (Static)

        function [py_bin, err] = resolve_python_bin(env_val, var_name)
        % RESOLVE_PYTHON_BIN  Resolve a Python path spec to an executable.
        %
        %   env_val can be:
        %     ''          — fall back to 'python3' on PATH
        %     '/path/to/python3'  — direct executable (file must exist)
        %     '/path/to/env'      — conda or venv root; tries bin/python3,
        %                           bin/python, then python.exe (Windows)
        %
        %   If the path cannot be resolved, the user is prompted interactively.
        %   var_name is the env-var name used in prompts and error messages.
        %   Returns py_bin (char) and err (empty on success, message on failure).
            if nargin < 2
                var_name = 'JWAVE_PYTHON';
            end
            err = '';
            if isempty(env_val)
                py_bin = 'python3';
                return
            end
            [py_bin, err] = test_jwave.try_resolve(env_val, var_name);
            if ~isempty(err)
                fprintf('[test_jwave] %s\n', err);
                py_bin = strtrim(input( ...
                    sprintf('Enter path to Python executable or conda/venv root for %s: ', var_name), 's'));
                if isempty(py_bin)
                    err = sprintf('No Python path provided for %s — cannot continue.', var_name);
                    return
                end
                [py_bin, err] = test_jwave.try_resolve(py_bin, var_name);
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

        function ref = get_ref_data(tc)
        % Return the k-Wave reference result, running it once and caching it.
            if isempty(tc.RefData)
                p = tc.make_params('kwave_ref', [48 48 80]);
                fprintf('[test_jwave] Running k-Wave reference simulation...\n');
                t0 = tic;
                out = prestus_pipeline(p);
                tc.RefData.elapsed     = toc(t0);
                tc.RefData.sensor_data = out.sensor_data;
                tc.RefData.kgrid       = out.kgrid;
                tc.RefData.parameters  = p;
                fprintf('[test_jwave] k-Wave reference done in %.1f s.\n', tc.RefData.elapsed);
            end
            ref = tc.RefData;
        end

        function p = make_params(tc, affix, grid_dims)
        % Build a minimal k-Wave water-medium parameter struct.
            sim_dir = fullfile(tc.SimPath, affix);
            mkdir(sim_dir);
            p = make_minimal_parameters(sim_dir);
            p.grid.default_dims   = grid_dims;
            p.grid.use_kWaveArray = 0;
            p.io.output_affix     = ['_' affix];
            p.io.overwrite_files  = 'always';
        end

        function p = make_jwave_params(tc, affix, grid_dims)
        % Build parameters for a jwave Helmholtz backend run.
            p = tc.make_params(affix, grid_dims);
            p.acoustic_backend = 'jwave';
            p.jwave.pml_size   = 20;
            p.jwave.precision  = 'single';
            p.jwave.dry_run    = false;

            % Pass resolved executable paths (not raw env-var values).
            p.jwave.python_env     = tc.PyBinCpu;
            p.jwave.python_env_gpu = tc.PyBinGpu;

            % use_gpu override (optional; [] = auto-detect via nvidia-smi)
            use_gpu_str = getenv('JWAVE_USE_GPU');
            if strcmpi(use_gpu_str, 'true')
                p.jwave.use_gpu = true;
            elseif strcmpi(use_gpu_str, 'false')
                p.jwave.use_gpu = false;
            else
                p.jwave.use_gpu = [];   % auto-detect
            end
        end

    end

end
