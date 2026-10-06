function sensor_data = jwave_wrapper(parameters, kgrid, medium, source, source_labels, sensor)
% JWAVE_WRAPPER  Run acoustic simulation via jwave (JAX-based Python backend)
%
% Translates PRESTUS kgrid/medium/source structs to HDF5, invokes
% run_jwave_sim.py, and reads back the pressure field in the same
% sensor_data format returned by acoustic_simulation.m.
%
% Uses the jwave Helmholtz solver: solves the heterogeneous Helmholtz equation
% with Fourier spectral collocation and GMRES, matching the solver used in the
% ITRUSST intercomparison benchmarks. Attenuation is applied via db2neper with
% fixed exponent y=2.0 inside jwave's wavevector operator.
%
% Both source paths are supported: direct geometry (use_kWaveArray == 0,
% makeBowl/makeArc voxels) and kWaveArray (use_kWaveArray == 1, band-limited
% off-grid source; recommended, as for k-Wave). A kWaveArray source arrives as
% a CWSourceMatrix and is written in factored form (per-voxel complex amplitude
% x [cos; -sin] basis). Axisymmetric grids are not supported. Thermal
% simulation is not supported.
%
% Use as:
%   sensor_data = jwave_wrapper(parameters, kgrid, medium, source, source_labels, sensor)
%
% Input:
%   parameters    - PRESTUS config; must contain parameters.jwave (see below)
%   kgrid         - kWaveGrid with time axis set (kgrid.t_array ready)
%   medium        - k-Wave medium struct (sound_speed, density, alpha_coeff,
%                   alpha_power); thermal fields are silently ignored
%   source        - struct with p_mask [grid_dims logical] and p [N_src x Nt]
%   source_labels - [grid_dims] int array (element IDs, 0 = inactive); used
%                   for provenance only, not required for the simulation
%   sensor        - (optional) k-Wave sensor struct from source_sensor_setup.m;
%                   sensor.record_start_index defines the steady-state window
%                   used both by the time-domain solvers' p_max_all/p_min_all
%                   and by the Helmholtz solver's source-amplitude DFT.
%                   Defaults to struct() (whole signal, 1-based start=1) when
%                   omitted, e.g. for standalone/benchmark callers.
%
% parameters.jwave fields (all optional with defaults):
%   .python_env     - path to CPU Python venv/conda env root or bin dir, or '' for system python3
%   .python_env_gpu - path to GPU (CUDA) Python env bin dir; leave empty to always use CPU
%   .use_gpu        - [] auto-detect via nvidia-smi (default) | true force GPU | false force CPU
%   .solver         - 'helmholtz' (default) | 'pseudo_spectral' | 'finite_differences'
%                     helmholtz:          frequency-domain GMRES, single-frequency
%                     pseudo_spectral:    PSTD time-domain with Stokes (y=2)
%                                         absorption; matches kspaceFirstOrder*
%                     finite_differences: same time-domain step as pseudo_spectral
%                                         but with local FD stencils (order fd_accuracy)
%   .fd_accuracy    - FD stencil accuracy order for finite_differences (default: 8)
%   .gmres_restart  - GMRES restart parameter (default: 10; helmholtz only;
%                     matches jwave's own helmholtz_solver default)
%   .gmres_maxiter  - Max GMRES outer iterations (default: 1000; helmholtz only;
%                     matches jwave's own helmholtz_solver default)
%   .helmholtz_operator - 'conservative' (default) | 'jwave'; helmholtz only.
%                     conservative: rho*div((1/rho)*grad p), source scaled by the
%                     local sound speed. jwave: jwave's built-in density term,
%                     which underestimates transmission through density jumps
%                     (~20 % through a skull plate; see run_jwave_sim.py).
%   .gmres_tol      - GMRES convergence tolerance (default: 1e-3; helmholtz only;
%                     matches jwave's own helmholtz_solver default). A tighter
%                     restart/maxiter/tol budget does not improve accuracy once
%                     the Helmholtz attenuation floor regularizes the lossless
%                     case (see run_jwave_sim.py's _apply_helmholtz_attenuation_floor)
%                     — it only costs more GMRES iterations/memory for no benefit.
%   .precision      - 'single' (default) | 'double'
%   .dry_run        - true to write inputs and exit without simulating (default false)
%   .convergence_threshold - relative Helmholtz spectral residual above which the
%                     result is rejected (default: 0.05; helmholtz only). Above
%                     threshold, jwave_wrapper ERRORS (not just warns) because an
%                     unconverged GMRES iterate is dominated by solver noise, not
%                     physics — it can simultaneously understate the true focal
%                     amplitude and fabricate spurious ringing near sharp
%                     impedance interfaces (skull/skin). The error message lists
%                     concrete fixes (GMRES restart/precision/PML-fraction).
%   .ignore_convergence_threshold - false (default) | true; if true, a residual
%                     above convergence_threshold only warns (previous behaviour)
%                     instead of erroring, and the unconverged result is returned.
%                     Only set this after having reviewed the field and diagnostics
%                     and having decided the residual is acceptable for your use case.
%
% PML is controlled by parameters.grid.pml_size (the field name is shared with
% k-Wave, but 'auto' resolves differently per backend).
%   grid.pml_size   - integer or 'auto'.  For jwave, 'auto' is resolved by
%                     resolve_jwave_pml_size to ceil(2.5 wavelengths in grid
%                     points), clipped to [10, 40] -- NOT via k-Wave's
%                     getOptimalPMLSize, which optimizes FFT prime factors for
%                     k-Wave's own backend and has no wavelength awareness
%                     (see resolve_jwave_pml_size.m for the rationale). The
%                     resulting scalar is broadcast to a per-axis vector before
%                     being written to the HDF5 input; run_jwave_sim.py's
%                     auto_pad_for_pml uses the per-axis values, and jwave's
%                     Medium (a single scalar pml_size) gets max() of them.
%
% Output:
%   sensor_data.p_max_all   [grid_dims] double [Pa]  peak positive pressure
%   sensor_data.p_min_all   [grid_dims] double [Pa]  peak negative pressure
%   sensor_data.p_final     [grid_dims] double [Pa]  pressure at final time step
%   sensor_data.jwave_input_h5   char  path to HDF5 input written
%   sensor_data.jwave_output_h5  char  path to HDF5 output read
%   sensor_data.jwave_pml_size   double  resolved PML thickness [grid points] (scalar or per-axis)
%
% See also: BABELBRAIN_WRAPPER, ACOUSTIC_SIMULATION, PRESTUS_PIPELINE

arguments
    parameters    (1,1) struct
    kgrid         (1,1)
    medium        (1,1) struct
    source        (1,1) struct
    source_labels = []
    sensor        = struct()
end

% sensor.record_start_index (from source_sensor_setup.m) marks the start of
% the "last N periods, after the wavefront has crossed the domain" steady-state
% window that k-Wave and jwave's time-domain solvers already use for
% p_max_all/p_min_all. The Helmholtz solver's source-amplitude DFT uses the
% same window (see run_jwave_sim.py:run_helmholtz) so all three jwave solvers
% share one definition of "steady state". Falls back to 1 (use the whole
% signal) for callers that don't pass sensor (e.g. water_baseline.m,
% benchmark_jwave.m dry-run calls), preserving prior behaviour for them.
record_start_index = 1;
if isfield(sensor, 'record_start_index') && ~isempty(sensor.record_start_index)
    record_start_index = sensor.record_start_index;
end

% -------------------------------------------------------------------------
% 0. Resolve parameters
% -------------------------------------------------------------------------

jw = struct();
if isfield(parameters, 'jwave')
    jw = parameters.jwave;
end

if isfield(jw, 'pml_size')
    warning('jwave_wrapper:ignoredField', ...
        ['parameters.jwave.pml_size is ignored. ' ...
         'Set parameters.grid.pml_size instead (shared with k-Wave).']);
end

python_env     = get_field(jw, 'python_env',     '');
python_env_gpu = get_field(jw, 'python_env_gpu', '');
use_gpu        = get_field(jw, 'use_gpu',        []);   % [] = auto-detect
% PML size is read from grid.pml_size (shared with k-Wave) so both backends
% use the same value without per-backend configuration.
pml_size       = get_field(parameters.grid, 'pml_size', 'auto');
solver         = get_field(jw, 'solver',          'helmholtz');
fd_accuracy    = get_field(jw, 'fd_accuracy',     8);
gmres_restart  = get_field(jw, 'gmres_restart',   10);
gmres_maxiter  = get_field(jw, 'gmres_maxiter',   1000);
gmres_tol      = get_field(jw, 'gmres_tol',       1e-3);
helmholtz_operator = get_field(jw, 'helmholtz_operator', 'conservative');
precision      = get_field(jw, 'precision',       'single');
dry_run        = get_field(jw, 'dry_run',         false);

valid_solvers = {'helmholtz', 'pseudo_spectral', 'finite_differences'};
if ~any(strcmp(solver, valid_solvers))
    error('jwave_wrapper: unknown solver ''%s''. Valid options: %s.', ...
          solver, strjoin(valid_solvers, ', '));
end
valid_operators = {'conservative', 'jwave'};
if ~any(strcmp(helmholtz_operator, valid_operators))
    error('jwave_wrapper: unknown helmholtz_operator ''%s''. Valid options: %s.', ...
          helmholtz_operator, strjoin(valid_operators, ', '));
end

% convergence_threshold — relative spectral residual above which the Helmholtz
% result is rejected (see step 5b below for the error/warning behaviour).
convergence_threshold = get_field(jw, 'convergence_threshold', 0.05);
% ignore_convergence_threshold — if true, a residual above convergence_threshold
% only warns (old behaviour) instead of erroring. Use only after having reviewed
% the field and accepted the residual; see the error message in step 5b.
ignore_convergence_threshold = get_field(jw, 'ignore_convergence_threshold', false);

% Resolve python executable, selecting GPU env when appropriate
[python_bin, using_gpu_env] = resolve_python_env(python_env, python_env_gpu, use_gpu);

% -------------------------------------------------------------------------
% 1. Ensure HDF5 plugin path is valid before any h5read/h5write call
% -------------------------------------------------------------------------
ensure_hdf5_plugin_path(python_bin);

% Locate the Python bridge script (same folder as this .m file)
this_dir   = fileparts(mfilename('fullpath'));
bridge_py  = fullfile(this_dir, 'run_jwave_sim.py');
if ~isfile(bridge_py)
    error('jwave_wrapper: run_jwave_sim.py not found at:\n  %s\n%s', ...
        bridge_py, 'This file must be in functions/jwave/ alongside jwave_wrapper.m.');
end

% -------------------------------------------------------------------------
% 2. Prepare output paths
% -------------------------------------------------------------------------

out_dir    = char(parameters.io.dir_cache);
affix      = char(parameters.io.output_affix);
subj_id    = parameters.subject_id;
medium_tag = parameters.simulation.medium;

input_h5  = fullfile(out_dir, sprintf('sub-%03d_%s_jwave_input%s.h5',  subj_id, medium_tag, affix));
output_h5 = fullfile(out_dir, sprintf('sub-%03d_%s_jwave_output%s.h5', subj_id, medium_tag, affix));

% -------------------------------------------------------------------------
% 3. Write HDF5 input
% -------------------------------------------------------------------------

% jwave's PML has none of k-Wave's FFT-prime-factorization sensitivity and
% only supports a single scalar thickness (jwave.geometry.Medium.pml_size),
% so 'auto' is resolved from a wavelength-based target instead of reusing
% k-Wave's getOptimalPMLSize. Sized from max(sound_speed): see
% resolve_jwave_pml_size.m.
pml_size = resolve_jwave_pml_size(pml_size, parameters.transducer(1).freq_hz, ...
    max(medium.sound_speed(:)), parameters.grid.resolution_mm / 1e3);

% jwave applies PML inside the declared grid (PMLInside is not configurable),
% so voxels within pml_size of any face are absorbing layers.  run_jwave_sim.py
% compensates via auto_pad_for_pml; this check exposes edge cases where the
% transducer sits unexpectedly close to the boundary.
jw_pml_max  = max(pml_size(:));
jw_grid_sz  = parameters.grid.dims(:)';
jw_nd       = numel(jw_grid_sz);
jw_src_idx  = find(source.p_mask(:));
if ~isempty(jw_src_idx)
    if jw_nd == 3
        [jw_ax{1}, jw_ax{2}, jw_ax{3}] = ind2sub(jw_grid_sz, jw_src_idx);
    else
        [jw_ax{1}, jw_ax{2}] = ind2sub(jw_grid_sz, jw_src_idx);
    end
    jw_src_in_pml = false;
    for jw_a = 1:jw_nd
        lo = min(jw_ax{jw_a});
        hi = max(jw_ax{jw_a});
        if lo <= jw_pml_max || hi > jw_grid_sz(jw_a) - jw_pml_max
            fprintf('[jwave_wrapper] Source spans voxels %d–%d on axis %d (PML boundary ±%d voxels). auto_pad will compensate.\n', ...
                lo, hi, jw_a, jw_pml_max);
            jw_src_in_pml = true;
        end
    end
    if jw_src_in_pml
        fprintf('[jwave_wrapper] run_jwave_sim.py auto_pad_for_pml will expand the grid and crop back; output stays in grid.dims space.\n');
    end
end
clear jw_pml_max jw_grid_sz jw_nd jw_src_idx jw_ax jw_src_in_pml jw_a

% --- Simulation summary ---
grid_dims   = parameters.grid.dims;
n_voxels    = prod(grid_dims);
n_src       = nnz(source.p_mask);
freq_kHz    = parameters.transducer(1).freq_hz / 1e3;
fprintf('[jwave_wrapper] ---- jwave simulation summary ----\n');
fprintf('[jwave_wrapper]   solver    : %s\n', solver);
fprintf('[jwave_wrapper]   precision : %s\n', precision);
fprintf('[jwave_wrapper]   backend   : %s\n', ternary(using_gpu_env, 'GPU (CUDA)', 'CPU'));
fprintf('[jwave_wrapper]   grid      : [%s] = %dk voxels\n', num2str(grid_dims), round(n_voxels/1e3));
fprintf('[jwave_wrapper]   grid res  : %.2f mm\n', parameters.grid.resolution_mm);
fprintf('[jwave_wrapper]   time      : Nt=%d  dt=%.3e s\n', kgrid.Nt, kgrid.dt);
fprintf('[jwave_wrapper]   source    : %d active voxels  freq=%.1f kHz\n', n_src, freq_kHz);
% jwave applies PML *inside* the specified grid (unlike k-Wave which adds it
% outside).  The scalar used by jwave is max over axes; report both.
pml_scalar   = max(pml_size(:));
pml_interior = grid_dims(:)' - 2 * pml_scalar;   % effective non-PML region
if isscalar(pml_size) || numel(unique(pml_size(:))) == 1
    pml_size_str = num2str(pml_scalar);
else
    pml_size_str = sprintf('[%s] → max=%d', num2str(pml_size(:)'), pml_scalar);
end
fprintf('[jwave_wrapper]   PML size  : %s voxels/side  location=inside-grid  interior=[%s]\n', ...
    pml_size_str, num2str(pml_interior));
if isscalar(medium.sound_speed)
    fprintf('[jwave_wrapper]   medium    : homogeneous  c=%.0f m/s  rho=%.0f kg/m3\n', ...
        medium.sound_speed, medium.density);
else
    fprintf('[jwave_wrapper]   medium    : heterogeneous  c=[%.0f–%.0f] m/s  rho=[%.0f–%.0f] kg/m3\n', ...
        min(medium.sound_speed(:)), max(medium.sound_speed(:)), ...
        min(medium.density(:)),     max(medium.density(:)));
end
fprintf('[jwave_wrapper] ------------------------------------\n');

% -------------------------------------------------------------------------
% 3b. Pre-flight: auto-reduce dt for Stokes absorption stability (PSTD/FD)
% -------------------------------------------------------------------------
% For time-domain solvers the explicit Stokes absorption step is numerically
% stable only when  tau_max * c_max^2 * (pi/dx)^2 * dt  <  ~0.30.
% High skull attenuation (>20 dB/cm/MHz) can push this into the unstable
% regime even at CFL=0.10.  We detect this and reduce dt before generating
% the HDF5 input so the source signals are consistent with the new timestep.
% Physics (alpha_coeff) are unchanged — only the timestep is finer.
kgrid_dt_write = kgrid.dt;
kgrid_Nt_write = kgrid.Nt;
if ismember(solver, {'pseudo_spectral', 'finite_differences'})
    % Absorption metric safety ceiling.  0.25 is the nominal edge of stability
    % for the explicit Stokes-absorption time integration; it is a *temporal*
    % effect (each step amplifies high-wavenumber noise) shared by both PSTD
    % and finite_differences, not specific to either solver's spatial
    % derivative operator.  Apply a safety factor so the auto-reduced dt lands
    % with headroom below the boundary instead of exactly on it.  Overridable
    % via parameters.jwave.absorb_safety.
    %
    % History note: both solvers were previously observed to produce
    % late-onset NaN while landing at exactly metric=0.2500 (zero margin) on
    % small, PML-dominated pFOV domains.  That NaN was NOT actually caused by
    % insufficient margin on this metric — the true cause was a missing
    % negative sign in run_jwave_sim.py's _absorb_tau_from_power_law, which
    % made the Stokes term energy-injecting instead of dissipative (fixed).
    % This safety factor is kept as a legitimate, independent margin on a
    % correctly-signed dissipative term (this metric bounds |tau|, so a sign
    % flip would not have changed its value or been caught by it) — not
    % because it was the actual fix for that bug.
    jw_absorb_target = 0.25;
    jw_absorb_safety = 0.80;   % default headroom for both time-domain solvers
    if isfield(parameters, 'jwave') && isfield(parameters.jwave, 'absorb_safety') ...
            && ~isempty(parameters.jwave.absorb_safety)
        jw_absorb_safety = parameters.jwave.absorb_safety;
    end
    jw_absorb_target = jw_absorb_target * jw_absorb_safety;
    jw_freq_hz = parameters.transducer(1).freq_hz;
    jw_f_MHz   = jw_freq_hz / 1e6;
    if isscalar(medium.alpha_coeff)
        jw_alpha_max = medium.alpha_coeff;
    else
        jw_alpha_max = max(medium.alpha_coeff(:));
    end
    if isfield(medium, 'alpha_power') && ~isempty(medium.alpha_power)
        if isscalar(medium.alpha_power)
            jw_alpha_power = medium.alpha_power;
        else
            jw_alpha_power = max(medium.alpha_power(:));
        end
    else
        jw_alpha_power = 1.0;
    end
    if isscalar(medium.sound_speed)
        jw_c_max = medium.sound_speed;
    else
        jw_c_max = max(medium.sound_speed(:));
    end
    % tau_max: Stokes relaxation time at transducer frequency
    jw_att_Np_m = jw_alpha_max * (jw_f_MHz ^ jw_alpha_power) * (log(10)/20) * 100;
    jw_omega    = 2 * pi * jw_freq_hz;
    jw_tau_max  = 2 * jw_att_Np_m * jw_c_max / jw_omega^2;
    jw_dx       = parameters.grid.resolution_mm / 1e3;
    jw_k_nyq    = pi / jw_dx;
    jw_absorb_metric = jw_tau_max * jw_c_max^2 * jw_k_nyq^2 * kgrid.dt;
    if jw_absorb_metric > jw_absorb_target && jw_tau_max > 0
        jw_dt_safe = jw_absorb_target / (jw_tau_max * jw_c_max^2 * jw_k_nyq^2);
        jw_t_end   = kgrid.Nt * kgrid.dt;
        jw_Nt_safe = ceil(jw_t_end / jw_dt_safe);
        fprintf('[jwave_wrapper] Stokes absorption metric %.4f > %.2f — ', ...
            jw_absorb_metric, jw_absorb_target);
        fprintf('auto-reducing dt: %.3e → %.3e s, Nt: %d → %d\n', ...
            kgrid.dt, jw_dt_safe, kgrid.Nt, jw_Nt_safe);
        fprintf('[jwave_wrapper] Physics (alpha_coeff) unchanged; timestep is finer.\n');
        % Resample source signals to new time axis.
        % pchip interpolation is accurate for smooth CW sinusoids.
        jw_t_old = kgrid.t_array;                              % [1 x Nt]
        jw_t_new = (0:jw_Nt_safe-1) * jw_dt_safe;             % [1 x Nt_safe]
        if isa(source.p, 'CWSourceMatrix')
            % kWaveArray CW source: rebuild on the new time axis (exact)
            jw_coeffs = factors(source.p);
            source.p = CWSourceMatrix(jw_coeffs(:, 1), jw_coeffs(:, 2), jw_omega, jw_dt_safe, jw_Nt_safe);
            clear jw_coeffs
        else
            source.p = interp1(jw_t_old', source.p', jw_t_new', 'pchip', 0)';
        end
        % Preserve the same physical recording window at the end of the run.
        jw_t_rec_window = (kgrid.Nt - record_start_index + 1) * kgrid.dt;
        record_start_index = max(1, ...
            jw_Nt_safe - round(jw_t_rec_window / jw_dt_safe));
        kgrid_dt_write = jw_dt_safe;
        kgrid_Nt_write = jw_Nt_safe;
    end
    clear jw_absorb_target jw_absorb_safety jw_freq_hz jw_f_MHz jw_alpha_max jw_alpha_power
    clear jw_c_max jw_att_Np_m jw_omega jw_tau_max jw_dx jw_k_nyq
    clear jw_absorb_metric jw_dt_safe jw_t_end jw_Nt_safe jw_t_old jw_t_new jw_t_rec_window
end

fprintf('[jwave_wrapper] Writing simulation input: %s\n', input_h5);
write_jwave_input(input_h5, kgrid, medium, source, parameters, pml_size, precision, kgrid_dt_write, kgrid_Nt_write, record_start_index);

if dry_run
    fprintf('[jwave_wrapper] Dry run — skipping simulation.\n');
    sensor_data.p_max_all      = [];
    sensor_data.p_min_all      = [];
    sensor_data.p_final        = [];
    sensor_data.jwave_input_h5 = input_h5;
    sensor_data.jwave_output_h5 = '';
    sensor_data.jwave_pml_size  = pml_size;
    return;
end

% -------------------------------------------------------------------------
% 4. Invoke Python bridge
% -------------------------------------------------------------------------

% run_command's Env option sets JAX_PLATFORMS natively on the child
% process -- no `VAR=val cmd` shell prefix syntax needed at all.
env_overrides = struct();
if using_gpu_env
    env_overrides.JAX_PLATFORMS = 'cuda,cpu';
else
    env_overrides.JAX_PLATFORMS = 'cpu';
end

% run_command invokes python directly via argv, no shell at all --
% input_h5/output_h5/parameters.io.log_file reach the process (or are
% used directly by MATLAB) as literal values regardless of their contents.
jwave_argv = {python_bin, bridge_py, '--input', input_h5, '--output', output_h5, ...
    '--solver', solver, ...
    '--gmres_restart', num2str(gmres_restart), ...
    '--gmres_maxiter', num2str(gmres_maxiter), ...
    '--gmres_tol', sprintf('%.1e', gmres_tol), ...
    '--fd_accuracy', num2str(fd_accuracy), ...
    '--helmholtz_operator', helmholtz_operator};
% parameters.simulation.debug drives run_jwave_sim.py's --debug: extra prints
% for PML/domain sizing in wavelengths, time-domain steady-state margin, and
% a Helmholtz compile-vs-solve time breakdown (see run_jwave_sim.py --help).
is_debug = isfield(parameters, 'simulation') && isfield(parameters.simulation, 'debug') ...
    && isequal(parameters.simulation.debug, 1);
if is_debug
    jwave_argv = [jwave_argv, {'--debug'}];
end

fprintf('[jwave_wrapper] Launching jwave simulation...\n');
fprintf('[jwave_wrapper] Command:\n  %s\n', strjoin(jwave_argv, ' '));

% Echo=true streams the Python subprocess's stdout/stderr to the MATLAB
% console live as it is produced (matching the cpp_cpu/cpp_gpu backends'
% behaviour, see kspaceFirstOrder3DC.m; matches the original system(cmd,
% '-echo')), rather than buffering it until the process exits. This is
% what makes run_jwave_sim.py's progress/ETA prints usable as real-time
% feedback instead of a wall of text dumped after the fact.
t_start = tic;
[status, cmd_out] = run_command(jwave_argv, 'Env', env_overrides, 'Echo', true);
elapsed = toc(t_start);

% Append Python stdout+stderr to the main PRESTUS log (no separate log
% file). Previously a shell `>> logfile 2>&1` redirect diverted the
% subprocess's output straight to the file (silently, with no console
% echo); now the full captured output is appended once the process
% finishes, in addition to the live console echo above.
if isfield(parameters, 'io') && isfield(parameters.io, 'log_file') && ~isempty(parameters.io.log_file)
    log_fid = fopen(parameters.io.log_file, 'a');
    if log_fid ~= -1
        fprintf(log_fid, '%s', cmd_out);
        fclose(log_fid);
    end
end

if status ~= 0
    error('jwave_wrapper: Python process exited with status %d.\n%s', status, cmd_out);
end

fprintf('[jwave_wrapper] jwave finished in %.1f s.\n', elapsed);

% -------------------------------------------------------------------------
% 5. Read HDF5 output
% -------------------------------------------------------------------------

if ~isfile(output_h5)
    error('jwave_wrapper: Output HDF5 not found.\n  Expected: %s', output_h5);
end

fprintf('[jwave_wrapper] Reading output: %s\n', output_h5);
sensor_data = read_jwave_output(output_h5, parameters);

sensor_data.jwave_input_h5  = input_h5;
sensor_data.jwave_output_h5 = output_h5;
sensor_data.jwave_pml_size  = pml_size;   % resolved grid.pml_size (reported as grid.pml_size_effective)

% -------------------------------------------------------------------------
% 5b. Read convergence diagnostics
% -------------------------------------------------------------------------

sensor_data.diagnostics = read_jwave_diagnostics(output_h5);
diag = sensor_data.diagnostics;

if isfield(diag, 'residual_norm') && ~isnan(diag.residual_norm)
    if diag.residual_norm > convergence_threshold
        msg = sprintf(['Helmholtz spectral residual %.4e exceeds threshold %.4e — GMRES did not converge.\n' ...
             'p_max_all/p_min_all/p_final are unreliable and must not be used as-is.\n' ...
             'Likely causes and suggested fixes, roughly in order of how often they are the actual cause:\n' ...
             '  1. gmres_restart/gmres_maxiter set above jwave''s own defaults (restart=10, maxiter=1000).\n' ...
             '     A larger restart in single precision can destabilize GMRES instead of helping it\n' ...
             '     converge (loss of Krylov-basis orthogonality). Currently restart=%d, maxiter=%d —\n' ...
             '     try the defaults first, do not just increase them further.\n' ...
             '  2. parameters.jwave.precision = ''double'' instead of ''%s'' — single precision can be\n' ...
             '     insufficient for ill-conditioned heterogeneous (e.g. skull) media.\n' ...
             '  3. PML consuming a large fraction of the domain (see the --debug PML-sizing print;\n' ...
             '     >25%% is a red flag). Reduce parameters.grid.pml_size, or enlarge\n' ...
             '     grid.acoustic_fov_length_mm/grid.acoustic_fov_diameter_mm so the PML is a smaller\n' ...
             '     share of the grid.\n' ...
             '  4. Decreasing parameters.jwave.gmres_tol (currently %.1e) — only likely to help once\n' ...
             '     1-3 have been ruled out.\n' ...
             '  5. If you have reviewed the field and accept this residual anyway, set\n' ...
             '     parameters.jwave.ignore_convergence_threshold = true to downgrade this to a warning\n' ...
             '     and proceed with the unconverged result.'], ...
            diag.residual_norm, convergence_threshold, gmres_restart, gmres_maxiter, precision, gmres_tol);
        if ignore_convergence_threshold
            warning('jwave_wrapper:helmholtzResidual', '%s', msg);
        else
            error('jwave_wrapper:helmholtzResidualExceeded', '%s', msg);
        end
    else
        fprintf('[jwave_wrapper] Helmholtz residual: %.4e (threshold %.4e) — OK\n', ...
            diag.residual_norm, convergence_threshold);
    end
end

% Verify output dimensions match PRESTUS grid
expected_dims = parameters.grid.dims(:)';
for chk_ = {'p_max_all', 'p_min_all', 'p_final'}
    arr = sensor_data.(chk_{1});
    if ~isempty(arr) && ~isequal(size(arr), expected_dims)
        error(['jwave_wrapper: output %s has size [%s] but expected PRESTUS grid [%s].\n' ...
            'Check dimension transpositions in run_jwave_sim.py.'], ...
            chk_{1}, num2str(size(arr)), num2str(expected_dims));
    end
end

end


% =========================================================================
% Local helper: write_jwave_input
% =========================================================================

function write_jwave_input(h5_path, kgrid, medium, source, parameters, pml_size, precision, dt_override, Nt_override, record_start_index)
% Write kgrid, medium, and source to HDF5 for run_jwave_sim.py.
% Writes to a local temp file first to avoid HDF5 file-locking errors on
% network/Docker filesystems where MATLAB's repeated h5create open/close
% cycles fail on the second H5F.open call.
%
% Optional dt_override / Nt_override: if provided, these replace kgrid.dt
% and kgrid.Nt in the HDF5 (used by the absorption-stability auto-dt-reduction
% block in jwave_wrapper to write a finer timestep without modifying kgrid).
%
% record_start_index (1-based, MATLAB convention): start of the steady-state
% recording window (see jwave_wrapper's top-level resolution of
% sensor.record_start_index). Defaults to 1 (whole signal) if omitted.

if nargin < 10 || isempty(record_start_index)
    record_start_index = 1;
end

if isfile(h5_path)
    delete(h5_path);
end

write_path = [tempname, '.h5'];
tmp_guard  = onCleanup(@() h5_tmp_cleanup(write_path));

cast_fn = @(x) cast(gather(x), precision);

% --- Grid ---
dims = numel(parameters.grid.dims);
h5create(write_path, '/grid/N',  [1 dims]);
h5write( write_path, '/grid/N',  int32(parameters.grid.dims(:)'));
h5create(write_path, '/grid/dx', [1 dims]);
h5write( write_path, '/grid/dx', cast_fn(parameters.grid.resolution_mm / 1e3 * ones(1, dims)));

% --- Time axis ---
Nt = kgrid.Nt;
dt = kgrid.dt;
if nargin >= 9 && ~isempty(Nt_override), Nt = Nt_override; end
if nargin >= 8 && ~isempty(dt_override), dt = dt_override; end
h5create(write_path, '/time/Nt', [1 1]);
h5write( write_path, '/time/Nt', int32(Nt));
h5create(write_path, '/time/dt', [1 1]);
h5write( write_path, '/time/dt', cast_fn(dt));

% --- Medium ---
thermal_fields = {'thermal_conductivity', 'specific_heat', 'perfusion_coeff'};
for f = thermal_fields
    if isfield(medium, f{1}), medium = rmfield(medium, f{1}); end
end

grid_sz = parameters.grid.dims(:)';
write_medium_field(write_path, 'sound_speed', medium.sound_speed, grid_sz, cast_fn);
write_medium_field(write_path, 'density',     medium.density,     grid_sz, cast_fn);
write_medium_field(write_path, 'alpha_coeff', medium.alpha_coeff, grid_sz, cast_fn);

alpha_power = get_alpha_power(medium, grid_sz);
write_medium_field(write_path, 'alpha_power', alpha_power, grid_sz, cast_fn);

% --- Source ---
% p_mask: store linear indices of active source voxels
src_indices = int32(find(source.p_mask(:))' - 1);  % 0-based for Python
h5create(write_path, '/source/mask_indices', size(src_indices));
h5write( write_path, '/source/mask_indices', src_indices);

if isa(source.p, 'CWSourceMatrix')
    % kWaveArray CW source: every voxel has its own BLI weight, so rows do not
    % deduplicate. Write the rank-2 factors instead: voxel i emits
    % coeffs(i,:) * basis with basis = [cos(wt); -sin(wt)] (2 x Nt).
    [coeffs, basis] = factors(source.p);
    h5create(write_path, '/source/p',      size(basis));
    h5write( write_path, '/source/p',      cast_fn(basis));
    h5create(write_path, '/source/coeffs', size(coeffs));
    h5write( write_path, '/source/coeffs', cast_fn(coeffs));
    fprintf('[jwave_wrapper] Source signals: %d voxels written as complex amplitudes (kWaveArray CW source)\n', ...
        size(coeffs, 1));
else
    % Source signal: deduplicate rows to reduce HDF5 size.
    % Voxels belonging to the same transducer element share the same waveform,
    % so storing only unique signals can cut the matrix from [N_src x Nt] to
    % [N_elements x Nt] — typically a 100-1000x reduction.
    src_p_cast = cast_fn(source.p);
    [unique_sigs, ~, sig_map] = unique(src_p_cast, 'rows', 'stable');
    n_unique = size(unique_sigs, 1);
    if n_unique < size(source.p, 1)
        fprintf('[jwave_wrapper] Source signals: %d voxels → %d unique waveforms (%.0f%% reduction in source data)\n', ...
            size(source.p,1), n_unique, (1 - n_unique/size(source.p,1))*100);
    end
    chunk_p = [min(256, n_unique), min(512, size(unique_sigs, 2))];
    h5create(write_path, '/source/p',       size(unique_sigs), 'Deflate', 4, 'ChunkSize', chunk_p);
    h5write( write_path, '/source/p',       unique_sigs);
    h5create(write_path, '/source/sig_map', [1, numel(sig_map)]);
    h5write( write_path, '/source/sig_map', int32(sig_map(:)' - 1));  % 0-based for Python
end

% Record start index (1-based MATLAB → 0-based Python)
h5create(write_path, '/source/record_start_index', [1 1]);
h5write( write_path, '/source/record_start_index', int32(record_start_index - 1));

% --- PML ---
% Write the full per-axis PML vector so run_jwave_sim.py can apply per-axis
% padding in auto_pad_for_pml.  getOptimalPMLSize returns a vector for 3D
% grids; broadcast a scalar to length ndims.  The Python bridge uses max()
% when constructing jwave's Medium (which accepts a scalar pml_size), but
% uses the per-axis values for the source-boundary padding step.
pml_size_vec = int32(pml_size(:)');
if numel(pml_size_vec) == 1
    pml_size_vec = repmat(pml_size_vec, 1, dims);
end
h5create(write_path, '/pml/size', [1 numel(pml_size_vec)]);
h5write( write_path, '/pml/size', pml_size_vec);

% --- Metadata ---
h5writeatt(write_path, '/', 'prestus_version', 'PRESTUS');
h5writeatt(write_path, '/', 'dims',            int32(dims));
h5writeatt(write_path, '/', 'freq_hz',         parameters.transducer(1).freq_hz);

movefile(write_path, h5_path);

end


function h5_tmp_cleanup(tmp_path)
if isfile(tmp_path), delete(tmp_path); end
end


function write_medium_field(h5_path, name, val, grid_sz, cast_fn)
% Write a medium field that may be scalar or a full spatial array.
dset = ['/medium/' name];
if isscalar(val)
    h5create(h5_path, dset, [1 1]);
    h5write( h5_path, dset, cast_fn(val));
else
    h5create(h5_path, dset, grid_sz);
    h5write( h5_path, dset, cast_fn(val));
end
end


function alpha_power = get_alpha_power(medium, grid_sz)
if isfield(medium, 'alpha_power')
    alpha_power = medium.alpha_power;
else
    alpha_power = 2.0;  % k-Wave fractional Laplacian default
end
% Broadcast scalar to match expectations on Python side
if isscalar(alpha_power)
    alpha_power = alpha_power * ones(grid_sz);
end
end


% =========================================================================
% Local helper: read_jwave_diagnostics
% =========================================================================

function diag = read_jwave_diagnostics(h5_path)
% Read scalar diagnostic fields from /diagnostics/ group in jwave output HDF5.
% Missing fields are returned as NaN so callers can use isfield + isnan.

diag = struct();
diag_fields = {'residual_norm'};
for f = diag_fields
    key = ['/diagnostics/' f{1}];
    try
        diag.(f{1}) = double(h5read(h5_path, key));
    catch
        diag.(f{1}) = NaN;
    end
end
end


% =========================================================================
% Local helper: read_jwave_output
% =========================================================================

function sensor_data = read_jwave_output(h5_path, parameters)
% Read p_max_all, p_min_all, p_final from jwave HDF5 output.

sensor_data = struct();

fields = {'p_max_all', 'p_min_all', 'p_final'};
for f = fields
    key = ['/pressure/' f{1}];
    try
        sensor_data.(f{1}) = double(h5read(h5_path, key));
    catch ME
        warning('jwave_wrapper: Could not read %s from %s:\n  %s', key, h5_path, ME.message);
        sensor_data.(f{1}) = [];
    end
end

end


% =========================================================================
% Utility
% =========================================================================

function [python_bin, is_gpu] = resolve_python_env(cpu_env, gpu_env, use_gpu)
% Select CPU or GPU Python environment and return the full path to the
% python3 executable, and whether the GPU env was chosen.
%
% use_gpu == true  : force GPU env (error if gpu_env is empty)
% use_gpu == false : force CPU env
% use_gpu == []    : auto-detect via nvidia-smi; use GPU env if a device is
%                    found and gpu_env is non-empty, otherwise fall back to
%                    CPU env with a warning

    if islogical(use_gpu) && use_gpu
        % Explicit GPU request
        if isempty(gpu_env)
            error('jwave_wrapper: parameters.jwave.use_gpu = true but python_env_gpu is not set.');
        end
        fprintf('[jwave_wrapper] GPU env selected (use_gpu = true).\n');
        chosen_env = gpu_env;
        is_gpu = true;

    elseif islogical(use_gpu) && ~use_gpu
        % Explicit CPU request
        chosen_env = cpu_env;
        is_gpu = false;

    else
        % Auto-detect: probe nvidia-smi, but also honour CUDA_VISIBLE_DEVICES.
        % nvidia-smi lists all node GPUs regardless of scheduler allocation;
        % CUDA_VISIBLE_DEVICES="" or "NoDevFiles" means no GPU was granted to
        % this job even if hardware is present.
        [smi_status, smi_out] = run_command({'nvidia-smi', '-L'});
        cvd = getenv('CUDA_VISIBLE_DEVICES');
        cvd_blocked = ~isempty(cvd) && (isempty(strtrim(cvd)) || strcmpi(strtrim(cvd), 'nodevfiles'));
        smi_out_trim = strtrim(smi_out);
        % Strip trailing punctuation so "No devices found." matches the check
        smi_out_norm = regexprep(smi_out_trim, '[.!]+$', '');
        smi_lists_device = (smi_status == 0) && ~isempty(smi_out_norm) && ...
            ~strcmpi(smi_out_norm, 'no devices found') && ...
            ~strcmpi(smi_out_norm, 'no devices were found') && ...
            ~contains(lower(smi_out_norm), 'no devices');
        gpu_found = smi_lists_device && ~cvd_blocked;

        if gpu_found && ~isempty(gpu_env)
            fprintf('[jwave_wrapper] GPU detected (%s) — using GPU env.\n', ...
                smi_out_trim(1:min(end,60)));
            chosen_env = gpu_env;
            is_gpu = true;
        else
            if gpu_found && isempty(gpu_env)
                warning('jwave_wrapper:noGpuEnv', ...
                    'GPU detected but python_env_gpu is not set — falling back to CPU env.');
            end
            chosen_env = cpu_env;
            is_gpu = false;
        end
    end

    % Resolve chosen_env path to a python3 executable
    if isempty(chosen_env)
        python_bin = 'python3';
    elseif isfile(chosen_env)
        % caller passed the binary directly
        python_bin = chosen_env;
    elseif isfile(fullfile(chosen_env, 'bin', 'python3'))
        % venv / conda root
        python_bin = fullfile(chosen_env, 'bin', 'python3');
    elseif isfile(fullfile(chosen_env, 'python3'))
        % bin directory passed directly
        python_bin = fullfile(chosen_env, 'python3');
    else
        warning('jwave_wrapper:pythonEnvNotFound', ...
            'Could not find python3 under ''%s'' — falling back to system python3.', chosen_env);
        python_bin = 'python3';
    end
end


function val = get_field(s, name, default)
if isfield(s, name)
    val = s.(name);
else
    val = default;
end
end


function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
