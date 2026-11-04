# jwave Acoustic Backend

!!! note "At a glance"
    - **What:** jwave as an alternative acoustic solver, with a frequency-domain Helmholtz solve for continuous-wave runs.
    - **When you need it:** To cross-check k-Wave with an independent implementation, or for a fast CW steady-state solve.
    - **Before you start:** A Python 3.11 environment with `jax`, `jwave`, `h5py` and `numpy` ([Python environment setup](#python-environment-setup)).

PRESTUS can run [jwave](https://github.com/ucl-bug/jwave) in place of k-Wave for the acoustic stage. Pre-processing (segmentation, grid, medium masks, placement) and post-processing (focal metrics, NIfTI export, HTML report) run unchanged. jwave (University College London) implements the k-Wave pseudo-spectral equations in [JAX](https://github.com/google/jax), which gives JIT compilation, automatic differentiation and NVIDIA GPU acceleration without compiled MATLAB/C++ binaries. It offers a time-domain solver equivalent to `kspaceFirstOrder*` and a frequency-domain Helmholtz solver that computes the CW steady state directly.

| | PRESTUS (k-Wave) | jwave |
|---|---|---|
| **Solver** | Pseudo-spectral FDTD | Pseudo-spectral (time-domain) or Helmholtz GMRES (frequency-domain) |
| **Physics** | Scalar linear acoustic | Scalar linear acoustic |
| **Skull model** | Scalar ρ, c, α | Scalar ρ, c, α |
| **Differentiable** | No | Yes (JAX automatic differentiation) |
| **GPU backends** | CUDA (via MATLAB `gpuArray` or C++ binary) | NVIDIA CUDA only |
| **Thermal simulation** | Yes (`kWaveDiffusion`) | Not supported |

jwave uses the same scalar acoustic model as k-Wave: no shear waves in bone and no Rayleigh pre-step. The `medium_properties` PRESTUS assigns apply to both equally, so results are directly comparable. See [k-Wave](doc_kwave.md) for the k-Wave backend options.

## Quick start

Set three things: the backend, the Python environment and the solver.

```yaml
acoustic_backend: jwave

jwave:
  python_env: '/home/user/.conda/envs/jwave_env/bin/python'   # CPU environment (see environment setup)
  solver: helmholtz          # 'helmholtz' (default) | 'pseudo_spectral' | 'finite_differences'
  precision: single          # 'single' | 'double'
  python_env_gpu: ''         # GPU (CUDA) environment; leave empty for CPU-only
  use_gpu: ~                 # null = auto-detect via nvidia-smi | true = force GPU | false = force CPU
  dry_run: false             # write HDF5 inputs only, skip simulation (for debugging)

# PML thickness uses the same grid.pml_size field as k-Wave (default: 'auto'),
# but jwave resolves 'auto' to ~2.5 wavelengths of PML (clipped to [10, 40]
# grid points), not k-Wave's FFT-prime-factorization heuristic.
# grid:
#   pml_size: 20             # override if needed
```

Or in MATLAB:

```matlab
parameters.acoustic_backend = 'jwave';

parameters.jwave.solver          = 'helmholtz';
parameters.jwave.precision       = 'single';
parameters.jwave.python_env      = '/home/user/.conda/envs/jwave_env/bin/python';
parameters.jwave.python_env_gpu  = '/home/user/.conda/envs/jwave_env_gpu/bin/python';
parameters.jwave.dry_run         = false;
% parameters.jwave.use_gpu = true;   % override auto-detection
% parameters.grid.pml_size = 20;     % override the wavelength-based 'auto' PML
```

`python_env` is the path of the CPU environment (its root, its `bin` directory or the Python executable). `python_env_gpu` points to the GPU (CUDA) environment. Leaving both empty falls back to the system `python3`, which is unlikely to have jwave installed. You get a `sensor_data` struct with the same fields as from k-Wave, so all downstream PRESTUS modules work unchanged.

??? note "Details: which environment runs"
    | `use_gpu` | GPU detected (`nvidia-smi`) | `python_env_gpu` set | Environment used |
    |---|---|---|---|
    | `[]` (auto) | yes | yes | GPU |
    | `[]` (auto) | yes | no | CPU (warning) |
    | `[]` (auto) | no | – | CPU |
    | `true` | – | yes | GPU |
    | `true` | – | no | **error** |
    | `false` | – | – | CPU |

## Python environment setup

jwave needs a dedicated Python environment. A separate environment avoids version conflicts with other PRESTUS Python tools (BabelBrain, SimNIBS).

**1. Create the environment**

```bash
# conda
conda create -n jwave_env python=3.11
conda activate jwave_env

# venv (no conda)
python3.11 -m venv ~/.venvs/jwave_env
source ~/.venvs/jwave_env/bin/activate
```

**2. Install jwave and its dependencies.** `jax[cpu]` and `jax[cuda12]` install different builds of `jaxlib` and are mutually exclusive: installing one overwrites the other. Pick the variant that matches your hardware.

```bash
# CPU only (macOS, Windows, or CPU-only Linux)
pip install "jax[cpu]" jwave h5py numpy

# NVIDIA GPU (CUDA 12, Linux x86-64 only)
pip install "jax[cuda12]" jwave h5py numpy
```

JAX CUDA wheels are published for Linux x86-64 only; the `jaxlib-cuda` wheel is large (~2 GB). A CUDA install falls back to CPU when no GPU is present at runtime. On macOS or Windows, always use `jax[cpu]`.

**3. Verify the installation**

```bash
python -c "import jwave, jax, h5py; print('jwave:', jwave.__version__); print('jax:', jax.__version__)"
python -c "import jax; print(jax.devices())"
# Expected on GPU host: [CudaDevice(id=0)]
# Expected on CPU host: [CpuDevice(id=0)]
```

**4. Find the Python executable path** and copy it into `jwave.python_env`:

```bash
which python   # inside the activated environment
# e.g. /home/user/.conda/envs/jwave_env/bin/python
#      /home/user/.venvs/jwave_env/bin/python
```

??? note "Details: parallel CPU and GPU environments"
    If you work on both CPU and GPU machines (a laptop and an HPC cluster with shared home directories), maintain two environments and point `jwave.python_env` and `jwave.python_env_gpu` to them:

    ```bash
    # CPU environment (default, works everywhere)
    conda create -n jwave_env python=3.11
    conda activate jwave_env
    pip install "jax[cpu]" jwave h5py numpy

    # GPU environment (Linux + CUDA 12 hosts only)
    conda create -n jwave_env_gpu python=3.11
    conda activate jwave_env_gpu
    pip install "jax[cuda12]" jwave h5py numpy
    ```

    ```yaml
    jwave:
      python_env: '/home/user/.conda/envs/jwave_env/bin/python'       # CPU
      # python_env: '/home/user/.conda/envs/jwave_env_gpu/bin/python'  # GPU
    ```

## Limitations

- **Acoustic only.** jwave has no bio-heat solver. Thermal simulation (`run_heating_sims`) uses k-Wave or BabelBrain.
- **No axisymmetric grids.** jwave requires a full 3-D domain.
- **GPU: NVIDIA CUDA only.** JAX does not officially support AMD ROCm or Apple Metal. CPU execution works on any platform.

## Solver options

### `helmholtz` (default)

Frequency-domain GMRES solver. It computes the complex pressure field `P(r)` at the drive frequency in a single linear solve instead of thousands of time steps, so it is much faster than `pseudo_spectral` for CW simulations, with equivalent accuracy for monochromatic excitation in lossy tissue. It returns `p_max_all` = |P| [Pa], `p_min_all` = −|P| and `p_final` = Re(P) [Pa]. GMRES convergence is controlled by `jwave.gmres_restart` and `jwave.gmres_maxiter`.

!!! note "Convergence is checked and enforced"
    After every Helmholtz solve, PRESTUS compares the relative GMRES residual (`/diagnostics/residual_norm`) against `jwave.convergence_threshold` (default `0.05`). If the residual exceeds the threshold, the pipeline errors, because an unconverged GMRES iterate is dominated by solver noise rather than physics. The error message gives ranked fixes: revert `gmres_restart`/`gmres_maxiter` to jwave's own defaults rather than enlarging them (a larger restart in single precision can destabilise GMRES); try `precision: double`; check the PML fraction of the domain via `--debug`; only then loosen `gmres_tol`. If you have reviewed the field and diagnostics and accept the residual anyway, set `jwave.ignore_convergence_threshold: true` to downgrade the error to a warning and keep the result.

??? note "Details: how the Helmholtz residual, water attenuation and source are handled"
    - **Residual.** The residual is computed exactly, by re-evaluating the operator and source scaling of the solve on the solved field (`conservative_residual`, or `compute_helmholtz_residual` for `helmholtz_operator: jwave`, in `run_jwave_sim.py`), so it reflects the true solver error rather than an approximation.
    - **Lossless water.** `alpha_coeff = 0` makes the Helmholtz operator near-singular (purely real wavenumber in the bulk), which stalls GMRES and underestimates peak pressure by ~18%. `run_jwave_sim.py` floors `alpha_coeff` to 0.002 dB/(cm·MHz²), the physical value for distilled water, before building the medium. The effect on the field is negligible (< 0.006 % attenuation over 100 mm at 500 kHz).
    - **Source-amplitude window.** The Helmholtz solver extracts the complex per-voxel drive amplitude via a single-frequency DFT of the CW source signal, restricted to `sensor.record_start_index:end`, the same steady-state window (last N periods, after the wavefront has crossed the domain) that k-Wave and jwave's time-domain solvers use for `p_max_all`/`p_min_all`. `createCWSignals()` applies a startup cosine up-ramp (default 4 periods) to the drive signal; a DFT over the full signal would average that lower-amplitude ramp in with the steady CW cycles and underestimate the drive amplitude. The steady-state window is built from an integer number of samples per period, which also removes most DFT scalloping loss.
    - **Sources.** Both source paths work: `grid.use_kWaveArray: 1` (band-limited off-grid source, as for k-Wave; recommended, and needed for flat elements) and `0` (`makeBowl`/`makeArc` voxels). The stair-stepped `makeBowl` bowl is ~9 % low in focal pressure in free water at 0.5 mm (ITRUSST BM1). The kWaveArray source reaches jwave as a `CWSourceMatrix` and is written to the HDF5 input in factored form (per-voxel complex amplitude × [cos; −sin] basis), since its BLI weights make every voxel's signal distinct.

#### Density term (`helmholtz_operator`)

The default `helmholtz_operator: conservative` solves ρ∇·((1/ρ)∇p) + k²p, with 1/ρ averaged onto the staggered points where ∇p is evaluated (as k-Wave staggers its density), and scales the source by the sound speed at each source voxel (as k-Wave does). `helmholtz_operator: jwave` uses jwave's built-in operator, which writes the density term as ∇²p − ∇ρ·∇p/ρ with a spectral gradient of ρ and scales the source by the minimum sound speed of the medium (`scale_source_helmholtz`). In homogeneous media the two agree.

??? note "Details: why the conservative form is the default"
    Across a sharp density jump (water/skull) jwave's built-in form does not enforce continuity of (1/ρ)∂p/∂n and loses energy at every interface. In 1D plane-wave tests against the analytic transmission it gives −19 % through a density-only plate and −26 % through the ITRUSST BM3 cortical plate (6.5 mm) at 0.5 mm, still −17 % at 0.125 mm. Attenuation and sound-speed contrasts alone are reproduced (≤ 3 %). In 3D this showed as PRESTUS jwave fields 23–29 % below the k-Wave reference on ITRUSST BM3–7. The conservative form gives −0.2 % (BM3 plate) and −3 % (BM4 skin + three-layer skull) at 0.5 mm in the same tests.

### `pseudo_spectral`

Pseudo-spectral symplectic time-domain integrator with the same equations as `kspaceFirstOrder*`. It runs `Nt` time steps with a PML absorbing boundary and records pressure over the last `Nt - record_start_index` steps (the steady-state window). It returns `p_max_all` (peak positive pressure), `p_min_all` (peak negative pressure) and `p_final` (pressure at the final recorded time step), all in Pa. Use it for direct equivalence with k-Wave output or when you need the transient field.

### `finite_differences`

Finite-difference time-domain solver. `jwave.fd_accuracy` sets the spatial stencil order (default 8).

## Configuration reference

| Parameter | Type | Default | Description |
|---|---|---|---|
| `solver` | string | `'helmholtz'` | Solver to use: `'helmholtz'` (default), `'pseudo_spectral'`, or `'finite_differences'` |
| `precision` | string | `'single'` | Floating-point precision: `'single'` (float32) or `'double'` (float64). Single is faster and sufficient for TUS. |
| `python_env` | string | `''` | Full path to the **bin directory** of the CPU jwave environment. Empty = system `python3`. |
| `python_env_gpu` | string | `''` | Full path to the **bin directory** of the GPU (CUDA) jwave environment. Empty = always use CPU env. |
| `use_gpu` | bool/empty | `[]` | GPU selection override. `[]` = auto-detect via `nvidia-smi`; `true` = force GPU env (errors if `python_env_gpu` unset); `false` = force CPU env. |
| `dry_run` | bool | `false` | If true, write HDF5 inputs but skip the simulation. The pipeline then stops after the acoustic stage (no analysis, thermal or report) and nothing is cached. Useful for inspecting the grid before running. |
| `cleanup_cache` | bool | `true` | Delete the HDF5 input (`*_jwave_input*.h5`, several GB for a head grid) from `io.dir_cache` after a successful run. It is rewritten on every run, so it is never reused. Dry runs and runs that error (e.g. Helmholtz residual above threshold) keep it. |
| `gmres_restart` | int | `10` | GMRES Krylov subspace size. Matches jwave's own `helmholtz_solver` default. Helmholtz solver only. |
| `gmres_maxiter` | int | `1000` | Maximum GMRES outer iterations. Matches jwave's own `helmholtz_solver` default. Helmholtz solver only. |
| `helmholtz_operator` | string | `'conservative'` | Density term: `'conservative'` or `'jwave'` (jwave's built-in). Helmholtz solver only. See [Density term](#density-term-helmholtz_operator). |
| `fd_accuracy` | int | `8` | Spatial stencil order. `finite_differences` solver only. |
| `convergence_threshold` | float | `0.05` | Relative spectral residual above which the result is rejected. Helmholtz solver only. By default this **errors**; see `ignore_convergence_threshold`. |
| `ignore_convergence_threshold` | bool | `false` | If true, exceeding `convergence_threshold` only warns (and returns the unconverged result) instead of erroring. Helmholtz solver only. |
| `grid.pml_size` | int/`'auto'` | `'auto'` | PML thickness in grid points (set under `grid:`, not `jwave:`, so both backends read the same field). `'auto'` resolves differently per backend: jwave uses a wavelength-based target of ~2.5λ, not k-Wave's FFT-prime-factorization heuristic (`resolve_jwave_pml_size.m`). |

??? note "Details: how `grid.pml_size: auto` is sized for jwave"
    With `grid.pml_size: auto`, jwave's PML is sized to about 2.5 wavelengths (clipped to [10, 40] grid points), in line with jwave's own example notebooks. The wavelength is taken at the **fastest** sound speed in the medium (`max(medium.sound_speed)`), not the slowest. jwave's Helmholtz PML (`complex_pml_on_grid`) is a fixed function of voxel position: it does not scale with frequency or with the local sound speed, so a given thickness absorbs fewer wavelengths where it borders a fast medium such as skull. On a head-model partial-FOV crop, a 20-voxel (15 mm) PML sized from water (1500 m/s) covered 3.0 wavelengths in water but only 1.8 in skull, and skull bordered most of the transverse boundary. Smoothing the skull/water interface did not reduce the Helmholtz-vs-k-Wave gap, which points to PML coverage rather than interface discontinuities. Sizing from the fastest medium gives every bordering tissue the full target thickness.

??? note "Details: output fields"
    `jwave_wrapper` returns a `sensor_data` struct with the same fields as `acoustic_simulation.m`:

    | Field | Size | Description |
    |---|---|---|
    | `p_max_all` | `[Nx Ny Nz]` double [Pa] | Peak positive pressure |
    | `p_min_all` | `[Nx Ny Nz]` double [Pa] | Peak negative (rarefaction) pressure |
    | `p_final` | `[Nx Ny Nz]` double [Pa] | Pressure at final time step (or Re(P) for Helmholtz) |
    | `jwave_input_h5` | string | Path to the HDF5 file written to `io.dir_cache`; empty once deleted by `cleanup_cache` |
    | `jwave_output_h5` | string | Path to the HDF5 file read from `io.dir_cache` |

??? note "Details: benchmark scripts"
    Two benchmark scripts are provided in `v05/code/`:

    - `benchmark_jwave.m`: staged debug script. Runs five progressive stages: Python environment check → HDF5 write → Python bridge → `jwave_wrapper` → full pipeline. Run stage by stage to isolate where a problem occurs.
    - `benchmark_jwave_pipeline.m`: single-run pipeline script, mirroring `benchmark_babelbrain_pipeline.m`. Sets `acoustic_backend = 'jwave'` and calls `prestus_pipeline_start` on the testdata subject.

??? note "Details: automatic differentiation"
    Because jwave runs in JAX, the simulation is fully differentiable. This is not exposed through the PRESTUS wrapper but is available by calling `run_jwave_sim.py` directly or extending it. All of the examples below use `jwave` and `jax` directly in Python; no MATLAB changes are needed.

    **Sensitivity map**: gradient of focal pressure with respect to every voxel's speed of sound. Cost: one reverse-mode pass, equivalent to a single forward simulation, regardless of grid size.

    ```python
    import jax
    from jwave import simulate_wave_propagation

    def focal_pressure(sound_speed):
        result = simulate_wave_propagation(domain, medium_with(sound_speed), source, ...)
        return result.pressure[focus_idx].real

    grad_fn = jax.grad(focal_pressure)
    sensitivity = grad_fn(medium.sound_speed)   # same shape as the grid
    ```

    **Phase optimisation**: optimise per-element drive phases to maximise focal pressure.

    ```python
    phases = jax.grad(lambda phi: -focal_pressure_from_phases(phi))(initial_phases)
    ```

    **Full-waveform inversion (FWI)**: estimate skull acoustic properties from measured transmissions.

    ```python
    loss = lambda sos: jnp.mean((simulate_pressure(sos) - measured_pressure) ** 2)
    sos_estimated = jax.scipy.optimize.minimize(loss, sos_init, method='BFGS')
    ```
