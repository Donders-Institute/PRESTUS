# jwave acoustic backend

!!! note "At a glance"
    - **What:** jwave as an alternative acoustic solver.
    - **When you need it:** Optional: to compare solvers.
    - **Before you start:** A Python environment with jwave and JAX.

PRESTUS can use [jwave](https://github.com/ucl-bug/jwave) as an alternative backend for acoustic simulation. By default PRESTUS uses k-Wave; the jwave backend can be enabled independently. All pre-processing stages (segmentation, grid setup, medium masks, transducer placement) and all downstream post-processing (focal metrics, NIfTI export, HTML report) run unchanged regardless of which backend is active.

**What is jwave?** jwave is an open-source Python library for differentiable acoustic wave simulation developed at University College London. It implements the same pseudo-spectral k-Wave equations but in [JAX](https://github.com/google/jax), giving access to JIT compilation, automatic differentiation, and NVIDIA GPU acceleration without requiring compiled MATLAB/C++ binaries. jwave provides two solvers: a time-domain pseudo-spectral integrator equivalent to `kspaceFirstOrder*`, and a frequency-domain Helmholtz solver that computes the CW steady-state directly without running a full transient simulation.

The main reasons to use jwave instead of (or alongside) k-Wave are:

- **Helmholtz solver for CW simulations.** For continuous-wave TUS the frequency-domain solve is exact and significantly faster than running a full time-domain simulation and extracting amplitude/phase at the end.
- **Automatic differentiation.** JAX's `jax.grad` can differentiate through the entire wave simulation. This enables gradient-based focal optimisation, sensitivity maps (`∂P_focus/∂SoS_voxel` for all voxels in a single reverse-mode pass), and full-waveform inversion (FWI) for skull property estimation — without deriving or implementing any adjoint manually.
- **Cross-validation.** Comparing results from an independent Python/JAX implementation against MATLAB k-Wave is a useful consistency check.

| | PRESTUS (k-Wave) | jwave |
|---|---|---|
| **Solver** | Pseudo-spectral FDTD | Pseudo-spectral (time-domain) or Helmholtz GMRES (frequency-domain) |
| **Physics** | Scalar linear acoustic | Scalar linear acoustic |
| **Skull model** | Scalar ρ, c, α | Scalar ρ, c, α |
| **Differentiable** | No | Yes (JAX automatic differentiation) |
| **GPU backends** | CUDA (via MATLAB `gpuArray` or C++ binary) | NVIDIA CUDA only |
| **Thermal simulation** | Yes (`kWaveDiffusion`) | Not supported |

See [Simulation backends](doc_backend.md) for the k-Wave-only backend options.

> **Note — jwave uses the same scalar acoustic model as k-Wave.** Unlike BabelBrain, jwave does not model shear waves in bone or include a Rayleigh pre-step. For skull simulations the material properties assigned by PRESTUS (`medium_properties`) apply to both k-Wave and jwave equally, so results are directly comparable.

---

## Limitations

- **Acoustic only** — jwave has no bio-heat solver. Thermal simulation (`run_heating_sims`) must still use k-Wave or BabelBrain.
- **No axisymmetric grids** — jwave requires a full 3-D domain.
- **GPU: NVIDIA CUDA only** — JAX does not officially support AMD ROCm or Apple Metal. CPU execution works on any platform.

---

## Python environment setup

jwave requires a dedicated Python environment. Creating a separate env avoids version conflicts with other PRESTUS Python tools (e.g. BabelBrain, SimNIBS) and makes the dependency set explicit.

### 1. Create the environment

**conda:**

```bash
conda create -n jwave_env python=3.11
conda activate jwave_env
```

**venv (no conda):**

```bash
python3.11 -m venv ~/.venvs/jwave_env
source ~/.venvs/jwave_env/bin/activate
```

### 2. Install jwave and its dependencies

`jax[cpu]` and `jax[cuda12]` install different builds of `jaxlib` and are **mutually exclusive** — installing one overwrites the other. Pick the variant that matches your hardware, or maintain two parallel environments (see below).

**CPU only (macOS, Windows, or CPU-only Linux):**

```bash
pip install "jax[cpu]" jwave h5py numpy
```

**NVIDIA GPU (CUDA 12, Linux x86-64 only):**

```bash
pip install "jax[cuda12]" jwave h5py numpy
```

A CUDA install falls back to CPU automatically when no GPU is present at runtime, so it works on CPU-only hosts too — but the `jaxlib-cuda` wheel is large (~2 GB) and Linux-only. On macOS or Windows, always use `jax[cpu]`.

> JAX CUDA wheels are published for Linux x86-64 only. GPU acceleration requires a Linux host with CUDA 12 drivers.

#### Parallel CPU and GPU environments

If you work on both CPU and GPU machines (e.g. a laptop and an HPC cluster with shared home directories), maintain two environments and point `parameters.jwave.python_env` to the appropriate one:

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

Then in your config, set `python_env` to whichever environment is appropriate for the current machine:

```yaml
jwave:
  python_env: '/home/user/.conda/envs/jwave_env/bin/python'       # CPU
  # python_env: '/home/user/.conda/envs/jwave_env_gpu/bin/python'  # GPU
```

### 3. Verify the installation

```bash
python -c "import jwave, jax, h5py; print('jwave:', jwave.__version__); print('jax:', jax.__version__)"
```

For GPU, confirm the device is visible:

```bash
python -c "import jax; print(jax.devices())"
# Expected on GPU host: [CudaDevice(id=0)]
# Expected on CPU host: [CpuDevice(id=0)]
```

### 4. Find the Python executable path

```bash
which python   # inside the activated environment
# e.g. /home/user/.conda/envs/jwave_env/bin/python
#      /home/user/.venvs/jwave_env/bin/python
```

Copy this path — it is what you pass to `parameters.jwave.python_env` (see [Configuration reference](#configuration-reference) below).

---

## Quick start

```yaml
acoustic_backend: jwave

jwave:
  solver: helmholtz          # 'helmholtz' (default) | 'pseudo_spectral' | 'finite_differences'
  precision: single          # 'single' | 'double'
  python_env: ''             # full path to CPU Python bin dir (see environment setup)
  python_env_gpu: ''         # full path to GPU (CUDA) Python bin dir; leave empty for CPU-only
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

`python_env` is the **full path to the Python bin directory** of the CPU environment. `python_env_gpu` points to the GPU (CUDA) environment. At runtime `jwave_wrapper` selects the environment as follows:

| `use_gpu` | GPU detected (`nvidia-smi`) | `python_env_gpu` set | Environment used |
|---|---|---|---|
| `[]` (auto) | yes | yes | GPU |
| `[]` (auto) | yes | no | CPU (warning) |
| `[]` (auto) | no | — | CPU |
| `true` | — | yes | GPU |
| `true` | — | no | **error** |
| `false` | — | — | CPU |

Leaving both fields empty falls back to the system `python3`, which is unlikely to have jwave installed.

---

## Solver options

### `helmholtz` (default)

Frequency-domain GMRES solver. Computes the complex pressure field `P(r)` at the drive frequency directly, without running a time-domain transient. Returns:

- `p_max_all` = |P| (pressure amplitude) [Pa]
- `p_min_all` = −|P|
- `p_final` = Re(P) (real part at phase zero) [Pa]

Significantly faster than `pseudo_spectral` for CW simulations because the steady-state is reached in a single linear solve rather than thousands of time steps. Equivalent accuracy for monochromatic excitation in lossy tissue media. GMRES convergence is controlled by `jwave.gmres_restart` and `jwave.gmres_maxiter`.

> **Convergence is checked and enforced.** After every Helmholtz solve, `jwave_wrapper` compares the relative GMRES residual (see `/diagnostics/residual_norm`) against `jwave.convergence_threshold` (default `0.05`). This residual is computed exactly, by re-evaluating the operator and source scaling of the solve on the solved field (`conservative_residual`, or `compute_helmholtz_residual` for `helmholtz_operator: jwave`, in `run_jwave_sim.py`), so it reflects the true solver error rather than an approximation. If the residual exceeds the threshold, **the pipeline errors**, since a genuinely unconverged GMRES iterate is dominated by solver noise rather than physics. The error message gives concrete, ranked fixes (revert `gmres_restart`/`gmres_maxiter` to jwave's own defaults rather than enlarging them — a larger restart in single precision can destabilize GMRES instead of helping; try `precision: double`; check the PML fraction of the domain via `--debug`; only then loosen `gmres_tol`). If you have reviewed the field and diagnostics and accept the residual anyway, set `jwave.ignore_convergence_threshold: true` to downgrade this back to a warning and keep the result.

> **Lossless water note.** Setting `alpha_coeff = 0` makes the Helmholtz operator near-singular (purely real wavenumber in the bulk), which causes GMRES to stall and underestimate peak pressure by ~18%. `run_jwave_sim.py` automatically floors `alpha_coeff` to 0.002 dB/(cm·MHz²) — the physical value for distilled water — before building the medium. This regularises GMRES convergence with negligible effect on the acoustic field (< 0.006 % attenuation over 100 mm at 500 kHz).

> **Source-amplitude window.** The Helmholtz solver extracts the complex per-voxel
> drive amplitude via a single-frequency DFT of the CW source signal. That DFT is
> restricted to `sensor.record_start_index:end` — the same steady-state window
> (last N periods, after the wavefront has crossed the domain) that k-Wave and
> jwave's own time-domain solvers already use for `p_max_all`/`p_min_all`. This
> matters because `createCWSignals()` applies a startup cosine up-ramp (default
> 4 periods) to the drive signal; DFT-ing the full signal averages that
> lower-amplitude ramp in with the steady CW cycles and systematically
> underestimates the drive amplitude (the time-domain solvers don't have this
> problem — they peak-track with `max`/`min` rather than average). Restricting
> to the steady-state window also removes most DFT scalloping loss "for free",
> since that window is built from an integer number of samples per period.

#### Density term (`helmholtz_operator`)

jwave's own Helmholtz operator writes the density term as ∇²p − ∇ρ·∇p/ρ, with a spectral gradient of ρ. Across a sharp density jump (water/skull) this form does not enforce continuity of (1/ρ)∂p/∂n and loses energy at every interface. In 1D plane-wave tests against the analytic transmission it gives −19 % through a density-only plate and −26 % through the ITRUSST BM3 cortical plate (6.5 mm) at 0.5 mm, still −17 % at 0.125 mm. Attenuation and sound-speed contrasts alone are reproduced (≤ 3 %). In 3D this showed as PRESTUS jwave fields 23–29 % below the k-Wave reference on ITRUSST BM3–7.

The default `helmholtz_operator: conservative` solves ρ∇·((1/ρ)∇p) + k²p instead, with 1/ρ averaged onto the staggered points where ∇p is evaluated (as k-Wave staggers its density), and scales the source by the sound speed at each source voxel (as k-Wave does; jwave's `scale_source_helmholtz` uses the minimum sound speed of the medium). Same tests: −0.2 % (BM3 plate) and −3 % (BM4 skin + three-layer skull) at 0.5 mm. In homogeneous media the two operators agree. `helmholtz_operator: jwave` reproduces earlier PRESTUS results.

> **Sources.** Both source paths work: `grid.use_kWaveArray: 1` (band-limited off-grid source, as for k-Wave; recommended, and needed for flat elements) and `0` (`makeBowl`/`makeArc` voxels). The stair-stepped `makeBowl` bowl is ~9 % low in focal pressure in free water at 0.5 mm (ITRUSST BM1). The kWaveArray source reaches jwave as a `CWSourceMatrix` and is written to the HDF5 input in factored form (per-voxel complex amplitude × [cos; −sin] basis), since its BLI weights make every voxel's signal distinct.

### `pseudo_spectral`

Pseudo-spectral symplectic time-domain integrator — the same equations as `kspaceFirstOrder*`. Runs for `Nt` time steps with a PML absorbing boundary, records pressure over the last `Nt - record_start_index` steps (the steady-state window), and returns:

- `p_max_all` — peak positive pressure [Pa]
- `p_min_all` — peak negative pressure [Pa]
- `p_final` — pressure field at the final recorded time step [Pa]

Use this solver when you want direct equivalence with k-Wave output or need the transient pressure field.

### `finite_differences`

Finite-difference time-domain solver. Use `jwave.fd_accuracy` to set the spatial stencil order (default 8).

---

## Configuration reference

| Parameter | Type | Default | Description |
|---|---|---|---|
| `solver` | string | `'helmholtz'` | Solver to use: `'helmholtz'` (default), `'pseudo_spectral'`, or `'finite_differences'` |
| `precision` | string | `'single'` | Floating-point precision: `'single'` (float32) or `'double'` (float64). Single is faster and sufficient for TUS. |
| `python_env` | string | `''` | Full path to the **bin directory** of the CPU jwave environment. Empty = system `python3`. |
| `python_env_gpu` | string | `''` | Full path to the **bin directory** of the GPU (CUDA) jwave environment. Empty = always use CPU env. |
| `use_gpu` | bool/empty | `[]` | GPU selection override. `[]` = auto-detect via `nvidia-smi`; `true` = force GPU env (errors if `python_env_gpu` unset); `false` = force CPU env. |
| `dry_run` | bool | `false` | If true, write HDF5 inputs but skip the simulation. The pipeline then stops after the acoustic stage (no analysis, thermal or report) and nothing is cached. Useful for inspecting the grid before running. |
| `gmres_restart` | int | `10` | GMRES Krylov subspace size. Matches jwave's own `helmholtz_solver` default. Helmholtz solver only. |
| `gmres_maxiter` | int | `1000` | Maximum GMRES outer iterations. Matches jwave's own `helmholtz_solver` default. Helmholtz solver only. |
| `helmholtz_operator` | string | `'conservative'` | Density term: `'conservative'` or `'jwave'` (jwave's built-in). Helmholtz solver only. See [Density term](#density-term-helmholtz_operator). |
| `fd_accuracy` | int | `8` | Spatial stencil order. `finite_differences` solver only. |
| `convergence_threshold` | float | `0.05` | Relative spectral residual above which the result is rejected. Helmholtz solver only. By default this **errors**; see `ignore_convergence_threshold`. |
| `ignore_convergence_threshold` | bool | `false` | If true, exceeding `convergence_threshold` only warns (and returns the unconverged result) instead of erroring. Helmholtz solver only. |
| `grid.pml_size` | int/`'auto'` | `'auto'` | PML thickness in grid points (set under `grid:`, not `jwave:`, so both backends read the same field — but `'auto'` resolves *differently* per backend: jwave uses a wavelength-based target of ~2.5λ, not k-Wave's FFT-prime-factorization heuristic; see `resolve_jwave_pml_size.m`). |

### PML thickness

With `grid.pml_size: auto`, jwave's PML is sized to about 2.5 wavelengths (clipped to
[10, 40] grid points), in line with jwave's own example notebooks. The wavelength is taken at
the **fastest** sound speed in the medium (`max(medium.sound_speed)`), not the slowest. jwave's
Helmholtz PML (`complex_pml_on_grid`) is a fixed function of voxel position: it does not scale
with frequency or with the local sound speed, so a given thickness absorbs fewer wavelengths
where it borders a fast medium such as skull. On a head-model partial-FOV crop, a 20-voxel
(15 mm) PML sized from water (1500 m/s) covered 3.0 wavelengths in water but only 1.8 in skull,
and skull bordered most of the transverse boundary. Smoothing the skull/water interface did not
reduce the Helmholtz-vs-k-Wave gap, which points to PML coverage rather than interface
discontinuities. Sizing from the fastest medium gives every bordering tissue the full target
thickness.

---

## Output fields

`jwave_wrapper` returns a `sensor_data` struct with the same fields as `acoustic_simulation.m`, so all downstream PRESTUS modules work unchanged:

| Field | Size | Description |
|---|---|---|
| `p_max_all` | `[Nx Ny Nz]` double [Pa] | Peak positive pressure |
| `p_min_all` | `[Nx Ny Nz]` double [Pa] | Peak negative (rarefaction) pressure |
| `p_final` | `[Nx Ny Nz]` double [Pa] | Pressure at final time step (or Re(P) for Helmholtz) |
| `jwave_input_h5` | string | Path to the HDF5 file written to `io.dir_cache` |
| `jwave_output_h5` | string | Path to the HDF5 file read from `io.dir_cache` |

---

## Benchmarking

Two benchmark scripts are provided in `v05/code/`:

- **`benchmark_jwave.m`** — staged debug script. Runs five progressive stages: Python environment check → HDF5 write → Python bridge → `jwave_wrapper` → full pipeline. Run stage by stage to isolate where a problem occurs.
- **`benchmark_jwave_pipeline.m`** — clean single-run pipeline script. Mirrors `benchmark_babelbrain_pipeline.m` in structure. Sets `acoustic_backend = 'jwave'` and calls `prestus_pipeline_start` on the testdata subject.

A `matlab.unittest.TestCase` class at `tests/benchmark_jwave.m` runs correctness checks (positive pressure, required output fields, spatial correlation with k-Wave ≥ 0.95) and a timing sweep across grid sizes.

---

## Advanced: automatic differentiation

Because jwave runs in JAX, the simulation is fully differentiable. This is not directly exposed through the PRESTUS wrapper but is available by calling `run_jwave_sim.py` directly or extending it.

**Sensitivity map** — gradient of focal pressure with respect to every voxel's speed of sound:

```python
import jax
from jwave import simulate_wave_propagation

def focal_pressure(sound_speed):
    result = simulate_wave_propagation(domain, medium_with(sound_speed), source, ...)
    return result.pressure[focus_idx].real

grad_fn = jax.grad(focal_pressure)
sensitivity = grad_fn(medium.sound_speed)   # same shape as the grid
```

Cost: one reverse-mode pass — equivalent to a single forward simulation, regardless of grid size.

**Phase optimisation** — optimise per-element drive phases to maximise focal pressure:

```python
phases = jax.grad(lambda phi: -focal_pressure_from_phases(phi))(initial_phases)
```

**Full-waveform inversion (FWI)** — estimate skull acoustic properties from measured transmissions:

```python
loss = lambda sos: jnp.mean((simulate_pressure(sos) - measured_pressure) ** 2)
sos_estimated = jax.scipy.optimize.minimize(loss, sos_init, method='BFGS')
```

All of the above require `jwave` and `jax` directly in Python; no MATLAB changes are needed.
