# jwave Configuration

jwave offers three solvers with different trade-offs between speed, memory, and accuracy. This page describes each solver and its convergence settings, and lists every `jwave` setting with its default.

On this page:

- [Solver Options](#solver-options): the Helmholtz, pseudo-spectral, and finite-difference solvers, convergence, and the density term.
- [Configuration Reference](#configuration-reference): every `jwave` setting with its default.

!!! note "At a glance"
    - **When you need it:** To choose a jwave solver or change its settings.
    - **Before you start:** jwave set up as in [jwave Acoustic Backend](doc_jwave.md).

## Solver Options

### Helmholtz Solver (Default)

Frequency-domain GMRES solver. It computes the complex pressure field `P(r)` at the drive frequency in a single linear solve instead of thousands of time steps, so it is much faster than `pseudo_spectral` for CW simulations, with equivalent accuracy for monochromatic excitation in lossy tissue. It returns `p_max_all` = |P| [Pa], `p_min_all` = −|P| and `p_final` = Re(P) [Pa]. GMRES convergence is controlled by `jwave.gmres_restart` and `jwave.gmres_maxiter`.

!!! note "Convergence is checked and enforced"
    After every Helmholtz solve, PRESTUS compares the relative GMRES residual (`/diagnostics/residual_norm`) against `jwave.convergence_threshold` (default `0.05`). If the residual exceeds the threshold, the pipeline errors, because an unconverged GMRES iterate is dominated by solver noise rather than physics. The error message gives ranked fixes: revert `gmres_restart`/`gmres_maxiter` to jwave's own defaults rather than enlarging them (jwave runs in single precision, where a larger restart can destabilise GMRES); check the PML fraction of the domain via `--debug`; only then loosen `gmres_tol`. If you have reviewed the field and diagnostics and accept the residual anyway, set `jwave.ignore_convergence_threshold: true` to downgrade the error to a warning and keep the result.

**How the Helmholtz residual, water attenuation, and source are handled**

- **Residual.** The residual is computed exactly, by re-evaluating the operator and source scaling of the solve on the solved field (`conservative_residual`, or `compute_helmholtz_residual` for `helmholtz_operator: jwave`, in `run_jwave_sim.py`), so it reflects the true solver error rather than an approximation.
- **Lossless water.** `alpha_coeff = 0` makes the Helmholtz operator near-singular (purely real wavenumber in the bulk), which stalls GMRES and underestimates peak pressure by ~18%. `run_jwave_sim.py` floors `alpha_coeff` to 0.002 dB/(cm·MHz²), the physical value for distilled water, before building the medium. The effect on the field is negligible (< 0.006 % attenuation over 100 mm at 500 kHz).
- **Source-amplitude window.** The Helmholtz solver extracts the complex per-voxel drive amplitude via a single-frequency DFT of the CW source signal, restricted to `sensor.record_start_index:end`, the same steady-state window (last N periods, after the wavefront has crossed the domain) that k-Wave and jwave's time-domain solvers use for `p_max_all`/`p_min_all`. `createCWSignals()` applies a startup cosine up-ramp (default 4 periods) to the drive signal; a DFT over the full signal would average that lower-amplitude ramp in with the steady CW cycles and underestimate the drive amplitude. The steady-state window is built from an integer number of samples per period, which also removes most DFT scalloping loss.
- **Sources.** Both source paths work: `grid.use_kWaveArray: 1` (band-limited off-grid source, as for k-Wave; recommended, and needed for flat elements) and `0` (`makeBowl`/`makeArc` voxels). The stair-stepped `makeBowl` bowl is ~9 % low in focal pressure in free water at 0.5 mm ([ITRUSST](https://itrusst.org/) BM1). The kWaveArray source reaches jwave as a `CWSourceMatrix` and is written to the HDF5 input in factored form (per-voxel complex amplitude × [cos; −sin] basis), since its BLI weights make every voxel's signal distinct.

#### Density Term (`helmholtz_operator`)

The default `helmholtz_operator: conservative` solves ρ∇·((1/ρ)∇p) + k²p, with 1/ρ averaged onto the staggered points where ∇p is evaluated (as k-Wave staggers its density), and scales the source by the sound speed at each source voxel (as k-Wave does). `helmholtz_operator: jwave` uses jwave's built-in operator, which writes the density term as ∇²p − ∇ρ·∇p/ρ with a spectral gradient of ρ and scales the source by the minimum sound speed of the medium (`scale_source_helmholtz`). In homogeneous media the two agree.

**Why the conservative form is the default**

Across a sharp density jump (water/skull) jwave's built-in form does not enforce continuity of (1/ρ)∂p/∂n and loses energy at every interface. In 1D plane-wave tests against the analytic transmission it gives −19 % through a density-only plate and −26 % through the ITRUSST BM3 cortical plate (6.5 mm) at 0.5 mm, still −17 % at 0.125 mm. Attenuation and sound-speed contrasts alone are reproduced (≤ 3 %). In 3D this showed as PRESTUS jwave fields 23–29 % below the k-Wave reference on ITRUSST BM3–7. The conservative form gives −0.2 % (BM3 plate) and −3 % (BM4 skin + three-layer skull) at 0.5 mm in the same tests.

### Pseudo-Spectral Solver

Pseudo-spectral symplectic time-domain integrator with the same equations as `kspaceFirstOrder*`. It runs `Nt` time steps with a [PML](doc_glossary.md#pml) absorbing boundary and records pressure over the last `Nt - record_start_index` steps (the steady-state window). It returns `p_max_all` (peak positive pressure), `p_min_all` (peak negative pressure), and `p_final` (pressure at the final recorded time step), all in Pa. Use it for direct equivalence with k-Wave output or when you need the transient field.

### Finite-Difference Solver

Finite-difference time-domain solver. `jwave.fd_accuracy` sets the spatial stencil order (default 8).

## Configuration Reference

The `jwave` keys (solver, environments, GPU selection, GMRES settings, convergence threshold, and cache cleanup), with their defaults, are listed under [jwave settings](doc_parameters_simulation.md#jwave-settings-jwave). The PML thickness is set by `grid.pml_size`, shared with k-Wave.

### PML Size

With `grid.pml_size: auto`, jwave's PML is sized to about 2.5 wavelengths (clipped to [10, 40] grid points), in line with jwave's own example notebooks. The wavelength is taken at the **fastest** sound speed in the medium (`max(medium.sound_speed)`), not the slowest. jwave's Helmholtz PML (`complex_pml_on_grid`) is a fixed function of voxel position: it does not scale with frequency or with the local sound speed, so a given thickness absorbs fewer wavelengths where it borders a fast medium such as skull. On a head-model partial-FOV crop, a 20-voxel (15 mm) PML sized from water (1500 m/s) covered 3.0 wavelengths in water but only 1.8 in skull, and skull bordered most of the transverse boundary. Smoothing the skull/water interface did not reduce the Helmholtz-vs-k-Wave gap, which points to PML coverage rather than interface discontinuities. Sizing from the fastest medium gives every bordering tissue the full target thickness.

**Output fields**

`jwave_wrapper` returns a `sensor_data` struct with the same fields as `acoustic_simulation.m`:

| Field | Size | Description |
|---|---|---|
| `p_max_all` | `[Nx Ny Nz]` double [Pa] | Peak positive pressure |
| `p_min_all` | `[Nx Ny Nz]` double [Pa] | Peak negative (rarefaction) pressure |
| `p_final` | `[Nx Ny Nz]` double [Pa] | Pressure at final time step (or Re(P) for Helmholtz) |
| `jwave_input_h5` | string | Path to the HDF5 file written to `io.dir_cache`; empty once deleted by `cleanup_cache` |
| `jwave_output_h5` | string | Path to the HDF5 file read from `io.dir_cache` |

**Checking the jwave environment**

`smoke_test_backends('jwave')` (in `tests/`) checks in under 30 s that jwave and its dependencies import and that the PRESTUS bridge scripts are present, without running a simulation. Run it after changing the Python environment.

**Automatic differentiation**

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
