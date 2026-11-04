# Acoustic simulations

!!! note "At a glance"
    - **What:** How the acoustic simulation is set up: grid sampling, intensity targeting, memory use and 2D axisymmetric runs.
    - **When you need it:** Read [ISPPA targeting](#isppa-targeting) for every study. The rest matters when you change grid settings or run out of memory.
    - **Before you start:** Nothing; the defaults are a sound starting point.

## What to set

The defaults run a 3D simulation with k-Wave's `kspaceFirstOrder3D` (see the available [simulation backends](doc_backend.md)). The settings most studies touch:

```yaml
grid:
  resolution_mm: 0.5        # spatial step; keep PPW >= 6 in the slowest medium
  source_cfl: 0.15          # time step relative to the stability limit
transducer:
  target_isppa_wcm2: 30     # free-water ISPPA the results are reported at
modules:
  run_water_baseline: auto  # free-water baseline for the ISPPA scaling
```

The simulation models compressional waves (P-waves) in isotropic media and computes pressure and velocity fields. Shear waves are not modelled explicitly; their scattering contributes to the bulk attenuation coefficient. This makes the run faster and less memory-intensive than an elastic simulation. An elastic skull model is available as a separate backend ([`kwave_elastic`](doc_kwave.md#elastic-skull-model-shear-waves), draft).

!!! note "Turning off acoustic simulations"
    With `modules.run_acoustic_sims: 0`, PRESTUS loads the existing acoustic results from the cache (subject to `io.overwrite_files`; see [Overwrite behaviour](doc_outputs.md#overwrite-behaviour)). This lets a follow-up thermal simulation reuse an earlier acoustic run. Other modules behave differently: turning them off may skip the step without loading existing data.

## Spatial and temporal sampling

Accuracy depends on resolving the shortest wavelength in the domain. Two parameters govern this.

**Grid resolution (`grid.resolution_mm`)** sets the spatial step `dx`. With the medium sound speed `c` and the fundamental frequency `f`, the number of spatial samples per wavelength is

```
PPW = c / (f × dx)
```

The worst case is the slowest medium (water or brain, ~1500 m/s). Keep **PPW ≥ 6**; PRESTUS warns below `grid.min_ppw`.

**CFL number (`grid.source_cfl`)** sets the time step `dt` relative to `dx` and the fastest medium (skull, ~2800 m/s), which sets the stability limit:

```
dt  = CFL × dx / c_max
PPP = PPW / CFL
```

`PPP` is the number of time steps per wave period. `resolution_mm` and `source_cfl` are independent, so you can refine or coarsen the time step without changing the grid. Both must stay within the stability limit.

Rule of thumb: the coarsest grid that satisfies a minimum PPW at frequency `f` in the slowest medium is `dx_max = c_min / (f × PPW_min)`. At 500 kHz in water (1500 m/s) with PPW = 6: `resolution_mm = 0.50 mm`.

??? note "Details: PPW across tissues, the stability recheck and the CFL default"
    **PPW varies across tissues.** The check against `grid.min_ppw` uses the slowest medium. At 500 kHz with 0.5 mm resolution, water/brain (1500 m/s) gives PPW ≈ 10 and skull (2800 m/s) PPW ≈ 37. The binding spatial constraint is always in the slow tissues.

    **The stability recheck can reduce `dt` below the CFL target.** After source setup, `prestus_pipeline` calls `checkStability`. If the CFL-derived `dt` exceeds the estimated stability limit for the heterogeneous medium, source setup is rerun with a smaller time step scaled by `grid.source_limit_fraction` (default 0.9). The final `dt` can therefore be smaller than `CFL × dx / c_max`. The PPW check runs only on the first source setup call; the adjusted call receives `dt` explicitly and bypasses it.

    **Overriding `source_ppw` decouples temporal from spatial resolution.** A manual `grid.source_ppw` only affects the time step (`PPP = source_ppw / CFL`). The spatial check against `grid.min_ppw` still uses the actual PPW from `c_min / (f × dx)`. A coarser `source_ppw` increases `dt` and reduces runtime without suppressing the spatial warning.

    **The default CFL of 0.15 is conservative.** k-Wave's own default is 0.3. PRESTUS uses 0.15 for extra stability margin in heterogeneous skull simulations, where the theoretical stability limit is an approximation ([`checkStability`](http://www.k-wave.org/documentation/checkStability.php)). This roughly doubles the number of time steps and the acoustic runtime. If simulations are stable, `source_cfl` can be raised towards 0.3; validate this against the `checkStability` output for your medium.

### Simulation duration (`t_end`)

With the default `grid.steady_state_transits: 1`, the simulation lasts long enough for the wave to cross the domain diagonal once, so the steady-state field is established everywhere:

```
t_end   = steady_state_transits × grid_diagonal / c_min
N_steps = round(t_end / dt)
```

where `grid_diagonal = sqrt(x_size² + y_size² + z_size²)` and `c_min` is the minimum sound speed across all media (water or brain, ~1500 m/s). Request more transits for compact, reflective domains, or set a fixed duration with `grid.sim_time_us`. The sensor records `p_max_all` and `p_final` over the last `grid.steady_state_num_periods` (default 3) wave periods.

## Amplitude ramping

A ramped pulse delivers less pulse-average intensity (ISPPA) than the temporal-peak intensity (ISPTP), because intensity scales with amplitude². When `timing.pulse_ramp_dur > 0`, `acoustic_analysis` scales the `acoustic_Ipa` map and every derived ISPPA value (per tissue, at target) by

```
f = ISPPA / ISPTP
```

before further analysis. With `pulse_ramp_dur` absent or zero, `f = 1`. The factor can also be computed directly with `compute_ramp_factor(pd, ramp_dur, shape)`:

| Parameter | Description |
|---|---|
| `pd` | Pulse duration [s] |
| `ramp_dur` | Ramp duration per edge [s]; applied at onset and offset, so the ramped interval per pulse is `2 × ramp_dur` |
| `shape` | `'linear'`: linear amplitude taper (intensity integral factor k = 1/3); `'tukey'`: raised-cosine taper (k = 3/8) |

Example: 20 ms pulse, 2 ms ramp per edge (10 % of `pd` per edge, 20 % total):

| Shape | k | f (ISPPA/ISPTP) |
|---|---|---|
| linear | 1/3 | ≈ 0.867 |
| tukey | 3/8 | ≈ 0.875 |

??? note "Details: ramp factor formula"
    ```
    t_flat = pd - 2 × ramp_dur
    f      = (k × ramp_dur + t_flat + k × ramp_dur) / pd
    ```

    Ramping reduces spectral splatter by smoothing the transducer drive at pulse onset and offset instead of switching abruptly.

## ISPPA targeting

The acoustic simulation runs at the drive `elem_p_amp` (from the study config or the transducer library). `transducer.target_isppa_wcm2` (default `30` W/cm²) sets the **free-water** ISPPA the results are reported at. Pressure is linear in drive, so PRESTUS scales the simulated field by

```
scale_p = sqrt(target_isppa / baseline_isppa)
```

where `baseline_isppa` is the free-water peak ISPPA of the same drive. All pressure fields used downstream (`p_max_all`, `p_min_all`, `p_final`) are scaled, so intensities, MI and heating are consistent. A scalar target scales one run; a vector runs [one thermal job per value](doc_multiintensity.md); an empty target (`[]`) reports the field at `elem_p_amp`.

The baseline comes from one of two sources (`isppa_baseline`):

1. a **free-water simulation** with the same grid, source and backend (`acoustic_water_baseline`; BabelBrain runs its own), recorded as `acoustic_provenance.freefield_isppa_wcm2` with its origin in `freefield_source`; or
2. the transducer's **calibration ISPPA** (`tx_calibration.isppa_wcm2`, from the library or a calibration file), used only when it matches the target to within 1 %. No scaling is then needed.

`modules.run_water_baseline` decides when the free-water simulation runs (`needs_water_baseline`):

| Value | Runs the baseline when |
|---|---|
| `auto` (default) | a target is set and the transducer has no calibration ISPPA, or the target differs from it by more than 1 % |
| `1` | always |
| `0` | only when the target differs from a calibration ISPPA |

A target is never reported unscaled. When a rescale is needed but no simulated baseline exists (for example `run_water_baseline: 0` with an uncalibrated transducer, or a cached result from a run without baseline and `run_source_setup: 0`), the run stops with `prestus:isppa:noBaseline` or `prestus:isppa:unsimulatedBaseline` (`isppa_scale_factor`). Scale factors above 4 or below 0.25 warn (`prestus:isppa:largeScale`): check that the target is in the linear regime of the simulated drive.

**Check what a run models.** The parameter summary at the start of the log and the HTML report (row *Modelled intensity*) state the free-water ISPPA the results report and how it relates to the calibration ISPPA of the drive. A target that differs from the calibration ISPPA, or the default target on a drive without calibration ISPPA, is marked ⚠️ and warns `prestus:isppa:targetVsCalibration`. A calibration file (`<equipment>-F<focus>mm-I<isppa>wpercm2.yaml`, [calibration](doc_calibration.md)) sets both `target_isppa_wcm2` and `tx_calibration.isppa_wcm2` to its calibrated intensity. If you copy only its `annular` drive into a config, also set `target_isppa_wcm2` to that intensity; otherwise the results are rescaled to the default 30 W/cm².

**Several transducers.** `target_isppa_wcm2` is set per transducer and refers to that transducer alone. Coherent runs set each transducer's drive to its own target before the simulation; async runs scale each field before combining ([Multi-transducer: ISPPA targets](doc_multitransducer.md#isppa-targets)). The older `calibration.target_isppa_wcm2` is read and moved to `transducer(1).target_isppa_wcm2` with a deprecation warning.

??? note "Details: axial profile plot"
    The baseline run also stores the free-water intensity along each transducer's beam axis (bowl apex to grid edge) as `acoustic_provenance.freefield_axial`. The axial profile plot (`sub-XXX_<medium>_axial_profile<affix>.png` / `.csv`, 3D grids) compares the simulated intensity with the free-water intensity of the simulated transducer settings: the k-Wave baseline profile when the baseline ran, and the analytical O'Neil profile for a single annular transducer (always), both at the drive level of the plotted field after scaling. The plot is omitted for combined async / FOV fields (scaled per source) and for the BabelBrain backend.

## Memory and source signal optimisation

Two components dominate acoustic memory:

- **Grid fields**: pressure, velocity and medium maps, proportional to the number of voxels. A 0.5 mm head grid takes roughly 7–8 GB; this scales with resolution, not transducer size.
- **Source signal matrix**: one waveform per active source voxel over all time steps. Negligible for a focused single-element transducer; for a large multi-element phased array it can exceed the grid fields.

PRESTUS compresses the source matrix automatically per backend (details below). A coarser grid or a [limited acoustic FOV](doc_limited_fov.md) reduces the grid fields; the `cpp_*` backends stream the source matrix from disk.

??? note "Details: how PRESTUS reduces source memory by backend"
    The source signal matrix has size `num_source_voxels × Nt`. For a kHz-range aperture at 0.5 mm resolution with the standard CFL of 0.15, `Nt` is on the order of several thousand, so the dense matrix reaches several gigabytes for large transducers.

    **`matlab_cpu` / `matlab_gpu`: lazy rank-2 factorisation for CW sources.** For continuous-wave (CW) stimulation every source voxel's time signal has the form

    ```
    p(voxel, t) = A_voxel · cos(ω t + φ_voxel)
    ```

    which two real numbers per voxel describe (amplitude and phase, or the real and imaginary part of a complex amplitude). PRESTUS stores only those two numbers per voxel (`CWSourceMatrix.m`) and evaluates each time column on demand during the simulation. Memory scales with `num_source_voxels` only, a reduction by a factor of ~2 000–5 000 for typical parameters. This path is used when `code_type` is `matlab_cpu` or `matlab_gpu` and all transducer elements share a single frequency; otherwise the dense matrix is used. k-Wave's own source scaling (`kspaceFirstOrder_scaleSourceTerms`), including the per-voxel loop for heterogeneous media, is applied to the stored amplitudes, so no dense matrix is built at any point.

    **`cpp_cpu` / `cpp_gpu`: HDF5 streaming from disk.** The source matrix is written to an HDF5 file before the simulation and the C++ binary reads it in blocks during execution. Peak memory is bounded by the block size rather than the full matrix, at the cost of I/O throughput.

    | Backend | Simulation device | Source signal representation | Peak additional RAM / VRAM |
    |---|---|---|---|
    | `matlab_gpu` | GPU | Lazy rank-2 (CW) or dense single | `2 × (num_src + Nt) × 4 B` (CW) |
    | `matlab_cpu` | CPU | Lazy rank-2 (CW) or dense single | `2 × (num_src + Nt) × 4 B` (CW) |
    | `cpp_cpu` / `cpp_gpu` | CPU / GPU | HDF5 file, streamed | Bounded by I/O block size |

    For a 10 000-voxel source at Nt = 9 000 time steps: dense = **360 MB**; CW lazy = **0.08 MB**.

## Axisymmetric simulations

For symmetric media (a 2D `phantom` or `water` medium), a 3D axisymmetric problem can be modelled on a 2D grid with [`kspaceFirstOrderAS`](http://www.k-wave.org/documentation/kspaceFirstOrderAS.php) (`grid.axisymmetric: 1`; see k-Wave's [focused bowl example](https://github.com/ucl-bug/k-wave/blob/main/k-Wave/examples/example_at_focused_bowl_AS.m)). This reduces the computational cost markedly compared with a full 3D simulation. [`kspaceFirstOrder2D`](http://www.k-wave.org/documentation/kspaceFirstOrder2D.php) is not a substitute: it models an infinite cylinder, not a focusing shell.

What applies in axisymmetric mode:

- Specify the grid as `[radial x 2, axial]`. The (2x) radial axis is assumed to be shorter than the axial axis.
- Only a **single transducer** is supported; several transducers with `axisymmetric == 1` raise an error.
- Only **annular transducers** explicitly support axisymmetric `kWaveArray` mode. Other transducer types are untested in this configuration.
- A 3D `phantom` is cut along its second dimension at the transducer's `trans_pos` when that is given in 3D, else at the centre slice. Output NIfTIs of `phantom` runs are written in the phantom file's voxel order and world space.
- The source keeps every grid point with a non-zero band-limited-interpolation weight, as the 3D path and k-Wave do.

!!! note "Axisymmetric 2D vs. 3D output"
    Axisymmetric simulations are 3D simulations performed on 2D input grids. When heating is requested in the pipeline, the acoustic outputs are provided in 3D, because PRESTUS assumes a 3D heating simulation follows. For 2D outputs, deactivate heating. A 2D heating simulation is possible by first running the acoustic 2D simulation without heating, then running heating with that output while the acoustic simulation is deactivated (or not overwritten).

### Backend compatibility

| Backend (`simulation.code_type`) | Axisymmetry support | Notes |
|---|---|---|
| `matlab_cpu` | ✓ | Uses `kspaceFirstOrderAS` |
| `matlab_gpu` | ✓ | Uses `kspaceFirstOrderAS` with `RadialSymmetry` option |
| `cpp_cpu` | ✓ | Uses the C++ OMP binary `kspaceFirstOrderASC` |
| `cpp_gpu` (CUDA) | ✗ | No CUDA binary supports axisymmetric geometry. PRESTUS falls back to the `matlab_gpu` path (`kspaceFirstOrderAS` with `gpuArray` data cast) with a warning. GPU acceleration is preserved; only the CUDA binary is bypassed. |

??? note "Details: axis mapping and PML in axisymmetric mode"
    `kspaceFirstOrder2D` uses a 2D Cartesian grid with X and Y axes, mapped onto the second and first dimensions of input matrices and coordinates ([y, x]). This mirrors 3D image voxel coordinates [i, j, k], which map onto [Y, X, Z]. `kspaceFirstOrderAS` uses axial and radial axes, mapped onto the first and second dimensions. Internally, `[radial x 2, axial]` is remapped onto `[axial, radial]` and only the columns from the beam axis (the transducer's lateral position) onward are retained. With the axis off centre, the expanded output repeats the outermost simulated value beyond the simulated radius. Medium, source and sensor are set up in [axial, radial] dimensions. After the simulation, values are either mirrored into a Cartesian 2D grid (`[radial x 2, axial]` again, as for `kspaceFirstOrder2D`) or radially expanded into a Cartesian 3D grid.

    **PML.** With `PMLInside=false` (the PRESTUS default) the PML is applied asymmetrically: the axial dimension receives `pml_size` voxels on both sides, the radial outer edge `pml_size` voxels on one side only. The inner radial boundary at `r = 0` is a symmetry axis and needs no PML. k-Wave handles this; no configuration change is needed.
