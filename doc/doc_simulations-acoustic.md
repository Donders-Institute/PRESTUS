# Acoustic simulations

PRESTUS supports 3D acoustic wave simulations using k-Wave’s pseudo-spectral method for time-domain 3D acoustic simulations. These simulations model compressional waves (P-waves), which involve pressure changes propagating through isotropic media *without* shear deformations. To describe wave propagation, they calculate pressure and velocity fields. They are faster, less memory-intensive, and require less assumptions than elastic wave simulations (e.g., using `kWaveElastic`; not implemented) because they ignore shear waves and complex material interactions. However, this also means that scattering via shear waves is not explicitly modeled (albeit contributing to the bulk attenuation coefficient resulting from absorption and scattering).

PRESTUS by default supports 3D simulations using `kspaceFirstOrder3D`.

See the available [simulation backends](doc_backend.md).

> **Turning off acoustic simulations.** PRESTUS uses the `run_acoustic_sims` flag to indicate whether to run simulations. If set to `false`, PRESTUS will first check for existing files and check whether `overwrite_files` is set to off to load the existing matrices. This is intended to make this more intuitive for workflows in which acoustic simulations are "turned off" to run follow-up thermal simulations. This behaviour is unique for this module (turning off other modules may not run the step at all incl. loading existing data).

## Spatial and temporal sampling

Acoustic simulation accuracy depends on resolving the shortest wavelength in the domain. Two parameters govern this:

- **Grid resolution (`grid.resolution_mm`)** sets the spatial step `dx`. Together with the medium sound speed and fundamental frequency, this determines the number of spatial samples per wavelength:

  ```
  PPW = c / (f × dx)
  ```

  The worst case (lowest PPW, highest risk of under-sampling) occurs in the slowest medium (typically water or brain at ~1500 m/s). A minimum of **PPW ≥ 6** is recommended; PRESTUS warns if this is violated (configurable via `grid.min_ppw`).

- **CFL number (`grid.source_cfl`)** sets the temporal resolution via the time step `dt`, relative to the spatial step `dx` and the maximum sound speed in the domain:

  ```
  dt = CFL × dx / c_max
  ```

  `c_max` is the fastest medium in the domain (typically skull at ~2800 m/s), which sets the stability limit. The number of temporal samples per wave period then follows from the spatial sampling (PPW) and the CFL:

  ```
  PPP = PPW / CFL
  ```

  Tools that only expose PPW implicitly couple spatial and temporal resolution. PRESTUS exposes them independently: `resolution_mm` controls spatial resolution, `source_cfl` controls temporal resolution given that spatial step. This means you can choose a coarser or finer time step independently of the grid — but both must be set consistently to remain within the stability limit.

> **Rule of thumb.** To find the coarsest grid that satisfies a minimum PPW at a given frequency and minimum sound speed:
> ```
> dx_max = c_min / (f × PPW_min)
> ```
> For example, at 500 kHz in water (1500 m/s) with PPW = 6: `resolution_mm = 0.50 mm`.

### Additional considerations

**PPW varies across tissues.** The minimum-PPW check (`grid.min_ppw`) uses the slowest medium sound speed, which is the worst case for spatial sampling. In practice, PPW varies considerably across the domain: at 500 kHz with 0.5 mm resolution, water/brain (1500 m/s) gives PPW ≈ 10 while skull (2800 m/s) gives PPW ≈ 37. The skull is always massively over-resolved spatially relative to the soft tissues — the binding spatial constraint is always in the slow tissues.

**The stability recheck can reduce `dt` below the CFL target.** After source setup, `prestus_pipeline` calls `checkStability` and, if the CFL-derived `dt` exceeds the estimated stability limit for the heterogeneous medium, reruns source setup with a smaller time step scaled by `grid.source_limit_fraction` (default 0.9). The final effective `dt` may therefore be smaller than `CFL × dx / c_max`. The PPW validation fires only on the first source setup call; the adjusted call bypasses it since `dt` is passed explicitly.

**Overriding `source_ppw` decouples temporal from spatial resolution.** If `grid.source_ppw` is set manually, it only affects the time step derivation (`PPP = source_ppw / CFL`). The spatial accuracy warning (`grid.min_ppw`) still computes the actual PPW from `c_min / (f × dx)` independently. It is therefore possible to set `source_ppw` to a coarser value to increase `dt` and reduce runtime without suppressing the spatial resolution warning.

**Default CFL of 0.15 is intentionally conservative.** k-Wave's own default CFL is 0.3; PRESTUS uses 0.15 to provide additional temporal stability margin in heterogeneous skull simulations, where the theoretical stability limit is only an approximation (see [`checkStability`](http://www.k-wave.org/documentation/checkStability.php)). This roughly doubles the number of time steps — and acoustic simulation runtime — compared to k-Wave's default. If runtime is a concern and simulations are stable, `source_cfl` can be increased toward 0.3, but this should be validated against `checkStability` output for the specific medium.

---

## Axisymmetric simulations

For debugging, it can be useful to specify and model symmetric 2D tissues (e.g, 2D `phantom` or `water` media). However, [`kspaceFirstOrder2D`](http://www.k-wave.org/documentation/kspaceFirstOrder2D.php) models an infinite cylinder, not a focusing shell. It uses a 2D Cartesian grid with **X and Y axes**. These are assumed to map onto the second and first dimensions of input matrices and coordinates (i.e., [y, x]), respectively. This may at first may seem unintuitive, but it is designed to facilitate the mapping to 3D image space voxel coordinates [i,j,k], which traditionally map onto [Y,X,Z]. In the future, this could also be resolved by assuming a XY mapping for k-Wave, and introducing a remapping of voxel space coordinates.

To model a 3D axisymmetric problem using a 2D computational grid (as recommended), [`kspaceFirstOrderAS`](http://www.k-wave.org/documentation/kspaceFirstOrderAS.php) can be used (see [this example benchmark](https://github.com/ucl-bug/k-wave/blob/main/k-Wave/examples/example_at_focused_bowl_AS.m). It significantly reduces computational cost compared to a full 3D simulation. The input structures (`kgrid`, `medium`, `source`, and `sensor`) are defined similarly to 2D simulations but interpreted in the axisymmetric coordinate system. The function `kspaceFirstOrderAS` accounts for wave propagation in such axisymmetric media and is functionally similar to `kspaceFirstOrder2D` but adapted for axisymmetry. kspaceFirstOrderAS uses a 2D axisymmetric coordinate system with **axial and radial axes**. These are assumed to map onto the first and second dimensions of input matrices and coordinates, respectively.

Multiple adjustments are made if axisymmetry is requested:

- The grid should be specified as `[radial x 2, axial]`.
- It is assumed that the (2x) radial axis is shorter than the axial axis.
- Internally, `[radial x 2, axial]` is remapped onto `[axial radial]` and only values from `[half+1:end]` are retained.
- Medium, source, and sensor are set up according to [axial, radial] dimensions.
- Following acoustic simulation, values are either mirrored into a cartesian 2D grid, or radially expanded into a cartesion 3D grid (see the note below). In cartesian 2D, the output will have [radial x 2, axial] dimensions again, in line with the default specification for kspaceFirstOrder2D.
- Only a **single transducer** is supported; specifying multiple transducers with `axisymmetric == 1` raises an error.
- Only **annular transducers** explicitly support axisymmetric `kWaveArray` mode. Other transducer types are not tested in this configuration.

> **Axisymmetric 2D vs. 3D output**
> Axisymmetric simulations are effectively 3D simulations that are performed efficiently on 2D input grids. By default, PRESTUS assumes that users want to follow-up an axisymmetric acoustic simulation with a 3D heating simulation. When heating is requested in the simulation pipeline, axisymmetric acoustic simulation outputs will be provided in 3D. If 2D outputs are desired, heating has to be deactivated. 2D heating simulations are possible by first running an acoustic 2D simulation without heating, and then running a heating simulation with the prior output, while deactivating (or not overwriting) the acoustic simulation.

### Backend compatibility

| Backend (`simulation.code_type`) | Axisymmetry support | Notes |
|---|---|---|
| `matlab_cpu` | ✓ | Uses `kspaceFirstOrderAS` |
| `matlab_gpu` | ✓ | Uses `kspaceFirstOrderAS` with `RadialSymmetry` option |
| `cpp_cpu` | ✓ | Uses the C++ OMP binary `kspaceFirstOrderASC` |
| `cpp_gpu` (CUDA) | ✗ | No CUDA binary supports axisymmetric geometry. PRESTUS automatically falls back to the `matlab_gpu` path (`kspaceFirstOrderAS` with `gpuArray` data cast) with a warning. GPU acceleration is preserved; only the CUDA binary is bypassed. |

### PML behaviour in axisymmetric mode

The PML is applied asymmetrically when `PMLInside=false` (the PRESTUS default): the **axial** dimension receives `pml_size` voxels on both sides, while the **radial** outer edge receives `pml_size` voxels on one side only. The inner radial boundary at `r = 0` is a symmetry axis and does not require a PML. k-Wave handles this automatically; no configuration change is needed.


