## Transducer Calibration

To emulate transducers, PRESTUS calibrates the phase and amplitude settings of virtual transducers so that their analytically predicted focal axis profiles match measurements made with real driving system–transducer setups in free water (measured in-house via hydrophones or provided by device manufacturers). Calibrated settings are stored in a transducer library and loaded at runtime.

Transducer calibration relies on an additional config (`calibration_config`) that should be loaded as `parameters.calibration` (see the example `calibration_standalone`).

---

## Calibration levels

The transducer library stores, for each measured depth, per-ring phase offsets, gains and the drive that reproduce the hydrophone profile (per-depth refinement, the default). Between measured depths they are interpolated; beyond them, geometric steering is used. See [doc_transducer_library.md](doc_transducer_library.md) for the model, the schema and the entry points.

| Level | Entry point | Phase source | Amplitude | k-Wave needed | Typical use |
|---|---|---|---|---|---|
| **Per-depth refinement** (default) | `build_transducer_library`, `refine_transducer_library`, `calibration_gui` | Analytical fit per measured depth (`fit_transducer_model`) | Analytical, per depth | No | Library entries from hydrophone profiles |
| **k-Wave check** (optional) | `kwave_check_refinement` (`KWave` option / GUI checkbox) | Unchanged | Corrected to free-water k-Wave ISPPA | Yes | Simulation-consistent drive at the refined depths |
| **Static joint model** (debugging) | `LIBRARY_MODE = 'static'`, `Mode = 'static'` | One depth-independent offset over all depths | Analytical fit per depth | No | Comparison with per-depth refinement |
| **Multi-depth joint fit** (legacy) | `calibration_setup` (`phase_mode: multi_depth`) | Joint optimisation across all depths | Analytical fit | No | Superseded by the static joint model |
| **Single-depth refinement** (legacy) | `calibrate_transducer` (`single_ref`) | Per-depth global search | Simulation-corrected | Optional | Older k-Wave calibration pipeline; results can be ingested with `update_transducer_library` |

Because the k-Wave check is optional, most entries stay at the "Analytical" amplitude level unless someone deliberately runs it. `load_transducer_from_library` does not assume that step happened: it resolves `tx_calibration.amp_kwave_verified` per requested depth and warns (`load_transducer_from_library:unverifiedAmplitude`) whenever the driven amplitude has not been checked — see [doc_transducer_library.md](doc_transducer_library.md#optional-free-water-k-wave-check).

### Focal depth convention (FLHM)

A "focal depth" quoted for a real transducer (e.g. a manufacturer steering table, or
`focal_depths_wrt_exit_plane` in a calibration config) is not the raw intensity-peak location —
it is the centre of the **Focal Length at Half Maximum** (FLHM), the midpoint of the contiguous
span around the peak where the axial profile stays at or above a threshold fraction of the peak.
Manufacturers differ on what "half maximum" means:

| `focal_reference` | Threshold | Convention |
|---|---|---|
| `flhm_3db` (default) | intensity ≥ peak / 2 | ITRUSST / Martin et al. (2025); Sonic Concepts CTX + TPO ("NeuroFUS") steering tables |
| `flhm_6db` | intensity ≥ peak / 4 (≡ pressure ≥ peak / 2) | Imasonic / IGT steering tables |
| `peak` | n/a — raw intensity-peak location | No FLHM centring |

Set `parameters.calibration.focal_reference` to match the convention the characterisation data
actually uses — mismatching it produces calibrated depths that look shifted by several mm/cm even
though the underlying phases/amplitudes are correct. The value is:

- Read by `compute_flhm_center` (via `focal_reference_to_frac`) wherever a calibration or
  post-hoc analysis needs to locate the FLHM centre of an axial profile.
- Echoed to the console during `calibration_setup` (`Focus from exit plane: 85.00 (FLHM -6dB)`).
- Recorded into the sidecar YAML's `meta.calibration.focal_reference` by `save_optimized_values`,
  so calibrated outputs stay traceable to the convention they were produced under.
- Used by `compute_focal_metrics` (and therefore `acoustic_analysis`) so post-hoc/in-pipeline
  focal-centre reporting stays consistent with how the transducer was calibrated. This is
  independent of the fixed ITRUSST −3 dB / −6 dB **volume** fields (`mask_3dB`/`mask_6dB`), which
  always use both thresholds regardless of `focal_reference`.

Defaults to `flhm_3db` when unset, preserving prior behaviour.

### Library build (per-depth refinement)

The library is built by `build_transducer_library.m` from the profiles in `local/profiles/`. By default (`LIBRARY_MODE = 'per_depth'`), every measured depth becomes a refinement fitted by `fit_transducer_model`: phase offsets and gains on top of geometric steering (`set_real_phases`, path-length geometry from the ring dimensions), plus the drive that reproduces the measured peak. The global model is geometric steering with the measured drive per depth. `meta.quality` records in-sample, held-out interpolation and geometric R². The result is written to `config/calibration/{combo_name}.yaml` (`LIBRARY_PATH`).

At runtime, load calibrated parameters with:
```matlab
parameters = load_transducer_from_library(combo_name, focal_depth_ep, isppa, equip_param);
```

### Single-depth refinement (legacy)

When the geometric baseline is insufficient at a specific depth — for example, due to systematic hardware delays — run `calibrate_transducer` in `single_ref` mode. This performs a global search at one reference depth to find optimised absolute phases and velocity, optionally followed by a k-Wave free-water simulation that corrects the analytical-to-simulation amplitude ratio.

The result is stored as a `calibration.focal_depths` entry with `is_refinement: true` in the library YAML. `load_transducer_from_library` reads it as a refinement knot (offset relative to geometric steering at that depth), interpolated and blended like any other refinement. For new entries, prefer per-depth refinement with the optional k-Wave check (see the levels table).

---

### Multi-depth joint fit (legacy)

Superseded by the static joint model of `fit_transducer_model` (`LIBRARY_MODE = 'static'`), which is itself kept only for debugging; per-depth refinement is the default.

When the hardware phase correction is expected to be **depth-independent** — i.e. the same per-element offset applies at every steering angle — the `multi_depth` mode finds a single correction vector by fitting all characterisation depths simultaneously.

Enable it by setting `calibration.phase_mode: multi_depth` in the calibration config. `calibration_setup` then calls `perform_joint_depth_fit` once instead of running `perform_global_search` independently per depth.

#### Algorithm

The optimisation finds a single vector `delta[1..n_elem]` (per-element phase offset in radians) and a scalar particle velocity `vel` that minimise the **mean of per-depth profile errors** across all `N` calibration depths:

$$\text{cost} = \frac{1}{N} \sum_{j=1}^{N} \mathrm{MSE}_j(\delta, v) + \lambda \cdot \overline{\delta^2}$$

where at each depth `j`:

$$\text{phases}_j = \text{geo\_phases}(j) + \delta$$

`geo_phases(j)` are the standard geometric steering phases recomputed fresh for each depth. The forward model (Rayleigh or O'Neil) is evaluated with `phases_j` and `vel` to produce a predicted axial intensity profile, which is compared against the empirical profile normalised to `desired_intensity` within the valid axial range.

The optional L2 term (weight `lambda`, default 0) regularises the magnitude of `delta` to prevent large hardware corrections without empirical support.

Optimisation starts from `delta = 0` (pure geometric steering) and is bounded to `±π` per element. The solver is selected by `opt_method` — `FEXminimize` (default, `popsize = 2000`) or MATLAB `GlobalSearch`.

#### Key difference from single-depth search

| | `multi_depth` | `single_depth` (per-depth global search) |
|---|---|---|
| Depths used | All simultaneously | One depth at a time |
| Shared parameters | One `delta` + one `vel` across all depths | Independent `delta` and `vel` per depth |
| Constraint | Joint consistency across all steering angles | Best fit at each individual depth |
| Typical use | Hardware offset is the same regardless of focus | Maximum per-depth accuracy |

The result is stored as `elem_phase_correction_deg` in the library YAML and as `calibration_type: multi_depth` in the `meta` block. At runtime, `load_transducer_from_library` adds these corrections on top of the interpolated geometric phases.

#### Configuration

```yaml
calibration:
  phase_mode: multi_depth       # enables joint fit
  opt_method: FEXminimize       # or GlobalSearch
  opt_upper_velocity: 0.2       # upper bound on velocity [m/s]
  opt_limits: [10, 80]          # axial fitting range [mm from bowl] (optional)
  opt_weights: 0                # 0 = uniform; ≥1 = Gaussian-weighted around focal peak
  opt_seed: 42                  # integer seed for reproducibility (optional)
  # lambda: 0                   # L2 regularisation weight (not yet a YAML key; pass via options struct)
```

A diagnostic figure `JointDepthFit.png` is saved to `outputs/img_calibration/` showing empirical (solid) vs model (dashed) profiles for every depth at the optimised parameters.

---

### Standalone calibration setup (single-depth refinement)

**Script: `examples/calibration_standalone`**

This script is the entry point for **single-depth refinements** — it runs `calibration_transducer` at specific depth × intensity combinations and writes per-depth override entries into the library. Use `build_transducer_library.m` to build the geometric baseline; use this script afterwards to refine specific depths.

#### Use cases

- Refine calibration at specific focal depths where the geometric baseline is insufficient
- Incorporate hardware-specific phase corrections
- Template to generate calibrated profiles for transducer equipment at the Donders Institute

#### Prerequisites

- Measured/estimated axial profiles (`path_input_axial`, usually `local/profiles/`; see [doc_transducer_library.md](doc_transducer_library.md#profile-folder-and-file-naming))
- Entries in the equipment configuration (one `<SERIAL>.yaml` per transducer in `PRESTUS/config/equipment/`)

Manufacturer phase tables are not needed: steering phases are always derived from ring geometry.

For the Donders, the empirical profiles (steering tables) can be found on the [FUSInitiative OneDrive](https://radbouduniversiteit.sharepoint.com/:f:/r/sites/FUSInitiative-SHAREDResearchinformation/Shared%20Documents/Software/PRESTUS/Acoustic%20profiling/)

Requested calibrations for unique TPO-transducer & focal depth & intensity combinations can be specified in `config_calibration.yaml` via `combinations`, `focal_depths_wrt_exit_plane`, and `desired_intensities`.

Each line in the configuration corresponds to a unique TPO-transducer setup. In the following example, an IGT transducer with ID `PCD15287_01001` would be emulated for two focal depths of 40 and 50 mm from the exit plane, each for free-water intensities of 30 and 60 W/cm² respectively. For a second transducer (`PCD15473_01001`), emulated phases and amplitudes would be provided for a depth of 40 mm at an intensity of 30 W/cm².

> This example fits a 32-channel transducer with 10 artificial channels. The number of emulated channels can impact the stability of the fitting solution. It is governed by the setup in `config_equipment.yaml`.

Transducer-TPO setups to be characterized:
```yaml
combinations:
  - IS_PCD15287_01001_IGT_32_ch_comb_10_ch
  - IS_PCD15473_01001_IGT_32_ch_comb_10_ch
```
List of focal depths (in mm) to be characterized:
```yaml
focal_depths_wrt_exit_plane:
  - [40, 50]
  - [40]
```
List of intensities (in free-water W/cm²) to be characterized:
```yaml
desired_intensities:
  - [30, 60]
  - [30]
```

#### Steps

- Define and initialize the simulation environment by setting paths and loading configuration files with equipment and user calibration data.
- Automatically load Donders-specific transducer information: Identify specific equipment combinations and extract associated transducer and driving system parameters, setting default initial amplitudes and phases.
- Set initial simulation to manufacturer data: Load measured characterization data of the actual transducer's acoustic field including axial intensity profiles and phase information provided by manufacturers.
- Interpolate requested distance from multiple empirically measured distances (`interpolate_axial_profiles`): the two neighbouring profiles are aligned on their focus (`calibration.profile_align`, default the column focal-depth labels) and blended linearly; outside the measured range only up to `calibration.max_extrapolation_mm` (default 0: error). The derived profile is saved as `Profile_<type>_at_F_<depth>_<equipment>.csv` (+ sidecar) and plotted to `path_output_profiles`.
- Add Focal Distance Offset (FDO; via `add_FDO`): Translate measured focal depths relative to the transducer exit plane into simulation-relevant coordinates centered on the transducer's mid-bowl. Missing distances are zero-interpolated.
- Select or interpolate the axial intensity profiles for the specified focal depth.
- Call `calibration_transducer`.

---

### Calibrate phase and amplitude settings (refinement pipeline)

**Function: `calibration_transducer`**

This function implements the single-depth refinement pipeline. It is called by `calibration_standalone` and produces single-depth calibrations that `update_transducer_library` ingests as refinements. For building the default library across all depths, use `build_transducer_library.m` instead.

#### Use cases

- Per-depth refinement with measured hydrophone profiles
- Dynamic integration into end-to-end simulation loops (e.g., iterating across an amplitude–distance parameter space)

#### Prerequisites

- `profile_empirical.axial_intensity`  
  Desired intensity profile along the focal beam axis [W/cm²]
- `profile_empirical.axial_distance_bowl`  
  Distance from the transducer bowl [mm]
- `desired_focal_distance_ep`  
  Requested focal distance from the transducer exit plane [mm]

#### Detailed pipeline

##### Step 1 — Scale the empirical profile to the desired intensity

`scale_real_intensity_profile` linearly rescales the input empirical profile so its peak equals `desired_intensity`. It simultaneously sets the initial `elem_amp` in `parameters.transducer.annular` to the pressure amplitude corresponding to `desired_intensity` via:

$$p = \sqrt{2 \cdot I_\mathrm{desired} \cdot 10^4 \cdot \rho \cdot c}$$

The profile is then truncated or NaN-padded to match the simulation axis length (derived from `grid.default_dims(end)`).

##### Step 2 — Compute analytical bootstrap profile

`extract_analytical_profile` evaluates the configured forward model using the nominal transducer parameters from the config (no k-Wave simulation is needed at this stage). This provides an initial velocity estimate and an intensity shape that phase optimisation can work with.

The forward model is selected by `calibration.forward_model`:

| Value | Description |
|---|---|
| `'oneil'` (default) | O'Neil closed-form solution via `focusedAnnulusONeil` |
| `'rayleigh'` | Rayleigh–Sommerfeld integral via `rayleigh_axial_intensity` |

> **Coordinate note:** `focusedAnnulusONeil` receives axial positions as `(axial_position - 0.5) × 10^{-3}` metres. The `−0.5` half-voxel shift converts from voxel-centre to voxel-edge coordinates.

##### Step 3 — Determine element phases

Two modes are available, selected by whether `calibration.elem_phase_correction_deg` is set:

**Global search (default):** `perform_global_search` minimises the weighted mean-squared error between the analytical forward model and the scaled empirical target profile over element phases [0, 2π] and particle velocity [0.001, `opt_upper_velocity`]. This runs purely analytically using the bootstrap profile from Step 2 — no k-Wave simulation is needed.

By default the global search is seeded from the **geometric steering phases** (computed analytically from the transducer geometry and target depth). Set `calibration.opt_use_initial_phases = true` to instead use the phases already in `transducer.annular.elem_phase_rad` (e.g. manufacturer-provided phases) as the starting point. Results can still vary across runs unless a seed is fixed via `opt_seed`.

Parameters that configure the search:

| Parameter | Default | Description |
|---|---|---|
| `opt_method` | `FEXminimize` | Optimization backend: `FEXminimize` (open-source, bundled) or `GlobalSearch` (MATLAB Global Optimization Toolbox) |
| `opt_weights` | `0` | Profile weighting: `0` = uniform across the full profile; `≥1` = Gaussian centred on the focal maximum with FWHM narrowing as weight increases (σ = focus_pos / weight) |
| `opt_limits` | full profile range | Distance range [mm] over which the error is evaluated; defaults to the non-NaN extent of the target profile |
| `opt_use_initial_phases` | `false` | Use phases from `transducer.annular.elem_phase_rad` as the starting point. Default (`false`) seeds from geometric steering phases. Set to `true` to use manufacturer-provided or previously calibrated phases instead. |
| `opt_seed` | *(none)* | Integer random seed for reproducibility; if unset, results may vary across runs |
| `opt_upper_velocity` | `0.2` m/s | Upper bound on particle velocity during search |
| `skip_front_peak_mm` | `0` | Excludes the first N mm from peak detection to avoid near-field artifacts. Applied in both the global search objective and the amplitude correction (Step 5). Does not restrict the fitting range — use `opt_limits` for that. |

The `FEXminimize` backend uses fixed internal settings: `popsize = 5000`, `FinDiffType = 'central'`, `TolCon = 1e-8`.

**Geometric correction mode:** If `calibration.elem_phase_correction_deg` is non-empty, the global search is skipped. The provided per-element offsets are added to the geometric steering phases from the config to form the combined element phases. This encodes systematic hardware delays (e.g. cable length mismatches) measured at a reference depth. Set `calibration.save_elem_correction = true` to persist global-search results as hardware corrections for reuse at other focal depths.

![calibration_fitting](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_fitting.png)
Example profile fit with uniform `opt_weights`.

##### Step 4 — Run free-water correction simulation (optional)

A k-Wave simulation is run in free water using the nominal source parameters. The simulated peak intensity is compared to the analytical prediction to compute `simulated_analytical_scaling` — the ratio by which the analytical model under- or over-estimates the simulation. This factor corrects the source amplitude in Step 6.

How simulations are run is determined by the main `parameters` (e.g., `simulation.code_type`). Additional settings in `parameters.calibration` can override default behaviour:

- `run_free_water_sim`  
  Run a free-water k-Wave simulation to compute the analytical-to-simulation scaling factor (`true`, default). Set to `false` for a fully simulation-free calibration — appropriate when the analytical model alone is considered sufficient (e.g. single-element transducers with no inter-element phase steering). The scaling factor defaults to 1 when skipped.
- `axisymmetric2D`  
  Override default 3D simulation to perform axisymmetric 2D water simulations (`1` = yes, `0` = no, default `0`).
- `force_kwavearray`  
  Force free-water simulations to use kWaveArray (`1`). If `0`, uses the setting in the default or study-specific config.
- `save_in_calibration_folder`  
  If `true` (default), all simulation outputs are redirected to `calibration.path_output` rather than the subject output folder.

##### Step 5 — Compute analytical reference profile with optimized phases

`compute_analytical_solution` evaluates the configured forward model with the final optimized phases to produce a reference analytical intensity profile and computes `simulated_analytical_scaling`:

$$\mathrm{simulated\_analytical\_scaling} = \frac{\max(I_\mathrm{simulated})}{\max(I_\mathrm{analytical})}$$

This ratio captures the systematic offset between the k-Wave simulation and the analytical model and is propagated through the amplitude calculation in Step 6. When `run_free_water_sim = false`, this scaling factor is 1 by construction.

#### Desired intensity loop

##### Step 5b — Amplitude correction to match target ISPPA

After phase determination optimises profile **shape**, the peak intensity of the analytical profile may not exactly equal `desired_intensity` because the search jointly optimises phases and velocity without a hard intensity constraint.

If only amplitude needs to be calibrated, this also does not require phase recalibration.

`fit_velocity_to_intensity` scales the velocity to match the target intensity. Since intensity scales with velocity squared ($I \propto v^2$), the scaled target velocity is:

$$v_\mathrm{target} = v_\mathrm{opt} \cdot \sqrt{\frac{I_\mathrm{target}}{I_\mathrm{peak}}}$$

`skip_front_peak_mm` is applied when finding $I_\mathrm{peak}$.

A warning is issued if `corrected_velocity` exceeds `opt_upper_velocity`.

##### Step 6 — Calculate the optimized source amplitude

The amplitude that yields the corrected velocity for use in the k-Wave simulation is:

$$
\mathrm{amp}_\text{optimized} = \mathrm{round}\!\left( \frac{v_\text{corrected}}{v_\text{original}} \cdot \sqrt{\frac{I_\text{analytical}}{I_\text{sim}}} \cdot \mathrm{amp}_\text{original} \right)
$$

where `v_original` and `amp_original` come from the initial simulation (Steps 4–5), `v_corrected` is the velocity that yields an analytic match. To correct for minor mismatches between analytical and simulated intensities, the velocity scaling is corrected via the square root of analytical/simulated intensity ratio (i.e., the sqrt of 1/`simulated_analytical_scaling` from Step 5).

##### Step 7 — Recalculate analytical solution with optimized parameters

`recompute_analytical_solution` evaluates the configured forward model with `opt_phases` and the corrected `opt_velocity` (from Step 5b) to produce an optimized analytical profile. This profile is plotted against the target and the reference analytical solution for visual inspection.

##### Step 8 — Rerun water simulation with optimized settings (Optional)

The pipeline reruns `prestus_pipeline_start` with `elem_amp = elem_amp_optimized`, `elem_phase_rad = opt_phases`, and `output_affix = '_optimized'`.

Whether a free-water simulation is run to verify the calibration results is governed by parameter `opt_amp_validation`.
By default, this is set to always, but it can also be set to 'never', 'initial' (only for the first target intensity), or 'final' (only for the final target intensity).

##### Step 9 — Extract simulated optimized intensity along the focal axis

Same procedure as Step 4, applied to the optimized simulation results.

| Initial simulation | Optimized simulation |
|--------------------|----------------------|
| ![calibration_initial_intensity](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_initial_intensity.png) | ![calibration_opt_intensity](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_opt_intensity.png) |

The black line indicates the **transducer bowl**, the red line the **transducer exit plane**, the white line the **maximum focal intensity** (see [distance definitions](doc_transducer.md#target-distance-parameters)).

##### Step 10 — Plot and save results

`plot_opt_sim_results` produces a comparison figure with the following legend entries: **Phase-optimised (Analytical)**, **Amplitude-calibrated (Analytical)**, and **Correction sim (Simulated)**. These labels reflect the role of each profile in the pipeline: the phase-optimised curve shows the analytical result after the global search; the amplitude-calibrated curve shows it after velocity scaling; the correction sim curve is the k-Wave free-water result used to compute `simulated_analytical_scaling`.

`save_optimized_values` writes two output files to `calibration.path_output_profiles`:

1. **CSV** (`calibration.filename_calibrated_CSV`) — a lookup table indexed by `[desired_intensity × focal_depth_ep]`. Each cell contains a string encoding `[opt_phases_deg, elem_amp]`. If the file already exists, the new entry is inserted and the table is re-sorted by intensity (rows) and focal depth (columns).

2. **YAML** (`<equipment_name>-F<focal_depth>mm-I<intensity>wpercm2.yaml`) — the full `transducer` parameter struct with optimized phases and amplitude, ready to be merged into a PRESTUS study config.

![calibration_optimized_analytical](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_optimized_analytical.png)

##### Step 11 — Cache cleanup

After all simulations complete, k-Wave `.mat` result files written to the cache folder during calibration are removed to free disk space. Set `parameters.io.save_acoustic_matrices = 1` (or `parameters.io.save_matrices = 1` while `save_acoustic_matrices` is left empty) to retain them for debugging.

---

## Library YAML format

The library YAML schema (`meta`, `global_model`, `refinement`, and the legacy forms that are still read) is documented in [doc_transducer_library.md](doc_transducer_library.md#library-yaml-schema).

---

## In-pipeline ISPPA targeting

Rather than scaling `elem_p_amp` before the acoustic simulation, PRESTUS uses a two-step approach: a **water baseline measurement** records the free-water ISPPA at the current drive level, and **post-hoc pressure scaling** is applied when loading the cached acoustic results for thermal analysis. This means the acoustic simulation never needs to be repeated when changing the target intensity.

| | `calibration_transducer` | ISPPA targeting |
|---|---|---|
| **When** | Once per transducer/depth/intensity, before any study | Per subject, at analysis time |
| **What is optimized** | Phases **and** amplitude | Neither — scaling is post-hoc |
| **Input** | Empirical hydrophone profile + target intensity [W/cm²] | Target free-water ISPPA [W/cm²] |
| **Output** | Calibrated YAML config ready for study use | Scaled pressure field fed into thermal sim |
| **Typical use** | Build a calibration library | Target a specific sonication intensity |

---

### Step 1 — Water baseline measurement

**Function: `water_baseline`**

Runs a fast homogeneous water simulation at the current `elem_p_amp` (after source setup, before the main acoustic simulation). The peak free-water ISPPA is computed from the sensor output and stored as `acoustic_provenance.freefield_isppa_wcm2` alongside the acoustic cache file. The free-water intensity along each transducer's beam axis (bowl apex → grid edge, distances from the apex) is stored as `acoustic_provenance.freefield_axial`. Combined async / FOV caches drop it, because the combined field is scaled per source.

This step is **disabled by default** (`modules.run_water_baseline = 0`). Enable it when ISPPA scaling provenance is needed. It adds a few minutes on a typical workstation (no tissue mapping, uses the configured `code_type`):

```yaml
modules:
  run_water_baseline: 1
```

> If `run_water_baseline` is disabled, `transducer.target_isppa_wcm2` will have no effect and a warning is issued.

---

### Step 2 — Post-hoc pressure scaling

When `transducer.target_isppa_wcm2` is set and the acoustic cache contains valid provenance, the pipeline scales `sensor_data.p_max_all` before acoustic analysis and thermal simulation:

$$p_\mathrm{scaled} = p_\mathrm{cached} \times \sqrt{\frac{I_\mathrm{target}}{I_\mathrm{baseline}}}$$

Intensity-derived fields (used in thermal simulation) then scale as $I \propto p^2$, so heat deposition is correctly proportional to the target intensity without rerunning the skull simulation.

A warning is issued if the scale factor exceeds 4× or is below 0.25×, indicating that the target is far from the simulated amplitude and the linear assumption may not hold.

### Configuration

`target_isppa_wcm2` is set **per transducer** in the transducer block:

```yaml
transducer:
  target_isppa_wcm2: 30   # Target free-water ISPPA [W/cm²]
  elem_p_amp: ...
  # ... other transducer fields
```

`modules.run_water_baseline` defaults to `0` and must be explicitly enabled to record free-water ISPPA provenance. The acoustic simulation runs at whatever `elem_p_amp` is configured; `target_isppa_wcm2` controls only the analysis and thermal stages.

> **Backwards compatibility:** The older form `calibration.target_isppa_wcm2` is still accepted and is automatically migrated to `transducer(1).target_isppa_wcm2` at load time with a deprecation warning. Update configs to use the per-transducer form.

### Multiple targets (multi-ISPPA mode)

Setting `target_isppa_wcm2` to a vector automatically triggers the multi-ISPPA pipeline, which runs one thermal job per target in parallel without repeating the acoustic simulation:

```yaml
transducer:
  target_isppa_wcm2: [10, 20, 30, 50]   # W/cm² — triggers multi-ISPPA mode
```

See [doc_multiintensity.md](doc_multiintensity.md) for full details.

### Provenance fields

The following fields are saved in the acoustic cache alongside `sensor_data`:

| Field | Description |
|---|---|
| `acoustic_provenance.freefield_isppa_wcm2` | Free-water ISPPA measured at the drive level used for the acoustic simulation [W/cm²] |

### Relationship to calibration

The typical workflow is:

1. **Once per transducer (library build)**: Run `build_transducer_library.m` to produce geometric-baseline phases and velocity fits for all measured depths. Results are stored in `config/transducer/{combo_name}.yaml`.
2. **Optionally, per depth (refinement)**: Run `calibration_transducer` at specific depths where greater accuracy is required. Refined entries override the library baseline at ±1 mm.
3. **Per subject simulation**: Load calibrated parameters via `load_transducer_from_library`. Set `transducer.target_isppa_wcm2` to the desired free-water intensity. The water baseline and post-hoc scaling handle the rest — no need to re-run the acoustic simulation when changing the target intensity.
