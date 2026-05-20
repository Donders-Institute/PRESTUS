## Transducer Calibration

To emulate transducers, PRESTUS optimises the velocity and phase settings of virtual transducers such that they produce focal axis profiles similar to those measured with real driving system–transducer setups in free water (measured either in-house via hydrophones or provided by device manufacturers). This calibrated emulation varies between different transducers, for different focal depth settings, and free-water pressures.

Transducer calibration relies on an additional config (`calibration_config`) that should be loaded as `parameters.calibration` (see the example `calibration_standalone`).

---

### Standalone calibration setup

**Script: `examples/calibration_standalone`**

#### Use cases

- Create a library of calibration settings for one (or multiple) transducer–depth settings
- Incorporate manufacturer information
- Template to set up and generate calibrated profiles for transducer equipment at the Donders Institute

#### Prerequisites

- Measured/estimated axial profiles (`path_input_axial`)
- Manufacturer-provided phase tables (`path_input_phase`)
- Entries in the equipment configuration (`PRESTUS/config/equipment/config_equipment.yaml`)

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
- Interpolate requested distance from multiple empirically measured distances.
- Add Focal Distance Offset (FDO; via `add_FDO`): Translate measured focal depths relative to the transducer exit plane into simulation-relevant coordinates centered on the transducer's mid-bowl. Missing distances are zero-interpolated.
- Select or interpolate the axial intensity profiles for the specified focal depth.
- Call `calibration_transducer`.

---

### Calibrate phase and amplitude settings

**Function: `calibration_transducer`**

#### Use cases

- Flexibility: only need to specify the desired axial profile
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

By default the initial phase guess is **random** (uniform over [0°, 360°]). Set `calibration.opt_use_initial_phases = true` to use the phases already in `transducer.annular.elem_phase_rad` (e.g. manufacturer-provided steering phases) as the starting point instead. This can improve convergence when good prior phases are available. Results can still vary across runs unless a seed is fixed via `opt_seed`.

Parameters that configure the search:

| Parameter | Default | Description |
|---|---|---|
| `opt_method` | `FEXminimize` | Optimization backend: `FEXminimize` (open-source, bundled) or `GlobalSearch` (MATLAB Global Optimization Toolbox) |
| `opt_weights` | `0` | Profile weighting: `0` = uniform across the full profile; `≥1` = Gaussian centred on the focal maximum with FWHM narrowing as weight increases (σ = focus_pos / weight) |
| `opt_limits` | full profile range | Distance range [mm] over which the error is evaluated; defaults to the non-NaN extent of the target profile |
| `opt_use_initial_phases` | `false` | Use phases from `transducer.annular.elem_phase_rad` as the starting point instead of random initialization. Useful when manufacturer-provided or previously calibrated phases are available. |
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

`plot_opt_sim_results` produces a comparison figure of the initial and optimized analytical and simulated profiles.

`save_optimized_values` writes two output files to `calibration.path_output_profiles`:

1. **CSV** (`calibration.filename_calibrated_CSV`) — a lookup table indexed by `[desired_intensity × focal_depth_ep]`. Each cell contains a string encoding `[opt_phases_deg, elem_amp]`. If the file already exists, the new entry is inserted and the table is re-sorted by intensity (rows) and focal depth (columns).

2. **YAML** (`<equipment_name>-F<focal_depth>mm-I<intensity>wpercm2.yaml`) — the full `transducer` parameter struct with optimized phases and amplitude, ready to be merged into a PRESTUS study config.

![calibration_optimized_analytical](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_optimized_analytical.png)

##### Step 11 — Cache cleanup

After all simulations complete, k-Wave `.mat` result files written to the cache folder during calibration are removed to free disk space. Set `parameters.io.save_acoustic_matrices = 1` or `parameters.io.save_matrices = 1` to retain them for debugging.

---

## In-pipeline ISPPA targeting

Rather than scaling `elem_amp` before the acoustic simulation, PRESTUS uses a two-step approach: a **water baseline measurement** records the free-water ISPPA at the current drive level, and **post-hoc pressure scaling** is applied when loading the cached acoustic results for thermal analysis. This means the acoustic simulation never needs to be repeated when changing the target intensity.

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

Runs a fast homogeneous water simulation at the current `elem_amp` (after source setup, before the main acoustic simulation). The peak free-water ISPPA is computed from the sensor output and stored as `acoustic_provenance.freefield_isppa_wcm2` alongside the acoustic cache file.

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
  elem_amp: ...
  # ... other transducer fields
```

`modules.run_water_baseline` defaults to `0` and must be explicitly enabled to record free-water ISPPA provenance. The acoustic simulation runs at whatever `elem_amp` is configured; `target_isppa_wcm2` controls only the analysis and thermal stages.

> **Backwards compatibility:** The older form `calibration.target_isppa_wcm2` is still accepted and is automatically migrated to `transducer(1).target_isppa_wcm2` at load time with a deprecation warning. Update configs to use the per-transducer form.

### Multiple targets (multi-ISPPA mode)

Setting `target_isppa_wcm2` to a vector automatically triggers the multi-ISPPA pipeline, which runs one thermal job per target in parallel without repeating the acoustic simulation:

```yaml
transducer:
  target_isppa_wcm2: [10, 20, 30, 50]   # W/cm² — triggers multi-ISPPA mode
```

See [doc_multi_isppa.md](doc_multi_isppa.md) for full details.

### Provenance fields

The following fields are saved in the acoustic cache alongside `sensor_data`:

| Field | Description |
|---|---|
| `acoustic_provenance.freefield_isppa_wcm2` | Free-water ISPPA measured at the drive level used for the acoustic simulation [W/cm²] |

### Relationship to `calibration_transducer`

The typical workflow is:

1. **Once per transducer/depth setup**: Run `calibration_transducer` to obtain optimized phases and an initial amplitude at a reference intensity. The output YAML encodes `elem_phase_deg` and `elem_amp`.
2. **Per subject simulation**: Load the calibration YAML into your study config. Set `transducer.target_isppa_wcm2` to the desired free-water intensity. The water baseline and post-hoc scaling handle the rest — no need to re-run the acoustic simulation when changing the target intensity.
