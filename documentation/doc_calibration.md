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
- Entries in the equipment configuration (`PRESTUS/configs/equipment/equipment_config.yaml`)

For the Donders, the empirical profiles (steering tables) can be found on the [FUSInitiative OneDrive](https://radbouduniversiteit.sharepoint.com/:f:/r/sites/FUSInitiative-SHAREDResearchinformation/Shared%20Documents/Software/PRESTUS/Acoustic%20profiling/)

Requested calibrations for unique TPO-transducer & focal depth & intensity combinations can be specified in `calibration_config.yaml` via `combinations`, `focal_depths_wrt_exit_plane`, and `desired_intensities`.

Each line in the configuration corresponds to a unique TPO-transducer setup. In the following example, an IGT transducer with ID `PCD15287_01001` would be emulated for two focal depths of 40 and 50 mm from the exit plane, each for free-water intensities of 30 and 60 W/cm² respectively. For a second transducer (`PCD15473_01001`), emulated phases and amplitudes would be provided for a depth of 40 mm at an intensity of 30 W/cm².

> This example fits a 32-channel transducer with 10 artificial channels. The number of emulated channels can impact the stability of the fitting solution. It is governed by the setup in `equipment_config.yaml`.

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

##### Step 2 — Run initial free-water simulation

How simulations are run is determined by the main `parameters` (e.g., `simulation.code_type`). Additional settings in `parameters.calibration` can override default behaviour:

- `axisymmetric2D`  
  Override default 3D simulation to perform axisymmetric 2D water simulations (`1` = yes, `0` = no, default `0`).
- `force_kwavearray`  
  Force free-water simulations to use kWaveArray (`1`). If `0`, uses the setting in the default or study-specific config.
- `save_in_calibration_folder`  
  If `true` (default), all simulation outputs are redirected to `calibration.path_output` rather than the subject output folder.

##### Step 3 — Extract simulated intensity along the focal axis

`extract_simulated_profile` retrieves `p_max_all` from the simulation results and:
- In 3D, takes the lateral slice at the transducer's lateral centre (`trans_pos(1:2)`)
- Extracts the 1D axial pressure profile from the transducer position onward
- Converts pressure to intensity: $I = p^2 / (2 \rho c) \times 10^{-4}$ [W/cm²]
- Back-computes particle velocity from `elem_amp` via the acoustic impedance relation: $v = \mathrm{elem\_amp} / (\rho \cdot c)$

The axial distance axis is expressed in mm from the transducer bowl.

##### Step 4 — Compute the analytical O'Neil solution

`compute_oneil_solution` calls `focusedAnnulusONeil` with the initial simulation velocity and phases to produce an analytical pressure profile, converted to intensity. This establishes a baseline for the optimization.

It also computes `simulated_analytical_scaling`:

$$\mathrm{simulated\_analytical\_scaling} = \frac{\max(I_\mathrm{simulated})}{\max(I_\mathrm{O'Neil})}$$

This ratio captures the systematic offset between the k-Wave simulation and the analytical model and is propagated through the amplitude calculation in Step 7.

> **Coordinate note:** `focusedAnnulusONeil` receives axial positions as `(axial_position - 0.5) × 10^{-3}` metres. The `−0.5` half-voxel shift converts from voxel-centre to voxel-edge coordinates.

##### Step 5 — Optimize element phases and velocity (global search)

`perform_global_search` minimises the weighted mean-squared error between the analytical O'Neil profile and the scaled empirical target profile over element phases [0, 2π] and particle velocity [0.001, `opt_upper_velocity`].

The initial phase guess is **random** (uniform over [0°, 360°]), not the manufacturer-provided phases. This means results can vary across runs unless a seed is fixed via `opt_seed`.

Parameters that configure the search:

| Parameter | Default | Description |
|---|---|---|
| `opt_method` | `FEXminimize` | Optimization backend: `FEXminimize` (open-source, bundled) or `GlobalSearch` (MATLAB Global Optimization Toolbox) |
| `opt_weights` | `0` | Profile weighting: `0` = uniform across the full profile; `≥1` = Gaussian centred on the focal maximum with FWHM narrowing as weight increases (σ = focus_pos / weight) |
| `opt_limits` | full profile range | Distance range [mm] over which the error is evaluated; defaults to the non-NaN extent of the target profile |
| `opt_seed` | *(none)* | Integer random seed for reproducibility; if unset, results may vary across runs |
| `opt_upper_velocity` | `0.2` m/s | Upper bound on particle velocity during search |
| `skip_front_peak_mm` | `0` | Excludes the first N mm from peak detection to avoid near-field artifacts. Applied in both the global search objective and the amplitude correction (Step 5b). Does not restrict the fitting range — use `opt_limits` for that. |

The `FEXminimize` backend uses fixed internal settings: `popsize = 5000`, `FinDiffType = 'central'`, `TolCon = 1e-8`.

![calibration_fitting](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/calibration_fitting.png)
Example profile fit with uniform `opt_weights`.

##### Step 5b — Amplitude correction to match desired peak intensity exactly (optional)

After the global search optimises profile **shape**, the peak intensity of the analytical profile may not exactly equal `desired_intensity` because the search jointly optimises phases and velocity without a hard intensity constraint.

`fit_velocity_to_intensity` corrects this analytically. Since intensity scales with velocity squared ($I \propto v^2$), the corrected velocity is:

$$v_\mathrm{corrected} = v_\mathrm{opt} \cdot \sqrt{\frac{I_\mathrm{desired} \cdot \mathrm{simulated\_analytical\_scaling}}{I_\mathrm{peak}}}$$

The `simulated_analytical_scaling` factor is included because `opt_source_amp` in Step 7 divides by it, so the analytical target must overshoot by that factor to yield the correct intensity in the final simulation.

`skip_front_peak_mm` is applied when finding $I_\mathrm{peak}$.

A warning is issued if `corrected_velocity` exceeds `opt_upper_velocity`.

This step is **enabled by default**. To disable:
```yaml
calibration:
  fit_velocity_to_intensity: false
```

When disabled, the velocity from the global search is used as-is, and the peak intensity of the optimized simulation may deviate from `desired_intensity`.

##### Step 6 — Recalculate the analytical O'Neil solution with optimized parameters

`recompute_oneil_solution` calls `focusedAnnulusONeil` with `opt_phases` and the corrected `opt_velocity` (from Step 5b) to produce the final analytical profile. This profile is plotted against the target and the original O'Neil solution for visual inspection.

##### Step 7 — Calculate the optimized source amplitude

The amplitude that yields the corrected velocity in the k-Wave simulation is:

$$\mathrm{elem\_amp\_optimized} = \mathrm{round}\!\left( \frac{v_\mathrm{corrected}}{v_\mathrm{original}} \cdot \frac{\mathrm{elem\_amp\_original}}{\mathrm{simulated\_analytical\_scaling}} \right)$$

where `v_original` and `elem_amp_original` come from the initial simulation (Steps 2–3), and `simulated_analytical_scaling` is from Step 4.

##### Step 8 — Rerun water simulation with optimized phases and amplitude

The pipeline reruns `prestus_pipeline_start` with `elem_amp = elem_amp_optimized`, `elem_phase_rad = opt_phases`, and `output_affix = '_optimized'`.

##### Step 9 — Extract simulated optimized intensity along the focal axis

Same procedure as Step 3, applied to the optimized simulation results.

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

---

## In-pipeline amplitude calibration

**Function: `calibration_amplitude_scaling`**

### Purpose and scope

This is a lightweight amplitude normalization step that runs *inside* the main simulation pipeline. It answers a different question than `calibration_transducer`:

| | `calibration_transducer` | `calibration_amplitude_scaling` |
|---|---|---|
| **When** | Once per transducer/depth/intensity, before any study | Every pipeline run (optional) |
| **What is optimized** | Phases **and** amplitude | Amplitude only |
| **Input** | Empirical hydrophone profile + target intensity [W/cm²] | Target peak free-water pressure [MPa] |
| **Output** | Calibrated YAML config ready for study use | `elem_amp` scaled in-place; pipeline continues |
| **Typical use** | Build a calibration library | Normalize a subject simulation to a specific free-water pressure |

Use `calibration_transducer` first to obtain calibrated phases and an initial amplitude. Then, if you want to run the full pipeline at a *specific free-water peak pressure* (e.g., 0.5 MPa), enable `calibration_amplitude_scaling` to fine-tune the amplitude automatically.

### How it works

1. A homogeneous water simulation is run using the current `kgrid`, `source`, and `sensor` (after source setup but before the main acoustic simulation). The backend is forced to `matlab_cpu` regardless of the main simulation's `code_type`.
2. The peak pressure across all sensor points is measured: `peak_pa = max(sensor_data.p_max_all)`.
3. A linear scale factor is computed: `scale = target_pa / peak_pa`.
4. `elem_amp` in all configured transducers is multiplied by `scale`.
5. The main acoustic simulation proceeds with the scaled amplitude.

The scaling is linear because k-Wave solves the linear wave equation: doubling `elem_amp` doubles peak pressure. No analytical model or phase fitting is involved — the peak is read directly from the simulation, so no `simulated_analytical_scaling` correction is required (unlike `calibration_transducer`, which derives amplitude indirectly via the O'Neil analytical solution).

### Prerequisites

- Phases must already be set correctly in the config (e.g., loaded from a `calibration_transducer` output YAML).
- `modules.run_source_setup` must be enabled (the calibration water sim reuses the source built in that stage).
- The sensor mask must cover the focal region so that `max(sensor_data.p_max_all)` captures the true focal peak. This is guaranteed when using the standard PRESTUS sensor setup.

### Configuration

In your study config (or default_config.yaml):

```yaml
calibration:
  target_isppa_wcm2: 30   # Target free-water ISPPA [W/cm²]

modules:
  run_amplitude_calibration: 1   # Enable in-pipeline amplitude calibration
```

`calibration.target_isppa_wcm2` is required when `run_amplitude_calibration = 1`. The pipeline will error if it is not set.

The target intensity is converted to a peak pressure internally via $p = \sqrt{2 \cdot I \cdot 10^4 \cdot \rho \cdot c}$ using water properties from `medium_properties.water`.

### Output

The calibration step writes two fields back into `parameters.calibration`:

| Field | Description |
|---|---|
| `calibration.measured_isppa_wcm2` | Free-water ISPPA measured at the current `elem_amp` before scaling [W/cm²] |
| `calibration.amplitude_scale_factor` | The multiplicative factor applied to `elem_amp` |

Both fields are available to the simulation report for traceability.

### Runtime cost

The calibration water simulation runs with `matlab_cpu` and uses a homogeneous medium (no tissue mapping), so it is substantially faster than the main simulation. On a typical workstation it adds a few minutes. On HPC, it runs locally on the submission node before the main job is dispatched — set `simulation.code_type` to a GPU or C++ backend for the main sim without concern; the calibration always uses CPU.

### Relationship to `calibration_transducer`

The typical workflow is:

1. **Once per transducer/depth setup**: Run `calibration_transducer` to obtain optimized phases and an initial amplitude at a reference intensity. The output YAML encodes `elem_phase_deg` and `elem_amp`.
2. **Per subject simulation**: Load the calibration YAML into your study config (overrides default `elem_amp` and `elem_phase_deg`). Enable `run_amplitude_calibration = 1` and set `calibration.target_isppa_wcm2` to fine-tune the amplitude to the exact free-water pressure you want for that simulation.
