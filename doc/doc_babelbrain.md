# BabelBrain backend

PRESTUS can use [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) as an alternative backend for acoustic simulation, thermal simulation, or both. By default PRESTUS uses k-Wave for both; the BabelBrain backends can be enabled independently. All pre-processing stages (segmentation, grid setup, medium masks, transducer placement) and all downstream post-processing (focal metrics, NIfTI export, HTML report) run unchanged regardless of which backend is active.

**What is BabelBrain?** BabelBrain is an open-source tool for transcranial focused ultrasound planning developed by Samuel Pichardo (University of Calgary). It uses a different acoustic solver called BabelViscoFDTD, which models the full viscoelastic wave equation — including shear waves in bone — rather than the scalar linear acoustic model used by k-Wave. BabelBrain also includes a Rayleigh integral step for computing the field through the water coupling layer before the skull, which can improve accuracy at the skull boundary. It supports GPU acceleration via CUDA, OpenCL, and Metal (Apple Silicon).

The main reasons to use BabelBrain instead of (or alongside) k-Wave are:

- **Independent solver, optional viscoelastic physics.** BabelViscoFDTD can model shear waves in bone. PRESTUS runs it in **fluid-bone mode by default** (shear off, PRESTUS tissue properties injected) so results are comparable to k-Wave/jwave; see [Tissue properties](#tissue-labels-and-material-properties).
- **Cross-validation.** Comparing results from two independent solvers is a common check for simulation-based studies.
- **Apple Silicon GPU acceleration.** On Apple Silicon, k-Wave's fast C++ binary (`kspaceFirstOrder-CUDA`) is unavailable because it requires NVIDIA CUDA; only the slower `matlab_gpu` path works. BabelBrain has a native Metal backend that makes full use of the M-series GPU.

| | PRESTUS (k-Wave) | BabelBrain |
|---|---|---|
| **Solver** | Pseudo-spectral FDTD | Explicit staggered-grid FDTD |
| **Physics** | Scalar linear acoustic (pressure/density) | Viscoelastic solver; run with shear = 0 (fluid bone) in PRESTUS |
| **Skull model** | Scalar ρ, c, α | Scalar ρ, c, α from PRESTUS `medium_properties` (or BabelBrain HU mapping via `ct_mapping`) |
| **Rayleigh pre-step** | None | Rayleigh integral for the water coupling layer |
| **Narrow-beam domain** | Optional (see [limited FOV](doc_limited_fov.md)) | Off by default (the pre-focal field would not be returned); enable with `babelbrain.tight_narrow_beam: true`. Always off when `grid.acoustic_fov_length_mm` is set (mask already cropped) |
| **GPU backends** | CUDA, MATLAB GPU arrays | CUDA, OpenCL, Metal (Apple Silicon) |

See [Simulation backends](doc_backend.md) for the k-Wave-only backend options.

> **Tissue properties and comparability.**
> On the default label-mask path, PRESTUS writes `medium_properties` to `<target>_medium_properties.json` and the Python bridge overwrites BabelBrain's built-in acoustic table (`MatFreq`) with those values: density, longitudinal sound speed, and attenuation converted from dB/(cm·MHz^y) to Np/m. Shear is set to 0. BabelBrain therefore solves with the **same tissue acoustics as k-Wave/jwave**. Remaining differences are the solver and BabelBrain's bone Q-correction, so validate against a k-Wave run before relying on absolute values.
> This does **not** apply to (a) the `ct_mapping` path, which uses BabelBrain's own HU models; (b) runs where writing the JSON failed, which print a warning and fall back to BabelBrain's built-in table (listed below); or (c) **thermal** properties. BabelBrain thermal uses its own table, which differs from PRESTUS in absorption fraction (see [Thermal simulation](#thermal-simulation)).

---

## Quick start

Both the acoustic and thermal simulations can be run in BabelBrain independently or together. The acoustic backend replaces k-Wave for pressure field computation; the thermal backend replaces k-Wave's BHTE solver and uses BabelBrain's own Pennes BHTE on the same simulation domain. All downstream PRESTUS post-processing (focal metrics, CEM43, NIfTI export, HTML report) runs unchanged in both cases.

**Acoustic only (default when using BabelBrain):**

```yaml
acoustic_backend: babelbrain

babelbrain:
  tx_system: R15473            # BabelBrain preset name; auto-detected if empty
  computing_backend: auto      # 'auto' | 'cuda' | 'opencl' | 'metal' | 'cpu'  (or 1-4)
  device_name: auto            # 'auto' selects device based on computing_backend; or set explicitly e.g. 'Apple'
  ppw: 12                      # points-per-wavelength floor (default 12); actual PPW = max(ppw, round(c_brain / (freq × res_mm × 1e-3))). See note below.
  do_refocusing: false         # skull back-propagation refocusing
  z_beyond_focus_mm: 40        # axial domain truncation past the focal point [mm]
  python_env: ''               # path to conda/venv environment; empty = system python3
  root: ''                     # path to BabelBrain source; auto-detected if empty
  tx_mech_adjust_z_m: 0.0      # ADDITIONAL transducer z offset into skin [m]; a placement
                               # correction (preset FL − |focus_pos − trans_pos|) is always added (see below)
  z_steering_m: ~              # ZSteering [m]; null = auto-derived (0 with explicit phases; see below)
  focal_calibration: true              # geometric fallback only: steer so the achieved FWHM-centre lands at focal_distance_ep;
                                       # locate-only (report) when explicit element phases exist (see below)
  focal_calibration_target: flhm_center # quantity placed at focal_distance_ep: 'flhm_center' | 'peak'
```

> **PPW and simulation resolution.** `ppw` sets a floor; PRESTUS resolves the actual PPW as `max(ppw, round(c_brain / (freq × res_mm × 1e-3)))` where `c_brain` is taken from `medium_properties.brain.sound_speed` (default 1500 m/s). BabelBrain internally sets its own spatial step to `SmallestSOS / (PPW × freq)`, where `SmallestSOS` is the minimum speed of sound across all tissue types in BabelBrain's own tables — including a built-in bone shear-speed reference (`DensityToSSoSPichardo`, ≈1102.5 m/s) that is normally the binding constraint under the fluid-bone model PRESTUS uses. This native spatial step is generally finer than `grid.resolution_mm`. BabelBrain reuses the input mask's voxel *array shape* directly as its own FDTD grid dimensions without resampling to this spatial step, so `write_babelbrain_input.m` resamples the mask (nearest-neighbor) from `grid.resolution_mm` to BabelBrain's own native resolution before writing it — this keeps the physical domain size BabelBrain simulates consistent with the PRESTUS grid, with no manual resolution matching required.

**Acoustic + thermal (both in BabelBrain):**

```yaml
acoustic_backend: babelbrain
thermal_backend:  babelbrain   # use BabelBrain Pennes BHTE instead of k-Wave

babelbrain:
  tx_system: R15473
  computing_backend: auto      # 'auto' | 'cuda' | 'opencl' | 'metal'  (thermal has no CPU backend)
  device_name: auto            # 'auto' selects device based on computing_backend
  # isppa_w_per_cm2: ~         # leave empty: derived from the scaled acoustic field (in-situ brain peak)
```

> **`thermal_backend: babelbrain`** requires `acoustic_backend: babelbrain`, because BabelBrain's thermal solver operates on the `*DataForSim.h5` file produced by the acoustic stage. Other combinations are rejected when the config is loaded.

PRESTUS searches for BabelBrain at `external/BabelBrain` (git submodule) automatically. Set `babelbrain.root` explicitly to override.

> **`tx_system`** is the name of the transducer preset inside BabelBrain (e.g., `R15473`, `CTX_500`). If left empty, PRESTUS attempts to match the transducer geometry to a known preset automatically. See [Transducer mapping](#transducer-mapping).

---

## Installation

BabelBrain is included as a git submodule at `external/BabelBrain`. Initialise it after cloning PRESTUS:

```bash
git submodule update --init external/BabelBrain
```

BabelBrain requires a dedicated Python environment — its compiled extension `BabelViscoFDTD` cannot be installed into a generic system Python. Use the conda environment file bundled with BabelBrain:

```bash
# Apple Silicon (Metal)
conda env create -f external/BabelBrain/environment_mac_arm64-314.yml
conda activate BabelBrain

# macOS Intel
conda env create -f external/BabelBrain/environment_mac_x64-314.yml
conda activate BabelBrain

# Linux (CUDA)
conda env create -f external/BabelBrain/environment_linux.yml
conda activate BabelBrain

# Windows
conda env create -f external/BabelBrain/environment_win-314.yml
conda activate BabelBrain
```

After activation, find the environment path and set it in your config:

```bash
conda activate BabelBrain
which python   # e.g. /opt/conda/envs/BabelBrain/bin/python
```

```yaml
babelbrain:
  python_env: /opt/conda/envs/BabelBrain   # directory containing bin/python
```

`babelbrain.python_env` **must** be set — leaving it empty falls back to system `python3`, which will not have `BabelViscoFDTD` and will fail with `No module named 'BabelViscoFDTD'`.

---

## BabelBrain output folder

When the BabelBrain backend is active, or when `io.export_babelbrain: true` is set, PRESTUS writes a dedicated output folder alongside the standard output directories:

```
<sim_output>/
  nii/          ← NIfTI outputs (unchanged)
  log/          ← logs and parameter YAMLs
  cache/        ← .mat cache files
  babelbrain/   ← BabelBrain input package
```

The `babelbrain/` folder contains everything needed to run or re-run BabelBrain independently:

| File | Description |
|---|---|
| `<target>_<freq>kHz_<PPW>PPW_BabelViscoInput.nii.gz` | Tissue mask in BabelBrain integer label format, including the focal voxel marker |
| `<target>_CalibWeights.h5` | Per-element calibration phases from PRESTUS calibration |
| `RunBabelBrain.py` | Ready-to-run Python script with all simulation parameters pre-filled |
| `babelbrain_config.yaml` | Parameter record for reproducibility |

> **Exporting inputs without switching the acoustic backend.** Set `io.export_babelbrain: true` with `acoustic_backend: kwave` (or no `acoustic_backend` key) to generate the BabelBrain input package while still using k-Wave for the simulation. This is useful for preparing a cross-validation or for running BabelBrain on a different machine.

---

## Tissue labels and material properties

### How tissue labels are passed to BabelBrain

BabelBrain identifies tissues from a single NIfTI file (`BabelViscoInput.nii.gz`) using fixed integer labels. PRESTUS writes this file automatically from its own medium mask, remapping layer names to BabelBrain's scheme:

| PRESTUS layer name | Default config indices | BabelBrain label | BabelBrain material |
|---|---|---|---|
| `water` | 0, 3, 6, 9, 10 | 0 | Water |
| `skin` | 5 | 1 | Skin |
| `skull_cortical` | 7 | 2 | Cortical bone |
| `skull_trabecular` | 8 | 3 | Trabecular bone |
| `skull` (undifferentiated) | 4 | 2 | Cortical bone (fallback) |
| `brain` | 1, 2 | 4 | Brain |
| *(focal voxel)* | — | 5 | — |

The focal voxel (label 5) is placed at the focus grid index. It is the only mechanism BabelBrain provides for specifying the target location.

### CT-based skull properties (`ct_mapping`)

By default PRESTUS passes a label mask only and BabelBrain applies its internal `MatFreq` literature values for each tissue type. Optionally, a pseudo-CT (HU) volume can be passed so that BabelBrain derives per-voxel skull properties from Hounsfield units using one of its fitted mapping models.

Set `parameters.babelbrain.ct_mapping` to a non-empty string to activate this path:

```yaml
babelbrain:
  ct_mapping: Webb-Marsac   # one of the values in the table below
  ct_use_petra: false        # true only for PETRA-sequence-derived pCT
  water_amp_calibration: true # scale output to the analytic Rayleigh free-water peak (see below)
```

**When to use this path.** Only when you want BabelBrain's own HU→acoustic-property conversion to apply to raw UTE or PETRA pseudo-CT Hounsfield values. This is the equivalent of BabelBrain's built-in `bUseCT` mode.

**File format.** BabelBrain's CT path does not accept a raw HU volume directly — it expects `CT.nii.gz` to already be a quantized index map (small integers indexing into a sorted table of unique HU levels present in the cortical/trabecular bone region), with that table saved alongside as `CT-cal.npz`. `write_babelbrain_input.m` writes the raw-HU `CT.nii.gz` (resampled to BabelBrain's own native resolution, identically to `BabelViscoInput.nii.gz`, so the two files share a grid); `run_babelbrain_case.py` then quantizes it (uniform round-to-level quantization, 10-bit, matching BabelBrain's own default) and writes the matching `CT-cal.npz` before invoking `RunCases()`. `bExtractAirRegions` is disabled for this path, since PRESTUS does not generate the separate air-region mask BabelBrain's own segmentation pipeline would normally supply.

> **Do not set `ct_mapping` when PRESTUS's pCT pipeline (`pct.enabled: 1`) has already converted HU to acoustic properties.** In that case, the mapping has been applied inside PRESTUS and k-Wave is already using voxel-wise values; BabelBrain receives only a label mask and uses its own `MatFreq` values for bone, which differ from the pCT-derived values k-Wave uses. There is currently no path to pass PRESTUS's pre-mapped per-voxel acoustic properties into BabelBrain (see [Limitations](#limitations-and-unsupported-prestus-modes)).

#### Supported `ct_mapping` values

The following HU→property mapping methods are available in BabelBrain. PRESTUS does not replicate all of them.

| `ct_mapping` value | Density | SoS | Attenuation | PRESTUS equivalent |
|---|---|---|---|---|
| `Webb-Marsac` (default) | Marsac linear HU→ρ | Webb 2018 tabulated (freq-dep.) | Webb 2020 tabulated (freq-dep.) | partial: `marsac` density + no Webb |
| `Aubry` | Porosity→ρ (Aubry 2003) | Porosity→SoS | Porosity→Att | `aubry` |
| `Pichardo` | HU→ρ air–tissue linear | Density→SoS (Pichardo 2D map, freq-dep.) | Density→Att (Pichardo map) | none |
| `McDannold` | HU→ρ air–tissue linear | Density→SoS (McDannold polynomial) | Density→Att (McDannold polynomial, freq-dep.) | none |
| `Marsac-Aubry` | Marsac linear HU→ρ | Porosity→SoS (Aubry) | Porosity→Att (Aubry) | partial |
| `Pichardo-Marsac` | Marsac linear HU→ρ | Pichardo 2D map (freq-dep.) | Pichardo map | none |
| `McDannold-Marsac` | Marsac linear HU→ρ | McDannold polynomial | McDannold polynomial (freq-dep.) | none |

The `Pichardo`, `McDannold`, and their `-Marsac` hybrids are **not replicated in PRESTUS** (`medium_pct_*` functions support `k-plan`, `k-wave`, `marsac`, and `aubry` only). If you want results comparable across backends, use `Webb-Marsac` (closest to the internal PRESTUS `marsac`/`mueller` combination) or `Aubry` (directly available in both).

#### PETRA pCT density formula

When `ct_use_petra: true`, BabelBrain uses Björn's formula for the density step inside `Webb-Marsac` mode (`SimNIBS_PETRApct_Density`). This is appropriate when the pseudo-CT HU values were generated from a PETRA sequence via SimNIBS's ZTE/PETRA pipeline. For UTE-derived pCT (the default PRESTUS pCT source), leave `ct_use_petra: false`.

#### UTE image as direct BabelBrain CT input

Providing the raw UTE NIfTI directly (rather than a pre-processed pseudo-CT) is technically supported by setting `ct_mapping` and passing a UTE volume as `pseudoCT`. The UTE HU-like intensities will then be fed into BabelBrain's mapping pipeline. However:

- UTE images do not use standard 120 kVp CT Hounsfield calibration, so `Webb-Marsac` and `McDannold` models (which are calibrated to conventional CT HU) will produce inaccurate results.
- `Aubry` (porosity-based) is less sensitive to the absolute HU scale and is the most defensible choice for raw UTE input.
- The recommended practice for UTE data is to run PRESTUS's pCT pipeline first (UTE → calibrated pseudo-CT HU → PRESTUS `marsac`/`aubry` mapping), which applies bias correction and skull-specific calibration before any HU→property conversion.

### BabelBrain material properties

BabelBrain's built-in database holds frequency-dependent material properties from literature fits, tabulated every 5 kHz from 100 to 1125 kHz. **On the default PRESTUS path this table is overwritten with `medium_properties`**. It is listed here for reference and applies only if material injection fails.

The tables below show example values at three frequencies commonly used in TUS (source: `TranscranialModeling/BabelIntegrationBASE.py`).

#### Density [kg/m³] — frequency-independent

| Tissue | Density |
|---|---|
| Water | 1000 |
| Skin | 1116 |
| Cortical bone | 1897 |
| Trabecular bone | 1738 |
| Brain | 1041 |
| White matter | 1041 |
| Gray matter | 1045 |
| CSF | 1007 |

#### Longitudinal speed of sound [m/s]

| Tissue | 300 kHz | 500 kHz | 650 kHz |
|---|---|---|---|
| Water | 1500 | 1500 | 1500 |
| Skin | 1537 | 1537 | 1537 |
| Brain | 1562 | 1562 | 1562 |
| White matter | 1537 | 1537 | 1537 |
| Gray matter | 1520 | 1520 | 1520 |
| CSF | 1507 | 1507 | 1507 |
| Cortical bone | **2452** | **2476** | **2494** |
| Trabecular bone | **2148** | **2205** | **2247** |

Bone longitudinal SoS is fitted as a linear function of frequency from two reference measurements at 270 kHz and 836 kHz (Pichardo 2011).

#### Shear speed of sound [m/s] — bone only; soft tissue = 0

| Tissue | 300 kHz | 500 kHz | 650 kHz |
|---|---|---|---|
| Cortical bone | 1473 | 1542 | 1594 |
| Trabecular bone | 1270 | 1313 | 1345 |
| All soft tissues | 0 | 0 | 0 |

Shear SoS is only non-zero for bone. This is what would let BabelBrain model mode conversion at the skull. PRESTUS sets bone shear to 0 (the `_patch_noncat_skull_attenuation` bridge patch), both to match the fluid solvers and because the non-CT shear path otherwise realises almost no skull attenuation.

#### Longitudinal attenuation [Np/m]

| Tissue | 300 kHz | 500 kHz | 650 kHz |
|---|---|---|---|
| Water | 0 | 0 | 0 |
| Skin | 1.4 | 2.3 | 3.0 |
| Brain | 2.1 | 3.5 | 4.5 |
| White matter | 3.1 | 5.1 | 6.6 |
| Gray matter | 1.3 | 2.2 | 2.9 |
| CSF | 0.03 | 0.05 | 0.06 |
| Cortical bone | 49 | 81 | 106 |
| Trabecular bone | 49 | 81 | 105 |

#### Shear attenuation [Np/m] — bone only

| Tissue | 300 kHz | 500 kHz | 650 kHz |
|---|---|---|---|
| Cortical bone | 99 | 164 | 214 |
| Trabecular bone | 99 | 164 | 214 |
| All soft tissues | 0 | 0 | 0 |

> **These built-in values are not used on the default path.** PRESTUS replaces them with `medium_properties` (see the note at the top of this page). They apply only when material injection fails, and they serve as a reference for how far the literature fits differ from PRESTUS defaults.
>
> **Comparison with PRESTUS defaults.** BabelBrain cortical bone longitudinal attenuation (81 Np/m at 500 kHz) is approximately 47% lower than the PRESTUS default (152 Np/m, derived from 13.3 dB/cm/MHz). BabelBrain cortical longitudinal SoS (2476 m/s at 500 kHz) is about 11% lower than the PRESTUS skull default (2800 m/s). Brain SoS: BabelBrain 1562 m/s vs. PRESTUS 1610 m/s (~3% difference). These differences would affect absolute pressure and heating values if the built-in table were used. See [Medium properties](doc_medium.md) for PRESTUS defaults.

<span style="font-size: 12px;">
  <strong>References for BabelBrain fitting functions:</strong><br>
  Pichardo, S. et al. Multi-frequency characterization of the speed of sound and attenuation coefficient for longitudinal transmission of freshly excised human skulls. <em>Phys Med Biol</em> <strong>56</strong>, 219–250 (2011). https://doi.org/10.1088/0031-9155/56/1/014<br>
  Pichardo, S. et al. In-vitro characterization of multi-frequency ultrasound transducers for skull-induced phase aberrations. <em>Phys Med Biol</em> <strong>62</strong>, 6938–6962 (2017). https://doi.org/10.1088/1361-6560/aa7ccc<br>
  Goss, S.A., Johnston, R.L. & Dunn, F. Comprehensive compilation of empirical ultrasonic properties of mammalian tissues. <em>J Acoust Soc Am</em> <strong>64</strong>, 423–457 (1978). https://doi.org/10.1121/1.382016<br>
  Labuda, C. et al. Speed of Sound and Attenuation in Ex-Vivo Human Brain Tissue Sections. <em>IEEE TUFFC</em> <strong>69</strong>, 1280–1291 (2022). https://doi.org/10.1109/TUFFC.2022.3161443
</span>

---

## Transducer mapping

BabelBrain only accepts named transducer presets. `babelbrain_transducer_map` resolves the PRESTUS transducer to a preset by:

1. **Serial-number lookup** — matches `tran.serial` from the equipment YAML against a built-in table.
2. **Geometry matching** — if no serial match, compares element count, radius of curvature, maximum outer diameter, and frequency against known preset specifications (tolerances: ±2 mm radius, ±3 mm OD, ±5% frequency).
3. **Custom CSV generation** — if no preset matches, writes an element CSV in BabelBrain annular-array format and prints step-by-step instructions for adding a new preset to BabelBrain.

Known mappings (verified against `BabelBrain/Babel_*/default.yaml`):

| PRESTUS serial / geometry | BabelBrain preset | Rings | R_curv [mm] | OD [mm] | Freq |
|---|---|---|---|---|---|
| IS_PCD15473_01001 | `R15473` | 10 | 100.0 | 65.0 | 300 kHz |
| IS_PCD15287_01001 | `R15287` | 4 | 62.9 | 64.0 | 500 kHz |
| IS_PCD15646_01001 | `R15646` | 10 | 72.0 | 72.0 | 500 kHz |
| CTX500 geometry | `CTX_500` | 7 | 62.9 | 64.0 | 500 kHz |
| CTX250 geometry | `CTX_250` | 7 | 62.9 | 64.0 | 250 kHz |
| Single-element bowl | `H317` / `H246` / `H301` | 1 | varies | varies | varies |
| IGT 64-channel phased array | `IGT64_500` | 64 | — | — | 500 kHz |

Set `babelbrain.tx_system` explicitly to bypass auto-detection.

---

## Calibration weights

When PRESTUS has explicit element phases for the transducer, these are automatically included in the BabelBrain input package as an `OptimizedWeightsFile`. BabelBrain applies these weights on top of ring phases it derives internally (Rayleigh back-propagation from its aim point `FocalLength + ZSteering`), so the file stores **deltas** relative to those internal phases:

```
w[n] = exp(i × (elem_phase[n] − geometric_phase_at_aim[n]))
```

The delta reference is the geometric phase at BabelBrain's actual aim. In the default explicit-phases mode (`ZSteering = 0`) the aim is the natural focus, where the reference is flat, so the delta carries the full drive and BabelBrain reproduces the config phases exactly (see "Transducer placement, ZSteering, and the drive phases" below). The Python bridge handles a format conversion at runtime (MATLAB writes split real/imag datasets; BabelBrain expects a native complex64 HDF5 dataset).

If no calibrated phases are available, BabelBrain falls back to geometric (Rayleigh) focusing, which is equivalent to the PRESTUS uncalibrated mode.

### Free-water amplitude calibration

BabelBrain's saved pressure amplitude is its raw FDTD field multiplied by a built-in empirical CFL/dispersion correction (`Correction = 100/(100 − ExpectedError(CFL))`, then ×√2, in `BabelIntegrationBASE`). At the operating CFL (~0.5 — the timestep is set by the slowest material in BabelBrain's table, bone, even in water-only mode) this polynomial carries a residual bias of a few % in pressure relative to the exact Rayleigh solution for the same drive (observed: ≈+6% pressure / +13% intensity for R15473 at 300 kHz).

With `babelbrain.water_amp_calibration: true` (default) the wrapper calibrates this out: it evaluates the in-repo Rayleigh–Sommerfeld model ([`rayleigh_axial_intensity.m`](../functions/calibration/rayleigh_axial_intensity.m)) for the exact per-ring drive, forms `k = p_Rayleigh_peak / p_BabelBrain_water_peak` from BabelBrain's own free-water field (always available — BabelBrain computes a water reference in every case), and scales `sensor_data.p_max_all` (and `freefield_isppa_wcm2` by `k²`). The factor is logged and stored in `sensor_data.amp_calibration`.

The water-derived scalar is valid for tissue runs because (a) the bias originates in the source-injection/FDTD-correction chain, which is identical in water-only and transcranial runs (same grid, same bone-limited timestep, same correction factor), and (b) the simulation is linear, so the entire field — including through skull — scales with the source. This is the numerical analogue of hydrophone free-water calibration of a physical transducer. A factor deviating >15% from 1 triggers a warning, since that magnitude indicates a setup mismatch (geometry/phases/source amplitude) rather than solver bias.

Caveats: requires explicit element phases (skipped otherwise). The raw `DataForSim.h5` consumed by `thermal_backend: babelbrain` is not rescaled, and does not need to be. BabelBrain's Pennes BHTE renormalizes the field so that the **peak ISPPA inside the brain** equals the `Isppa` it receives (`PressureRatio = sqrt(Isppa·2ρc)/max(p_amp[brain])` in `CalculateTemperatureEffects.py`), so any global amplitude scalar cancels there. PRESTUS passes as `Isppa` the in-brain peak `p²/(2ρc)` of the calibrated, ISPPA-scaled `sensor_data.p_max_all` (`babelbrain_insitu_isppa.m`), which is the same field k-Wave thermal consumes. The two thermal backends therefore deposit heat from the same acoustic amplitude. (Before this fix, the free-water calibration ISPPA was passed, which dropped the transcranial loss and overestimated BabelBrain heating.)

---

## Thermal simulation

### BabelBrain thermal backend

When `thermal_backend: babelbrain` is set, PRESTUS calls BabelBrain's `CalculateThermalProcess` on the `*DataForSim.h5` acoustic output. BabelBrain solves the Pennes bio-heat equation on the same simulation domain using its own GPU-accelerated relaxation solver.

The timing parameters from `parameters.timing` are converted automatically:

| PRESTUS parameter | BabelBrain parameter | Conversion | Notes |
|---|---|---|---|
| `timing.pd` / `timing.pri` | `DutyCycle` | `pd / pri` | |
| `timing.pri` | `PRF` | `1 / pri` [Hz] | |
| `timing.ptd` | `DurationUS` | direct [s] | ON duration of one pulse train |
| `timing.ptri - timing.ptd` | `DurationOff` | direct [s] | Inter-train cooling interval; 0 if `ptri ≤ ptd` |
| `timing.n_ptri_reps` | `Repetitions` | rounded | Number of (on + off) cycles |
| `timing.post_ptri_dur` | `PauseBetweenGroupedSonications` | direct [s]; 0 if unset/NaN | Cooldown after the last repetition (`NumberGroupedSonications = 1`) |
| `thermal.temp_0.brain` | `BaselineTemperature` | direct [°C] | |
| in-brain peak of `sensor_data.p_max_all` | `Isppa` | `max p²/(2ρc)` over brain [W/cm²] | Override with `babelbrain.isppa_w_per_cm2`; see ISPPA note below |

> **`DurationUS` and multi-PTRI protocols.** `DurationUS` is the ON duration of one pulse train (`timing.ptd`), not the total `ptrd`. BabelBrain runs `Repetitions = n_ptri_reps` cycles of `DurationUS` on + `DurationOff` off, so cooling between trains is modelled. With a single sonication group, BabelBrain then runs `PauseBetweenGroupedSonications` once more with no heat source, starting at `ptrd` — the same point where the k-Wave thermal path starts `post_ptri_dur`. `TempEndFUS`/`DoseEndFUS` (end of the last ON period) and `FinalTemp`/`FinalDose` (after the cooldown) are read into the end-of-sonication and final maps respectively. BabelBrain's thermal step is fixed at 10 ms, so run time scales with the total simulated duration including the cooldown.

> **ISPPA normalisation.** BabelBrain always rescales the entire pressure field so that the peak ISPPA in brain equals the requested `Isppa` (`PressureAdjust = sqrt(Isppa × 2 × ρ × c)`). This is an **in-situ** quantity, not the free-water calibration ISPPA. PRESTUS derives it automatically from the scaled acoustic field and logs it next to the free-water value and the implied transmission. Set `babelbrain.isppa_w_per_cm2` only to force a specific in-situ brain ISPPA.

> **`sel_p` and refocusing.** When `babelbrain.do_refocusing: true`, the thermal simulation automatically uses the `/p_amp_refocus` pressure field (the skull-corrected map) rather than `/p_amp`. No manual configuration is needed.

### Thermal tissue properties

!!! warning "Known limitation: BabelBrain thermal uses its own tissue thermal properties"
    `thermal_backend: babelbrain` ignores the thermal fields of PRESTUS `medium_properties` (`thermal_conductivity`, `specific_heat_capacity`, `perfusion`, `absorption_fraction`), including uncertainty-mode variants. BabelBrain has no input for them, so PRESTUS uses BabelBrain's values unchanged. The acoustic field amplitude is matched between the thermal backends; the tissue thermal properties are not. Use `thermal_backend: kwave` when PRESTUS thermal properties must apply.

BabelBrain thermal does not read PRESTUS `medium_properties`. On the label-mask path it uses a built-in table (`CalculateTemperatureEffects.py`):

| Tissue | Conductivity [W/m/°C] | Specific heat [J/kg/°C] | Perfusion | Absorption fraction | PRESTUS absorption fraction |
|---|---|---|---|---|---|
| Water | 0.60 | 4178 | 0 | 0 | 1 |
| Skin | 0.37 | 3391 | 106 | 0.85 | 1 |
| Cortical bone | 0.32 | 1313 | 10 | 0.16 | 0.28 |
| Trabecular bone | 0.31 | 2274 | 30 | 0.15 | 0.28 |
| Brain | 0.51 | 3630 | 559 | 0.85 | 1 |

Conductivity, specific heat and perfusion match the PRESTUS defaults, apart from trabecular conductivity (0.31 vs 0.32). The absorption fraction scales heat deposition in both solvers, so BabelBrain thermal deposits about 15% less heat in skin and brain and roughly 45% less in skull than k-Wave thermal for the same acoustic field. Take this into account when comparing temperature rise or CEM43 across thermal backends. BabelBrain also does not compute the ISO CEM43 fields (returned as NaN) or a focal-plane time series.

On the CT path (`ct_mapping`), BabelBrain uses the same soft-tissue values and assigns all bone voxels the mean of the cortical and trabecular values (specific heat 1793.5, conductivity 0.315, perfusion 20, absorption fraction 0.155).

### Post-processing with PRESTUS

After BabelBrain thermal runs, `read_babelbrain_thermal` reads the output, converts temperature rise and CEM43 fields to PRESTUS conventions, and resamples the volumes from BabelBrain's simulation grid back into the PRESTUS acoustic grid. The returned `results_heating` struct is identical in structure to the output of `thermal_simulation`, so:

- `thermal_analysis` computes focal temperature, CEM43 metrics, and safety checks unchanged
- `nifti_thermal` exports temperature rise, CEM43, and heat-rise NIfTIs to the standard `nii/` folder
- The HTML report includes BabelBrain thermal results alongside or instead of k-Wave results

The resampling uses the `planimg` affine to map BabelBrain's RAS-mm coordinate system to the PRESTUS acoustic grid.

### Combination modes

The acoustic and thermal backends can be mixed independently:

| `acoustic_backend` | `thermal_backend` | Scenario |
|---|---|---|
| `kwave` (default) | `kwave` (default) | Standard PRESTUS, no BabelBrain |
| `babelbrain` | `kwave` | BabelBrain acoustics, PRESTUS thermal |
| `babelbrain` | `babelbrain` | Full BabelBrain pipeline, PRESTUS post-processing |

---

## Limitations and unsupported PRESTUS modes

Not all PRESTUS pipeline features are currently compatible with the BabelBrain backend. The table below summarises known limitations.

| PRESTUS feature | BabelBrain support | Details |
|---|---|---|
| **Single transducer** | ✓ Full | Standard use case. |
| **Multi-transducer — async/incoherent** (`simulation.transducer_coupling: async`) | ✓ Supported | `async_transducer_pipeline` runs one independent acoustic simulation per transducer (each a single-transducer BabelBrain call), then combines intensity fields incoherently via `combine_async_intensity`. `acoustic_backend: babelbrain` is inherited by each per-transducer sub-run automatically. |
| **Multi-transducer — coherent** (`simulation.transducer_coupling: coherent`) | ✓ Supported | `babelbrain_wrapper` detects `numel(parameters.transducer) > 1` and recurses once per transducer with a unique `output_affix` (`_t01`, `_t02`, …). Each sub-run reads `/p_complex` (the complex Fourier pressure field stored by BabelBrain alongside `/p_amp`) from its `DataForSim.h5`. The complex fields are summed coherently; `p_max_all = abs(Σ p_complex_i)`. Unlike k-Wave's native coherent mode (all transducers in one simulation), this treats transducer interactions as independent and ignores mutual acoustic coupling between transducers. |
| **Limited acoustic FOV** (`grid.acoustic_fov_*`) | ✓ Supported | The NIfTI passed to BabelBrain is the same beam-centred, FOV-cropped box used for k-Wave/jwave (`acoustic_fov_crop.m`/`acoustic_grid_fov.m`). BabelBrain's own narrow-beam domain cropping (`bTightNarrowBeamDomain`) is disabled in this case — it would crop the already-cropped box a second time, and its ray-traced lateral footprint estimate has been found to under-shoot the true beam width near the skull. BabelBrain simulates the full FOV-cropped box with no further internal shrink, so its domain matches the other backends' voxel-for-voxel. |
| **Axisymmetric 2D** (`simulation.axisymmetric: 1`) | ✗ Not supported | BabelBrain is always a 3D solver. Axisymmetric mode must use the k-Wave backend. |
| **`multi_isppa` mode** | ⚠ Untested | `multi_isppa_pipeline` derives each acoustic/thermal job from the full parameter struct, so `acoustic_backend`/`thermal_backend` carry over. Not yet validated end-to-end with BabelBrain. |
| **`phantom` / `water` test medium** | ⚠ Untested | BabelBrain expects a realistic segmented NIfTI. Purely water or phantom media are geometrically valid but BabelBrain's material property tables will apply skull properties to any non-background label, which may produce unexpected results. |
| **Post-hoc water simulation** | ⚠ Untested | The post-hoc water loop re-enters the pipeline with `simulation.medium: water` and inherits `acoustic_backend`, so with BabelBrain it runs a BabelBrain water-only simulation. |
| **Thermal tissue properties** | ⚠ Known limitation | BabelBrain thermal uses its built-in conductivity, specific heat, perfusion and absorption fraction; PRESTUS `medium_properties` thermal fields (and their uncertainty variants) are not applied. The absorption fraction differs most (skin/brain 0.85 vs 1.0, skull ≈0.15 vs 0.28), so BabelBrain heating is lower than k-Wave's for the same field, mostly in skull. See [Thermal tissue properties](#thermal-tissue-properties). |
| **BabelBrain thermal without BabelBrain acoustic** | ✗ Not supported | `thermal_backend: babelbrain` requires the `*DataForSim.h5` file produced by `acoustic_backend: babelbrain`. It cannot consume a k-Wave acoustic result. |
| **CT-based skull properties into BabelBrain** | ✓ Supported (via `ct_mapping`) | Set `babelbrain.ct_mapping` to a non-empty mapping method (e.g. `'Webb-Marsac'`) and pass the pseudoCT HU volume as the 8th argument to `babelbrain_wrapper`. `babelbrain_wrapper` writes a `CT.nii.gz` alongside `BabelViscoInput.nii.gz` and passes `--use_ct`, `--ct_mapping`, and (optionally) `--use_petra` flags to `run_babelbrain_case.py`, which forwards them to BabelBrain's `RUN_SIM().RunCases()`. BabelBrain then performs its own HU→property conversion internally. **Note:** use this path only when passing raw HU values from a UTE/PETRA scan. If PRESTUS's pCT pipeline (`pct.enabled: 1`) has already converted HU to acoustic properties, use the default label-based path (`ct_mapping` empty) — the per-voxel acoustic properties k-Wave uses cannot be passed to BabelBrain. See [CT-based skull properties](#ct-based-skull-properties-ct_mapping) for available mapping methods. |
| **Uncertainty mode** | ⚠ Untested | The liberal/conservative configs vary `medium_properties`. These are forwarded to BabelBrain through the material-injection JSON (density, sound speed, attenuation), so variants should reach the acoustic solve on the label-mask path. They have no effect with `ct_mapping`, and never on BabelBrain **thermal** properties. Not yet validated end-to-end. |

---

## Developer notes: PRESTUS integration architecture

> This section describes how the BabelBrain integration is wired inside PRESTUS. It is intended for developers working on the backend or extending GPU platform support.

### Call chain

```
prestus_pipeline.m
  └─ acoustic_stage.m            ← dispatches on parameters.acoustic_backend (run_acoustic_backend)
       └─ babelbrain_wrapper.m   ← MATLAB orchestrator
            ├─ babelbrain_transducer_map.m   (preset resolution)
            ├─ write_babelbrain_input.m       (NIfTI label mask)
            ├─ export_babelbrain_weights.m    (calibration HDF5)
            └─ system() → run_babelbrain_case.py   ← Python bridge
                 └─ BabelBrain RUN_SIM().RunCases()
```

`babelbrain_wrapper.m` never imports Python directly. It constructs a command string and calls `system()`, which invokes `run_babelbrain_case.py` as a subprocess. All BabelBrain Python dependencies are therefore isolated in the conda environment pointed to by `babelbrain.python_env`; the MATLAB process only needs to call `system()`.

`run_babelbrain_case.py` adds the BabelBrain source tree (`external/BabelBrain`) to `sys.path` at runtime, imports `TranscranialModeling.BabelIntegrationBASE` dynamically, and calls `RUN_SIM().RunCases()`. The output `*DataForSim.h5` is written to `<dir_cache>/<subject_id>/`; PRESTUS reads it back with `read_babelbrain_output.m`.

### GPU backend selection and platform support

`computing_backend` maps directly to BabelBrain's `ComputingBackend` integer:

| Value | Backend | Supported hardware | Typical environment |
|---|---|---|---|
| `1` | CUDA | NVIDIA GPUs | Linux HPC / workstation with CUDA toolkit |
| `2` | OpenCL | AMD GPUs, Intel GPUs, some NVIDIA | Linux, any OpenCL-capable GPU |
| `3` | Metal | Apple Silicon (M1/M2/M3/M4) | macOS on ARM |
| `4` | CPU | Any | Fallback; much slower |

**Mac / Apple Silicon.** Set `computing_backend: 3` (Metal). This is the recommended path on M-series Macs; `device_name` can be left at `'auto'`. k-Wave's CUDA binary is unavailable on Apple Silicon, so BabelBrain + Metal is the only GPU-accelerated acoustic solver available on that platform.

**Linux / NVIDIA.** Set `computing_backend: 1`. With `device_name: auto` the Python bridge auto-detects the first CUDA GPU via `nvidia-smi`, so you normally do not need to set `device_name` explicitly on SLURM clusters. Ensure the SLURM job requests a GPU (`--gres=gpu:1`).

**Linux / AMD or other OpenCL GPU.** Set `computing_backend: 2`. The bridge uses `pyopencl` to auto-detect the first available GPU device. Install `pyopencl` in the BabelBrain conda environment if not already present.

**macOS Intel.** Metal is available but performance is limited on Intel integrated graphics. Use `computing_backend: 3` with `environment_mac_x64-314.yml`; for CPU-only runs use `computing_backend: 4`.

**Windows.** PRESTUS has not been tested on Windows. BabelBrain ships an `environment_win-314.yml`; the MATLAB `system()` call uses Unix-style path separators which would need adjustment.

### Device name auto-detection

`run_babelbrain_case.py::resolve_device_name()` handles `device_name: auto` by detecting the appropriate device for the selected backend:

- **CUDA (backend 1):** calls `nvidia-smi --query-gpu=name` and uses the first result.
- **OpenCL (backend 2):** iterates `pyopencl.get_platforms()` and picks the first GPU device.
- **Metal (backend 3):** resolves to `'Apple'`; no detection needed.

Set `device_name` explicitly to target a specific GPU when multiple are present.

### Transducer placement, ZSteering, and the drive phases

BabelBrain builds its simulation domain with the transducer apex at `z_mm ≈ 0` and the **FocalSpot (target) voxel at `z_mm = ` preset `FocalLength`** (read from `Babel_{tx_system}/default.yaml`) — i.e. it assumes the transducer sits at its *natural* focal distance from the target. PRESTUS's geometry instead places the apex at `trans_pos` (bowl rear centre), i.e. `|focus_pos − trans_pos|` from the target. That distance equals `focal_distance_bowl` only when the target sits at the nominal focal depth; planned/pFOV placements generally differ (e.g. 71 mm vs 92 mm for R15473 at 85 mm EP with a 64 mm-deep target). On the way back, [`resample_babelbrain_vol.m`](../functions/babelbrain/resample_babelbrain_vol.m) anchors the `z_mm = FocalLength` plane onto PRESTUS's `focus_pos` (`prestus_dim3 = (z_mm − focal_length_mm)/res + focus_pos`), so any difference between the two distances **rigidly translates the pasted field by `|focus_pos − trans_pos| − FocalLength`** along the beam axis (−7.7 mm for R15473 with the target at 85 mm EP; −28.75 mm with the target at 64 mm EP). Electronic steering cannot compensate for that translation: a translated field of a differently-steered transducer is not the true field — the near-field lobe structure shifts with the translation while the drive physics differ.

The wrapper therefore corrects the geometry **mechanically**:

```
TxMechanicalAdjustmentZ = preset_FocalLength − |focus_pos − trans_pos| · res_mm   (+ any user standoff)
```

It is computed by [`babelbrain_placement_correction.m`](../functions/babelbrain/babelbrain_placement_correction.m) (shared with `export_babelbrain_inputs`). Only when `trans_pos` is unavailable does it fall back to `preset_FocalLength − focal_distance_bowl`, which assumes the target sits at the nominal focal depth. Using `focal_distance_bowl` for a target at a different depth places BabelBrain's transducer too far from (or too close to) the head, shifting the whole field — skull incidence, near-field lobes and focus — relative to the anatomy.

This shifts the Tx *inside BabelBrain's domain* (`BabelIntegrationANNULAR_ARRAY`, `CalculateRayleighFieldsForward`) so the apex sits at its true distance from the target; the focus-anchored resample is then also apex-aligned (distance-preserving). BabelBrain's internal Rayleigh phase reference aims at `center_z = ZDim[FocalSpot] + ZSteering + TxMechanicalAdjustmentZ`, so the mechanical adjustment shifts Tx and aim together — with `ZSteering = 0` the aim is always the **natural focus relative to the apex**.

On top of the corrected geometry, the drive is set in one of two modes:

- **Explicit element phases** (`tr.annular.elem_phase_rad` present — calibrated/analytical configs): `ZSteering = 0` and the weights file carries the **full drive** as a delta from the natural-focus reference (where BabelBrain's internal ring phases are flat, all rings being equidistant from the centre of curvature). BabelBrain then reproduces *exactly* the per-ring phases that k-Wave/jwave/the Rayleigh model use — solver-benchmark parity by construction. `focal_calibration` runs **locate-only** in this mode: it reports the predicted focus of the fixed drive but never replaces it.
- **No explicit phases** (geometric fallback): BabelBrain steers electronically, with the aim measured physically from the apex:

  ```
  ZSteering = aim_bowl − preset_FocalLength
  ```

  where `aim_bowl` is the calibrated steer `D*` when `focal_calibration: true` (default; [`calibrate_focal_steering.m`](../functions/calibration/calibrate_focal_steering.m) root-finds the geometric steer whose achieved `flhm_center`/`peak` lands on `focal_distance_ep`, correcting the diffraction focal shift), or `focal_distance_bowl` (uncalibrated).

Pin `babelbrain.z_steering_m` to override the auto-derivation entirely (also bypasses calibration; the weight-delta reference follows the pinned aim `FocalLength + ZSteering`).

**Free-water cross-check.** With explicit phases, the BabelBrain axial profile must overlay the analytic Rayleigh profile of the same phases in *both* the near-field lobes and the focal lobe. A rigid axial offset of ≈ `FocalLength − |focus_pos − trans_pos|` (near-field lobes shifted by the same amount as the peak) indicates the mechanical placement correction is missing; a *shape* mismatch of the near-field lobes indicates the effective drive phases differ (e.g. weights delta referenced to the wrong aim).

**Bowl rim vs. the source plane.** BabelBrain injects the Rayleigh field on the first domain slice on the transducer side and requires every source element to lie behind it; if the bowl rim protrudes into the domain it snaps the transducer back in whole voxels (`Warning: Need to apply correction to reposition Tx for …` in the BabelBrain log, and a non-zero `source_residual_mm` in the `_prestus_anchor.json` sidecar). This happens whenever any part of the bowl lies inside the grid: the rim sits one bowl depth (5.4 mm for R15473) on the target side of the apex, so the bowl spans `d_apex … d_apex + bowl_depth` from the transducer-side face (`d_apex` = apex distance into the grid, negative if the apex is outside). Examples: a pFOV crop that starts at `trans_pos` (source pushed 5.5 mm too far back), or the bone-slab phantom with the apex 7.5 mm inside the grid (13.5 mm, which shifted the BabelBrain focus ~14 mm proximally). [`write_babelbrain_input.m`](../functions/babelbrain/write_babelbrain_input.m) therefore trims `ceil((d_apex + bowl_depth)/dx) + 1` transducer-side slices from the BabelBrain mask (it is water inside the bowl concavity; a warning is raised if it contains tissue). The trimmed slab is zero-filled when the field is pasted back.

### Output alignment: tissue-marker anchoring and the ZDomainStart fix

**The problem.** BabelBrain's own `MaterialMap` — its internal record of the tissue mask it actually simulated on — is byte-identical to the mask PRESTUS wrote (`write_babelbrain_input.m`); BabelBrain never alters tissue. So if the returned field and material map are pasted back onto the PRESTUS grid at the wrong position, *every* downstream consumer that assumes k-Wave/jwave's "one shared grid, no internal resampling" invariant — `acoustic_analysis.m`'s per-tissue metrics (`Isppa_brain`/`Isppa_skull`), `thermal_analysis.m`'s per-tissue thermal/CEM43 metrics, the axial-profile and intensity-overlay plots, and any exported NIfTI compared against `medium_masks` in an external viewer — is silently wrong at tissue boundaries, most visibly at the skull-brain interface where BabelBrain's elastic-skull solver produces a real, physically-expected reflection hotspot.

Two independent effects combined to create a multi-millimetre misalignment (observed: ~5.5 mm on a benchmark pFOV run):

1. **`ZDomainStart`.** `CalculateRayleighFieldsForward` places the transducer at `Tx_z += TxMechanicalAdjustmentZ + ZDomainStart`, where `ZDomainStart = CalculateDomainZReference()` is a skin-surface heuristic (`GlobalStartSkin − LineOfSightStartSkin`) meant for BabelBrain's own GUI workflow with oblique trajectories through a full-head domain. PRESTUS already computes the correct `TxMechanicalAdjustmentZ` itself (see above) and hands BabelBrain a beam-aligned mask, so this extra term silently moves the source off the requested position.
2. The (now-removed) resample anchor compensated for source mis-placement by shifting the paste-back to match the *source*, which correctly re-aligned the field to the transducer but dragged the (unchanged) tissue away from the original `medium_masks`.

**The fix.**

- `run_babelbrain_case.py` monkey-patches `SimulationConditionsBASE.CalculateDomainZReference` to return `0.0` (`_patch_zdomain_reference`, same policy as the existing skull-attenuation patch — no edits under `external/`), so the source is placed at exactly the requested `TxMechanicalAdjustmentZ`. Disable with the Python bridge's `--no_zdomain_fix` flag (not currently exposed as a `babelbrain.*` config option — same as `--no_skull_atten_fix`) to reproduce the un-patched behaviour.
- `read_babelbrain_output.m` / `read_babelbrain_thermal.m` anchor the resample on **`/TargetLocation`** — the label-5 focal-marker voxel PRESTUS itself wrote at `focus_pos` — rather than on the source position. This is the exact inverse of the forward mask transform, so the tissue map lands back on `medium_masks` by construction (verified on a benchmark run: overall label agreement 83% → 99%, bone 20% → 98%). With the source correctly placed (previous point), one anchor now aligns tissue *and* field simultaneously — the marker-anchored paste is source-consistent too.
- The passed-in `focal_length_mm` (preset `FocalLength`) is used only as a **fallback** when `/TargetLocation` is missing from the HDF5 (older files).
- A per-run sidecar (`<target>_prestus_anchor.json`) now records the **realized-minus-requested source position** as a diagnostic (not an anchor correction): with the `ZDomainStart` fix active, the only remaining residual should be BabelBrain's whole-voxel source snap-back (≤ 1 native voxel). `read_babelbrain_output.m` warns when the residual exceeds 1.5 native voxels, which points at a missing/disabled `ZDomainStart` patch rather than a coordinate bug.

This is fully general — no `medium_masks`/`segmentation` special-casing is needed anywhere downstream. `thermal_analysis.m` (and any other consumer that trusts `medium_masks`) automatically inherits correctly-aligned data because the BabelBrain acoustic step now returns output in the same PRESTUS-space convention k-Wave/jwave always provided.

### BabelBrain output format

The bridge writes `*DataForSim.h5` under `<dir_cache>/<subject_id>/`. `read_babelbrain_output.m` reads:

- `/p_amp` — peak positive pressure [Pa], shape `[Nx, Ny, Nz]`
- `/p_amp_refocus` — skull-corrected pressure (only when `do_refocusing: true`)

`p_max_all` in the returned `sensor_data` struct is set to `/p_amp` (or `/p_amp_refocus` when refocusing is active), matching the convention used by `acoustic_wrapper`.

### FOV/domain-coverage and focal-anchor debugging

`read_babelbrain_output.m` permutes and axis-flip-corrects the raw HDF5 field so that BabelBrain's
`[Nx, Ny, Nz]` axes line up with PRESTUS's `[lateral X, lateral Y, axial]` convention (dim 3 = beam
axis, matching `acoustic_grid_fov.m`'s `beam_dim = 3`). This axis handling is verified against the
shared tissue-layer ordering in benchmark axial profiles (e.g. `plot_axial_profile.m` output), which
lines up identically between the k-Wave and BabelBrain backends.

A BabelBrain axial profile with no reflections and exactly-zero pressure over a stretch nearest the
transducer/skull has two distinct possible causes that look identical in the final profile. **Inspect
the raw `*DataForSim.h5` directly before assuming either one** (`h5py`/`h5disp` on `/MaterialMap`,
`/p_amp`, `/TargetLocation`, `/z_vec`) — the BabelBrain log lines `ZShrink_L` / `ZShrink_R` (printed
by `BabelIntegrationBASE.py` during domain sizing) tell you immediately whether narrow-beam cropping
actually touched the proximal side for a given run; it is easy to assume it did when it didn't.

1. **Domain-coverage limit.** BabelBrain's own domain sizing (`--tight_narrow_beam` cropping when no
   PRESTUS-side FOV crop is active, or its baseline PML/padding-offset sizing when the grid is already
   FOV-cropped) determines how far its returned volume actually extends. `resample_babelbrain_vol.m`
   maps that volume onto the PRESTUS grid with **zero-fill for any voxel outside the BabelBrain domain**
   (`affine_resample_3d(..., method, 0)`), so a BabelBrain domain that genuinely stops short of the
   skull reads back as "no reflections" when the data there is simply absent rather than physically
   zero. Check the `ZShrink_L`/`ZShrink_R`/`XShrink_*`/`YShrink_*` log lines — nonzero values there
   only occur when narrow-beam mode is actually active (no PRESTUS-side FOV crop). When it is active,
   its ray-traced footprint estimate has been found to under-shoot the true beam width laterally near
   the skull (the aperture is wide near the transducer, narrowing only toward the focus) — cross-check
   against `BabelViscoInput.nii.gz`'s own voxel counts (does it actually contain skin/cortical/
   trabecular voxels?) vs. the exported `DataForSim.h5` `MaterialMap`'s voxel counts to see whether
   narrow-beam cropped away tissue that was genuinely present in the input. Also check the "BabelBrain
   native resolution" log line for a mismatch against `grid.resolution_mm` — BabelBrain reuses the
   input mask's voxel array shape directly as its own grid dimensions, so a resolution mismatch there
   silently rescales the whole domain independent of any shrink settings.

2. **Focal-anchor error.** The z-axis resample anchor is BabelBrain's own `/TargetLocation` — the
   index where it placed the label-5 focal marker `write_babelbrain_input.m` wrote into the input
   mask, i.e. exactly PRESTUS's `focus_pos` voxel (see
   [Output alignment](#output-alignment-tissue-marker-anchoring-and-the-zdomainstart-fix)). The
   passed-in `focal_length_mm` (preset `FocalLength` / `TxMechanicalAdjustmentZ` / `ZSteering` chain
   — see [Transducer placement](#transducer-placement-zsteering-and-the-drive-phases)) is used only
   as a **fallback** when `/TargetLocation` is absent from the file (older BabelBrain HDF5s) — it is
   never used to override or cross-check the marker when present. `/TargetLocation` is computed by
   BabelBrain *before* it flips `p_amp`/`MaterialMap` along the depth axis, so reading it correctly
   means indexing the **raw, unflipped** `/z_vec`, not the flip-corrected coordinate array
   (`zdim`/`z_mm`) used for the data arrays — get this backwards and the anchor lands on the wrong
   voxel. When the fallback *is* in use and disagrees with what the marker would have given (only
   possible on files without `/TargetLocation`), the entire field is rigidly shifted along the beam
   axis when resampled — real, non-zero BabelBrain data (including skull reflections) landing at the
   wrong depth, some of it outside the PRESTUS grid (reads as the same "zero near the skull" symptom
   as case 1).

When `simulation.debug: 1` and a BabelBrain→PRESTUS resample occurs, `read_babelbrain_output.m` writes
diagnostics to `<dir_output>/debug/babelbrain/`:

- `sub-<id>_babelbrain_fov_debug<output_affix>.png` — a stacked figure:
  1. BabelBrain-native central-axis **pressure** profile (its own grid, before resampling), with the
     `/TargetLocation` tissue-marker anchor at `z=0` and, only when the passed fallback anchor
     differs from the marker actually used, a second line showing where that (unused) fallback would
     have placed it.
  2. The same axis, but for BabelBrain's own `/MaterialMap` (tissue index: 0=water, 1=skin,
     2=cortical, 3=trabecular, 4=brain) — shows directly where BabelBrain itself thinks the skull is,
     independent of any PRESTUS-side coordinate assumption.
  3. The resampled PRESTUS-grid pressure profile, with the region outside the BabelBrain domain
     shaded (zero-filled). The depth axis is beam-sign-corrected, so the shaded (truncated) side
     always corresponds to the actual proximal/distal side of the domain regardless of which way the
     beam runs along grid dimension 3.
- `sub-<id>_babelbrain_fov_debug<output_affix>.csv` — coverage + anchor numbers as a table:
  `bb_depth_rel_focus_min_mm` / `_max_mm`, `prestus_depth_rel_focus_min_mm` / `_max_mm`,
  `zero_filled_proximal_mm` / `_distal_mm`, and `anchor_correction_applied_mm` (how far the geometric
  estimate differs from `/TargetLocation`; normally ~0).
- `sub-<id>_babelbrain_material_map<output_affix>.nii.gz` — BabelBrain's `/MaterialMap`, resampled
  onto the PRESTUS grid with **nearest-neighbor** interpolation (never linear — this is integer label
  data) so it can be overlaid directly against `medium_masks` / the T1 in an external viewer.

Run-time console lines summarize both checks even without inspecting the files:

```
[read_babelbrain_output] FOV coverage check: BabelBrain domain covers depth-from-focus [...] mm; PRESTUS grid requires [...] mm.
[read_babelbrain_output] -> zero-filled (no BabelBrain data): X mm proximal (transducer side) / Y mm distal (far side) of the PRESTUS grid.
[read_babelbrain_output] Anchor check: focal_length_mm=... mm (estimate) vs /TargetLocation=... mm (BabelBrain ground truth, mismatch ... mm) — using /TargetLocation for the resample anchor.
```

### BabelBrain pipeline stages and what to check at each one

```
write_babelbrain_input.m        PRESTUS medium_masks → BabelBrain label NIfTI (0-5)
  Cortical/trabecular skull labels are derived from a single combined
  (skull+cortical+trabecular) smoothed mask in preproc_medium_mask.m, with the
  trabecular assignment constrained to that combined shell — trabecular cannot
  extend past the cortical envelope (e.g. touch skin directly). A central-line
  tissue-transition dump of the input NIfTI should never show skin→trabecular
  with no cortical between them; if it does, check the raw segmentation itself.
  The pre-flip (line ~132) must exactly cancel BabelBrain's own np.flip(axis=2)
  on load — verify both sides together if either changes.
  The NIfTI header (voxel -> RAS affine) is built by grid_to_ras_affine.m,
  composing planimg.inv_transf (head-alignment rotation), planimg.t1_header
  (T1's own affine), any acoustic-FOV-crop offset (planimg.fov_offset_ac,
  set in prestus_pipeline.m right after acoustic_fov_crop), and the dim-3
  flip above — in that order. This is what makes the written NIfTI overlay
  correctly onto the T1/medium_masks; a purely diagonal/origin-only affine
  (no rotation) will misalign whenever grid.mode rotates the head, which is
  the default.

run_babelbrain_case.py          Invokes RUN_SIM().RunCases(); patches material DB
  _inject_prestus_materials only runs when --medium_properties_json is passed AND
  --no_skull_atten_fix is NOT set — without the skull-attenuation fix, injected bone
  attenuation is never realised (QCorr stays 1.0), reproducing near-zero skull
  attenuation even though the Material table looks correctly populated.
  BabelBrain computes its own SpatialStep from PPW and reuses the input mask's voxel
  array shape directly as its own grid dimensions, without resampling to that
  SpatialStep — write_babelbrain_input.m resamples the mask to BabelBrain's own
  native resolution before writing so voxel spacing actually matches (see
  bb_native_res_mm there); check the log's "BabelBrain native resolution" line if the
  effective grid resolution is in question.
  _patch_zdomain_reference (default on; --no_zdomain_fix disables) forces
  CalculateDomainZReference to 0 so the source lands at exactly the requested
  TxMechanicalAdjustmentZ — see "Output alignment" above. Check the run log's
  "ZDomainStart neutralized: ... mm" line to see the magnitude BabelBrain would
  otherwise have added.

BabelIntegrationBASE.py         Domain sizing, Rayleigh forward-prop, FDTD, flip-before-save
  --tight_narrow_beam domain shrink (XShrink/YShrink/ZShrink_L/R): only active when
  no PRESTUS-side FOV crop is configured (babelbrain_wrapper.m disables it otherwise,
  since the mask is already cropped to a beam-centred box matching k-Wave/jwave).
  When active, check the run log's own `ZShrink_L <n>` / `ZShrink_R <n>` lines;
  ZShrink_L>0 means the proximal (transducer/skull) side was actually cropped. It is
  often 0 even with narrow-beam active, since only the far/beyond-focus side is
  targeted by z_beyond_focus_mm — don't assume cropping without checking the log.
  The lateral (XShrink/YShrink) crop has been found to under-estimate the beam's
  true width near the skull, cropping away skull tissue that is genuinely present
  in the input mask — this was the reason narrow-beam mode is disabled by default
  whenever a PRESTUS-side FOV crop already exists.
  /TargetLocation is BabelBrain's own record of where it placed the focal-spot
  marker PRESTUS wrote into the input mask — this IS the resample anchor (see
  "Output alignment" above), not merely a cross-check.
  p_amp/MaterialMap are np.flip(axis=2)'d before being added to DataForSim; z_vec is
  added AFTER that flip loop and is therefore NOT flipped — read_babelbrain_output.m's
  zdim = z_vec_bb(end:-1:1) compensates for exactly this asymmetry. Do not change one
  side without the other.

read_babelbrain_output.m / read_babelbrain_thermal.m   h5 read, flip-correction, resample
  Pressure/complex fields resample with 'linear' (resample_babelbrain_vol.m default);
  MaterialMap (an INTEGER label field) is resampled with method='mode' (majority-vote
  over the finer BabelBrain-native voxels landing in each PRESTUS cell) — 'nearest'
  would drop thin structures (e.g. sub-voxel cortical bone) that miss the single
  sample point.

babelbrain_wrapper.m             Free-water ISPPA baseline (separate --water_only run)
  The water-only baseline run is only regenerated if its Water_DataForSim.h5 doesn't
  already exist in the cache dir (is_water_medium branch) — a stale Water_DataForSim.h5
  left over from an earlier run (different transducer config, different geometry) will
  be silently reused for ISPPA calibration (unless io.overwrite_files is 'always'/'force',
  which deletes it first). If the ISPPA-normalized field looks wrong, check the
  Water_DataForSim.h5's mtime against the main DataForSim.h5's.
  The water-only run must receive the same --source_amp_pa as the transcranial run;
  without it BabelBrain drives at its ~50 kPa default and freefield_isppa_wcm2 comes out
  ~(156.6/50)^2 ≈ 10x too low for R15473 (the wrapper passes it since this fix — water
  baselines cached from before it are wrong and must be regenerated).
```

---

## Version tracking

The BabelBrain version is recorded automatically at runtime and stored with the simulation results:

```matlab
sensor_data.babelbrain_version.version_txt   % e.g. '0.8.1'
sensor_data.babelbrain_version.git_describe  % e.g. '0.8.0-19-g5237cd2'
sensor_data.babelbrain_version.git_commit    % full SHA
sensor_data.babelbrain_version.bb_root       % path used at runtime
```

This information is written into the cached `.mat` file via `acoustic_provenance.babelbrain_version` and therefore survives pipeline restarts.


---

<span style="font-size: 12px;">
  <strong>BabelBrain citation:</strong><br>
  Pichardo, S. BabelBrain: An Open-Source Application for Prospective Modeling of Transcranial Focused Ultrasound for Neuromodulation Applications. <em>IEEE Trans Ultrason Ferroelectr Freq Control</em> <strong>70</strong>, 587–599 (2023). https://doi.org/10.1109/TUFFC.2023.3274046
</span>
