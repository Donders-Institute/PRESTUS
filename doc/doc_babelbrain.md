# BabelBrain backend

!!! note "At a glance"
    - **What:** BabelBrain as an alternative acoustic and thermal solver.
    - **When you need it:** Optional: to compare with or run BabelBrain from within PRESTUS.
    - **Before you start:** A BabelBrain Python environment.

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
> On the default label-mask path, PRESTUS writes `medium_properties` to `<target>_medium_properties.json` and the Python bridge overwrites BabelBrain's built-in acoustic table (`MatFreq`) with those values: density, longitudinal sound speed, and attenuation converted from dB/(cm·MHz^y) to Np/m. Shear is set to 0. This covers skin, cortical bone, trabecular bone and brain (an undifferentiated `skull` fills both bone slots), so BabelBrain solves with the **same tissue acoustics as k-Wave/jwave** for those tissues. Remaining differences are the solver and the bone Q-correction (1.0 on this path; see *Realised bone absorption* below), so validate against a k-Wave run before relying on absolute values.
> This does **not** apply to (a) **water**: BabelBrain builds its water (coupling) material from the module-level `Material` left over from the `MatFreq` construction loop (the 1120 kHz entry), not from `MatFreq[frequency]`, so water stays at BabelBrain's 1000 kg/m³, 1500 m/s, no attenuation (PRESTUS default 994 kg/m³, 1500 m/s, 2.17e-3 dB/(cm·MHz²)); (b) the `ct_mapping` path, where skull properties come from BabelBrain's own HU models (skin and brain are still injected); (c) runs where writing the JSON failed, or where the drive frequency has no `MatFreq` entry (5 kHz steps), which print a warning and fall back to BabelBrain's built-in table (listed below); (d) voxel-wise pseudo-CT properties computed by PRESTUS; or (e) **thermal** properties. BabelBrain thermal uses its own table, which differs from PRESTUS in absorption fraction (see [Thermal simulation](#thermal-simulation)).

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
  q_correction: []             # bone Q-factor correction of the relaxation attenuation model; [] = 1.0 non-CT / 3.0 CT
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

By default PRESTUS passes a label mask, and BabelBrain's internal `MatFreq` acoustic table is overwritten with the PRESTUS `medium_properties` for skin, bone and brain (not water; see the note at the top of this page), so each tissue type gets one set of properties. Optionally, a pseudo-CT (HU) volume can be passed so that BabelBrain derives per-voxel skull properties from Hounsfield units using one of its fitted mapping models.

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

#### Internal air (`air_regions`)

On the pseudo-CT path, BabelBrain can treat internal air (frontal and sphenoid sinuses, mastoid cells) as a reflector. Air is not simulated as a material, because its low sound speed would force a tiny timestep. Instead, BabelBrain's solver never updates stress or particle velocity inside the air mask (`ReflectorMask` in BabelViscoFDTD), so the region acts as a pressure-release boundary with reflection close to −1. BabelBrain's own GUI enables this by default for CT/ZTE/PETRA input (*Generate air regions*). PRESTUS enables it with:

```yaml
babelbrain:
  ct_mapping: Webb-Marsac
  air_regions: true
  air_hu_window: [-1200, -400]   # default, as BabelBrain
  air_keep_boundary: false
```

`run_babelbrain_case.py` writes `<target>_<f>kHz_<ppw>PPW_AirRegions.nii.gz` next to the CT input, on the same grid, using BabelBrain's recipe:
1. Take the pseudo-CT voxels within `air_hu_window`.
2. Apply a 3-voxel median filter and find connected components (26-connectivity).
3. Remove the largest component (the air around the head).

PRESTUS additionally removes every component that touches the domain boundary (`air_keep_boundary: false`). PRESTUS hands BabelBrain an already-cropped grid, so the exterior air can be split into several pieces. That includes the gap between transducer and scalp, which is −1000 HU in the pseudo-CT but water in the simulation and would otherwise block the beam. It then sets `bExtractAirRegions`. The log reports the voxel count and volume of the internal air that was kept.

Limitations:

- **Pseudo-CT only.** It needs `ct_mapping`. On the label path, PRESTUS has no air label (SimNIBS background is water), so the option is ignored with a warning.
- **Free-water runs never get the mask.** The water baseline and `water_amp_calibration` are unaffected.
- **Air depends on the pseudo-CT.** In PRESTUS's UTE pseudo-CT, label-0 voxels are −1000 HU and skull voxels are mapped from UTE intensity, so air cells inside the skull mask are captured only where their mapped HU falls in the window. Check the written `AirRegions.nii.gz` against the anatomy.
- **Air touching the domain boundary is dropped.** Internal air cut by the domain boundary is removed together with the exterior air. Set `air_keep_boundary: true` if a sinus is cut by the domain edge and no exterior air touches the boundary.
- **No k-Wave equivalent.** k-Wave cannot freeze fields voxel-wise (pseudospectral), so k-Wave and BabelBrain runs differ wherever the beam meets internal air.

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
> **Comparison with PRESTUS defaults.** BabelBrain cortical bone longitudinal attenuation (81 Np/m at 500 kHz) is about 6% higher than the PRESTUS default (≈77 Np/m at 500 kHz, from 13.3 dB/(cm·MHz) with y = 1). BabelBrain cortical longitudinal SoS (2476 m/s at 500 kHz) is about 12% lower than the PRESTUS skull default (2800 m/s). Brain SoS: BabelBrain 1562 m/s vs. PRESTUS 1546 m/s (~1% difference). These differences would affect absolute pressure and heating values if the built-in table were used. See [Medium properties](doc_medium.md) for PRESTUS defaults.

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

With `babelbrain.water_amp_calibration: true` (default) the wrapper calibrates this out: it evaluates the in-repo Rayleigh–Sommerfeld model (`rayleigh_axial_intensity.m`) for the exact per-ring drive, forms `k = p_Rayleigh_peak / p_BabelBrain_water_peak` from BabelBrain's own free-water field (always available — BabelBrain computes a water reference in every case), and scales `sensor_data.p_max_all` (and `freefield_isppa_wcm2` by `k²`). The factor is logged and stored in `sensor_data.amp_calibration`.

The water-derived scalar is valid for tissue runs because (a) the bias originates in the source-injection/FDTD-correction chain, which is shared by water-only and transcranial runs (same grid, same source) — note however that the timestep, and hence the `Correction` factor, is **not** identical: with the fluid-bone model the transcranial timestep is limited by the cortical sound speed while the water-only run uses a larger one (sub-602 pFOV, R15473 300 kHz 6 PPW: CFL/water 0.300 → Correction 0.290 vs 0.471 → 0.465). The polynomial tracks the resulting ~1/CFL amplitude scaling to within ~2% in pressure (0.290/0.465 = 0.624 vs CFL ratio 0.637), so the water-derived scalar still transfers to ~4% in intensity; and (b) the simulation is linear, so the entire field — including through skull — scales with the source. This is the numerical analogue of hydrophone free-water calibration of a physical transducer. A factor deviating >15% from 1 triggers a warning, since that magnitude indicates a setup mismatch (geometry/phases/source amplitude) rather than solver bias.

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

After BabelBrain thermal runs, `read_babelbrain_thermal` reads the output, converts temperature rise and CEM43 fields to PRESTUS conventions, and resamples the volumes from BabelBrain's simulation grid onto the PRESTUS thermal grid (`parameters_th` from `thermal_stage`: the full grid when an acoustic FOV crop was used, possibly offset / at `grid.thermal_resolution_mm` with a separate thermal grid). The returned `results_heating` struct is identical in structure to the output of `thermal_simulation`, so:

- `thermal_analysis` computes focal temperature, CEM43 metrics, and safety checks unchanged
- `nifti_thermal` exports temperature rise, CEM43, and heat-rise NIfTIs to the standard `nii/` folder
- The HTML report includes BabelBrain thermal results alongside or instead of k-Wave results

The resampling (`resample_babelbrain_vol`) is a diagonal scale + shift anchored on the `/TargetLocation` focal marker, mapped onto the thermal grid's focus position (`babelbrain_thermal_wrapper` takes it, and the step, from the thermal grid: `grid.thermal_focus_pos_exact` / `thermal_trans_pos_exact` / `thermal_resolution_mm` set by `thermal_grid_setup`, else `transducer(1).focus_pos`; `sensor_data.focus_pos` is on the acoustic grid and only a fallback). Two orientation details differ from the acoustic readback: the thermal arrays are stored in `z_vec` order (BabelBrain's `CalculateTemperatureEffects` flips `p_amp`/`MaterialMap` on load and saves without flipping back), so `z_vec` is **not** reversed, unlike for `DataForSim.h5`; and temperatures are resampled as rises over `BaselineTemperature`, so voxels outside BabelBrain's thermal domain read the baseline rather than 0 °C.

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
| **Single transducer** | ✓ Full | Standard use case, with `grid.mode: focal_axis` (default) or `transducer_axis`, which put the beam along the grid z axis. See [Geometry limits](#geometry-limits). |
| **Multi-transducer — async/incoherent** (`simulation.transducer_coupling: async`) | ✗ Not supported | Errors at config load and in `async_transducer_pipeline`. See [Geometry limits](#geometry-limits). Use `acoustic_backend: kwave` or `jwave` with `thermal_backend: kwave`. |
| **Multi-transducer — coherent** (`simulation.transducer_coupling: coherent`) | ✗ Not supported | Errors at config load and in `babelbrain_wrapper`. See [Geometry limits](#geometry-limits). |
| **`grid.mode: ras_plus`** | ✗ Not supported | The beam is generally off the grid z axis; `babelbrain_wrapper` errors. See [Geometry limits](#geometry-limits). |
| **Limited acoustic FOV** (`grid.acoustic_fov_*`) | ✓ Supported | The NIfTI passed to BabelBrain is the same beam-centred, FOV-cropped box used for k-Wave/jwave (`acoustic_fov_crop.m`/`acoustic_grid_fov.m`). BabelBrain's own narrow-beam domain cropping (`bTightNarrowBeamDomain`) is disabled in this case — it would crop the already-cropped box a second time, and its ray-traced lateral footprint estimate has been found to under-shoot the true beam width near the skull. BabelBrain simulates the full FOV-cropped box with no further internal shrink, so its domain matches the other backends' voxel-for-voxel. |
| **Axisymmetric 2D** (`simulation.axisymmetric: 1`) | ✗ Not supported | BabelBrain is always a 3D solver. Axisymmetric mode must use the k-Wave backend. |
| **`multi_isppa` mode** | ⚠ Untested | `multi_isppa_pipeline` derives each acoustic/thermal job from the full parameter struct, so `acoustic_backend`/`thermal_backend` carry over. Not yet validated end-to-end with BabelBrain. |
| **`phantom` / `water` test medium** | ⚠ Untested | BabelBrain expects a realistic segmented NIfTI. Purely water or phantom media are geometrically valid but BabelBrain's material property tables will apply skull properties to any non-background label, which may produce unexpected results. |
| **Post-hoc water simulation** | ⚠ Untested | The post-hoc water loop re-enters the pipeline with `simulation.medium: water` and inherits `acoustic_backend`, so with BabelBrain it runs a BabelBrain water-only simulation. |
| **Thermal tissue properties** | ⚠ Known limitation | BabelBrain thermal uses its built-in conductivity, specific heat, perfusion and absorption fraction; PRESTUS `medium_properties` thermal fields (and their uncertainty variants) are not applied. The absorption fraction differs most (skin/brain 0.85 vs 1.0, skull ≈0.15 vs 0.28), so BabelBrain heating is lower than k-Wave's for the same field, mostly in skull. See [Thermal tissue properties](#thermal-tissue-properties). |
| **BabelBrain thermal without BabelBrain acoustic** | ✗ Not supported | `thermal_backend: babelbrain` requires the `*DataForSim.h5` file produced by `acoustic_backend: babelbrain`. It cannot consume a k-Wave acoustic result. |
| **CT-based skull properties into BabelBrain** | ✓ Supported (via `ct_mapping`) | Set `babelbrain.ct_mapping` to a non-empty mapping method (e.g. `'Webb-Marsac'`) and pass the pseudoCT HU volume as the 8th argument to `babelbrain_wrapper`. `babelbrain_wrapper` writes a `CT.nii.gz` alongside `BabelViscoInput.nii.gz` and passes `--use_ct`, `--ct_mapping`, and (optionally) `--use_petra` flags to `run_babelbrain_case.py`, which forwards them to BabelBrain's `RUN_SIM().RunCases()`. BabelBrain then performs its own HU→property conversion internally. **Note:** use this path only when passing raw HU values from a UTE/PETRA scan. If PRESTUS's pCT pipeline (`pct.enabled: 1`) has already converted HU to acoustic properties, use the default label-based path (`ct_mapping` empty) — the per-voxel acoustic properties k-Wave uses cannot be passed to BabelBrain. See [CT-based skull properties](#ct-based-skull-properties-ct_mapping) for available mapping methods. |
| **Water (coupling medium) properties** | ⚠ Known limitation, by design | PRESTUS `medium_properties.water` is **not** applied: BabelBrain builds its water material from the module-level `Material` (the last, 1120 kHz, entry of the `MatFreq` construction loop), so water is always 1000 kg/m³, 1500 m/s, no attenuation (PRESTUS default: 994 kg/m³, 1500 m/s, 2.17e-3 dB/(cm·MHz²)). The `Water` line in the bridge's "Injected PRESTUS materials" log only updates `MatFreq[frequency]['Water']`, which BabelBrain reads solely in `GetSmallestSOS` (grid spacing), not as the simulated water. Water is deliberately not injected: the water density and sound speed also set BabelBrain's source pressure-to-velocity factor (`FactorConvPtoU`), so changing them would shift the source amplitude and require re-checking the free-water calibration, for a ~0.6% impedance difference and negligible attenuation. Studies that need identical water across backends should set `medium_properties.water` to 1000 kg/m³ / 1500 m/s. |
| **Uncertainty mode** | ⚠ Untested | The liberal/conservative configs vary `medium_properties`. These are forwarded to BabelBrain through the material-injection JSON (density, sound speed, attenuation), so variants should reach the acoustic solve on the label-mask path. They have no effect with `ct_mapping`, and never on BabelBrain **thermal** properties. Not yet validated end-to-end. |

### Geometry limits

BabelBrain simulates **one transducer per run, with its beam along the z axis of its domain**. In its own workflow it guarantees this by resampling the head into a grid aligned to the trajectory (`BabelDatasetPreps.GetSkullMaskFromSimbNIBSSTL`: the domain is rotated by the trajectory matrix). PRESTUS instead hands BabelBrain the PRESTUS simulation grid as it is: `write_babelbrain_input` resamples it to BabelBrain's resolution and orients the depth axis, but does not rotate it. BabelBrain then places the transducer on that grid's z axis, through the focus marker, at the apex distance |`focus_pos` − `trans_pos`|. A beam that is not along z would be simulated in the wrong direction while still aiming at the right focus.

PRESTUS therefore refuses these cases (`babelbrain_check_geometry`, error IDs `PRESTUS:babelbrain:multiTransducer` and `PRESTUS:babelbrain:obliqueBeam`):

- **More than one transducer**, with `acoustic_backend` or `thermal_backend` `babelbrain`, coherent or async. Any grid shared by several transducers puts at least one beam off z (`grid.mode: focal_axis` tilts each beam by about half the angle between them), and a per-transducer grid would need the fields resampled between grids, which is not accurate for the complex fields a coherent sum needs. BabelBrain's thermal step runs on its own single-transducer acoustic output, so it cannot take a combined heat source either. Checked at config load (`load_parameters`), at the start of `async_transducer_pipeline`, and in `babelbrain_wrapper`. The Nexus flags it before launch.
- **A beam more than 2° off the grid z axis**, checked in `babelbrain_wrapper` from `trans_pos` and `focus_pos`. With one transducer, `focal_axis` and `transducer_axis` put the beam on z; `ras_plus` generally does not.

Multi-transducer support for BabelBrain would need each transducer simulated in its own beam-aligned domain and the fields combined on one grid afterwards (as `simulation.async_grid_mode: per_transducer` does for k-Wave), plus a thermal step that accepts a combined heat source. Until then, use `acoustic_backend: kwave` or `jwave` with `thermal_backend: kwave` for multi-transducer runs.

---

## Developer notes

The call chain, GPU backend selection, placement and ZSteering, output alignment and stage-by-stage debugging are described in [BabelBrain integration internals](doc_babelbrain_internals.md).

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
