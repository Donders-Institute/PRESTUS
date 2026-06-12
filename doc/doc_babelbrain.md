# BabelBrain acoustic backend

PRESTUS can use [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) as an alternative backend for acoustic simulation, thermal simulation, or both. By default PRESTUS uses k-Wave for both; the BabelBrain backends can be enabled independently. All pre-processing stages (segmentation, grid setup, medium masks, transducer placement) and all downstream post-processing (focal metrics, NIfTI export, HTML report) run unchanged regardless of which backend is active.

**What is BabelBrain?** BabelBrain is an open-source tool for transcranial focused ultrasound planning developed by Samuel Pichardo (University of Calgary). It uses a different acoustic solver called BabelViscoFDTD, which models the full viscoelastic wave equation — including shear waves in bone — rather than the scalar linear acoustic model used by k-Wave. BabelBrain also includes a Rayleigh integral step for computing the field through the water coupling layer before the skull, which can improve accuracy at the skull boundary. It supports GPU acceleration via CUDA, OpenCL, and Metal (Apple Silicon).

The main reasons to use BabelBrain instead of (or alongside) k-Wave are:

- **Shear wave physics in bone.** Mode conversion at the skull is a real effect; BabelBrain models it explicitly. k-Wave's scalar model absorbs this into an effective attenuation coefficient.
- **Cross-validation.** Comparing results from two independent solvers is a common check for simulation-based studies.
- **Apple Silicon GPU acceleration.** On Apple Silicon, k-Wave's fast C++ binary (`kspaceFirstOrder-CUDA`) is unavailable because it requires NVIDIA CUDA; only the slower `matlab_gpu` path works. BabelBrain has a native Metal backend that makes full use of the M-series GPU.

| | PRESTUS (k-Wave) | BabelBrain |
|---|---|---|
| **Solver** | Pseudo-spectral FDTD | Explicit staggered-grid FDTD |
| **Physics** | Scalar linear acoustic (pressure/density) | Full viscoelastic (longitudinal + shear stress/strain) |
| **Skull model** | Scalar ρ, c, α | Viscoelastic: separate longitudinal and shear SoS and attenuation |
| **Rayleigh pre-step** | None | Rayleigh integral for the water coupling layer |
| **Narrow-beam domain** | Optional (see [limited FOV](doc_limited_fov.md)) | Default (`bTightNarrowBeamDomain`) |
| **GPU backends** | CUDA, MATLAB GPU arrays | CUDA, OpenCL, Metal (Apple Silicon) |

See [Simulation backends](doc_backend.md) for the k-Wave-only backend options.

> **Important — acoustic material properties are not shared with PRESTUS.**
> BabelBrain assigns tissue properties (speed of sound, density, attenuation, shear parameters) from its own internal frequency-fitted literature tables. It does not read `medium_properties` from the PRESTUS config. As a result, BabelBrain and k-Wave use different baseline property values (see [Tissue labels and material properties](#tissue-labels-and-material-properties) for a direct comparison), and **the liberal/conservative property variants used in uncertainty mode have no effect on BabelBrain acoustic simulations.** Results from the two backends should not be compared directly without accounting for these differences.

---

## Quick start

Both the acoustic and thermal simulations can be run in BabelBrain independently or together. The acoustic backend replaces k-Wave for pressure field computation; the thermal backend replaces k-Wave's BHTE solver and uses BabelBrain's own Pennes BHTE on the same simulation domain. All downstream PRESTUS post-processing (focal metrics, CEM43, NIfTI export, HTML report) runs unchanged in both cases.

**Acoustic only (default when using BabelBrain):**

```yaml
acoustic_backend: babelbrain

babelbrain:
  tx_system: R15473            # BabelBrain preset name; auto-detected if empty
  computing_backend: 3         # 1=CUDA  2=OpenCL  3=Metal  4=CPU
  device_name: Apple           # GPU device identifier string
  ppw: 6                       # points-per-wavelength floor; actual PPW = max(ppw, round(c_brain / (freq × res_mm × 1e-3))). See note below.
  do_refocusing: false         # skull back-propagation refocusing
  z_beyond_focus_mm: 40        # axial domain truncation past the focal point [mm]
  python_env: ''               # path to conda/venv environment; empty = system python3
  root: ''                     # path to BabelBrain source; auto-detected if empty
  tx_mech_adjust_z_m: 0.0      # transducer z offset into skin [m] (0 = geometry in NIfTI)
  z_steering_m: ~              # ZSteering [m]; null = auto-derived from focal distance
```

> **PPW and simulation resolution.** `ppw` sets a floor; PRESTUS resolves the actual PPW as `max(ppw, round(c_brain / (freq × res_mm × 1e-3)))` where `c_brain` is taken from `medium_properties.brain.sound_speed` (default 1500 m/s). BabelBrain then internally sets its spatial step to `SmallestSOS / (PPW × freq)`, where `SmallestSOS` is the minimum speed of sound across all tissue types in BabelBrain's own tables — including bone shear waves (~1473 m/s at 300 kHz), which can be slightly lower than water or brain longitudinal SoS. In practice the NIfTI voxel size written by PRESTUS (`grid.resolution_mm`) must match BabelBrain's expected `SpatialStep`; BabelBrain will print a warning if they differ. To ensure consistency, `grid.resolution_mm` should equal `c_min / (freq × PPW × 1e-3)` using the same `c_min` as BabelBrain's `SmallestSOS`. Because shear SoS is slightly below 1500 m/s, using the default `c_brain = 1500 m/s` for `ppw_auto` is conservative and will not underestimate the required PPW.

**Acoustic + thermal (both in BabelBrain):**

```yaml
acoustic_backend: babelbrain
thermal_backend:  babelbrain   # use BabelBrain Pennes BHTE instead of k-Wave

babelbrain:
  tx_system: R15473
  computing_backend: 3
  device_name: Apple
  isppa_w_per_cm2: 5.0         # ISPPA at focus [W/cm²]; default: from transducer calibration
```

> **`thermal_backend: babelbrain`** requires `acoustic_backend: babelbrain` to have run first, because BabelBrain's thermal solver operates on the `*DataForSim.h5` file produced by the acoustic stage.

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
```

**When to use this path.** Only when you want BabelBrain's own HU→acoustic-property conversion to apply to raw UTE or PETRA pseudo-CT Hounsfield values. This is the equivalent of BabelBrain's built-in `bUseCT` mode.

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

BabelBrain assigns frequency-dependent material properties internally from literature fits computed at the simulation frequency. They are tabulated every 5 kHz from 100 to 1125 kHz. You do not need to configure these; BabelBrain applies them automatically based on the tissue labels above.

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

Shear SoS is only non-zero for bone. This is what enables BabelBrain to model mode conversion at the skull that k-Wave cannot capture.

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

> **Comparison with PRESTUS defaults.** BabelBrain cortical bone longitudinal attenuation (81 Np/m at 500 kHz) is approximately 47% lower than the PRESTUS default (152 Np/m, derived from 13.3 dB/cm/MHz). BabelBrain cortical longitudinal SoS (2476 m/s at 500 kHz) is about 11% lower than the PRESTUS skull default (2800 m/s). Brain SoS: BabelBrain 1562 m/s vs. PRESTUS 1610 m/s (~3% difference). These differences affect absolute pressure and heating values; results from the two solvers should not be compared directly without accounting for property differences. See [Medium properties](doc_medium.md) for PRESTUS defaults.

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

When PRESTUS has calibrated element phases for the transducer, these are automatically included in the BabelBrain input package as an `OptimizedWeightsFile`. BabelBrain applies these weights on top of its own geometric (Rayleigh-steered) phases, so the file stores **delta corrections** relative to geometric phases rather than absolute phases:

```
w[n] = exp(i × (calibrated_phase[n] − geometric_phase[n]))
```

This is computed automatically from the PRESTUS calibration result. The Python bridge handles a format conversion at runtime (MATLAB writes split real/imag datasets; BabelBrain expects a native complex64 HDF5 dataset).

If no calibrated phases are available, BabelBrain falls back to geometric (Rayleigh) focusing, which is equivalent to the PRESTUS uncalibrated mode.

---

## Thermal simulation

### BabelBrain thermal backend

When `thermal_backend: babelbrain` is set, PRESTUS calls BabelBrain's `CalculateThermalProcess` on the `*DataForSim.h5` acoustic output. BabelBrain solves the Pennes bio-heat equation on the same simulation domain using its own GPU-accelerated relaxation solver.

The timing parameters from `parameters.timing` are converted automatically:

| PRESTUS parameter | BabelBrain parameter | Conversion | Notes |
|---|---|---|---|
| `timing.pd` / `timing.pri` | `DutyCycle` | `pd / pri` | |
| `timing.pri` | `PRF` | `1 / pri` [Hz] | |
| `timing.ptrd` | `DurationUS` | direct [s] | One pulse-train block; see note below |
| `timing.ptri - timing.ptrd` | `DurationOff` | direct [s] | Inter-train cooling interval; 0 if `ptri ≤ ptrd` |
| `thermal.temp_0.brain` | `BaselineTemperature` | direct [°C] | |
| `transducer.set_intensity_w_per_cm2` | `Isppa` | direct [W/cm²] | See ISPPA note below |

> **`DurationUS` and multi-PTRI protocols.** `DurationUS` maps to `timing.ptrd` — the duration of one pulse-train repetition block. For protocols with `n_ptri_reps > 1`, BabelBrain's `Repetitions` parameter controls the number of blocks. PRESTUS emits a warning when `n_ptri_reps > 1` so you can verify the intended protocol.

> **ISPPA normalisation.** BabelBrain always rescales the entire pressure field so that ISPPA at the focal voxel equals the requested `Isppa` value (`PressureAdjust = sqrt(Isppa × 2 × ρ × c)`). When calibrated phases are provided via `OptimizedWeightsFile`, this rescaling overrides the calibrated amplitude. Set `transducer.set_intensity_w_per_cm2` (or `babelbrain.isppa_w_per_cm2`) to the intended in-situ ISPPA at the focus.

> **`sel_p` and refocusing.** When `babelbrain.do_refocusing: true`, the thermal simulation automatically uses the `/p_amp_refocus` pressure field (the skull-corrected map) rather than `/p_amp`. No manual configuration is needed.

`babelbrain.isppa_w_per_cm2` can be set explicitly to override the transducer calibration value.

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
| `kwave` | `babelbrain` | PRESTUS acoustics, BabelBrain thermal (via `external_acoustic_nifti` handoff) |
| `babelbrain` | `babelbrain` | Full BabelBrain pipeline, PRESTUS post-processing |

---

## Limitations and unsupported PRESTUS modes

Not all PRESTUS pipeline features are currently compatible with the BabelBrain backend. The table below summarises known limitations.

| PRESTUS feature | BabelBrain support | Details |
|---|---|---|
| **Single transducer** | ✓ Full | Standard use case. |
| **Multi-transducer — async/incoherent** (`simulation.transducer_coupling: async`) | ✓ Supported | `async_transducer_pipeline` runs one independent acoustic simulation per transducer (each a single-transducer BabelBrain call), then combines intensity fields incoherently via `combine_async_intensity`. `acoustic_backend: babelbrain` is inherited by each per-transducer sub-run automatically. |
| **Multi-transducer — coherent** (`simulation.transducer_coupling: coherent`) | ✓ Supported | `babelbrain_wrapper` detects `numel(parameters.transducer) > 1` and recurses once per transducer with a unique `output_affix` (`_t01`, `_t02`, …). Each sub-run reads `/p_complex` (the complex Fourier pressure field stored by BabelBrain alongside `/p_amp`) from its `DataForSim.h5`. The complex fields are summed coherently; `p_max_all = abs(Σ p_complex_i)`. Unlike k-Wave's native coherent mode (all transducers in one simulation), this treats transducer interactions as independent and ignores mutual acoustic coupling between transducers. |
| **Limited acoustic FOV** (`grid.acoustic_fov_*`) | ✗ Not supported | The NIfTI passed to BabelBrain always covers the full grid. BabelBrain performs its own narrow-beam domain cropping internally (`bTightNarrowBeamDomain`). Setting `acoustic_fov_*` alongside `acoustic_backend: babelbrain` will cause a shape mismatch during the back-projection step and an error. |
| **Axisymmetric 2D** (`simulation.axisymmetric: 1`) | ✗ Not supported | BabelBrain is always a 3D solver. Axisymmetric mode must use the k-Wave backend. |
| **`multi_isppa` mode** | ✗ Not supported | `multi_isppa_pipeline` dispatches internally without propagating `acoustic_backend`, so all ISPPA sweep runs silently use k-Wave. |
| **`phantom` / `water` test medium** | ⚠ Untested | BabelBrain expects a realistic segmented NIfTI. Purely water or phantom media are geometrically valid but BabelBrain's material property tables will apply skull properties to any non-background label, which may produce unexpected results. |
| **Post-hoc water simulation** | ⚠ Untested | The free-field water reference loop at the end of the pipeline always calls the k-Wave path and bypasses `use_babelbrain`. The reference field will be from k-Wave regardless of the acoustic backend setting. |
| **BabelBrain thermal without BabelBrain acoustic** | ✗ Not supported | `thermal_backend: babelbrain` requires the `*DataForSim.h5` file produced by `acoustic_backend: babelbrain`. It cannot consume a k-Wave acoustic result. |
| **CT-based skull properties into BabelBrain** | ✓ Supported (via `ct_mapping`) | Set `babelbrain.ct_mapping` to a non-empty mapping method (e.g. `'Webb-Marsac'`) and pass the pseudoCT HU volume as the 8th argument to `babelbrain_wrapper`. `babelbrain_wrapper` writes a `CT.nii.gz` alongside `BabelViscoInput.nii.gz` and passes `--use_ct`, `--ct_mapping`, and (optionally) `--use_petra` flags to `run_babelbrain_case.py`, which forwards them to BabelBrain's `RUN_SIM().RunCases()`. BabelBrain then performs its own HU→property conversion internally. **Note:** use this path only when passing raw HU values from a UTE/PETRA scan. If PRESTUS's pCT pipeline (`pct.enabled: 1`) has already converted HU to acoustic properties, use the default label-based path (`ct_mapping` empty) — the per-voxel acoustic properties k-Wave uses cannot be passed to BabelBrain. See [CT-based skull properties](#ct-based-skull-properties-ct_mapping) for available mapping methods. |
| **Uncertainty mode** | ✗ Not supported | BabelBrain does not accept different acoustic medium property variations. The liberal and conservative configs vary `medium_properties` (skull/brain sound speed, density, attenuation coefficient). k-Wave reads these directly when building the simulation medium. BabelBrain does not: it assigns material properties internally from frequency-fitted literature tables (see [Tissue labels and material properties](#tissue-labels-and-material-properties)), and `babelbrain_wrapper` never forwards `medium_properties` to the Python bridge. The only use of `medium_properties.brain.sound_speed` inside `babelbrain_wrapper` is to compute the minimum required PPW — it does not affect the BabelBrain simulation itself. |

---

## Developer notes: PRESTUS integration architecture

> This section describes how the BabelBrain integration is wired inside PRESTUS. It is intended for developers working on the backend or extending GPU platform support.

### Call chain

```
single_subject_pipeline.m
  └─ acoustic_wrapper.m          ← dispatches on parameters.acoustic_backend
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

**Mac / Apple Silicon.** The default config (`computing_backend: 3`, `device_name: Apple`) targets Metal. This is the recommended path on M-series Macs. k-Wave's CUDA binary is unavailable on Apple Silicon, so BabelBrain + Metal is the only GPU-accelerated acoustic solver available on that platform.

**Linux / NVIDIA.** Set `computing_backend: 1`. The Python bridge auto-detects the first CUDA GPU via `nvidia-smi` when `device_name` is left at the default sentinel `'Apple'`, so you normally do not need to set `device_name` explicitly on SLURM clusters. Ensure the SLURM job requests a GPU (`--gres=gpu:1`).

**Linux / AMD or other OpenCL GPU.** Set `computing_backend: 2`. The bridge uses `pyopencl` to auto-detect the first available GPU device. Install `pyopencl` in the BabelBrain conda environment if not already present.

**macOS Intel.** Metal is available but performance is limited on Intel integrated graphics. Use `computing_backend: 3` with `environment_mac_x64-314.yml`; for CPU-only runs use `computing_backend: 4`.

**Windows.** PRESTUS has not been tested on Windows. BabelBrain ships an `environment_win-314.yml`; the MATLAB `system()` call uses Unix-style path separators which would need adjustment.

### Device name auto-detection

`run_babelbrain_case.py::resolve_device_name()` handles the common case where `device_name` is left at the default `'Apple'` but a non-Metal backend is selected:

- **CUDA (backend 1):** calls `nvidia-smi --query-gpu=name` and uses the first result.
- **OpenCL (backend 2):** iterates `pyopencl.get_platforms()` and picks the first GPU device.
- **Metal (backend 3):** `'Apple'` is already the correct device name; no detection needed.

Set `device_name` explicitly to target a specific GPU when multiple are present.

### ZSteering derivation

BabelBrain's Rayleigh pre-step always focuses at `FocalLength + ZSteering`. PRESTUS derives `ZSteering` automatically:

- **With calibrated phases** (`tr.annular.elem_phase_rad` populated): the phases already encode the steering offset, so `ZSteering = 0`.
- **Without calibrated phases**: `ZSteering = focal_distance_mm/1000 − preset FocalLength`, read from `Babel_{tx_system}/default.yaml` in the BabelBrain source tree.

Override with `babelbrain.z_steering_m` if the auto-derivation is wrong for a custom preset.

### BabelBrain output format

The bridge writes `*DataForSim.h5` under `<dir_cache>/<subject_id>/`. `read_babelbrain_output.m` reads:

- `/p_amp` — peak positive pressure [Pa], shape `[Nx, Ny, Nz]`
- `/p_amp_refocus` — skull-corrected pressure (only when `do_refocusing: true`)

`p_max_all` in the returned `sensor_data` struct is set to `/p_amp` (or `/p_amp_refocus` when refocusing is active), matching the convention used by `acoustic_wrapper`.

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
