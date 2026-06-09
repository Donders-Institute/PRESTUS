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
  ppw: 6                       # points-per-wavelength (never less than the auto-computed minimum)
  do_refocusing: false         # skull back-propagation refocusing
  z_beyond_focus_mm: 40        # axial domain truncation past the focal point [mm]
  python_env: ''               # path to conda/venv environment; empty = system python3
  root: ''                     # path to BabelBrain source; auto-detected if empty
  tx_mech_adjust_z_m: 0.0      # transducer z offset into skin [m] (0 = geometry in NIfTI)
  z_steering_m: ~              # ZSteering [m]; null = auto-derived from focal distance
```

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

BabelBrain requires Python 3.9+ with BabelViscoFDTD and GPU compute libraries. The recommended setup uses the conda environment file bundled with BabelBrain:

```bash
# Apple Silicon
conda env create -f external/BabelBrain/environment_mac_arm64-311.yml
conda activate BabelBrain

# Linux (CUDA)
conda env create -f external/BabelBrain/environment_linux.yml
conda activate BabelBrain
```

Set `babelbrain.python_env` to the full path of the environment, or leave empty to use the active `python3`.

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
| **Multi-transducer** (`parameters.transducer(2:N)`) | ✗ Not supported | `babelbrain_wrapper` reads `transducer(1)` only. The multi-transducer workflow requires per-transducer acoustic runs followed by field combination (`assemble_limited_fov_fields`); BabelBrain thermal cannot accept a combined multi-source field. |
| **Limited acoustic FOV** (`grid.acoustic_fov_*`) | ✗ Not supported | The NIfTI passed to BabelBrain always covers the full grid. BabelBrain performs its own narrow-beam domain cropping internally (`bTightNarrowBeamDomain`). Setting `acoustic_fov_*` alongside `acoustic_backend: babelbrain` will cause a shape mismatch during the back-projection step and an error. |
| **Axisymmetric 2D** (`simulation.axisymmetric: 1`) | ✗ Not supported | BabelBrain is always a 3D solver. Axisymmetric mode must use the k-Wave backend. |
| **`multi_isppa` mode** | ✗ Not supported | `multi_isppa_pipeline` dispatches internally without propagating `acoustic_backend`, so all ISPPA sweep runs silently use k-Wave. |
| **`phantom` / `water` test medium** | ⚠ Untested | BabelBrain expects a realistic segmented NIfTI. Purely water or phantom media are geometrically valid but BabelBrain's material property tables will apply skull properties to any non-background label, which may produce unexpected results. |
| **Post-hoc water simulation** | ⚠ Untested | The free-field water reference loop at the end of the pipeline always calls the k-Wave path and bypasses `use_babelbrain`. The reference field will be from k-Wave regardless of the acoustic backend setting. |
| **BabelBrain thermal without BabelBrain acoustic** | ✗ Not supported | `thermal_backend: babelbrain` requires the `*DataForSim.h5` file produced by `acoustic_backend: babelbrain`. It cannot consume a k-Wave acoustic result. |

> **Recommended workaround for multi-transducer + BabelBrain.** Run each transducer independently with `acoustic_backend: babelbrain`; combine the resulting pressure fields with `assemble_limited_fov_fields` (treating each BabelBrain output as a full-grid result); then run thermal with `thermal_backend: kwave` on the combined field.

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
