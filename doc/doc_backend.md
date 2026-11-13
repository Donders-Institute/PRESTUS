# Simulation Backends

[k-Wave](https://github.com/ucl-bug/k-wave) is the default backend and needs no extra installation. Another backend serves to cross-check a k-Wave result with an independent solver, or to model what the fluid skull leaves out: `kwave_elastic` adds shear waves in the skull, which lowered the intracranial pressure in an ex vivo comparison with a fluid skull model. The fluid-skull backends are designed to give comparable acoustics from the same tissue properties, whereas [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain) thermal uses its own absorption fractions and gives lower heating.

On this page:

- [Comparing Results Across Backends](#comparing-results-across-backends): what agrees between backends, and what differs by design or as a known limitation.

!!! note "At a glance"
    - **When you need it:** You need nothing but k-Wave to start. Read this page to choose another backend, and [k-Wave](doc_kwave.md) to choose between CPU, GPU, [MATLAB](https://www.mathworks.com/products/matlab.html) and C++.
    - **Before you start:** k-Wave (installed with PRESTUS).

PRESTUS runs the acoustic and the thermal simulation through interchangeable backends. Pre-processing (segmentation, medium mapping, grid, transducer placement) and post-processing (focal metrics, CEM43, [NIfTI](https://nifti.nimh.nih.gov/) export, HTML report) are shared, so switching backend is a config change:

```yaml
acoustic_backend: kwave    # 'kwave' (default) | 'kwave_elastic' | 'babelbrain' | 'jwave'
thermal_backend:  kwave    # 'kwave' (default) | 'babelbrain'
```

| | [k-Wave](doc_kwave.md) (default) | [BabelBrain](doc_babelbrain.md) | [jwave](doc_jwave.md) |
|---|---|---|---|
| **Acoustic solver** | Pseudo-spectral time domain | Staggered-grid FDTD ([BabelViscoFDTD](https://github.com/ProteusMRIgHIFU/BabelViscoFDTD)) | Helmholtz (frequency domain) or pseudo-spectral, in [JAX](https://github.com/jax-ml/jax) |
| **Thermal solver** | k-Wave `kWaveDiffusion` (Pennes) | BabelBrain Pennes BHTE | none (use k-Wave thermal) |
| **Tissue acoustics** | PRESTUS `medium_properties` | PRESTUS `medium_properties` (injected; fluid bone), or BabelBrain's own HU mapping via `babelbrain.ct_mapping` | PRESTUS `medium_properties` |
| **Tissue thermal properties** | PRESTUS `medium_properties` | BabelBrain's built-in table (differs in absorption fraction, see below) | none |
| **GPU** | [CUDA](https://developer.nvidia.com/cuda-toolkit) (C++ binary), MATLAB GPU arrays | CUDA, [OpenCL](https://www.khronos.org/opencl/), [Metal](https://developer.apple.com/metal/) | CUDA (JAX) |
| **Not supported** | none | Axisymmetric grids | Axisymmetric grids |
| **Extra install** | k-Wave submodule (+ binaries for `cpp_*`) | BabelBrain submodule + [conda](https://docs.conda.io/) env | [Python](https://www.python.org/) env with jwave/JAX |

Valid combinations:

| `acoustic_backend` | `thermal_backend` | |
|---|---|---|
| `kwave` / `kwave_elastic` / `jwave` | `kwave` | Supported |
| `babelbrain` | `kwave` | Supported |
| `babelbrain` | `babelbrain` | Supported |
| `kwave` / `kwave_elastic` / `jwave` / external NIfTI | `babelbrain` | Rejected at config load (BabelBrain thermal needs BabelBrain's `DataForSim.h5`) |

After setting up or updating a BabelBrain or jwave environment, run `smoke_test_backends` (in `tests/`) for a quick check of imports, dependencies and bridge scripts.

## Comparing Results Across Backends

- **Acoustics: designed to be comparable.** The three fluid-skull backends (k-Wave, BabelBrain, jwave) use the same drive phases and the same PRESTUS tissue properties. BabelBrain models bone as a fluid, like k-Wave and jwave (shear waves off), and its output is calibrated against the analytic free-water field (`babelbrain.water_amp_calibration`). What remains is solver discretisation and BabelBrain's bone-attenuation Q-correction. Compare against a k-Wave run before relying on absolute BabelBrain values. Setting `babelbrain.ct_mapping` switches BabelBrain to its own HU-to-property models, and the backends then no longer share tissue acoustics.
- **Elastic skull: intentionally different.** [`kwave_elastic`](doc_kwave.md#elastic-skull-model-shear-waves) adds shear waves in the skull. In an ex vivo comparison ([Gao et al.](https://doi.org/10.1109/TBME.2025.3603836)), fluid models overestimated the intracranial pressure by ~40 % at shallow cortical targets. Use `elastic.shear_mapping: off` for a fluid reference run through the same solver. The elastic model is a draft and not yet validated in PRESTUS.
- **Thermal: same acoustic amplitude, different tissue table (known limitation).** BabelBrain thermal receives the in-brain peak ISPPA of the same scaled field that k-Wave thermal consumes, so the heat source amplitude matches. BabelBrain's built-in thermal table matches the PRESTUS defaults for conductivity, specific heat and perfusion, except trabecular conductivity (0.31 vs 0.32 W/m/°C), but uses a different **absorption fraction**: skin/brain 0.85 (PRESTUS 1.0), skull 0.15–0.16 (PRESTUS 0.28). Expect lower BabelBrain heating, particularly in skull. PRESTUS thermal `medium_properties` (and uncertainty variants) cannot be passed to BabelBrain thermal; use `thermal_backend: kwave` where they must apply. `timing.post_ptri_dur` is applied as a post-sequence cooldown (`PauseBetweenGroupedSonications`); ISO CEM43 is not computed. Details: [BabelBrain thermal](doc_babelbrain.md#thermal-simulation).
