# Simulation backends

!!! note "At a glance"
    - **What:** The simulation engines PRESTUS can use: k-Wave (default), an elastic k-Wave skull model, BabelBrain and jwave.
    - **When you need it:** You do not need anything but k-Wave to start. Read the k-Wave section to choose between CPU, GPU, MATLAB and C++.
    - **Before you start:** k-Wave (installed with PRESTUS).

PRESTUS runs the acoustic and the thermal simulation through interchangeable backends. Pre-processing (segmentation, medium mapping, grid, transducer placement) and post-processing (focal metrics, CEM43, NIfTI export, HTML report) are shared, so switching backend is a config change:

```yaml
acoustic_backend: kwave    # 'kwave' (default) | 'kwave_elastic' | 'babelbrain' | 'jwave'
thermal_backend:  kwave    # 'kwave' (default) | 'babelbrain'
```

| | k-Wave (default) | [BabelBrain](doc_babelbrain.md) | [jwave](doc_jwave.md) |
|---|---|---|---|
| **Acoustic solver** | Pseudo-spectral time domain | Staggered-grid FDTD (BabelViscoFDTD) | Helmholtz (frequency domain) or pseudo-spectral, in JAX |
| **Thermal solver** | k-Wave `kWaveDiffusion` (Pennes) | BabelBrain Pennes BHTE | — (use k-Wave thermal) |
| **Tissue acoustics** | PRESTUS `medium_properties` | PRESTUS `medium_properties` (injected; fluid bone) — or BabelBrain's own HU mapping via `babelbrain.ct_mapping` | PRESTUS `medium_properties` |
| **Tissue thermal properties** | PRESTUS `medium_properties` | BabelBrain's built-in table (differs in absorption fraction, see below) | — |
| **GPU** | CUDA (C++ binary), MATLAB GPU arrays | CUDA, OpenCL, Metal | CUDA (JAX) |
| **Not supported** | — | Axisymmetric grids | Axisymmetric grids, `kWaveArray` sources |
| **Extra install** | k-Wave submodule (+ binaries for `cpp_*`) | BabelBrain submodule + conda env | Python env with jwave/JAX |

Valid combinations:

| `acoustic_backend` | `thermal_backend` | |
|---|---|---|
| `kwave` / `kwave_elastic` / `jwave` | `kwave` | Supported |
| `babelbrain` | `kwave` | Supported |
| `babelbrain` | `babelbrain` | Supported |
| `kwave` / `kwave_elastic` / `jwave` / external NIfTI | `babelbrain` | Rejected at config load (BabelBrain thermal needs BabelBrain's `DataForSim.h5`) |

### Are results comparable across backends?

- **Acoustics — designed to be comparable.** The three fluid-skull backends (k-Wave, BabelBrain, jwave) use the same drive phases and the same PRESTUS tissue properties. BabelBrain models bone as a fluid, like k-Wave and jwave (shear waves off), and its output is calibrated against the analytic free-water field (`babelbrain.water_amp_calibration`). What remains is solver discretisation and BabelBrain's bone-attenuation Q-correction. Compare against a k-Wave run before relying on absolute BabelBrain values. Setting `babelbrain.ct_mapping` switches BabelBrain to its own HU→property models, and the backends then no longer share tissue acoustics.
- **Elastic skull — intentionally different.** `kwave_elastic` adds shear waves in the skull, which typically *lowers* the transcranial pressure (Gao et al. report ~40 % overestimation by fluid models at shallow cortical targets). Use `elastic.shear_mapping: off` for a fluid reference run through the same solver.
- **Thermal — same acoustic amplitude, different tissue table (known limitation).** BabelBrain thermal receives the in-brain peak ISPPA of the same scaled field that k-Wave thermal consumes, so the heat source amplitude matches. BabelBrain's built-in thermal table matches the PRESTUS defaults for conductivity, specific heat and perfusion, but uses a different **absorption fraction**: skin/brain 0.85 (PRESTUS 1.0), skull 0.15–0.16 (PRESTUS 0.28). Expect lower BabelBrain heating, particularly in skull. PRESTUS thermal `medium_properties` (and uncertainty variants) cannot currently be passed to BabelBrain thermal; use `thermal_backend: kwave` where they must apply. `timing.post_ptri_dur` is applied as a post-sequence cooldown (`PauseBetweenGroupedSonications`); ISO CEM43 is not computed. Details: [BabelBrain thermal](doc_babelbrain.md#thermal-simulation).

After setting up or updating a BabelBrain or jwave environment, run `smoke_test_backends` (in `tests/`) for a quick check of imports, dependencies and bridge scripts.

---

## k-Wave backend

PRESTUS allows k-Wave deployment using different computing setups (`parameters.simulation.code_type`). This will primarily affect acoustic simulations. With either of the GPU variants, thermal simulations will also use GPU acceleration (if k-Wave 1.4.1 is provided).

- `matlab_cpu`
    - `kspaceFirstOrder3D` | `kspaceFirstOrder2D` | `kspaceFirstOrderAS` 
- `matlab_gpu`
    - MATLAB with GPU acceleration via `DataCast`.
    - Requires Parallel Computing Toolbox
- `cpp_cpu` (C++)
    - `kspaceFirstOrder3DC` warps the C++ binary `kspaceFirstOrder-OMP`. 
    - Requires the specific binary from [k-Wave downloads](http://www.k-wave.org/download.php). 
- `cpp_gpu`
    - `kspaceFirstOrder3DG` wraps the C++/CUDA binary `kspaceFirstOrder-CUDA`, offloading computation to NVIDIA GPUs.
    - Requires a CUDA-capable Nvidia GPU and the specific binary from [k-Wave downloads](http://www.k-wave.org/download.php).

To use C++ with GPU support you may need to recompile the binary for your GPU version (e.g., for the A100 GPUs on the Donders HPC). 

- Download the [LINUX source files for k-Wave 1.3](http://www.k-wave.org/download.php). This includes the Makefile to compile the binary.
- To compile the binaries, you can follow the instruction by [TU Delft](https://qiweb.tudelft.nl/sysman/kwave_hpc.html). Note that you need to include your GPU type in the Makefile.
- Example scripts for compiling binaries on the Donders HPC can be found in `PRESTUS/examples/hpc_compile...`
- For the Donders HPC, you can download a precompiled binary [here](https://github.com/jkosciessa/PRESTUS_bin/raw/refs/heads/main/donders_kwave_bin/kspaceFirstOrder-CUDA).
- Place the updated binary into `PRESTUS/external/k-Wave/binaries/`

---

## k-Wave elastic backend (shear waves)

`acoustic_backend: kwave_elastic` runs k-Wave's `pstdElastic3D` / `pstdElastic2D` instead of `kspaceFirstOrder*`, so the skull supports shear as well as compressional waves. The motivation is Gao et al. (2026), who compared fluid and elastic k-Wave skull models against hydrophone scans through an ex vivo skull at normal and 40° incidence for shallow (cortical) foci. In those scans the fluid skull model overestimated the intracranial pressure by ~40 % on average and the −6 dB focal area by up to 170 %. The elastic model was also less sensitive to small transducer position errors.

**Status: draft.** Not yet validated in PRESTUS against measurements or the ITRUSST benchmarks.

### Configuration

```yaml
acoustic_backend: kwave_elastic
simulation:
  code_type: matlab_gpu      # cpp_* is overridden to matlab_* (no C++/CUDA elastic solver)
grid:
  axisymmetric: 0            # not supported
elastic:                     # defaults shown; see doc_parameters.md
  shear_mapping: pichardo    # pichardo | off
  shear_density_min: 1200
  poisson_min: 0.1
  alpha_coeff_shear: 43.63   # dB/(cm MHz^y)
  alpha_power_shear: 1.662
  alpha_coeff_long: 27.31    # dB/(cm MHz^y)
  alpha_power_long: 1.722
  long_attenuation: elastic
```

### Skull elastic properties

Compressional speed and density come from the usual medium setup (layer values or pseudo-CT mapping). The elastic settings follow Pichardo et al. (2017), the viscoelastic model behind BabelBrain, who fitted shear speed and attenuation to transcranial scans at 40° incidence (10 locations, 3 formalin-fixed skulls, 270 and 836 kHz). Every value is evaluated at the transducer's fundamental frequency f0:

| Property | Pichardo et al. 2017 | PRESTUS default at f0 |
|---|---|---|
| Shear speed | cS = 0.458·ρ + 545.52 (270 kHz); cS = 0.386·ρ + 815.51 (836 kHz) | slope and intercept interpolated linearly in f. At 500 kHz: 0.429·ρ + 655.2, i.e. 1448 m/s at ρ = 1850 (cortical layer), 1384 m/s at ρ = 1700 (trabecular layer). Extrapolated, with a warning, outside 270–836 kHz |
| Shear attenuation | 57 Np/m (270 kHz), 373 Np/m (836 kHz), one value for the whole skull | power law through both points: 43.63 dB/(cm·MHz^1.662); 0.68 Np/cm at 300 kHz, 1.59 Np/cm at 500 kHz |
| Longitudinal attenuation | 33 Np/m (270 kHz), 231 Np/m (836 kHz), from Pichardo et al. 2011; the shear values were fitted with these | power law through both points: 27.31 dB/(cm·MHz^1.722); 0.40 Np/cm at 300 kHz, 0.95 Np/cm at 500 kHz |

The shear-speed fits are BabelBrain's `DensityToSSoSPichardo` before averaging: BabelBrain uses their mean over the two frequencies (0.422·ρ + 680.5).

**Why the skull longitudinal attenuation is replaced:** the shear attenuation was fitted together with this longitudinal attenuation, so the pair is used as a set. The fluid/pCT skull values are fits for a fluid model, which have to absorb mode-conversion losses into the longitudinal attenuation. Set `long_attenuation: medium` to keep them anyway. The sources disagree by about 2×: Gao et al. (2026) used 0.5 Np/cm at 500 kHz with the same k-Wave solver (Pichardo pair: 0.95; PRESTUS fluid default 13.3 dB/(cm·MHz): 0.77).

**Validity:** Pichardo et al. found good agreement at 270 kHz (delay error 5 % of a wavelength) but only global agreement at 836 kHz (17 %), where trabecular scattering degraded the spatial prediction of an isotropic model. Expect the elastic model to be most reliable at low frequencies.

**Poisson-ratio floor:** combining a density-based shear speed with a density-based compressional speed can give unphysical bone. For example, the Gao et al. shear mapping with the k-Plan cL mapping gives ν ≈ 0 at ρ ≈ 1300 kg/m³; BabelBrain disables shear on its CT path for this reason. The Pichardo mapping is better behaved, but at ρ = 1200 kg/m³ it still gives ν from 0.17 (250 kHz) down to ≈0 (836 kHz), because cS rises with frequency while cL does not. `poisson_min` caps cS in those voxels and the run log reports how many.

**Resolution:** the slowest shear wave sets the grid requirement. A warning is issued below `min_ppw_shear` (6; Gao et al. used 6.6 with k-Wave, Pichardo et al. 8 with 4th-order FDTD). At 0.5 mm, 300 kHz gives ≈7.4 points per wavelength for cS = 1105 m/s (ρ = 1200); 500 kHz needs ≈0.39 mm.

### How it works

| Step | Fluid `kwave` | `kwave_elastic` |
|---|---|---|
| Medium | `sound_speed`, `density`, `alpha_coeff`/`alpha_power` | the same maps, plus `sound_speed_shear` and `alpha_coeff_shear` in the skull and (by default) `elastic.alpha_coeff_long` as skull longitudinal attenuation (`medium_elastic_properties`, called from `acoustic_wrapper` with the final density map and skull masks) |
| Absorption | power law, `alpha_power` (remapped to y = 2 by `fit_alpha_power`) | Kelvin-Voigt, fixed y = 2; every coefficient is rescaled so the attenuation at the source frequency is exact |
| Source | pressure source `p_mask`, `p` | equivalent normal-stress source `s_mask = p_mask`, `sxx = syy = szz = −p`. k-Wave scales both by 2·dt·c0/dx and records p = −(sxx+syy+szz)/3, so the injected pressure is the same |
| Solver | `kspaceFirstOrder3D/2D/AS`, MATLAB or C++ | `pstdElastic3D/2D`, MATLAB CPU or GPU (`DataCast`) |
| Free-water baseline | fluid solver | same elastic solver (cS = 0), so ISPPA scaling uses one forward model |

The fluid `kwave_medium` returned to the pipeline is unchanged; thermal simulation, NIfTI export and the report use it as before. With `simulation.debug: 1`, the shear speed and shear absorption grids are written to `debug/medium/matrix_*_shear*.nii.gz`.

### Limitations and cost

- **Compute and memory:** `pstdElastic3D` holds ≈65–75 single-precision arrays of the padded grid (split stress and velocity fields, staggered material maps, derivatives), i.e. ≈0.27–0.30 GB per million voxels — ≈26–30 GB for a 512×512×384 padded grid. Runtime is several times that of the C++/CUDA fluid solver (MATLAB time loop, ~3× the FFT work); Gao et al. report ~15 min for 320×256×384 on an RTX 3090 at CFL 0.3. Keep `grid.source_cfl` ≤ 0.3; the PRESTUS default 0.15 is conservative.
- **Not supported:** axisymmetric grids and C++/CUDA binaries. Nonlinearity (`BonA`) is not passed to the elastic solver.
- **Property uncertainty:** the shear properties come from few ex vivo skulls (Pichardo et al. 2017: 3 formalin-fixed skulls plus 1 for validation; Gao et al. 2026: 1 skull). Gao et al. note that no single property set fits all incidence angles, so treat absolute values with caution.
- **Caches:** the acoustic result cache is not keyed on the backend. Use a distinct `io.output_affix` when comparing `kwave` and `kwave_elastic` runs.

<span style="font-size: 12px;">
**References:**
Gao, Y., et al. Influence of shear waves on transcranial ultrasound propagation in cortical brain regions. *IEEE Trans. Biomed. Eng.* **73**(3):1298–1306 (2026). https://doi.org/10.1109/TBME.2025.3603836 <br>
Pichardo, S., et al. A viscoelastic model for the prediction of transcranial ultrasound propagation: application for the estimation of shear acoustic properties in the human skull. *Phys. Med. Biol.* **62**(17):6938 (2017). https://doi.org/10.1088/1361-6560/aa7ccc <br>
White, P. J., Clement, G. T., & Hynynen, K. Longitudinal and shear mode ultrasound propagation in human skull bone. *Ultrasound Med. Biol.* **32**(7):1085–1096 (2006). https://doi.org/10.1016/j.ultrasmedbio.2006.03.015
</span>
