# Elastic Skull (Shear Waves)

Bone carries shear waves as well as pressure waves, and at oblique incidence part of the energy crosses the skull as shear. The elastic k-Wave model simulates both. This page covers its settings, the elastic skull properties, and its limits and cost.

On this page:

- [Elastic Skull Model (Shear Waves)](#elastic-skull-model-shear-waves): settings, limitations, and cost.

!!! note "At a glance"
    - **When you need it:** When the beam meets the skull obliquely, or to test how much shear waves change a result.
    - **Before you start:** A working k-Wave setup ([k-Wave](doc_kwave.md)).

## Elastic Skull Model (Shear Waves)

`acoustic_backend: kwave_elastic` runs k-Wave's `pstdElastic3D` / `pstdElastic2D` instead of `kspaceFirstOrder*`, so the skull supports shear as well as compressional waves. The motivation is [Gao et al. (2026)](https://doi.org/10.1109/TBME.2025.3603836), who compared fluid and elastic k-Wave skull models against hydrophone scans through an ex vivo skull at normal and 40° incidence for shallow (cortical) foci. In those scans the fluid skull model overestimated the intracranial pressure by ~40 % on average and the −6 dB focal area by up to 170 %. The elastic model was also less sensitive to small transducer position errors.

**Status: draft.** Not yet validated in PRESTUS against measurements or the [ITRUSST](https://itrusst.org/) benchmarks.

### Elastic Model Settings

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

Compressional speed and density come from the usual medium setup (layer values or pseudo-CT mapping). The shear speed and the two attenuations follow [Pichardo et al. (2017)](https://doi.org/10.1088/1361-6560/aa7ccc), the viscoelastic model behind [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain), evaluated at the transducer's fundamental frequency. Expect the elastic model to be most reliable at low frequencies: Pichardo et al. found good agreement at 270 kHz (delay error 5 % of a wavelength) but only global agreement at 836 kHz (17 %), where trabecular scattering degraded the spatial prediction of an isotropic model.

??? note "Skull elastic properties"
    Pichardo et al. (2017) fitted shear speed and attenuation to transcranial scans at 40° incidence (10 locations, 3 formalin-fixed skulls, 270 and 836 kHz). Every value is evaluated at the transducer's fundamental frequency f0:

    | Property | Pichardo et al. 2017 | PRESTUS default at f0 |
    |---|---|---|
    | Shear speed | cS = 0.458·ρ + 545.52 (270 kHz); cS = 0.386·ρ + 815.51 (836 kHz) | slope and intercept interpolated linearly in f. At 500 kHz: 0.429·ρ + 655.2, i.e. 1448 m/s at ρ = 1850 (cortical layer), 1384 m/s at ρ = 1700 (trabecular layer). Extrapolated, with a warning, outside 270–836 kHz |
    | Shear attenuation | 57 Np/m (270 kHz), 373 Np/m (836 kHz), one value for the whole skull | power law through both points: 43.63 dB/(cm·MHz^1.662); 0.68 Np/cm at 300 kHz, 1.59 Np/cm at 500 kHz |
    | Longitudinal attenuation | 33 Np/m (270 kHz), 231 Np/m (836 kHz), from [Pichardo et al. 2011](https://doi.org/10.1088/0031-9155/56/1/014); the shear values were fitted with these | power law through both points: 27.31 dB/(cm·MHz^1.722); 0.40 Np/cm at 300 kHz, 0.95 Np/cm at 500 kHz |

    The shear-speed fits are BabelBrain's `DensityToSSoSPichardo` before averaging: BabelBrain uses their mean over the two frequencies (0.422·ρ + 680.5).

    **Why the skull longitudinal attenuation is replaced.** The shear attenuation was fitted together with this longitudinal attenuation, so the pair is used as a set. The fluid/pCT skull values are fits for a fluid model, which have to absorb mode-conversion losses into the longitudinal attenuation. Set `long_attenuation: medium` to keep them anyway. The sources disagree by about 2×: Gao et al. (2026) used 0.5 Np/cm at 500 kHz with the same k-Wave solver (Pichardo pair: 0.95; PRESTUS fluid default 13.3 dB/(cm·MHz): 0.77).

    **Poisson-ratio floor.** Combining a density-based shear speed with a density-based compressional speed can give unphysical bone. For example, the Gao et al. shear mapping with the [k-Plan](https://k-plan.io/) cL mapping gives ν ≈ 0 at ρ ≈ 1300 kg/m³; BabelBrain disables shear on its CT path for this reason. The Pichardo mapping is better behaved, but at ρ = 1200 kg/m³ it still gives ν from 0.17 (250 kHz) down to ≈0 (836 kHz), because cS rises with frequency while cL does not. `poisson_min` caps cS in those voxels and the run log reports how many.

    **Resolution.** The slowest shear wave sets the grid requirement. A warning is issued below `min_ppw_shear` (6; Gao et al. used 6.6 with k-Wave, Pichardo et al. 8 with 4th-order FDTD). At 0.5 mm, 300 kHz gives ≈7.4 points per wavelength for cS = 1105 m/s (ρ = 1200); 500 kHz needs ≈0.39 mm.

??? note "How the elastic run differs from the fluid run"
    | Step | Fluid `kwave` | `kwave_elastic` |
    |---|---|---|
    | Medium | `sound_speed`, `density`, `alpha_coeff`/`alpha_power` | the same maps, plus `sound_speed_shear` and `alpha_coeff_shear` in the skull and (by default) `elastic.alpha_coeff_long` as skull longitudinal attenuation (`medium_elastic_properties`, called from `acoustic_wrapper` with the final density map and skull masks) |
    | Absorption | power law, `alpha_power` (remapped to y = 2 by `fit_alpha_power`) | Kelvin-Voigt, fixed y = 2; every coefficient is rescaled so the attenuation at the source frequency is exact |
    | Source | pressure source `p_mask`, `p` | equivalent normal-stress source `s_mask = p_mask`, `sxx = syy = szz = −p`. k-Wave scales both by 2·dt·c0/dx and records p = −(sxx+syy+szz)/3, so the injected pressure is the same |
    | Solver | `kspaceFirstOrder3D/2D/AS`, MATLAB or C++ | `pstdElastic3D/2D`, MATLAB CPU or GPU (`DataCast`) |
    | Free-water baseline | fluid solver | same elastic solver (cS = 0), so I<sub>SPPA</sub> scaling uses one forward model |

    The fluid `kwave_medium` returned to the pipeline is unchanged; thermal simulation, [NIfTI](https://nifti.nimh.nih.gov/) export and the report use it as for a fluid run. With `simulation.debug: 1`, the shear speed and shear absorption grids are written to `debug/medium/matrix_*_shear*.nii.gz`.

### Elastic Model Limitations and Cost

- **Compute and memory:** `pstdElastic3D` holds ≈65–75 single-precision arrays of the padded grid (split stress and velocity fields, staggered material maps, derivatives), i.e. ≈0.27–0.30 GB per million voxels, or ≈26–30 GB for a 512×512×384 padded grid. Runtime is several times that of the C++/CUDA fluid solver (MATLAB time loop, ~3× the FFT work); Gao et al. report ~15 min for 320×256×384 on an RTX 3090 at CFL 0.3. Keep `grid.source_cfl` ≤ 0.3; the PRESTUS default 0.15 is conservative.
- **Not supported:** axisymmetric grids and C++/CUDA binaries. Nonlinearity (`BonA`) is not passed to the elastic solver.
- **Property uncertainty:** the shear properties come from few ex vivo skulls (Pichardo et al. 2017: 3 formalin-fixed skulls plus 1 for validation; Gao et al. 2026: 1 skull). Gao et al. note that no single property set fits all incidence angles, so treat absolute values with caution.
- **Caches:** the acoustic result cache is not keyed on the backend. Use a distinct `io.output_affix` when comparing `kwave` and `kwave_elastic` runs.

??? note "References"
    - Gao, Y., et al. Influence of shear waves on transcranial ultrasound propagation in cortical brain regions. *IEEE Trans. Biomed. Eng.* **73**(3):1298–1306 (2026). https://doi.org/10.1109/TBME.2025.3603836
    - Pichardo, S., et al. A viscoelastic model for the prediction of transcranial ultrasound propagation: application for the estimation of shear acoustic properties in the human skull. *Phys. Med. Biol.* **62**(17):6938 (2017). https://doi.org/10.1088/1361-6560/aa7ccc
    - White, P. J., Clement, G. T., & Hynynen, K. Longitudinal and shear mode ultrasound propagation in human skull bone. *Ultrasound Med. Biol.* **32**(7):1085–1096 (2006). https://doi.org/10.1016/j.ultrasmedbio.2006.03.015
