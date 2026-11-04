# (pseudo-)CT

!!! note "At a glance"
    - **What:** Optional subject-specific skull properties: a pseudo-CT estimated from a UTE/PETRA MRI replaces the single skull value of each layer.
    - **When you need it:** When you have a UTE/PETRA scan and want skull properties that follow the individual skull.
    - **Before you start:** A UTE/PETRA image per subject, FSL and ANTs.

![pct_ute_density](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/ute_pct_density.png)

PRESTUS turns a UTE-based image into pseudo-Hounsfield units (pHU) and maps them to density, sound speed and attenuation within the skull layer of a `layered` simulation.

!!! warning "Experimental"
    (pseudo-)HU mapping is an experimental feature in active development. The most suitable mapping remains an active area of research; treat all mappings as experimental. See [this issue](https://github.com/Donders-Institute/PRESTUS/issues/43).

## Setting up a pseudo-CT run

1. Acquire an MR sequence with detailed bone contrast (PETRA UTE; example parameters below).
2. Run SimNIBS with the T1w and the PETRA UTE scan (instead of the T2w) as inputs, either in PRESTUS (**default**) or in the SimNIBS GUI. The UTE image is the one in `path.t2_pattern`, which SimNIBS registers to the T1 as `m2m_sub-XXX/T2_reg.nii.gz` (`UTE_reg.nii.gz` is also read).
3. Set `pct.enabled: 1`. PRESTUS runs `pct_create_pseudoCT.sh` after segmentation if `pseudoCT.nii.gz` does not exist yet. This needs **FSL** (`fslmaths`), **ANTs** (`N4BiasFieldCorrection`, `SmoothImage`) and **MATLAB** on `PATH`; set `startup.fsl_bin_path` and `startup.ants_bin_path` so PRESTUS can prepend them. On HPC, also set `startup.matlab_bin_path`.

Minimum config:

```yaml
path:
  t2_pattern: 'sub-%1$03d_UTE.nii.gz'  # UTE image (used for both SimNIBS and pseudoCT)

startup:
  fsl_bin_path: '/opt/fsl/bin'           # prepended to PATH for fslmaths
  ants_bin_path: '/opt/ants/bin'          # prepended to PATH for N4BiasFieldCorrection, SmoothImage
  matlab_bin_path: '/opt/matlab/R2024a/bin/matlab'  # required on HPC; local falls back to matlabroot

pct:
  enabled: 1
  skull_mapping: 'kosciessa' # UTE→HU linear mapping algorithm
  debug: 1                   # keep intermediate images in the debug_pct/ subfolder
  mapping_density: 'k-plan'
  density_calibration: ''    # HU->density table for k-plan; '' = k-Plan default (see below)
  mapping_soundspeed: 'k-plan'
  mapping_attenuation: 'k-plan'
```

Outputs go to `<path.pct>/sub-XXX/`, or to `<parent of path.seg>/pct/sub-XXX/` when `path.pct` is empty:

| File / folder | Content |
|---|---|
| `pseudoCT.nii.gz`, `tissues_mask.nii.gz` | pseudo-CT and tissue mask |
| `pseudoCT.json` | sidecar: skull mapping, segmentation folder, PRESTUS version |
| `debug_pct/` | intermediate images and the log |
| `log_hpc/` | cluster job logs |

**Reuse.** An existing `pseudoCT.nii.gz` is used and not regenerated. Its `pseudoCT.json` must record the current `pct.skull_mapping`; otherwise the run stops with `prestus:cache:signatureMismatch` (see `io.cache_signature_check`). To regenerate, move the folder aside or set `path.pct` to another folder. A pseudo-CT without a sidecar is used with the warning `prestus:cache:unverified`. A pseudo-CT inside `m2m_sub-XXX/` is used with the warning `PRESTUS:pct:legacyLocation`; move it to the pseudo-CT folder. See [Starting from earlier stage outputs](doc_outputs.md#starting-from-earlier-stage-outputs).

**Manual.** You can also run `pct_create_pseudoCT.sh` in bash before starting PRESTUS. An example wrapper is at `examples/createPseudoCT.sh`. Required modules: SimNIBS, FSL, ANTs.

??? note "Details: example PETRA UTE parameters"
    From Carpino et al. (2023):

    - TR: 3.32 ms
    - TE: 0.07 ms
    - voxel size: 0.8 mm3
    - 352 sagittal slices
    - flip angle: 2°
    - FoV: 294 mm

## From UTE to pseudo-HU

### pseudoCT generation

`pct.skull_mapping` sets the linear UTE→pHU mapping of the skull. Soft tissue is set to 42 HU and air to -1000 HU for every mapping.

| `skull_mapping` | pHU |
|---|---|
| `kosciessa` (**default**) | individualised: a line fitted per subject through pHU_trabecular=300 & pHU_trabecular_cortical = 700 |
| `miscouridou` | pHU = −2085 UTE + 2329 |
| `carpino` | pHU = -2194 UTE + 2236 |
| `wiesinger` | pHU = -2000 (UTE-1) + 42 |
| `treeby` | pHU = -2929.6 UTE + 3274.9 (ucl-bug/petra-to-ct; fitted to 80 kVp CT: use with `density_calibration: ucl_lowdose_80kvp`, see [below](#pairing-a-utehu-mapping-with-its-ct-calibration)) |
| `custom` | pHU = `pct.custom_slope` UTE + `pct.custom_intercept` (any other linear map) |

??? note "Details: the generation steps"
    - UTE: threshold at 0, bias field correction, then normalisation by the soft tissue peak (normalised UTE intensity = 1). `pct_soft_tissue_peak` finds the peak.
    - pHU soft tissue (normalised UTE intensity = 1): 42 HU (Wiesinger et al., 2018).
    - pHU air: -1000 HU (Wiesinger et al., 2018; Miscouridou et al., 2022).
    - pHU skull: the linear mapping set by `pct.skull_mapping` (table above).
    - `kosciessa` (`pct_skullmapping_fit`) fits a least-squares line through three points from the UTE histograms of charm cortical and trabecular bone: the lower quarter-max crossing of the cortical histogram → 1900 HU, the crossover of the two histograms → 700 HU, and the upper quarter-max crossing of the trabecular histogram around its x-weighted peak → 300 HU. Bins 0.9–1.1 are replaced by a spline against soft-tissue partial volume. The anchors, width level and interpolation window are options of `pct_skullmapping_fit`.
    - 3D Gaussian smoothing (σ=0.8 voxels via ANTs' SmoothImage) prevents abrupt tissue transitions in the simulation grid. Within an eroded skull mask, unsmoothed values are kept to limit partial volume effects.

### Pairing a UTE→HU mapping with its CT calibration

A UTE→HU line returns HU on the scale of the CT it was fitted to. That scale depends on scanner, tube voltage (kVp) and reconstruction kernel: at 80 kVp, bone has about 1.35 times the HU it has at about 120 kVp. The simulation needs density, and turning HU into density needs the calibration of that same CT. **Set `pct.skull_mapping`, `pct.mapping_density: k-plan` and `pct.density_calibration` together:**

| `skull_mapping` | CT the line was fitted to | `mapping_density` | `density_calibration` |
|---|---|---|---|
| `kosciessa`, `miscouridou`, `carpino`, `wiesinger` | clinical CT, about 120 kVp (assumed; for `kosciessa` the anchors 1900/700/300 HU are assumed on this scale) | any | `''` (k-Plan default) |
| `treeby` | low-dose CT, 80 kVp, BONEPLUS kernel (UCL) | `k-plan` | `ucl_lowdose_80kvp` |
| `custom` with petra2density's line (−2080.02, 2133.23) | the cph2025 phantom calibration | `k-plan` | `cph2025` |
| `custom` with your own line | your CT (any scanner and kVp) | `k-plan` | path to that CT's phantom table |

```yaml
pct:
  skull_mapping: 'treeby'
  mapping_density: 'k-plan'
  density_calibration: 'ucl_lowdose_80kvp'   # the 80 kVp CT treeby was fitted to
  mapping_soundspeed: 'k-plan'               # follows from density
  mapping_attenuation: 'k-plan'
```

`pct.density_calibration` values:

| Value | Table |
|---|---|
| `''` (default) | k-Plan default (`config/ct_calibration/kplan_default.csv`) |
| `'ucl_lowdose_80kvp'` | UCL low-dose CT (GE Revolution, 80 kVp, BONEPLUS; CIRS 062M phantom): the CT behind `skull_mapping: treeby`. Bone has about 1.35× the default's HU for the same density |
| `'cph2025'` | petra2density's phantom calibration, the CT behind its line −2080.02·UTE + 2133.23 (use `skull_mapping: custom`). Bone has about 0.85× the default's HU |
| a path | any CSV with two columns, HU and density [kg/m³] (lines starting with `#` are comments), e.g. a phantom scan of your own scanner or a table that goes with a study-specific `custom` line |

Tables are applied piecewise linearly and extrapolated linearly beyond their end points. For a line fitted at another kVp, scan a density phantom with the same protocol and pass its table; no HU conversion between kVp is needed. `medium_setup` prints the resolved table and records it in the medium statistics.

Only `k-plan` reads the calibration; with another density mapping it is ignored, with a warning. The other mappings assume a fixed, clinical HU scale. Use them with lines on the clinical scale only; otherwise use `k-plan` with the matching calibration. PRESTUS warns when `skull_mapping: treeby` is used without `mapping_density: k-plan` and `density_calibration: ucl_lowdose_80kvp`.

??? note "Details: HU scale assumed by the other mappings"
    - `k-wave`: Schneider et al.; `marsac`: 300–2000 HU; `air-tissue` (needed by `pichardo` and `mcdannold`): about −1000 to 1000 HU.
    - `webb` sound speed and attenuation: GE 120 kVp, B kernel by default; scanner-specific values via `pct.webb_*`.
    - `aubry` and `mueller` rescale by the skull's own minimum and maximum HU. They are insensitive to the HU scale but remove between-subject differences in density level.

## From pseudo-HU to acoustic properties

With `pct.enabled: 1`, choose a mapping for each of density, sound speed and attenuation. The pseudo-CT is read from the pseudo-CT folder above. This is a subtype of a `layered` medium. Mappings apply only within the uniform skull mask. A property set to `"none"`, or with its field missing, keeps the uniform `skull` value. Defaults below are those in `config_default.yaml`.

### Mapping skull density

Set by `pct.mapping_density`.

| Value | Mapping |
|---|---|
| `k-plan` (**default**) | 4-part piecewise linear fit based on k-Plan defaults; `pct.density_calibration` replaces the table ([see above](#pairing-a-utehu-mapping-with-its-ct-calibration)) |
| `k-wave` | 4-part piecewise linear fit based on Schneider et al. (1996) |
| `marsac` | linear between `HU_min = 300` and `HU_max = 2000`, Marsac et al. (2017) |
| `aubry` | porosity-based mixture of bone and water, Aubry et al. (2003) |
| `air-tissue` | linear air-to-bone calibration (BabelBrain-compatible); **use with `pichardo` or `mcdannold` sound speed or attenuation** |

??? note "Details: density equations"
    Several mappings use the porosity of Aubry et al. (2003), with the max. skull HU denominator from Guo et al. (2019):

    $$ \phi_\text{skull} = 1 - \frac{\text{HU}}{\text{HU}_\text{max}} $$

    **`k-plan`.** Piece-wise linear mapping between HU and mass density in kg/m^3 ([k-Plan documentation](https://dispatch.k-plan.io/static/docs/planning-images.html#ct-calibration)):

    $$ \rho_\text{skull} = \text{fit_pairwiselinear}(\text{HU}, \text{density}) $$

    ```
    HU = [-990, 60, 1000, 1950]
    density = [1.2, 1060, 1530, 2150]
    ```

    The minimum is regularised to the configured water density. Sound speed and attenuation of the k-plan pipeline follow from density. k-Plan's bone segmentation starts at density values of 1150 kg/m3, which internally regularizes data ranges; PRESTUS does not enforce this. **Thermal simulations** use a fixed bone density:

    $$ \rho_\text{skull} = 1850 \, \text{kg/m}^3 $$

    **`k-wave`.** Pseudo-CT values are shifted by +1000 and thresholded at 300 to align with [hounsfield2density](http://www.k-wave.org/documentation/hounsfield2density.php):

    $$ \rho_\text{skull} = \text{hounsfield2density}(\text{HU} + 1000) $$

    Density is regularised at the minimum to the configured water density and at a max. skull density `rho_bone = 2100` kg/m3.

    **`marsac`.** Skull pHU are clipped to `HU_min = 300` and `HU_max = 2000`. Unlike Marsac et al., 2017, PRESTUS clips pHU < HU_min instead of excluding them from the skull mask, which would create skull holes.

    $$ \rho_\text{skull} = \rho_\text{water} + (\rho_\text{bone} - \rho_\text{water}) \cdot \frac{\text{HU} - \text{HU}_\text{min}}{\text{HU}_\text{max} - \text{HU}_\text{min}} $$

    with max. skull density `rho_bone = 2100` kg/m3.

    **`aubry`.** Mixture of bone and water density, both from the configuration:

    $$ \rho_\text{skull} = \rho_\text{water} \cdot \phi_\text{skull} + \rho_\text{bone} \cdot (1 - \phi_\text{skull}) $$

    **`air-tissue`.** Linear HU→density between air (HU ≈ −1000) and cortical bone (HU ≈ 1000), as `HUtoDensityAirTissue` in BabelBrain's Pichardo and McDannold workflows:

    $$ \rho_\text{skull} = 1.01237293 \cdot \text{HU} + 1013.67 $$

    The minimum is regularised to water density. `pichardo` and `mcdannold` need density from this calibration.

### Mapping skull sound speed

Set by `pct.mapping_soundspeed`.

| Value | Mapping |
|---|---|
| `k-plan` (**default**) | density-based, [k-Plan](https://dispatch.k-plan.io/static/docs/simulation-pipeline.html#converting-the-primary-planning-image-to-material-properties) |
| `marsac` | density-based, Marsac et al. (2017) |
| `aubry` | porosity-based, Aubry et al. (2003) |
| `webb` | linear in HU, Webb et al. (2018); scanner-specific slope and intercept via `pct.webb_sos_slope` and `pct.webb_sos_intercept` |
| `pichardo` | frequency-dependent ρ→c, Pichardo et al. 2D interpolation map; needs `pct.mapping_density: air-tissue` |
| `mcdannold` | polynomial ρ→c, McDannold et al.; needs density first (`air-tissue` or `marsac`) |

For `k-plan`, `webb`, `pichardo` and `mcdannold`, the minimum sound speed is regularised to water sound speed.

??? note "Details: sound speed equations"
    **`k-plan`.**

    $$ c_\text{skull} = 1.33 \cdot \rho_\text{skull} + 167 $$

    **`marsac`.**

    $$ c_\text{skull} = c_\text{water} + (c_\text{bone} - c_\text{water}) \cdot \frac{\rho_\text{skull} - \rho_\text{water}}{\rho_\text{bone} - \rho_\text{water}} $$

    with max. speed of sound in skull `c_skull = 3360` m/s and max. skull density `rho_bone = 2100` kg/m3. Water properties come from the configuration.

    **`aubry`.**

    $$ c_\text{skull} = c_\text{water} \cdot \phi_\text{skull} + c_\text{bone} \cdot (1 - \phi_\text{skull}) $$

    Tissue sound speeds come from the configuration.

    **`webb`.**

    $$ c_\text{skull} = \text{slope} \cdot \text{HU} + \text{intercept} \times 1000 $$

    Defaults (GE 120 kVp, B kernel, 0.5 × 0.6 mm resolution): slope = 0.75, intercept = 1.32 (c = 0.75 · HU + 1320 m/s). Set scanner-specific values from Webb et al. Tables I and II with `pct.webb_sos_slope` and `pct.webb_sos_intercept`.

    **`pichardo`.** Reads `MapPichardo.h5` from the BabelBrain submodule: longitudinal speed of sound [m/s] as a function of density [kg/m³] and frequency [MHz]. Interpolation is bilinear with nearest-neighbour extrapolation. MATLAB cannot read the Blosc-compressed H5 without the filter plugin, so PRESTUS uses the same table in `functions/medium/MapPichardo.mat`. Requires `pct.mapping_density: air-tissue` (or another density algorithm run first), the BabelBrain submodule initialised at `external/BabelBrain`, and the simulation frequency in `parameters.transducer(1).freq_hz`.

    **`mcdannold`.**

    $$ c_\text{skull} = \frac{1}{p(\rho_\text{skull})} $$

    where $p(\rho) = -2.79 \times 10^{-18}\,\rho^4 + 5.31 \times 10^{-16}\,\rho^3 + 1.69 \times 10^{-10}\,\rho^2 - 7.63 \times 10^{-7}\,\rho + 1.24 \times 10^{-3}$. Requires density first (`pct.mapping_density: air-tissue` or `marsac`).

### Mapping skull attenuation

Set by `pct.mapping_attenuation`.

| Value | Mapping |
|---|---|
| `k-plan` (**default**) | uniform: the skull `alpha_coeff` and `alpha_power` of `medium_properties`; to replicate k-Plan set `alpha_coeff = 13.3` and `alpha_power = 1` |
| `mueller` | HU-based, Mueller et al. (2017), bounds from Yaakub et al. (2023) |
| `aubry` | porosity-based, Aubry et al. (2003) |
| `webb` | exponential HU-based model, Webb et al. (2020); scanner-specific values via `pct.webb_att_alpha0`, `pct.webb_att_beta` and `pct.webb_att_c` |
| `pichardo` | frequency-dependent ρ→Att from Pichardo et al. 2D interpolation map; prerequisites as for `pichardo` sound speed |
| `mcdannold` | polynomial ρ→Att from McDannold et al. with linear frequency scaling; needs density first |

??? note "Details: attenuation equations"
    **`k-plan`.**

    $$ \alpha_\text{skull} = \alpha_\text{coeff} $$

    Reads `alpha_coeff` and `alpha_power` of `medium_properties.skull` from the config, and warns when they differ from k-Plan's fixed `13.3` and `1`.

    **`mueller`.**

    $$ \alpha_\text{skull} = \alpha_\text{min} + (\alpha_\text{max} - \alpha_\text{min}) \left(1 - \frac{\text{HU} - \text{HU}_\text{min}}{\text{HU}_\text{max} - \text{HU}_\text{min}}\right)^{0.5} $$

    with `α_min` = 4 and `α_max` = 8.7 (see Yaakub et al., 2023). These are 𝛼(f) at 500 kHz (see Aubry, J.-F., 2022 for benchmark simulations). PRESTUS converts them to `alpha_0` in `𝛼(f) = alpha_0 x f[MHz] ^ y` as `alpha_0 = 𝛼(f)/0.5^y`, with `y` the skull `alpha_power`.

    **`aubry`.**

    $$ \alpha_\text{skull} = \alpha_\text{min} + (\alpha_\text{max} - \alpha_\text{min}) \cdot \phi_\text{skull}^{0.5} $$

    with `α_min = 0.2` and `α_max = 8` (Aubry et al., 2003).

    **`webb`.**

    $$ \alpha_\text{skull}(f) = \alpha_0 \cdot \left(\frac{f}{1\,\text{MHz}}\right)^{\beta} \cdot \exp(c \cdot \text{HU}) $$

    Output is k-Wave `alpha_coeff` [dB/cm/MHz$^\beta$] with `alpha_power = β`. Defaults (GE 120 kVp, B kernel, 0.5 × 0.6 mm): α₀ = 26, β = 1.3, c = −0.0016. Set scanner-specific values from Webb et al. Table IV with `pct.webb_att_alpha0`, `pct.webb_att_beta` and `pct.webb_att_c`.

    **`pichardo`.** Reads longitudinal attenuation [Np/m] from `MapPichardo.h5` as a function of density [kg/m³] and frequency [MHz]. The value at the simulation frequency is converted to k-Wave `alpha_coeff` [dB/cm/MHz] with `alpha_power = 1` (linear frequency scaling assumed). Prerequisites as for `pichardo` sound speed; density must be filled first.

    **`mcdannold`.**

    $$ \alpha_\text{skull}(f) = p(\rho_\text{skull}) \cdot \frac{f}{660\,\text{kHz}} $$

    where $p(\rho) = 1.36 \times 10^{-10}\,\rho^4 - 1.41 \times 10^{-6}\,\rho^3 + 5.40 \times 10^{-3}\,\rho^2 - 9.02\,\rho + 5.71 \times 10^{3}$ [Np/m]. Converted to k-Wave `alpha_coeff` [dB/cm/MHz] with `alpha_power = 1`. Requires density first.

### Predicting the modelled skull properties

To see the skull properties a run will model before simulating, call `pct_predict_skull_properties(parameters)` or use **Predict skull properties** on the Nexus Pseudo-CT sub-tab. It applies the configured mappings to the pseudo-CT skull voxels and writes their distributions to `cache/sub-NNN_<medium>_pct_prediction<affix>.json`, alongside the skull values of the uncertainty bound variants. See [Pseudo-CT runs](doc_uncertainty.md#pseudo-ct-runs).

??? note "References"
    - Aubry, J.-F., et al. Experimental demonstration of noninvasive transskull adaptive focusing based on prior computed tomography scans. *J Acoust Soc Am* **113**, 84–93 (2003).
    - Aubry, J.-F. et al. Benchmark problems for transcranial ultrasound simulation: Intercomparison of compressional wave modelsa). J. Acoust. Soc. Am. 152, 1003–1019 (2022).
    - Carpino et al. (2024). Transcranial ultrasonic stimulation of the human amygdala to modulate threat learning. MSc thesis.
    - Fat, D. L. et al. The Hounsfield value for cortical bone geometry in the proximal humerus: an in vitro study. Skelet. Radiol. 41, 557–568 (2012).
    - Guo, S., et al. Feasibility of ultrashort echo time images using full-wave acoustic and thermal modeling for transcranial MRI-guided focused ultrasound (tcMRgFUS) planning. *Phys. Med. Biol.* **64**, 095008 (2019).
    - Marsac, L. et al. Ex Vivo Optimisation of a Heterogeneous Speed of Sound Model of the Human Skull for Non-Invasive Transcranial Focused Ultrasound at 1 MHz. International Journal of Hyperthermia. 33. 635–645 (2017).
    - McDannold, N., et al. Transcranial magnetic resonance imaging-guided focused ultrasound surgery of brain tumors: initial findings in 3 patients. *Neurosurgery* **66**(2):323–332 (2010).
    - Miscouridou, M., Pineda-Pardo, J. A., Stagg, C. J., Treeby, B. E. & Stanziola, A. Classical and Learned MR to Pseudo-CT Mappings for Accurate Transcranial Ultrasound Simulation. IEEE Trans. Ultrason., Ferroelectr., Freq. Control 69, 2896–2905 (2022).
    - Mueller, J. K., Ai, L., Bansal, P. & Legon, W. Numerical Evaluation of the Skull for Human Neuromodulation with Transcranial Focused Ultrasound. Journal of Neural Engineering. 14. 066012 (2017).
    - Pichardo, S., et al. Multi-frequency characterization of the speed of sound and attenuation coefficient for longitudinal transmission of freshly excised human skulls. *Phys. Med. Biol.* **56**(1):219–250 (2011). DOI: 10.1088/0031-9155/56/1/014
    - Schneider, U., Pedroni, E., and Lomax A., "The calibration of CT Hounsfield units for radiotherapy treatment planning," Phys. Med. Biol., 41, pp. 111-124 (1996).
    - Webb, T. D., et al. Measurements of the Relationship Between CT Hounsfield Units and Acoustic Velocity and How It Changes With Photon Energy and Reconstruction Kernel. *IEEE Trans Ultrason Ferroelectr Freq Control* **65**(7):1111–1124 (2018). DOI: 10.1109/TUFFC.2018.2827899
    - Webb, T. D., et al. Measurements of the Relationship Between CT Hounsfield Units and Acoustic Attenuation and How It Changes With Photon Energy and Reconstruction Kernel. *IEEE Trans Ultrason Ferroelectr Freq Control* **68**(5):1532–1545 (2020). DOI: 10.1109/TUFFC.2020.3039743
    - Wiesinger, F. et al. Zero TE-based pseudo-CT image conversion in the head and its application in PET/MR attenuation correction and MR-guided radiation therapy planning. Magn. Reson. Med. 80, 1440–1451 (2018).
    - Yaakub, S. N. et al. Pseudo-CTs from T1-Weighted MRI for Planning of Low-Intensity Transcranial Focused Ultrasound Neuromodulation: An Open-Source Tool. Brain Stimulation. 16. 75–78 (2023).
