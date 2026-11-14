# Skull Property Mappings

![Pipeline steps; this page belongs to 3 Grid & medium](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/fig_stage-strip-medium_diagram.svg){ .stage-strip }

The pseudo-CT holds pseudo-Hounsfield units, while the simulation needs density, sound speed and attenuation in each skull voxel. Several published mappings make this conversion, and they assume different HU scales and give different skull properties. This page lists the mappings PRESTUS offers for each property, their equations and prerequisites, and how to check the resulting skull properties before simulating.

On this page:

- [From Pseudo-HU to Acoustic Properties](#from-pseudo-hu-to-acoustic-properties): how the mappings are applied.
- [Mapping Skull Density](#mapping-skull-density): the `pct.mapping_density` options.
- [Mapping Skull Sound Speed](#mapping-skull-sound-speed): the `pct.mapping_soundspeed` options.
- [Mapping Skull Attenuation](#mapping-skull-attenuation): the `pct.mapping_attenuation` options.
- [Predicting the Modelled Skull Properties](#predicting-the-modelled-skull-properties): the skull property distributions a run will model, and the references.

!!! note "At a glance"
    - **When you need it:** When you change `pct.mapping_density`, `pct.mapping_soundspeed` or `pct.mapping_attenuation` from the `k-plan` default, or want the equations behind a mapping.
    - **Before you start:** A pseudo-CT ([Setting Up a Pseudo-CT Run](doc_pseudoCT.md#setting-up-a-pseudo-ct-run)) and a density calibration that matches its UTE-to-HU mapping ([Matching the CT Calibration](doc_pseudoCT_generation.md#matching-the-ct-calibration)).

## From Pseudo-HU to Acoustic Properties

With `pct.enabled: 1`, choose a mapping for each of density, sound speed and attenuation. The pseudo-CT is read from the pseudo-CT folder (see [Setting Up a Pseudo-CT Run](doc_pseudoCT.md#setting-up-a-pseudo-ct-run)). This is a subtype of a `layered` medium. Mappings apply only within the uniform skull mask. A property set to `"none"`, or with its field missing, keeps the uniform `skull` value. Defaults below are those in `config_default.yaml`.

## Mapping Skull Density

Set by `pct.mapping_density`.

| Value | Mapping |
|---|---|
| `k-plan` (**default**) | 4-part piecewise linear fit based on k-Plan defaults; `pct.density_calibration` replaces the table (see [Matching the CT Calibration](doc_pseudoCT_generation.md#matching-the-ct-calibration)) |
| `k-wave` | 4-part piecewise linear fit based on Schneider et al. (1996) |
| `marsac` | linear between `HU_min = 300` and `HU_max = 2000`, [Marsac et al. (2017)](https://doi.org/10.1080/02656736.2017.1295322) |
| `aubry` | porosity-based mixture of bone and water, [Aubry et al. (2003)](https://doi.org/10.1121/1.1529663) |
| `air-tissue` | linear air-to-bone calibration ([BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain)-compatible); **use with `pichardo` or `mcdannold` sound speed or attenuation** |

??? note "Density equations"
    Several mappings use the porosity of Aubry et al. (2003), with the max. skull HU denominator from [Guo et al. (2019)](https://doi.org/10.1088/1361-6560/ab12f7):

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

## Mapping Skull Sound Speed

Set by `pct.mapping_soundspeed`.

| Value | Mapping |
|---|---|
| `k-plan` (**default**) | density-based, [k-Plan](https://dispatch.k-plan.io/static/docs/simulation-pipeline.html#converting-the-primary-planning-image-to-material-properties) |
| `marsac` | density-based, Marsac et al. (2017) |
| `aubry` | porosity-based, Aubry et al. (2003) |
| `webb` | linear in HU, [Webb et al. (2018)](https://doi.org/10.1109/TUFFC.2018.2827899); scanner-specific slope and intercept via `pct.webb_sos_slope` and `pct.webb_sos_intercept` |
| `pichardo` | frequency-dependent ρ→c, [Pichardo et al.](https://doi.org/10.1088/0031-9155/56/1/014) 2D interpolation map; needs `pct.mapping_density: air-tissue` |
| `mcdannold` | polynomial ρ→c, [McDannold et al.](https://doi.org/10.1227/01.NEU.0000360379.95800.2F); needs density first (`air-tissue` or `marsac`) |

For `k-plan`, `webb`, `pichardo` and `mcdannold`, the minimum sound speed is regularised to water sound speed.

??? note "Sound speed equations"
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

## Mapping Skull Attenuation

Set by `pct.mapping_attenuation`.

| Value | Mapping |
|---|---|
| `k-plan` (**default**) | uniform: the skull `alpha_coeff` and `alpha_power` of `medium_properties`; to replicate k-Plan set `alpha_coeff = 13.3` and `alpha_power = 1` |
| `mueller` | HU-based, [Mueller et al. (2017)](https://doi.org/10.1088/1741-2552/aa843e), bounds from [Yaakub et al. (2023)](https://doi.org/10.1016/j.brs.2023.01.838) |
| `aubry` | porosity-based, Aubry et al. (2003) |
| `webb` | exponential HU-based model, [Webb et al. (2021)](https://doi.org/10.1109/TUFFC.2020.3039743); scanner-specific values via `pct.webb_att_alpha0`, `pct.webb_att_beta` and `pct.webb_att_c` |
| `pichardo` | frequency-dependent ρ→Att from Pichardo et al. 2D interpolation map; prerequisites as for `pichardo` sound speed |
| `mcdannold` | polynomial ρ→Att from McDannold et al. with linear frequency scaling; needs density first |

??? note "Attenuation equations"
    **`k-plan`.**

    $$ \alpha_\text{skull} = \alpha_\text{coeff} $$

    Reads `alpha_coeff` and `alpha_power` of `medium_properties.skull` from the config, and warns when they differ from k-Plan's fixed `13.3` and `1`.

    **`mueller`.**

    $$ \alpha_\text{skull} = \alpha_\text{min} + (\alpha_\text{max} - \alpha_\text{min}) \left(1 - \frac{\text{HU} - \text{HU}_\text{min}}{\text{HU}_\text{max} - \text{HU}_\text{min}}\right)^{0.5} $$

    with `α_min` = 4 and `α_max` = 8.7 (see Yaakub et al., 2023). These are 𝛼(f) at 500 kHz (see [Aubry, J.-F., 2022](https://doi.org/10.1121/10.0013426) for benchmark simulations). PRESTUS converts them to `alpha_0` in `𝛼(f) = alpha_0 x f[MHz] ^ y` as `alpha_0 = 𝛼(f)/0.5^y`, with `y` the skull `alpha_power`.

    **`aubry`.**

    $$ \alpha_\text{skull} = \alpha_\text{min} + (\alpha_\text{max} - \alpha_\text{min}) \cdot \phi_\text{skull}^{0.5} $$

    with `α_min = 0.2` and `α_max = 8` (Aubry et al., 2003).

    **`webb`.**

    $$ \alpha_\text{skull}(f) = \alpha_0 \cdot \left(\frac{f}{1\,\text{MHz}}\right)^{\beta} \cdot \exp(c \cdot \text{HU}) $$

    Output is [k-Wave](https://github.com/ucl-bug/k-wave) `alpha_coeff` [dB/cm/MHz$^\beta$] with `alpha_power = β`. Defaults (GE 120 kVp, B kernel, 0.5 × 0.6 mm): α₀ = 26, β = 1.3, c = −0.0016. Set scanner-specific values from Webb et al. Table IV with `pct.webb_att_alpha0`, `pct.webb_att_beta` and `pct.webb_att_c`.

    **`pichardo`.** Reads longitudinal attenuation [Np/m] from `MapPichardo.h5` as a function of density [kg/m³] and frequency [MHz]. The value at the simulation frequency is converted to k-Wave `alpha_coeff` [dB/cm/MHz] with `alpha_power = 1` (linear frequency scaling assumed). Prerequisites as for `pichardo` sound speed; density must be filled first.

    **`mcdannold`.**

    $$ \alpha_\text{skull}(f) = p(\rho_\text{skull}) \cdot \frac{f}{660\,\text{kHz}} $$

    where $p(\rho) = 1.36 \times 10^{-10}\,\rho^4 - 1.41 \times 10^{-6}\,\rho^3 + 5.40 \times 10^{-3}\,\rho^2 - 9.02\,\rho + 5.71 \times 10^{3}$ [Np/m]. Converted to k-Wave `alpha_coeff` [dB/cm/MHz] with `alpha_power = 1`. Requires density first.

## Predicting the Modelled Skull Properties

To see the skull properties a run will model before simulating, call `pct_predict_skull_properties(parameters)` or use **Predict skull properties** on the Nexus Pseudo-CT sub-tab. It applies the configured mappings to the pseudo-CT skull voxels and writes their distributions to `cache/sub-NNN_<medium>_pct_prediction<affix>.json`, alongside the skull values of the uncertainty bound variants. See [Pseudo-CT runs](doc_uncertainty.md#pseudo-ct-runs).

??? note "References"
    - Aubry, J.-F., et al. Experimental demonstration of noninvasive transskull adaptive focusing based on prior computed tomography scans. *J Acoust Soc Am* **113**, 84–93 (2003).
    - Aubry, J.-F. et al. Benchmark problems for transcranial ultrasound simulation: Intercomparison of compressional wave models. J. Acoust. Soc. Am. 152, 1003–1019 (2022).
    - Carpino et al. (2024). Transcranial ultrasonic stimulation of the human amygdala to modulate threat learning. MSc thesis.
    - Guo, S., et al. Feasibility of ultrashort echo time images using full-wave acoustic and thermal modeling for transcranial MRI-guided focused ultrasound (tcMRgFUS) planning. *Phys. Med. Biol.* **64**, 095008 (2019).
    - Lim Fat, D. et al. The Hounsfield value for cortical bone geometry in the proximal humerus: an in vitro study. Skelet. Radiol. 41, 557–568 (2012). https://doi.org/10.1007/s00256-011-1255-7
    - Marsac, L. et al. Ex Vivo Optimisation of a Heterogeneous Speed of Sound Model of the Human Skull for Non-Invasive Transcranial Focused Ultrasound at 1 MHz. International Journal of Hyperthermia. 33. 635–645 (2017).
    - McDannold, N., et al. Transcranial magnetic resonance imaging-guided focused ultrasound surgery of brain tumors: initial findings in 3 patients. *Neurosurgery* **66**(2):323–332 (2010).
    - Miscouridou, M., Pineda-Pardo, J. A., Stagg, C. J., Treeby, B. E. & Stanziola, A. Classical and Learned MR to Pseudo-CT Mappings for Accurate Transcranial Ultrasound Simulation. IEEE Trans. Ultrason., Ferroelectr., Freq. Control 69, 2896–2905 (2022).
    - Mueller, J. K., Ai, L., Bansal, P. & Legon, W. Numerical Evaluation of the Skull for Human Neuromodulation with Transcranial Focused Ultrasound. Journal of Neural Engineering. 14. 066012 (2017).
    - Pichardo, S., et al. Multi-frequency characterization of the speed of sound and attenuation coefficient for longitudinal transmission of freshly excised human skulls. *Phys. Med. Biol.* **56**(1):219–250 (2011). DOI: 10.1088/0031-9155/56/1/014
    - Schneider, U., Pedroni, E., and Lomax A., "The calibration of CT Hounsfield units for radiotherapy treatment planning," Phys. Med. Biol., 41, pp. 111-124 (1996).
    - Webb, T. D., et al. Measurements of the Relationship Between CT Hounsfield Units and Acoustic Velocity and How It Changes With Photon Energy and Reconstruction Method. *IEEE Trans Ultrason Ferroelectr Freq Control* **65**(7):1111–1124 (2018). DOI: 10.1109/TUFFC.2018.2827899
    - Webb, T. D., et al. Acoustic Attenuation: Multifrequency Measurement and Relationship to CT and MR Imaging. *IEEE Trans Ultrason Ferroelectr Freq Control* **68**(5):1532–1545 (2021). DOI: 10.1109/TUFFC.2020.3039743
    - Wiesinger, F. et al. Zero TE-based pseudo-CT image conversion in the head and its application in PET/MR attenuation correction and MR-guided radiation therapy planning. Magn. Reson. Med. 80, 1440–1451 (2018).
    - Yaakub, S. N. et al. Pseudo-CTs from T1-Weighted MRI for Planning of Low-Intensity Transcranial Focused Ultrasound Neuromodulation: An Open-Source Tool. Brain Stimulation. 16. 75–78 (2023).
