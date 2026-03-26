# Simulation outputs

PRESTUS provides multiple outputs. Some of these outputs are optional (or can be deactivated upon request to save space.)

---

## Summary table

Filename: `sub-XXX_<simulation_medium>_output_table<affix>.csv`

All quantitative metrics are written to a single CSV file. Acoustic metrics are written first; thermal metrics are appended if heating simulations were run. Tissue-specific columns are only present for layered simulations.

### Intensity nomenclature

All intensity values are **pulse-average intensity (IPA)**, computed from the steady-state peak pressure amplitude as:

```
I = p² / (2ρc)   [W/cm²]
```

where `p` is the peak pressure amplitude, `ρ` is density, and `c` is sound speed. The cycle-averaging factor of 2 reflects that `p` is the amplitude of a sinusoid — this yields the time-average intensity *during* the pulse, i.e. IPA, not the temporal peak (ITP, which would lack the factor of 2) and not the temporal average over the full repetition interval (ISPTA, which would additionally apply the duty cycle PD/PRI).

The **spatial peak** of IPA is ISPPA. When IPA is reported at a specific location (target point, radius average) rather than at the spatial maximum, it is labelled IPA.

The **mechanical index** at each voxel is defined as MI = p⁻ / √f, where p⁻ is the peak negative pressure in MPa and f is the centre frequency in MHz. Note: PRESTUS uses `p_max_all` (peak positive pressure) as a proxy for p⁻; a direct negative-pressure map would require additional sensor configuration.

An overview of the nomenclature is provided in [Darmani et al. (2022, Clinical Neurophysiology)](https://doi.org/10.1016/j.clinph.2021.12.010).

![intensity_nomenclature](https://ars.els-cdn.com/content/image/1-s2.0-S1388245721008920-gr2.jpg)

<br>
<span style="font-size: 12px;">
  <strong>References:</strong><br>
  Darmani, G., Bergmann, T. O., Butts Pauly, K., Caskey, C. F., de Lecea, L., Fomenko, A., Fouragnan, E., Legon, W., Murphy, K. R., Nandi, T., Phipps, M. A., Pinton, G., Ramezanpour, H., Sallet, J., Yaakub, S. N., Yoo, S. S., & Chen, R. (2022). Non‑invasive transcranial ultrasound stimulation for neuromodulation. <em>Clinical Neurophysiology, 133</em>(7), 167–177. https://doi.org/10.1016/j.clinph.2022.03.011
</span>
<br>

### Acoustic metrics

#### Global

| Column | Unit | Description |
|---|---|---|
| `subject_id` | — | Subject identifier |
| `freq_Hz` | Hz | Transducer centre frequency |
| `Isppa` | W/cm² | ISPPA — maximum pulse-average intensity across the entire simulation grid |
| `Isppa_after_exitplane` | W/cm² | Maximum IPA beyond the transducer exit plane (or within brain if modelled) |
| `real_focal_distance_mm` | mm | Distance from the transducer face to the location of maximum IPA |
| `Ipa_target` | W/cm² | IPA at the specified target coordinate |
| `Ipa_target_radius` | W/cm² | Mean IPA within a sphere of radius `analysis.focus_area_radius` mm centred on the target |
| `trans_pos_vox` | voxels | Transducer position in simulation grid coordinates |
| `focus_pos_vox` | voxels | Target focus position in simulation grid coordinates |

Water/free-field simulations additionally include `Psptp` (spatial peak temporal peak pressure) and `Ptp_target` (temporal peak pressure at the target).

#### Tissue-specific *(layered simulations only)*

| Column | Unit | Description |
|---|---|---|
| `Isppa_brain` | W/cm² | Maximum IPA within brain tissue (GM + WM) |
| `Isppa_skull` | W/cm² | Maximum IPA within skull bone |
| `Isppa_skin` | W/cm² | Maximum IPA within skin |
| `Psptp_brain` | Pa | Spatial peak temporal peak pressure within brain tissue |
| `Psptp_skull` | Pa | Spatial peak temporal peak pressure within skull bone |
| `Psptp_skin` | Pa | Spatial peak temporal peak pressure within skin |
| `Ptp_target` | Pa | Temporal peak pressure at the target coordinate |
| `MI_brain` | — | Maximum mechanical index within brain tissue (GM + WM) |
| `MI_skull` | — | Maximum mechanical index within skull bone |
| `MI_skin` | — | Maximum mechanical index within skin |
| `MI_tc` | — | **Transcranial MI** — maximum mechanical index across all intracranial voxels (WM, GM, CSF, blood; skull bone excluded). Primary safety metric for cavitation risk at the sonication target. ITRUSST limit: 1.9 |
| `Ix_brain_vox`, `Iy_brain_vox`, `Iz_brain_vox` | voxels | Grid coordinates of peak IPA within brain |
| `halfmax_ISPPA_volume_brain_mm3` | mm³ | Volume of the −6 dB focal region within brain (voxels where IPA ≥ 50 % of `Isppa_brain`) |

---

### Thermal metrics *(only written when heating simulations are run)*

Temperature is simulated over the full sonication protocol (on/off cycles) using the Pennes bioheat equation. All values are spatial maxima within the respective tissue.

#### Global

| Column | Unit | Description |
|---|---|---|
| `maxT` | °C | Global maximum temperature reached at any point during sonication |
| `endT` | °C | Global maximum temperature at the end of the last sonication pulse |
| `maxCEM43` | min | Global maximum CEM43 over the full protocol |
| `maxCEM43end` | min | Global maximum CEM43 at end of last pulse |

#### Tissue-specific

| Column | Unit | Description |
|---|---|---|
| `maxT_brain/skull/skin` | °C | Maximum temperature in each tissue during sonication. ITRUSST limit: 39 °C |
| `endT_brain/skull/skin` | °C | Maximum temperature in each tissue at end of last pulse |
| `riseT_brain/skull/skin` | °C | Maximum temperature rise above baseline (maxT − baseline). ITRUSST limit: 2 °C |
| `rise_endT_brain/skull/skin` | °C | Temperature rise above baseline at end of last pulse |
| `CEM43_brain/skull/skin` | min | Maximum CEM43 in each tissue. ITRUSST limits: brain 2 min, skull 16 min, skin 21 min |
| `CEM43_end_brain/skull/skin` | min | CEM43 in each tissue at end of last pulse |

---

### Non-significant risk limits

The HTML report flags each metric against the ITRUSST consensus safety limits for transcranial ultrasound (Aubry et al., 2025, Brain Stimulation):

| Metric | Limit |
|---|---|
| MI transcranial | 1.9 |
| Temperature rise (a) | 2 °C |
| Maximum temperature (b) | 39 °C |
| CEM43 (brain) (c) | 2 min |
| CEM43 (skull) (c) | 16 min |
| CEM43 (skin) (c) | 21 min |

For thermal estimates, non-significant risk can be characterized by a, b, or c. NSR limits above and in the HTML are provided for quick reference. Always check the original reference for details and check for possible consensus updates. PRESTUS does not make independent recommendations.

---

## NIfTI images

PRESTUS exports 3D volumes (or 2D for axisymmetric simulations) in NIfTI format, provided in subject-native space (`_orig_coord`) and MNI-152 space (`_MNI`).

Filename pattern: `sub-XXX_final_<type>_MNI<affix>.nii.gz`

| Type | Description |
|---|---|
| `medium_masks` | Tissue label mask (integer indices) |
| `intensity` | Pulse-average intensity (IPA) [W/cm²] |
| `MI` | Mechanical index |
| `pressure` | Peak pressure [Pa] |
| `heating` | Maximum temperature during sonication [°C] |
| `heating_end` | Temperature at end of last pulse [°C] |
| `heatrise` | Temperature rise above baseline [°C] |
| `heatrise_end` | Temperature rise above baseline at end of last pulse [°C] |
| `CEM43` | CEM43 thermal dose [min] |
| `CEM43_end` | CEM43 at end of last pulse [min] |

---

## MATLAB structures

PRESTUS saves simulation parameters and, by default, intermediate matrices. If matrices are detected in the output folder and `io.overwrite_files` is set to `never`, they are loaded instead of recomputing. Matrix saving can be disabled globally via `io.save_matrices = 0`.

| File | Description |
|---|---|
| `sub-XXX_<medium>_parameters<affix>.mat` | Full parameters struct used for the simulation |
| `sub-XXX_after_rotating_and_scaling<affix>.mat` | Head volume after grid scaling |
| `sub-XXX_after_cropping_and_smoothing<affix>.mat` | Cropped head including medium masks |
| `sub-XXX_<medium>_kwave_source<affix>.mat` | k-Wave source definition |
| `sub-XXX_<medium>_results<affix>.mat` | Acoustic simulation sensor data |
| `sub-XXX_<medium>_heating_res<affix>.mat` | Thermal simulation results |

---

## Figures

PRESTUS saves PNG overlays of intensity and temperature on the segmentation for quick visual inspection.

![PRESTUS_fig_example_thermal](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/thermal_fig_examples.png)

---

## HTML report

A self-contained HTML report is written to the output directory (`sub-XXX_<medium>_report<affix>.html`). It includes:

- **Safety dashboard** — color-coded cards for all ITRUSST-limited metrics (green / amber / red)
- **Simulation summary** — key acoustic and thermal values at a glance
- **Configuration summary** — parameters used for the run
- **Medium properties** — per-tissue acoustic and thermal properties (layered only)
- **Acoustic results** — intensity overlay images and the full acoustic metrics table
- **Thermal results** — temperature and CEM43 plots (if heating simulations were run)
- **Post-hoc water simulation** — free-field reference results (if run)
- **Log** — full pipeline execution log

![PRESTUS_html](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/html_report.png)

A full example report on a phantom can be found here: [🔗 Download Example HTML](https://raw.githubusercontent.com/jkosciessa/PRESTUS_bin/main/examples/example_report.html)
