# Metrics Table Columns

The metrics CSV holds every number a run reports: intensities, pressures, mechanical indices, focal volumes and, with heating, temperatures and thermal doses. Each column has a precise definition and unit, and several are checked against the ITRUSST non-significant risk limits, so the column you quote decides what your reported value means.

On this page:

- [Metrics table](#metrics-table): the file and where targeting accuracy is written.
    - [Intensity nomenclature](#intensity-nomenclature): how IPA, ISPPA and the mechanical index are computed.
    - [Acoustic metrics](#acoustic-metrics): global and tissue-specific acoustic columns.
    - [Thermal metrics](#thermal-metrics): global and tissue-specific thermal columns.
    - [Non-significant risk limits](#non-significant-risk-limits): the reference levels the report flags against.

!!! note "At a glance"
    - **When you need it:** To report a value from the CSV, or to check what a column measures.
    - **Before you start:** A finished run; see [Simulation Outputs](doc_outputs.md) for where the CSV is written.

---

## Metrics Table

Filename: `sub-XXX_<simulation_medium><affix>.csv`

One CSV holds all quantitative metrics: acoustic metrics first, then thermal metrics if heating was simulated. Tissue-specific columns exist for layered simulations only.

Targeting accuracy (`modules.run_targeting_analysis`, layered or phantom media) goes to `sub-XXX_<simulation_medium><affix>_targeting.csv`; see [Targeting accuracy](doc_group.md#targeting-accuracy).

### Intensity Nomenclature

All intensity values are **pulse-average intensity (IPA)**, computed from the steady-state peak pressure amplitude:

```
I = p² / (2ρc)   [W/cm²]
```

`p` is the peak pressure amplitude, `ρ` density and `c` sound speed. The factor 2 averages over the sinusoid's cycle. The result is the temporal-peak intensity (ITP). IPA, the time-average intensity *during* the pulse, is ITP times the ramp factor of `timing.pulse_ramp_dur` and `timing.pulse_ramp_shape`; for pulses without a ramp IPA = ITP. The average over the repetition interval (ISPTA) also applies the duty cycle PD/PRI.

The **spatial peak** of IPA is ISPPA. IPA at a specific location (target point, radius average) is labelled IPA.

The **mechanical index** per voxel is MI = p⁻ / √f, with p⁻ the peak rarefaction (negative) pressure in MPa and f the centre frequency in MHz (IEC 62359 / FDA). PRESTUS uses the absolute value of `p_min_all` from [k-Wave](https://github.com/ucl-bug/k-wave). For multi-transducer async combinations it uses the per-voxel most-negative pressure across transducers (`p_min_async`).

[Darmani et al. (2022, Clinical Neurophysiology)](https://doi.org/10.1016/j.clinph.2021.12.010) give an overview of the nomenclature.

![intensity_nomenclature](https://ars.els-cdn.com/content/image/1-s2.0-S1388245721008920-gr2.jpg)

??? note "References"
    Darmani, G., Bergmann, T. O., Butts Pauly, K., Caskey, C. F., de Lecea, L., Fomenko, A., Fouragnan, E., Legon, W., Murphy, K. R., Nandi, T., Phipps, M. A., Pinton, G., Ramezanpour, H., Sallet, J., Yaakub, S. N., Yoo, S. S., & Chen, R. (2022). Non‑invasive transcranial ultrasound stimulation for neuromodulation. *Clinical Neurophysiology, 135*, 51–73. https://doi.org/10.1016/j.clinph.2021.12.010

### Acoustic Metrics

Focal distances come in two references: columns ending in `_ep_mm` count from the exit plane, as `focal_distance_ep` does; columns ending in `_bowl_mm` from the bowl apex (`trans_pos`). Group and comparison reports read the earlier `real_focal_distance_mm` as `real_focal_distance_bowl_mm`. See [Transducer-Target Distance](doc_transducer.md#transducer-target-distance).

??? note "Global acoustic columns"
    | Column | Unit | Description |
    |---|---|---|
    | `subject_id` | none | Subject identifier |
    | `freq_Hz` | Hz | Transducer centre frequency |
    | `Isppa` | W/cm² | ISPPA: maximum pulse-average intensity across the entire simulation grid |
    | `Isppa_after_exitplane` | W/cm² | Maximum IPA beyond the transducer exit plane (or within brain if modelled) |
    | `real_focal_distance_bowl_mm` | mm | Distance from the bowl apex (`trans_pos`) to the location of maximum IPA |
    | `real_focal_distance_ep_mm` | mm | Distance from the exit plane to the location of maximum IPA; compare with `focal_distance_ep` |
    | `Ipa_target` | W/cm² | IPA at the specified target coordinate |
    | `Ipa_target_radius` | W/cm² | Mean IPA within a sphere of radius `analysis.focus_area_radius` mm centred on the target |
    | `trans_pos_vox` | voxels | Transducer position in simulation grid coordinates |
    | `focus_pos_vox` | voxels | Target focus position in simulation grid coordinates |

    Water/free-field simulations also include `Psptp` (spatial peak temporal peak pressure), `Ptp_target` (temporal peak pressure at the target) and `MI` (global mechanical index from peak rarefaction pressure, unitless; [ITRUSST](https://itrusst.org/) limit 1.9, labelled "MI (free water)" in the report). They also include the ITRUSST focal metrics (`focal_vol_3dB_mm3`, `focal_vol_6dB_mm3`, `focal_centroid_3dB_bowl_mm`, `focal_centroid_3dB_ep_mm`, `focal_flhm_center_bowl_mm`, `focal_flhm_center_ep_mm`, `focal_fwhm_axial_mm`, `focal_fwhm_lateral_min_mm`, `focal_fwhm_lateral_max_mm`), computed beyond the transducer exit plane.

??? note "Tissue-specific acoustic columns (layered simulations only)"
    | Column | Unit | Description |
    |---|---|---|
    | `Isppa_brain` | W/cm² | Maximum IPA within brain tissue (GM + WM) |
    | `Isppa_skull` | W/cm² | Maximum IPA within skull bone |
    | `Isppa_skin` | W/cm² | Maximum IPA within skin |
    | `Psptp_brain` | Pa | Spatial peak temporal peak pressure within brain tissue |
    | `Psptp_skull` | Pa | Spatial peak temporal peak pressure within skull bone |
    | `Psptp_skin` | Pa | Spatial peak temporal peak pressure within skin |
    | `Psptp` | Pa | Global spatial peak temporal peak pressure across the entire grid (also written for water/free-field runs). Max pressure safety limit: 2 MPa |
    | `Ptp_target` | Pa | Temporal peak pressure at the target coordinate |
    | `MI_brain` | none | Maximum mechanical index within brain tissue (GM + WM) |
    | `MI_skull` | none | Maximum mechanical index within skull bone |
    | `MI_skin` | none | Maximum mechanical index within skin |
    | `MI_tc` | none | **Transcranial MI**: maximum mechanical index across all intracranial voxels (WM, GM, CSF, blood; skull bone excluded). Primary safety metric for cavitation risk at the sonication target. ITRUSST limit: 1.9 |
    | `Ix_brain_vox`, `Iy_brain_vox`, `Iz_brain_vox` | voxels | Grid coordinates of peak IPA within brain |
    | `halfmax_ISPPA_volume_brain_mm3` | mm³ | Volume of the −3 dB intensity region within brain (IPA ≥ 50 % of `Isppa_brain`); equals `focal_vol_3dB_mm3` |
    | `focal_vol_3dB_mm3` | mm³ | Volume of the −3 dB focal zone within brain (IPA ≥ 50 % of peak; 10·log₁₀(0.5) = −3 dB) |
    | `focal_vol_6dB_mm3` | mm³ | Volume of the −6 dB focal zone within brain (IPA ≥ 25 % of peak = pressure ≥ 50 % of peak; FWHM) |
    | `focal_centroid_3dB_bowl_mm` | mm | Distance from the bowl apex to the geometric centroid of the −3 dB focal region |
    | `focal_centroid_3dB_ep_mm` | mm | Distance from the exit plane to the geometric centroid of the −3 dB focal region |
    | `focal_flhm_center_bowl_mm` | mm | Distance from the bowl apex to the FLHM centre on the beam axis (`calibration.focal_reference`) |
    | `focal_flhm_center_ep_mm` | mm | Distance from the exit plane to the FLHM centre; the manufacturer-style focal distance |
    | `focal_fwhm_axial_mm` | mm | −3 dB focal width (FWHM of IPA) through the brain peak along the transducer axis (`trans_pos` to `focus_pos`); crossings of half the peak interpolated between grid points. Empty when the focal lobe reaches the grid edge |
    | `focal_fwhm_lateral_min_mm`, `focal_fwhm_lateral_max_mm` | mm | Narrowest and widest −3 dB focal width through the brain peak in the plane perpendicular to the transducer axis, measured every 5°; equal for a round focus |

### Thermal Metrics

Thermal columns are only written when heating simulations are run. Temperature is simulated over the full sonication protocol (on/off cycles) with the [Pennes bioheat equation](https://doi.org/10.1152/jappl.1948.1.2.93). All values are spatial maxima within the respective tissue.

??? note "Global thermal columns"
    | Column | Unit | Description |
    |---|---|---|
    | `maxT` | °C | Global maximum temperature reached at any point during sonication |
    | `endT` | °C | Global maximum temperature at the end of the last sonication pulse |
    | `maxCEM43` | min | Global maximum CEM43 (kWave form) over the full protocol |
    | `maxCEM43end` | min | Global maximum CEM43 (kWave form) at end of last pulse |
    | `maxCEM43iso` | min | Global maximum CEM43 (ISO form) over the full protocol |
    | `maxCEM43isoend` | min | Global maximum CEM43 (ISO form) at end of last pulse |

??? note "Tissue-specific thermal columns"
    | Column | Unit | Description |
    |---|---|---|
    | `maxT_brain/skull/skin` | °C | Maximum temperature in each tissue during sonication. ITRUSST limit: 39 °C |
    | `endT_brain/skull/skin` | °C | Maximum temperature in each tissue at end of last pulse |
    | `riseT_brain/skull/skin` | °C | Maximum temperature rise above baseline (maxT − baseline). ITRUSST limit: 2 °C |
    | `rise_endT_brain/skull/skin` | °C | Temperature rise above baseline at end of last pulse |
    | `CEM43_brain/skull/skin` | min | Maximum CEM43 (kWave form) in each tissue. ITRUSST limits: brain 2 min, skull 16 min, skin 21 min |
    | `CEM43_end_brain/skull/skin` | min | CEM43 (kWave form) in each tissue at end of last pulse |
    | `CEM43iso_brain/skull/skin` | min | Maximum CEM43 (ISO form) in each tissue |
    | `CEM43iso_end_brain/skull/skin` | min | CEM43 (ISO form) in each tissue at end of last pulse |

Both CEM43 formulations are always computed; see [Thermal simulations](doc_simulations-thermal.md) for the formulas. The HTML report shows the one selected by `parameters.thermal.cem43_iso`.

### Non-Significant Risk Limits

The HTML report flags each metric against the ITRUSST consensus non-significant risk limits for transcranial ultrasound ([Aubry et al., 2025](https://doi.org/10.1016/j.brs.2025.10.007), Brain Stimulation):

| Metric | Limit |
|---|---|
| MI transcranial | 1.9 |
| MI (free water) | 1.9 |
| Max pressure | 2 MPa |
| Temperature rise (a) | 2 °C |
| Maximum temperature (b) | 39 °C |
| CEM43 (brain) (c) | 2 min |
| CEM43 (skull) (c) | 16 min |
| CEM43 (skin) (c) | 21 min |

For thermal estimates, non-significant risk can be characterised by a, b, or c. These limits are for quick reference. Always check the original reference and possible consensus updates. PRESTUS does not make independent recommendations.
