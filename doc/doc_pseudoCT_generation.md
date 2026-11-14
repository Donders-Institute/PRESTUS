# From UTE to Pseudo-HU

![Pipeline steps; this page belongs to 3 Grid & medium](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/fig_stage-strip-medium_diagram.svg){ .stage-strip }

The pseudo-CT is only as good as the line that turns UTE intensity into pseudo-Hounsfield units (pHU), and that line returns HU on the scale of the CT it was fitted to. This page describes the available UTE-to-HU mappings, how PRESTUS generates the pseudo-CT, and how to pair each mapping with the HU-to-density calibration of the same CT.

On this page:

- [Choosing the UTE-to-HU Mapping](#choosing-the-ute-to-hu-mapping): the `pct.skull_mapping` options and the generation steps.
- [Matching the CT Calibration](#matching-the-ct-calibration): pairing `pct.skull_mapping` with `pct.density_calibration`.

!!! note "At a glance"
    - **When you need it:** When you change `pct.skull_mapping` from the default, use your own UTE-to-HU line, or set `pct.density_calibration`.
    - **Before you start:** A pseudo-CT run set up as in [Setting Up a Pseudo-CT Run](doc_pseudoCT.md#setting-up-a-pseudo-ct-run).

## Choosing the UTE-to-HU Mapping

`pct.skull_mapping` sets the linear UTE→pHU mapping of the skull. Soft tissue is set to 42 HU and air to -1000 HU for every mapping.

| `skull_mapping` | pHU |
|---|---|
| `kosciessa` (**default**) | individualised: a line fitted per subject through three UTE histogram points set to 1900, 700 and 300 HU (see below) |
| `miscouridou` | pHU = −2085 UTE + 2329 |
| `carpino` (deprecated) | pHU = -2194 UTE + 2236. Warns `prestus:pct:deprecatedSkullMapping`; the same line is `custom` with `custom_slope: -2194`, `custom_intercept: 2236` |
| `wiesinger` | pHU = -2000 (UTE-1) + 42 |
| `treeby` | pHU = -2929.6 UTE + 3274.9 ([ucl-bug/petra-to-ct](https://github.com/ucl-bug/petra-to-ct); fitted to 80 kVp CT: use with `density_calibration: ucl_lowdose_80kvp`, see [below](#matching-the-ct-calibration)) |
| `custom` | pHU = `pct.custom_slope` UTE + `pct.custom_intercept` (any other linear map) |

??? note "The generation steps"
    - UTE: threshold at 0, bias field correction, then normalisation by the soft tissue peak (normalised UTE intensity = 1). `pct_soft_tissue_peak` finds the peak.
    - pHU soft tissue (normalised UTE intensity = 1): 42 HU ([Wiesinger et al., 2018](https://doi.org/10.1002/mrm.27134)).
    - pHU air: -1000 HU (Wiesinger et al., 2018; [Miscouridou et al., 2022](https://doi.org/10.1109/TUFFC.2022.3198522)).
    - pHU skull: the linear mapping set by `pct.skull_mapping` (table above).
    - `kosciessa` (`pct_skullmapping_fit`) fits a least-squares line through three points from the UTE histograms of charm cortical and trabecular bone: the lower quarter-max crossing of the cortical histogram → 1900 HU, the crossover of the two histograms → 700 HU, and the upper quarter-max crossing of the trabecular histogram around its x-weighted peak → 300 HU. Bins 0.9–1.1 are replaced by a spline against soft-tissue partial volume. The anchors, width level and interpolation window are options of `pct_skullmapping_fit`.
    - 3D Gaussian smoothing (σ=0.8 voxels via ANTs' SmoothImage) prevents abrupt tissue transitions in the simulation grid. Within an eroded skull mask, unsmoothed values are kept to limit partial volume effects.

## Matching the CT Calibration

A UTE→HU line returns HU on the scale of the CT it was fitted to. That scale depends on scanner, tube voltage (kVp) and reconstruction kernel: at 80 kVp, bone has about 1.35 times the HU it has at about 120 kVp. The simulation needs density, and turning HU into density needs the calibration of that same CT. **Set `pct.skull_mapping`, `pct.mapping_density: k-plan` and `pct.density_calibration` together:**

| `skull_mapping` | CT the line was fitted to | `mapping_density` | `density_calibration` |
|---|---|---|---|
| `kosciessa`, `miscouridou`, `carpino` (deprecated), `wiesinger` | clinical CT, about 120 kVp (assumed; for `kosciessa` the anchors 1900/700/300 HU are assumed on this scale) | any | `''` ([k-Plan](https://k-plan.io/) default) |
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

??? note "HU scale assumed by the other mappings"
    - `k-wave`: [Schneider et al.](https://doi.org/10.1088/0031-9155/41/1/009); `marsac`: 300–2000 HU; `air-tissue` (needed by `pichardo` and `mcdannold`): about −1000 to 1000 HU.
    - `webb` sound speed and attenuation: GE 120 kVp, B kernel by default; scanner-specific values via `pct.webb_*`.
    - `aubry` and `mueller` rescale by the skull's own minimum and maximum HU. They are insensitive to the HU scale but remove between-subject differences in density level.
