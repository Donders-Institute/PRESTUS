# Calibration Parameters

Calibration fits a library entry to measured hydrophone profiles. These runs read their own config, `config_calibration.yaml`, loaded as `parameters.calibration`; this page lists its keys. How the fits work is described in [Transducer Calibration](doc_calibration.md).

On this page:

- [Transducer Calibration (`calibration`)](#transducer-calibration-calibration): keys of `config_calibration.yaml`.

## Transducer Calibration (`calibration`)

| **Parameter** | **Description** | **Comments** |
|---|---|---|
| `path_input_axial` | Directory containing axial intensity profiles; `char_data_path` in the equipment YAML is relative to it. | Usually `<PRESTUS>/local/profiles/` (see doc_transducer_library.md, Appendix). |
| `path_output` | Directory for saving free-water simulation results. | |
| `path_output_profiles` | Directory for saving optimised profile data. | |
| `filename_calibrated_CSV` | Filename of calibrated CSV data. | Mandatory only when not generated within standalone script. |
| `save_in_calibration_folder` | Save in `path_output` (`TRUE`) or `sim_path` (`FALSE`). | If `TRUE`, results are appended to existing calibration data. |
| `combinations` | Equipment combinations (must refer to combos in the transducer library, `config/equipment/`, or your own library). | Multiple combinations can be specified. |
| `focal_depths_wrt_exit_plane` | List of focal depths to characterise [mm]. | |
| `desired_intensities` | Desired free-water intensities [W/cm²]. | |
| `add_FDO` | Append Focal Distance Offset (bowl-to-exit-plane distance)? | Set to `1` if zero point in profiles reflects exit plane rather than bowl. |
| `axisymmetric2D` | Use axisymmetric 2D instead of default 3D free-water simulations? | `1 = yes`, `0 = no` |
| `force_kwavearray` | Force use of kWaveArray for free-water simulations? | If `0`, uses the setting in the default/study config. |
| `opt_method` | Optimisation method. | [`FEXminimize`](https://github.com/rodyo/FEX-minimize) (open source) / `GlobalSearch` (MATLAB [Global Optimization Toolbox](https://www.mathworks.com/products/global-optimization.html)) |
| `opt_limits` | Distance limits for optimisation [mm]. | |
| `opt_weights` | Weighting of original profile during fitting. | `1` = equal; `>1` = Gaussian (narrower with larger values) |
| `opt_seed` | Random seed for optimisation. | |
| `opt_use_initial_phases` | Use phases from `transducer.annular.elem_phase_deg` as the initial guess instead of random initialization. | `false` (default); set to `true` when manufacturer-provided or previously calibrated phases are available. |
| `opt_upper_velocity` | Upper velocity in global search | |
| `opt_phase_precession` | Constrain element phases during optimisation. | `false` (default) = independent; `'linear'` = strict ramp (`phase_start + i * phase_step`, 2 free params); `'monotonic'` = ordered phases via cumulative non-negative increments (N params, monotonically constrained). |
| `opt_amp_validation` | Run free-water simulation to validate phases and amplitudes?| 'always' (default) / 'never' / 'initial' (only for the first target intensity) / 'final' (only for the final target intensity) | |
| `skip_front_peak_mm` | Distance from profile start to ignore [mm]. | Avoids near-field artefacts in peak/FWHM calculations. |
| `profile_align` | Focus alignment when interpolating measured profiles to an unmeasured focal depth. | `'label'` (default) = column focal-depth headers; `'peak'` = intensity peak beyond `skip_front_peak_mm`. See `interpolate_axial_profiles`. |
| `max_extrapolation_mm` | How far outside the measured depth range a desired depth may lie [mm]. | `0` (default) = no extrapolation (error). Otherwise the nearest measured profile, shifted to the desired depth. |
| `focal_reference` | Convention for what a quoted focal depth means. | ![](assets/icons/icon_focalref-flhm-3db_glyph.svg){ .choice-icon } `flhm_3db` (default, ITRUSST, intensity ≥ peak/2) / ![](assets/icons/icon_focalref-flhm-6db_glyph.svg){ .choice-icon } `flhm_6db` ([Imasonic](https://www.imasonic.com/)/IGT, intensity ≥ peak/4) / ![](assets/icons/icon_focalref-peak_glyph.svg){ .choice-icon } `peak` (raw intensity peak). Recorded into the sidecar YAML's `meta.calibration.focal_reference`. See [Transducer Calibration](doc_calibration.md). |
