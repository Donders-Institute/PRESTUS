# Matrix Arrays: Steering Correction

The focal peak a matrix array delivers changes with the steering target in ways the geometric model does not predict. The steering correction scales the drive by a factor measured at a set of targets, so that a simulated condition matches what the device delivers. Only users of matrix (2D) arrays need this page.

On this page:

- [Steering Correction of Matrix Arrays](#steering-correction-of-matrix-arrays): how the correction is applied, building an entry from hydrophone peaks, and its limits.

!!! note "At a glance"
    - **When you need it:** You simulate a matrix (2D) array steered to targets away from its geometric focus.
    - **Before you start:** Hydrophone focal-peak measurements at 30–40 steering targets at one drive voltage, and the array's equipment YAML.

## Steering Correction of Matrix Arrays

A matrix (2D) array's delivered focal peak can differ from the geometric model by a factor that changes with the steering target (element directivity, crosstalk, manufacturing spread). The steering correction measures that factor at a set of targets and scales the drive by it.

**How it is applied.** `create_matrix_karray` multiplies every element amplitude by k(target) and records it in `transducer.matrix.steering_k`. `elem_p_amp` stays the nominal drive setting, and the simulated field is what the device delivers at that setting. With [ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting) the results are reported at the target ISPPA either way; the correction changes the drive that achieves them: `elem_p_amp × isppa_scale_p`, or in volts `drive_ref_v × elem_p_amp × isppa_scale_p / elem_p_amp_ref`. With `target_isppa_wcm2: []` the results are reported at the corrected drive. Switch it off with `transducer.matrix.use_steering_correction: false`.

**Building an entry.** Measure the focal peak at 30–40 targets spread through the steering range you will use, including its edges and the geometric focus, plus a few kept back to check the fit. Use one drive voltage and the geometric phases PRESTUS computes for each target. Write one CSV row per target, then build the entry with `build_steering_correction`.

??? note "CSV columns and build call"
    | Column | Unit | Required | Meaning |
    |---|---|---|---|
    | `target_x_mm`, `target_y_mm`, `target_z_mm` | mm | yes | Steering target in the transducer frame: origin at the apex (curved) or array centre (flat), z along the beam |
    | `peak_pressure_mpa` | MPa | yes | Largest pressure amplitude in the plane through the target, normal to the beam axis |
    | `drive_v` | V | no | Drive voltage; rows at other voltages are scaled to the median voltage |
    | `peak_x_mm`, `peak_y_mm` | mm | no | Where the peak was found |
    | `model_pressure_mpa` | MPa | no | A modelled peak at the same drive, used instead of the Rayleigh model |

    Then, with the array's geometry in `tr` (from its equipment YAML) and `tr.matrix.elem_p_amp` the element amplitude that stands for the measurement drive:

    ```matlab
    % run from the PRESTUS folder; writes the library entry next to the shipped ones
    sc = build_steering_correction(parameters, tr, 'steering_peaks.csv', ...
        struct('out_file', fullfile('config', 'calibration', '<serial>.yaml'), 'serial', '<serial>'));
    ```

`resolve_transducer_from_serial` attaches the entry to a matrix transducer given by `serial` (the entry for `<serial>` with `combo.ds_serial` first, then the driving-system-agnostic one, or `transducer.library_path`; see [Library keys](doc_transducer_library.md#using-a-library-entry)). Without an entry the array is simulated uncorrected. The entry stores the fit, its cross-validated errors, `elem_p_amp_ref` and `drive_ref_v`. An inline `transducer.matrix.steering_correction` takes precedence.

**Limits.** Only the focal peak is corrected, not side lobes or field shape. The model assumes the beam along +z of the transducer frame and does not support Clover setups. The Rayleigh model treats each element as a baffled piston in free water; check it against k-Wave with a free-water run at a few targets.

??? note "What is fitted"
    At each measured target x, `k(x) = measured / modelled focal peak`, both at the same drive. The modelled peak is the largest free-field pressure in the plane through the target, normal to the beam axis, from a Rayleigh integral over the array (`matrix_focal_peak`, `matrix_rayleigh_pressure`). The measured peak is the same quantity from a hydrophone scan of that plane. `fit_steering_correction` fits k as one constant (`scale`) or a full quadratic surface in x, y, z (`quadratic`, 10 coefficients, linear least squares), and scores both by leave-one-out cross-validation. The default `auto` keeps the quadratic surface only with at least 15 targets and a median cross-validated error at most 0.8 times that of the constant. A quadratic map is not extrapolated: targets outside the measured range use the value at its edge, with a warning.

    [Dadgar-Kiani et al. (2026)](https://doi.org/10.1016/j.crmeth.2026.101568) fit such a correction with a Gaussian of the difference in MPa; PRESTUS fits the ratio, which applies to the drive at any level.

??? note "References"
    Dadgar-Kiani et al. (2026). *Calibrated simulations for dynamic focusing of ultrasound through the temporal window.* Cell Reports Methods 6:101568.
