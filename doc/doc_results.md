# Reading Your Results

A finished run leaves a folder of files. This page tells you which ones to open first, what to check, and how to read the safety summary. Every file and column is described in [Outputs & Reports](doc_outputs.md); terms are defined in the [Glossary](doc_glossary.md).

## Where to look

Results are in `<path.sim>/sub-XXX/`. Open these first:

| File | What it is |
|---|---|
| `sub-XXX_<medium>_report<affix>.html` | The simulation report: safety summary, key numbers, figures and the run log. Start here. |
| `sub-XXX_placement_<mode><affix>.html` | The placement report: where the transducer and the target ended up, drawn on the T1 (head simulations). |
| `sub-XXX_<medium><affix>.csv` | All metrics in one table, for your own analysis. |
| `nii/` | 3D maps (intensity, pressure, MI, temperature) in the participant's T1 space and in MNI space, to view in any NIfTI viewer. |

`<medium>` is `water` or `layered`, and `<affix>` is your `io.output_affix` (empty by default).

## Check, in this order

### 1. Did it simulate what you meant?

Open the placement report and the report's figures:

- **The transducer sits on the scalp**, not inside the head or far from it, and points at the target.
- **The focus lies at the target.** The report's intensity overlay should show the hot spot in the target region. In the table, compare `real_focal_distance_mm` (where the peak actually is) with your focal depth setting, and `Ipa_target` (intensity at the target) with `Isppa_brain` (the peak in the brain). A peak far from the target means the placement or the focal depth setting needs a second look.
- **The configuration summary** in the report lists the transducer, focal depth, target intensity and medium that were used. Check that they are the ones you intended.

For many runs, the targeting analysis (`sub-XXX_<medium><affix>_targeting.csv`) measures the distance between the simulated and the planned focus ([Targeting accuracy](doc_group.md#targeting-accuracy)).

### 2. Is the intensity what you expect?

- `transducer.target_isppa_wcm2` is the **free-water** intensity (30 W/cm² by default). The intensity **in the brain** (`Isppa_brain`, `Ipa_target`) is lower, mostly because the skull reflects and absorbs ultrasound. How much lower depends on the person, the skull thickness and the angle of incidence.
- If the run log shows `prestus:isppa:largeScale`, the field was scaled by more than 4× (or less than 0.25×) from the simulated drive. The result is then an extrapolation; check that it is still realistic for your device.

### 3. Is it within the safety reference levels?

The top of the report shows a **safety verdict** and one card per metric. Each card is coloured by how close the value is to its ITRUSST non-significant-risk (NSR) reference level (Aubry et al., 2025):

| Colour | Meaning |
|---|---|
| Green | below 50 % of the reference level |
| Amber | 50–100 % of the reference level |
| Red | above the reference level |
| Grey | not computed for this run (e.g. no thermal simulation, or a metric that does not apply to this medium) |

| Metric | Reference level |
|---|---|
| Mechanical index, transcranial (`MI_tc`) or in free water (`MI`) | 1.9 |
| Peak pressure (`Psptp`) | 2 MPa |
| Temperature rise (`riseT_*`) | 2 °C |
| Maximum temperature (`maxT_*`) | 39 °C |
| Thermal dose (`CEM43_*`) | brain 2, skull 16, skin 21 min |

A red card does **not** by itself mean the protocol is unsafe, and the verdict then says *inspect* rather than *fail*. The ITRUSST thermal criteria are alternatives: an exposure qualifies as non-significant risk through the temperature rise, the maximum temperature **or** the thermal dose. Check the flagged values against the criteria that apply to your protocol, and read the original consensus for details. PRESTUS reports the comparison; it makes no safety recommendation of its own.

Thermal metrics only appear when a thermal simulation ran (`modules.run_heating_sims: 1`, with a sonication protocol in `timing.*`). The report shows the CEM43 form selected by `thermal.cem43_iso` ([Thermal](doc_simulations-thermal.md#thermal-dose-cem43-calculation)).

### 4. Anything unusual?

| What you see | What it may mean |
|---|---|
| The highest intensity is in the skin or skull rather than the brain | Expected to some degree near the transducer; a hot spot at the scalp can also mean poor contact or a transducer placed into the head. Check the placement report. |
| The focus is much shallower or deeper than the setting | Strong skull refraction, or a target depth outside the transducer's range (`min_foc`–`max_foc` of the device). |
| Warnings about points per wavelength (PPW) in the log | The grid is too coarse for the frequency; results near the focus may be inaccurate ([Spatial and temporal sampling](doc_simulations-acoustic.md#spatial-and-temporal-sampling)). |
| `load_transducer_from_library:uncalibratedAmplitude` in the log | The transducer's drive at this depth was not checked in k-Wave; absolute intensities carry that uncertainty ([Transducer Library](doc_calibration.md#optional-free-water-k-wave-drive-calibration)). |

More known messages and their fixes are in [Troubleshooting](doc_troubleshooting.md).

## Uncertainty of the estimate

Simulated intensities and temperatures depend on assumed tissue properties, which are uncertain, especially for the skull. To see how much your results could change, run the [uncertainty mode](doc_uncertainty.md) (default tissue properties and the ITRUSST bounds that maximise and minimise heating and intensity), or use subject-specific skull properties from a [pseudo-CT](doc_pseudoCT.md).

## Methods text

The report contains an automatically written methods paragraph (with references) describing the software, transducer, grid, tissue properties, protocol and resulting exposure. Copy it with *Copy methods text*, and check it against your study before using it in a manuscript.
