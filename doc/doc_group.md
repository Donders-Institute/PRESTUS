# Group analysis

To facilitate group analysis, PRESTUS automatically outputs key results as Nifti images in MNI space. Projection to MNI space can be done either using SimNIBS’ subject2mni command, or in MATLAB by reading in MNI2conform_12DOF.txt. 

PRESTUS provides two complementary group-level workflows:

- **MNI-space group plots** — `create_group_MNI_plots` (via the `prestus_group_start` entry point) render per-metric maps in MNI space with a shared colour scale and optional ROI overlays. Documented immediately below.
- **Group HTML report** — `prestus_group_report_start` aggregates the per-subject CSV/PNG outputs into a single self-contained, interactive HTML file. See [Group HTML report](#group-html-report).

Neither workflow runs simulations; both read existing per-subject outputs.

## MNI-space group plots

### create_group_MNI_plots

#### Rationale
- This will ensure that all figures will use the same scale
- Allows one to add a ROI MNI mask to visualise targeting accuracy
- Adds additional columns to each csv with information about intensity in the ROI mask, percentage of fwhm voxels within ROI etc.
- When using structural MRI's, it also allows one to normalise the brightness so the discrepancy in average brightness between individual T1's is reduced.

#### Necessary input:
- subject_list: list of all subjects that you want to include in the figures.
- parameters: load the '.yaml' file using 'load_parameters.m' just as you would do in a pipeline.
- options: a structure similar to 'parameters'. Most of these are (as the name would suggest) optional. only one of the two following parameters needs to be selected:
	- options.slice_to_plot = 0 (give the number of the slice)
	- options.plot_max_intensity = 0 (turn on by changing to 1)

#### Plotting a ROI mask
- Use a mask saved in MNI space
- Load it into the function under the option 'ROI_MNI_mask'

##### Optional parameters:
    options.ROI_MNI_mask (:,:,:) # Loads in a 3d matrix of your ROI mask (in MNI space)
    options.slice_label = 'y' # Selects the axis along which your slice is made
    options.rotation = 90;  # Rotates your structural background figure
    options.plot_heating = 1 # Binary option to enable or disable heating plots
    options.outputs_suffix = '' # Allows you to add a string at the end of the name of your output to say, differentiate between plots of maximum intensity and a given slice
    options.isppa_thresholds = [] # Manually set the Ipa range that is to be plotted
    options.add_FWHM_boundary = 0 # Binary option to add a dotted boundary around the FWHM of intensity
    options.add_ROI_boundary = 1 # Binary option to add a red boundary around your ROI mask
    options.skip_missing = 0 # Binary option to skip subjects with missing files
    options.brightness_correction = 0 # Binary option to normalise the brightness between structural T1's
    options.average_target_brightness = 100 # Average targetted brightness, only used when brightness correction is enabled

## Group HTML report

`prestus_group_report_start` aggregates the per-subject outputs (CSV tables and PNG figures) of every subject sharing the same `simulation.medium` and `io.output_affix` into one **self-contained, interactive HTML report**. It runs no simulations — it only reads outputs that already exist on disk.

The report is written to:

```
<path.sim>/group_<medium>_report<affix>.html
```

with its box-plot images saved alongside in `<path.sim>/group_plots/`. See [Outputs](doc_outputs.md#group-outputs) for the full file listing.

### Running the report

**Manually**, point it at your project config (a YAML path or a loaded `parameters` struct). Subjects are auto-discovered, or you can pass an explicit list:

```matlab
% Auto-discover all completed subjects under path.sim
prestus_group_report_start('config/config_study.yaml');

% Restrict to an explicit subject list
prestus_group_report_start('config/config_study.yaml', [1 2 3 5 7]);
```

**Automatically**, set the module flag so the report is regenerated at the end of every per-subject pipeline run:

```yaml
modules:
  generate_group_report: 1   # regenerate the group report after each subject
```

When enabled, `prestus_pipeline` calls `discover_group_subjects` and rebuilds the report over all subjects completed so far. The call is best-effort (wrapped in `try`/`catch`): a group-report problem never aborts the subject's own pipeline. See [Modules](doc_modules.md).

### Requirements

Each included subject must already have its per-subject output table on disk:

```
<path.sim>/sub-NNN/sub-NNN_<medium><affix>.csv
```

Per-subject PNGs (positioning, intensity, and max-temperature overlays) are embedded when present. The report reads `path.sim`, `simulation.medium`, `io.output_affix`, and `modules.run_heating_sims` (the latter gates the thermal plots).

### Report contents

- **Header** — subject count, medium, affix, and simulation path.
- **Subject filter** — checkboxes to toggle subjects in / out.
- **Exposure dashboard** — per-metric mean ± SD with min–max range, colour-coded against ITRUSST non-significant-risk limits. The Mechanical sub-grid always shows both **MI (free water)** and **MItc** (transcranial) tiles, plus a **Max pressure** tile (≤ 2 MPa). Only one of the two MI tiles has real data for any given medium (free water for water/free-field runs, MItc for layered runs) — the other renders as an empty, grayed-out "N/A" tile, so MI (free water) is effectively grayed out whenever MItc is available, and vice versa. Per-tissue `MI_brain`/`MI_skull`/`MI_skin` are not shown on this board (they remain in the per-subject CSV and the Mechanical Index box plot below).
- **Subject roster** — one row per subject with key metrics and a link to that subject's own report.
- **Group acoustic summary** — two box plots: **Intensity** (Isppa / Ipa_target, W/cm²) and **Mechanical Index** (per tissue, unitless).
- **Group thermal summary** *(layered medium with heating simulations only)* — three box plots: **maximum temperature**, **temperature rise**, and **CEM43 thermal dose**, each on its own scale.
- **Per-subject cards** — collapsible panels embedding each subject's figures.
- **Configuration summary** — the parameters used for the run.

All group box plots (`render_group_boxplot` in `generate_group_report.m`) auto-scale their figure width to the number of columns being plotted (`max(600, 140 × columns)` px, fixed 420 px height) and are exported as 150 dpi PNGs into `group_plots/` alongside the HTML report.

> **Note:** The subject filter updates the dashboard and per-subject cards live, but the group box plots are static images rendered over the full subject set — they do not change when subjects are toggled.
