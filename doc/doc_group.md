# Group Analysis

A study simulates many subjects or conditions, and comparing them needs the same settings applied to every run and the results collected in one place. PRESTUS runs every row of a study table with its own settings, measures per run how well the focus hits the target, and summarises the runs in group MNI maps and a group report. The report and the maps never run simulations; they read the per-run outputs that already exist on disk.

On this page:

- [Targeting Accuracy](#targeting-accuracy): how far each simulated focus lies from the planned focus, an MNI target, and an MNI ROI.
- [Group MNI Maps](#group-mni-maps): focal-volume overlap, mean and SD intensity, and mean heating per condition.
- [Group HTML Report](#group-html-report): one interactive report across runs, with a one-row-per-run CSV.

!!! note "At a glance"
    - **When you need it:** For studies with more than a handful of runs.
    - **Before you start:** A study config that works for one subject.

Running many subjects or conditions from a study manifest is described in [Study Setup](doc_group_study-setup.md).

## Targeting Accuracy

With `modules.run_targeting_analysis: 1` (default, layered or phantom media, `run_nifti_creation` on), each run writes `sub-XXX_<medium><affix>_targeting.csv` next to its per-subject CSV. The analysis reads the focal-mask [NIfTIs](https://nifti.nimh.nih.gov/) and never fails the pipeline; unavailable measures are `NaN`. Targets and ROIs can differ per run: set `analysis.target_mni_mm` or `analysis.target_roi_mni` as manifest columns. With several transducers, the planned focus and beam axis of `transducer(1)` are used.

??? note "Targeting CSV columns"
    | Columns | Meaning |
    |---|---|
    | `condition`, `placement_mode`, `target_label`, `transducer_name`, `n_transducers`, `focal_distance_ep_planned_mm`, `target_isppa_wcm2` | The run's settings, shown in the report's study design table. |
    | `focal_centroid_ras_*`, `peak_ras_*` | Centroid of the −3 dB focal volume (brain-masked) and the intensity peak inside it, in T1 RAS mm. |
    | `planned_focus_ras_*`, `planned_trans_ras_*` | Focus and transducer positions resolved by placement, recorded before the grid stage, in RAS mm. |
    | `focus_error_mm` | Distance between the focal centroid and the planned focus. |
    | `focus_error_axial_mm`, `focus_error_lateral_mm` | That error along and across the planned beam axis (transducer → planned focus); axial `> 0` means the simulated focus is deeper than planned. |
    | `peak_error_mm` | Distance between the peak and the planned focus. |
    | `focal_centroid_mni_*`, `peak_mni_*` | The same positions in MNI mm (needs `io.save_MNI`). |
    | `target_mni_*`, `target_mni_source` | The intended MNI target: `analysis.target_mni_mm`, else the MNI target of the placement mode (`placement.heuristic` / `placement.plantus` `.mni_target_mm`, `placement.mni.focus_pos_mm`). |
    | `target_error_mni_mm`, `peak_target_error_mni_mm` | Distance from the MNI focal centroid (or peak) to the MNI target. |
    | `roi_volume_mm3`, `roi_overlap_mm3`, `roi_focal_fraction`, `roi_coverage`, `roi_mean_intensity`, `roi_max_intensity`, `roi_centroid_dist_mm` | With `analysis.target_roi_mni` (an MNI NIfTI mask, any grid): ROI volume, overlap with the focal volume, fraction of the focal volume inside the ROI, fraction of the ROI covered, intensity in the ROI (`intensity_quantity` = Ipa, else Itp), and the distance between the centroids. |

The same analysis runs post hoc with `targeting_analysis(parameters, 'PlannedFocus', f, 'PlannedTrans', t)` (T1 voxel positions).

## Group MNI Maps

`prestus_group_report_start` builds these by default (`'Maps', false` to skip), for layered media only (`simulation.medium` containing `layered` in the study or subject config); other media get the report without maps. They use each run's MNI NIfTIs and are computed per condition, or across all runs when there are no conditions. `<out>` is `<path.sim>/group/<study.name>/` for a study and `<path.sim>/group/<medium><affix>/` otherwise.

| File in `<out>/maps/` | Content |
|---|---|
| `group_<cond>_focal_overlap_3dB.nii.gz` | Fraction of runs whose −3 dB focal volume covers each voxel (0–1). |
| `group_<cond>_<Ipa\|Itp>_mean.nii.gz`, `…_sd.nii.gz` | Voxelwise mean and SD intensity (Ipa when every run has it, else Itp when every run has that; no intensity maps otherwise). |
| `group_<cond>_heatrise_mean.nii.gz` | Mean temperature rise (when every run has heating maps). |

`<out>/plots/group_<cond>_focal_overlap.png` shows three orthogonal MNI slices through the overlap peak, with the ROI outline, the per-run focal centroids, and the targets. The background is the mean MNI-space T1 of the subjects (`<parent of path.seg>/mni/sub-XXX/T1_to_MNI_post-hoc.nii.gz`, or the copy inside `m2m_sub-XXX/toMNI/` of earlier outputs, when present). Runs whose MNI grid differs from the first run's are skipped with a warning.

## Group HTML Report

```matlab
prestus_group_report_start('config/config_pilot.yaml');           % study: runs from the manifest
prestus_group_report_start('config/config_single.yaml');            % no manifest: auto-discover subjects
prestus_group_report_start('config/config_single.yaml', [1 2 3 5]); % explicit subject list
prestus_group_report_start(cfg, [], 'Maps', false);                 % report only
```

| Mode | Runs | Written to |
|---|---|---|
| Study (`study.manifest` set) | Included manifest rows | `<path.sim>/group/<name>/group_<name>_report.html`, `group_<name>_runs.csv`, `plots/`, `maps/` |
| Subjects | All subjects in `path.sim` with `sub-XXX_<medium><affix>.csv` (or the given list) | `<path.sim>/group_<medium>_report<affix>.html`, `group_<medium>_runs<affix>.csv`, `group_plots/` |

With `modules.generate_group_report: 1`, `prestus_pipeline` rebuilds the subject-mode report after each subject (best-effort, without group maps). Runs started by `prestus_study_start` skip this and build one study report at the end.

The report reads each run's `sub-XXX_<medium><affix>.csv` (last row) and, when present, the targeting sidecar and the per-subject PNGs. Runs without a CSV are listed as missing. It opens with a safety banner for runs above any [ITRUSST](https://itrusst.org/) non-significant-risk reference level, a group safety dashboard, and a per-run filter; condition comparison and targeting sections appear when the study has conditions or targeting sidecars.

??? note "Report sections"
    - **Header**: study, manifest, runs, subjects, conditions, missing runs, medium, and a link to the one-row-per-run CSV.
    - **Safety banner**: runs above any ITRUSST non-significant-risk reference level, with links to their rows.
    - **Filter**: checkboxes per run and per-condition buttons. The dashboard, roster, and cards update live; box plots and figures show the full set.
    - **Group safety dashboard and other metrics**: mean ± SD, range, and N per ITRUSST metric, coloured by the worst run against the NSR limit. The Mechanical sub-grid shows both MI (free water) and MItc tiles; the tile that does not apply to the medium is greyed out.
    - **Study design** (study, or more than one condition): runs per condition, the manifest keys and how many values each takes, and each run's placement mode, target, transducer, planned depth, target intensity, and manifest settings.
    - **Condition comparison** (more than one condition): mean ± SD (n) of the key metrics per condition, plus box plots. It is descriptive only: runs of one subject under several conditions are not independent.
    - **Targeting accuracy** (targeting sidecars present): focus and MNI-target errors, ROI measures per condition and per run, and a figure of the MNI centroids and the axial/lateral error.
    - **Group MNI maps** (when built): overlap figures and links to the NIfTIs.
    - **Outliers**: runs whose value lies far from their condition's median (modified z-score 0.6745·(x − median)/MAD, |z| > 3.5; Iglewicz & Hoaglin, 1993), for conditions with at least 5 runs. They are also flagged in the roster.
    - **Data analysis**: client-side distribution, correlation, and by-group views (group by subject, condition, placement mode, target, or transducer when they vary).
    - **Roster**: one sortable row per run with key metrics, NSR and outlier flags, and a link to the run's own report.
    - **Group acoustic and thermal summaries**: box plots (I<sub>SPPA</sub> per tissue, MI per tissue; maximum temperature, temperature rise, and CEM43 when thermal results exist).
    - **Per-run cards**: the run's settings and figures.
    - **Configuration summary**: of the loaded (study) config. Per-run settings are in the study design table.

??? note "Older entry point (`prestus_group_start`)"
    `prestus_group_start(subject_list, config)` runs `prestus_group_report_start` (group MNI maps + report). Plot options of earlier versions are ignored; for an ROI, set `analysis.target_roi_mni`. ROI statistics (intensity in the ROI, overlap, centroid distance) are part of the [targeting analysis](#targeting-accuracy).
