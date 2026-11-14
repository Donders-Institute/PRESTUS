# Nexus Calibration and Localite Pages

The Calibration and Localite pages bring `calibration_gui` and `neuronav_gui` into the Nexus. Both compute in a MATLAB session that the Nexus starts when an action first needs it, so their first step waits for MATLAB to start. This page describes that session, the two pages, and the standalone MATLAB GUIs that the **Tools** menu opens.

On this page:

- [MATLAB Session](#matlab-session): when the shared MATLAB session starts and ends.
- [Calibration Page](#calibration-page): calibrating the study's transducers to hydrophone profiles.
- [Localite Page](#localite-page): ingesting Localite neuronavigation recordings.
- [Standalone MATLAB GUIs](#standalone-matlab-guis): opening `neuronav_gui` and `calibration_gui` from the Nexus.

!!! note "At a glance"
    - **When you need it:** You calibrate a study's transducer to hydrophone profiles, or ingest Localite neuronavigation recordings.
    - **Before you start:** MATLAB (the binary set on the [Run tab](doc_nexus_tabs.md#run-tab)); for calibration, a hydrophone profile CSV; for Localite, the recordings under `path.localite_raw`.

## MATLAB Session

The Calibration and Localite pages compute in MATLAB, with the same functions as `calibration_gui` and `neuronav_gui`.

??? note "When the session starts and ends"
    The Nexus keeps one MATLAB session open, started the first time an action needs it, e.g. *Load*, *Predict* or *Fit* on the Calibration page or *Detect series* on the Localite page; showing a page does not start it. The MATLAB binary is the one set on the Run tab. The first call takes as long as MATLAB needs to start (15–60 s); later ones take only the computation. The pages show the session's state (*starting*, *ready*, *working*), and **Stop** ends it.

    - The session ends by itself after 30 min without use, since it holds a MATLAB licence, and when the Nexus quits or stops answering. The next call starts a new one, and the pages open their transducer or session again.
    - Long work does not use the session: batch ingestion and the k-Wave drive check are written as scripts and run on the Run tab, where the log streams and *Stop* works.

??? developer "How the session exchanges data"
    The session is `nexus_bridge`, started with `matlab -batch`. The Nexus and MATLAB exchange JSON files in a temporary folder. MATLAB runs only `nexus_calibration`, `nexus_neuronav` and `nexus_ping` from it.

## Calibration Page

*Transducer ▸ Calibration* is the transducer-calibration GUI for the study's own transducers. It does what `calibration_gui` does ([Transducer library](doc_transducer_library.md)), computed by the same functions (`calibration_desired_profiles`, `calibration_evaluate`, `calibration_fit`, `calibration_save_entry`).

The controls are grouped in three steps: **1 Data** (profile CSV, measured and desired depths), **2 Fit** (where to calibrate, reference ISPPA, *Fit*, *Evaluate*) and **3 Save** (save, *Export*, k-Wave drive check). The fitter's settings (ring profile, free parameters, forward model, near-field skip, loss, Huber δ, starts) are under **Fit options**, collapsed by default, because the defaults suit most fits.

??? note "The plots"
    Until MATLAB has answered, the right-hand side shows the current library entry from its YAML alone:

    - **Quality:** the global model's R² per calibrated depth, with the geometric and per-depth-fit medians as reference lines.
    - **Drive:** the drive amplitude per depth.
    - **Phases:** each ring's phase correction per depth, unwrapped; steps over 30° between neighbouring depths are marked.

    The calibrated range is shaded and this study's depths are drawn as lines, coloured by whether the entry covers them. The same plot sits under the transducer drawing on the Settings sub-tab. After *Predict*, *Evaluate* or *Fit*, the *Profiles* view shows the profiles with a bar per depth coloured by R² (green ≥ 0.9, amber ≥ 0.8, red below); clicking a bar shows that depth.

??? note "The controls"
    - **Transducer** lists the config's transducers that use a device from `config/equipment`; the driving system and `transducer.library_path` are taken from the config. The current library entry is the one PRESTUS would use (the user folder first, then `config/calibration`).
    - **Profile CSV** defaults to the combination's `char_data_path` under `local/profiles`. Loading it shows the geometric model and the library entry against the measured profiles (R² per depth).
    - **Desired depths** and **Predict** model any depths through the runtime library path, with the measured, interpolated or (within *extrapolate ≤*) extrapolated profile at each depth. The desired depths are prefilled with the depths this study uses with the device: the focal depth settings (`focal_distance_ep` or `_bowl` of every matching transducer, and the manifest's `transducer[(i)].focal_distance_ep` columns) and each matching transducer's target depth in the placement viewer (exit plane → `focus_pos`, 0.1 mm). They follow the config and placement each time the page is shown, until you type your own; **Study depths** fills them in again. Without either, the depths of the entry's `meta.fit`, then the device's defaults, are used. The plot then shows focus against target per depth, and the table the phases and drive each depth resolves to.
    - **Ring profile** (*keep the entry's* / *fit* / *uniform rings*) selects the within-ring velocity profile ([Within-ring velocity profile](doc_calibration_refinement.md#within-ring-velocity-profile)); *fit* fits it on the measured profiles before the per-depth fits. **Options…** sets the profile fit and a plane-scan file (off-axis data) fitted with the axial profiles.
    - **Fit** calibrates per-depth refinements at the ticked measured depths or at the desired depths. **Reference ISPPA** defaults to `transducer.target_isppa_wcm2`. *Advanced…* holds the fitter's further options (as in `calibration_gui`; see [Calibration](doc_calibration_refinement.md#refining-a-library-entry)), the desired-profile tolerance and the k-Wave drive check's grid and correction, and loads a preset.
    - **Save entry** writes `<serial>[_<ds>].yaml` to `transducer.library_path`. If that is not set, it asks for a folder and sets it in the config. The entry records how it was fitted (`meta.fit`), and opening it later restores those settings. The device panel on the Settings sub-tab then shows the new entry as the one PRESTUS will use.
    - **k-Wave drive check…** (per-depth fits) writes the unsaved fit and `calibration/kwave_drive_<combo>.m` next to the config. The script runs `kwave_calibrate_drive` with the study config's grid and platform, then saves the drive-corrected entry. Run it on the Run tab, or yourself (e.g. on the cluster); it takes minutes per depth.
    - *Export profiles…* and *Export settings…* write the same CSVs as `calibration_gui`.

Matrix transducers and new equipment entries (`calibration_gui`'s *View geometry* and *New…*) are only available in `calibration_gui`.

## Localite Page

*Placement ▸ Localite* ingests Localite recordings with the study config's settings (`path.localite_raw`/`_post`, `placement.localite.*`, `io.overwrite_files`, `neuronav.*`). These are edited on the Settings sub-tab (*Edit settings*), not on the page. The config is saved before each step, because MATLAB reads it with `load_parameters` as a run does. The steps are those of `neuronav_gui` ([doc_placement_neuronav.md](doc_placement_neuronav.md)).

??? note "Detect series and write outputs"
    1. **Detect series** (`neuronav_detect_series`, `neuronav_name_series`) lists the series of the chosen subject and session with their QC (pulses, position and rotation deviation, tracking flags). Each series is named after the coil table, the planned target nearest its beam axis, or `series_<k>`. Tabs show the pulse positions in 3-D (with the beam direction and the planned targets), the deviation per pulse, the marker files found and the log.
    2. **Write outputs** (`neuronav_check_target_map`, `neuronav_ingest_markers`) writes the ticked series under their (editable) names and colours to `path.localite_post/sub-XXX/ses-XX/`. When the inputs or the settings change after detection, the page asks you to detect again.

??? note "Using the written outputs, batch ingestion and MNI view"
    - **Use in this study** sets `placement.localite.session` and `.target` to the written session and target(s) (one per transducer).
    - **Import into study…** builds manifest rows from all ingested sessions (Study tab).
    - **Ingest all sessions…** scans every session (`neuronav_scan_sessions`) and lists the planned-target labels found, whose target names you can edit (*Group labels by stem*, *Use labels as found*). *Save label map and run batch* writes the map to `neuronav.target_labels` and the write mode to `neuronav.batch_write`, then runs `ingest_<study>.m` (`neuronav_ingest_batch` on the ticked sessions) on the Run tab.
    - **MNI view…** maps `neuronav_ingest_summary.csv` to MNI space with each subject's SimNIBS warp (`neuronav_summary_to_mni`, needs `path.seg`). It plots transducer, focus and beam axis over the MNI152 head outline, with filters per subject, session and target.

## Standalone MATLAB GUIs

The **Tools** menu opens the Localite and Calibration pages, and also opens PRESTUS's standalone MATLAB GUIs in an interactive MATLAB session (they need the MATLAB desktop):

- **neuronav_gui in MATLAB…** ([Localite ingestion](doc_placement_neuronav_workflows.md#gui-workflow)). The study config is saved first and passed to the GUI.
- **calibration_gui in MATLAB…** saves entries to `transducer.library_path` when that is set. It also offers *View geometry* and *New…* (adding a transducer to `config/equipment`), which the Calibration page does not.
