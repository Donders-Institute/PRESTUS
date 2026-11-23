# Nexus Tabs and Pages

Most Nexus controls explain themselves: each field shows its name, its help from `config_default.yaml` (hover or **?**), its config key, and its default. This page covers what the window does not show: which config keys a control writes, what a button does beyond the obvious, where displayed data comes from, and what a feature needs.

On this page:

- [Tab Overview](#tab-overview): what each tab holds, and how fields appear and hide.
- [Menus and Shortcuts](#menus-and-shortcuts): undo, search, offline docs, and errors.
- [Segmentation Step](#segmentation-step): running SimNIBS charm on its own.
- [QuickStart Tab](#quickstart-tab): the four steps and the config keys they write.
- [Advanced Options](#advanced-options): hidden settings and locked tables.
- [Transducer Tab](#transducer-tab): devices, calibration, and several transducers.
- [Calibration Page](#calibration-page): calibrating the study's transducers to hydrophone profiles.
- [Placement Viewer](#placement-viewer): what the viewer writes and what its checks measure.
- [Localite Page](#localite-page): ingesting Localite neuronavigation recordings.
- [Issues Panel](#issues-panel): what is checked before a run.
- [Run Tab](#run-tab): what each run button writes and runs, and submitting over SSH.
- [Results Tab](#results-tab): where the shown values come from.
- [Study Tab](#study-tab): the run table and the study commands.
- [Timing Diagram and TUS Calculator](#timing-diagram-and-tus-calculator): editing the protocol and exchanging it with the TUS calculator.
- [MATLAB Session](#matlab-session): when the session behind the Calibration and Localite pages starts and ends.
- [Standalone MATLAB GUIs](#standalone-matlab-guis): opening `neuronav_gui` and `calibration_gui` from the Nexus.

!!! note "At a glance"
    - **When you need it:** To find out what a control writes or runs, where a value comes from, or why a field or mode is not offered.
    - **Before you start:** The Nexus installed and started ([PRESTUS Nexus (GUI)](doc_nexus.md)); [A First Simulation](doc_nexus.md#a-first-simulation) walks through the main tabs once.

## Tab Overview

| Tab | Holds |
|---|---|
| QuickStart | A first simulation in four steps (below). |
| Setup | Subject, paths, environment, input and output. |
| Transducer | Device and driving system, geometry, drive and steering; *Calibration* and *Equipment* sub-tabs. |
| Placement | Placement method and positions; *Viewer*, *Beam axis*, *3D*, and the method pages (Localite, PlanTUS, TUS entry). |
| Medium | Medium type; head model, tissue properties, uncertainty, and skull CT sub-tabs. |
| Sonication | Pulse timing, multi-transducer coupling, thermal model. |
| Solver | Acoustic and thermal solvers, their compute device, grid, time reversal. |
| Study | The run table for several subjects or conditions. |
| Run | Pipeline stages, where to run, the run plan, the YAML, launch, and jobs. |
| Results | The runs in the results folder, and follow-up runs. |

Fields appear only when they apply to the current choices, e.g. the Localite settings only with `placement.mode: localite`. In the *Placement method* row, a mode is **unavailable** (greyed; the tooltip names the missing input) when PRESTUS cannot produce an input it needs: no T1, no SimNIBS segmentation and no `startup.simnibs_bin_path` to run charm (PlanTUS also needs the mesh `sub-XXX.msh`), or no recorded Localite poses. **Needs setup** means charm will run first or the mode's tool paths are still empty. The mode in use always stays selectable.

Figures export as SVG, PDF, or PNG (3D views as PNG only), in the light theme unless **Keep current theme** is ticked.

## Menus and Shortcuts

- **Undo / Redo** (Ctrl+Z / Ctrl+Shift+Z) cover edits on every tab in one history, including viewer clicks and the Study tab's table. A button that writes several values is one step. *Reload from disk* and opening another study clear the history; saving does not.
- **Find setting…** (Ctrl+F) searches names, config keys, choices, and help text, and opens the setting on its tab. Settings that do not apply to the current choices are listed greyed.
- **Help ▸ Use local documentation files** opens the Markdown files in `doc/` instead of the published pages, e.g. offline. F1 opens the page of the shown tab.
- **Tools ▸ Usage statistics…** records which [usage statistics](doc_telemetry.md#choosing-what-to-share) levels you share, in `~/.prestus/telemetry.json`. This is needed because the Nexus runs MATLAB with `matlab -batch`, where PRESTUS never asks. The Nexus itself sends nothing. What you share shows at the right of the status bar; click it to change it. The Nexus asks at launch until you decide, unless `PRESTUS_TELEMETRY` is set. If you allowed an earlier version of the notice, it asks once more.
- Ctrl+1 to Ctrl+9 switch tabs; Ctrl+Shift+V validates in MATLAB; Ctrl+R runs locally.

??? note "When something goes wrong"
    - An unexpected error shows a dialog instead of closing the Nexus; the traceback is written to `nexus_errors.log` in the [cache folder](doc_nexus.md#installation-options).
    - Saving writes to a temporary file first, so a full or read-only share leaves the old config intact.
    - The MATLAB session must start within 180 s, and a request may take up to 600 s; after a time-out the session restarts on the next request.

## Segmentation Step

charm takes an hour or more, and its output (`m2m_sub-XXX`) feeds the layered head model, heuristic, MNI, PlanTUS, and TUS_entry placement, and the viewer's scalp checks. It can therefore run on its own first.

- **Run segmentation** (Medium tab, *Segmentation* card) runs only charm for this subject (`modules.segmentation_only: 1`, `simulation.medium: layered`): as a cluster job where the Run tab submits there, else in MATLAB on this computer.
- **Segment again…** sets `io.overwrite_simnibs: 1` for that run only: charm deletes the m2m folder and segments again. Placements, head models, and MNI conversions made from the old segmentation then no longer match it.
- *Incomplete* means tissue labels or MNI warps are missing: charm did not finish or failed.

## QuickStart Tab

QuickStart writes the same config keys as the other tabs, and every setting it does not show keeps its PRESTUS default.

1. **Head image:** the T1 sets `path.anat`, `path.t1_pattern` (as `sub-%1$03d_…` when the name holds `sub-<n>`), and `subject_id`. A T2 next to it (`T1w`→`T2w`) and an `m2m_sub-XXX` folder in, under, or next to the anatomy folder (`path.seg`) are found automatically. Without a segmentation, it looks for SimNIBS 4 (`charm` on the `PATH`, then conda environments and the installers' default folders) to set `startup.simnibs_bin_path`. The results folder defaults to `prestus_sim` next to the study config (or next to the anatomy folder while the config is unsaved), never inside the anatomy folder.
2. **Transducer:** device, driving system, [focal depth setting](doc_glossary.md#focal-depth-setting), and free-water [I<sub>SPPA</sub>](doc_glossary.md#isppa). A device with one calibrated driving system selects it. Without a calibration entry of its own, a line names the generic entry PRESTUS falls back to, or what to do instead.
3. **Target:** the coordinate space decides the placement mode. An MNI target alone gives `placement.mode: heuristic` (which needs a target name); an MNI target with a scalp point gives `mni` (the transducer is moved onto the scalp, [exit plane](doc_glossary.md#exit-plane) `placement.skin_gap_mm` outside it). A target and bowl apex in the participant's T1 give `manual` with `placement.manual.coordinate_space: ras` or `voxel`, used as given. For a participant-space target without a transducer position, use PlanTUS on the Placement tab.
4. **Check and run:** sets `platform` and lists the errors left. *Continue to Run* saves the study.

## Advanced Options

Each tab hides the settings most studies leave at their defaults; **Advanced** shows them, on the whole tab or per card. An advanced setting that holds a value is always shown. Settings the chosen setup does not read stay hidden too, e.g. anatomy, SimNIBS, and targeting settings with a `water` medium.

??? note "Locked tables: tissue layers, tissue properties, uncertainty values"
    - *Tissue layers* and the *Tissue properties* table (`medium_properties`) are read-only until **🔒 Unlock** (with advanced options shown); the unlock lasts for the session. Emptying a cell resets it. Values from a loaded config apply whether locked or not.
    - With a `layered` medium and a segmentation, a tissue none of whose labels occur in `final_tissues.nii.gz` is greyed as not modelled, as PRESTUS drops its layer. For a charm segmentation this is *Skull (single layer)* (label 4); with `pct.enabled`, the cortical and trabecular columns, which merge into the single skull layer. Volumes are for the whole head, not the simulation grid.
    - A tissue sets either the absorption fraction or the absorption coefficient; the table shows the other in brackets as implied. Typing into an implied cell switches the tissue to that field.
    - On the *Uncertainty* sub-tab, an edited variant value sets [`simulation.uncertainty_overrides`](doc_uncertainty.md#customising-medium-property-ranges) and raises a warning.

??? note "Skull CT sub-tab"
    For a pseudo-CT, the UTE→pseudo-HU line is read from `pct_create_pseudoCT.sh` or the subject's fit. The HU→property curves are evaluated by MATLAB (`pct_mapping_curves`). The in-skull property distributions come from the default run when it exists; otherwise **Predict skull properties** predicts them without simulating. With the BabelBrain backend, the sub-tab also holds BabelBrain's own HU mapping (`babelbrain.ct_mapping`, `babelbrain.air_*`). See [Pseudo-CT runs](doc_uncertainty.md#pseudo-ct-runs).

??? developer "Where locked cards are defined"
    Locked cards are set by `locked: true` under `sections:` in `nexus/prestus_nexus/overlay.yaml`. The `tuning:` list there (which keys affect speed, accuracy, memory, or disk) is kept for a later *Tune for* control and is not shown.

## Transducer Tab

Devices and driving systems come from `config/equipment`. Fields filled from the library are marked *from library*; *Geometry* starts collapsed when the library supplies it.

- With several transducers, transducer 2 and later inherit transducer 1's values unless set.
- *Sweep…* next to the target I<sub>SPPA</sub> sets a list of targets, which runs the [multi-I<sub>SPPA</sub> pipeline](doc_multiintensity.md).
- The *Calibration* sub-tab compares, fits, and saves library entries ([Calibration page](#calibration-page)). With **Advanced**, the *Equipment* sub-tab adds transducers and driving systems to your library ([New transducer or driving system in the Nexus](doc_transducer_library_maintenance.md#new-transducer-or-driving-system-in-the-nexus)).

## Calibration Page

*Transducer ▸ Calibration* is the transducer-calibration GUI. It computes in the Nexus's [MATLAB session](#matlab-session) and needs a hydrophone profile CSV. It does what `calibration_gui` does ([Transducer library](doc_transducer_library.md)), computed by the same functions (`calibration_desired_profiles`, `calibration_evaluate`, `calibration_fit`, `calibration_save_entry`). It works on the study's transducers, or, with **Any library transducer**, on any transducer and driving system of the equipment library: entries are then saved to a folder of your own library (*Folder…*) without changing the config, and the k-Wave drive check asks for a simulation config (grid, platform, HPC; or PRESTUS defaults). **View geometry** shows the selected transducer on the *Equipment* sub-tab.

The plots fill the top of the page and the controls sit below them, grouped in three steps: **1 Data** (profile CSV, measured and desired depths), **2 Fit** (where to calibrate, reference [I<sub>SPPA</sub>](doc_glossary.md#isppa), *Fit*, *Evaluate*), and **3 Save** (save, *Export*, *Report*, k-Wave drive check). The fitter's settings (ring profile, free parameters, forward model, near-field skip, loss, Huber δ, starts) are under **Fit options**, collapsed by default, because the defaults suit most fits.

??? note "The plots"
    Until MATLAB has answered, the plot area shows the current library entry from its YAML alone:

    - **Quality:** the global model's R² per calibrated depth, with the geometric and per-depth-fit medians as reference lines.
    - **Drive:** the drive amplitude per depth.
    - **Phases:** each ring's phase correction per depth, unwrapped; steps over 30° between neighbouring depths are marked.

    The calibrated range is shaded and this study's depths are drawn as lines, coloured by whether the entry covers them. The same plot sits under the transducer drawing on the Settings sub-tab. After *Predict*, *Evaluate*, or *Fit*, the *Profiles* view shows the profiles with a bar per depth coloured by R² (green ≥ 0.9, amber ≥ 0.8, red below); clicking a bar shows that depth.

??? note "The controls"
    - **Transducer** lists the config's transducers that use a device from `config/equipment`; the driving system and `transducer.library_path` are taken from the config. The current library entry is the one PRESTUS would use (the user folder first, then `config/calibration`).
    - **Profile CSV** defaults to the combination's `char_data_path` under `local/profiles`. Loading it shows the geometric model and the library entry against the measured profiles (R² per depth).
    - **Desired depths** and **Predict** model any depths through the runtime library path, with the measured, interpolated, or (within *extrapolate ≤*) extrapolated profile at each depth. The desired depths are prefilled with the depths this study uses with the device: the focal depth settings (`focal_distance_ep` or `_bowl` of every matching transducer, and the manifest's `transducer[(i)].focal_distance_ep` columns) and each matching transducer's target depth in the placement viewer (exit plane → `focus_pos`, 0.1 mm). They follow the config and placement each time the page is shown, until you type your own; **Study depths** fills them in again. Without either, the depths of the entry's `meta.fit`, then the device's defaults, are used. The plot then shows focus against target per depth, and the table the phases and drive each depth resolves to.
    - **Ring profile** (*keep the entry's* / *fit* / *uniform rings*) selects the within-ring velocity profile ([Within-ring velocity profile](doc_calibration_hydrophone-profiles.md#within-ring-velocity-profile)); *fit* fits it on the measured profiles before the per-depth fits. **Options…** sets the profile fit and a plane-scan file (off-axis data) fitted with the axial profiles.
    - **Fit** calibrates per-depth refinements at the ticked measured depths or at the desired depths. **Reference I<sub>SPPA</sub>** defaults to `transducer.target_isppa_wcm2`. *Advanced…* holds the fitter's further options (as in `calibration_gui`; see [Calibration](doc_calibration_hydrophone-profiles.md#refining-a-library-entry)), the desired-profile tolerance and the k-Wave drive check's grid and correction, and loads a preset.
    - **Save entry** writes `<serial>[_<ds>].yaml` to `transducer.library_path`. If that is not set, it asks for a folder and sets it in the config. The entry records how it was fitted (`meta.fit`), and opening it later restores those settings. The device panel on the Settings sub-tab then shows the new entry as the one PRESTUS will use.
    - **k-Wave drive check…** (per-depth fits) writes the unsaved fit and `calibration/kwave_drive_<combo>.m` next to the config. The script runs `kwave_calibrate_drive` with the study config's grid and platform, then saves the drive-corrected entry; it takes minutes per depth. **Submit as a job** sends it to the cluster as a CPU job (`kwave_drive_submit`: 8 GB, wall time = depths × `hpc.timelimit` + 1 h), which submits one simulation job per depth and waits for them on a compute node; over SSH (Run tab) the Nexus writes a copy with the cluster's paths (`kwave_drive_<combo>_cluster.m`). The job appears in the Run tab's job list. *Run on the Run tab* runs it in the Nexus's MATLAB instead, which on a cluster login node keeps MATLAB waiting there; or run the script yourself.
    - *Export profiles…* and *Export settings…* write the same CSVs as `calibration_gui`.
    - **Report** opens the HTML report of the saved library entry (`generate_calibration_report`: provenance, fit quality, the settings a study gets at the fitted depths, and the amplitude lookup), written to a temporary folder.

## Placement Viewer

The viewer sets positions for manual placement and shows the result of the other modes read-only. Positions are written in the space of `placement.manual.coordinate_space` (`voxel` or `ras`; [Manual placement](doc_placement_manual.md#entering-coordinates)). A clicked transducer position is the bowl apex.

- **MNI coordinates** are converted with SimNIBS (`mni2subject_coords`) when you click, which needs the subject's `m2m_sub-XXX` folder; `startup.simnibs_bin_path` alone is not enough. **Set focus** uses the nonlinear warp and **Set transducer** the affine, as in [MNI placement](doc_placement_mni.md#how-a-position-is-resolved). If the warp is missing, **Set focus** falls back to the affine and says so. A study with `coordinate_space: mni` is not shown in the viewer, because PRESTUS converts its values only at run time.
- **Two depths:** the *target depth* (exit plane → `focus_pos`) comes from the placement; the *focal depth setting* (`transducer.focal_distance_ep`) is set on the equipment. For annular arrays the setting selects the calibrated phases, so the simulated focus lands at the setting even when the target is deeper or shallower. A difference above 2 mm is flagged.
- **Beam checks** (with a segmentation): skull on the beam (mean ± SD over 37 rays from the aperture to the focus), mean skull incidence (flagged above 20°), the exit-plane-to-scalp gap (negative means the transducer cuts into the head; a large SD means uneven contact), and the share of rays through an eye. These are the measures of `placement.beam_qc`. With manual placement, the red flags also appear in the Issues panel.
- *Snap to scalp* keeps the beam axis and puts the exit plane a chosen gap outside the scalp.
- Trajectories import from Brainsight, 3D Slicer, and k-Plan; see [Importing a Trajectory](doc_placement.md#importing-a-trajectory).
- Other modes show the positions the last run's placement stage resolved (`sub-XXX_placement<affix>.json`). To get them without simulating, run once with only `modules.run_transducer_placement` enabled.

??? note "Beam axis, 3D, and several transducers"
    - **Beam axis** shows the last layered run's `medium_masks`, or before a run the SimNIBS labels mapped through `layers` without PRESTUS's smoothing. The free-water profile is read from the run's axial profile (`img/sub-XXX_<medium>_axial_profile<affix>.csv`); before a run, it is predicted in MATLAB from the device's library model at the focal depth setting (annular library devices only). *Absolute* needs a run, since the prediction has no absolute scale.
    - **3D** draws each transducer in a cylindrical housing sized from `housing:` in its equipment entry; without one, the aperture's diameter, reaching 15 mm behind the bowl apex. The face is removed by default so the view does not identify the participant. It needs graphics support (Qt Quick 3D).
    - **Several transducers:** *Grid orientation* sets `grid.mode` ([Grid Alignment](doc_multitransducer.md#grid-alignment)). It warns when `focal_axis` falls back to transducer 1's axis (`grid.focal_axis_min_norm`), and when a beam tilts 20° or more on a grid with fewer than 6 points per wavelength in water. Each FOV box is resolved from its own transducer: coherent coupling simulates their union, async coupling one box each ([Partial Acoustic FOV](doc_multitransducer.md#partial-acoustic-fov)).

??? note "PlanTUS, heuristic, and TUS entry pages"
    - **PlanTUS** shows the skin maps of the last run (`plantus_metrics.mat` in the placement run folder). **Place here** sets `placement.plantus.vertex` and places the transducer on the same maps without recomputing them. *Choose in Workbench…* needs a display and always runs on this computer ([Running PlanTUS](doc_placement_plantus.md#running-plantus)).
    - **Heuristic** shows the candidates of the last run (`tpars_sub-XXX_<target>.csv`). **Place here** sets `placement.heuristic.candidate`; *Use the heuristic's choice* clears it. Distances are in voxels, as the heuristic computes them.
    - **TUS entry** shows the result from `tus_entry_report.csv` ([TUS_entry placement](doc_placement_tus_entry.md)).
    - Computing maps, candidates, or entries runs as a placement-only job: on the cluster where the Run tab submits there, else in MATLAB on this computer.

## Localite Page

*Placement ▸ Localite* ingests the Localite recordings under `path.localite_raw` in the Nexus's [MATLAB session](#matlab-session), with the study config's settings (`path.localite_raw`/`_post`, `placement.localite.*`, `io.overwrite_files`, `neuronav.*`). These are edited on the Settings sub-tab (*Edit settings*), not on the page. The config is saved before each step, because MATLAB reads it with `load_parameters` as a run does. The steps are those of `neuronav_gui` ([Neuronavigation Read-In](doc_placement_neuronav.md)).

??? note "Detect series and write outputs"
    1. **Detect series** (`neuronav_detect_series`, `neuronav_name_series`) lists the series of the chosen subject and session with their QC (pulses, position and rotation deviation, tracking flags). Each series is named after the coil table, the planned target nearest its beam axis, or `series_<k>`. Tabs show the pulse positions in 3-D (with the beam direction and the planned targets), the deviation per pulse, the marker files found, and the log.
    2. **Write outputs** (`neuronav_check_target_map`, `neuronav_ingest_markers`) writes the ticked series under their (editable) names and colours to `path.localite_post/sub-XXX/ses-XX/`. When the inputs or the settings change after detection, the page asks you to detect again.

??? note "Using the written outputs, batch ingestion, and MNI view"
    - **Use in this study** sets `placement.localite.session` and `.target` to the written session and target(s) (one per transducer).
    - **Import into study…** builds manifest rows from all ingested sessions (Study tab).
    - **Ingest all sessions…** scans every session (`neuronav_scan_sessions`) and lists the planned-target labels found, whose target names you can edit (*Group labels by stem*, *Use labels as found*). *Save label map and run batch* writes the map to `neuronav.target_labels` and the write mode to `neuronav.batch_write`, then runs `ingest_<study>.m` (`neuronav_ingest_batch` on the ticked sessions) on the Run tab.
    - **MNI view…** maps `neuronav_ingest_summary.csv` to MNI space with each subject's SimNIBS warp (`neuronav_summary_to_mni`, needs `path.seg`). It plots transducer, focus, and beam axis over the MNI152 head outline, with filters per subject, session, and target.

## Issues Panel

The panel lists what PRESTUS would otherwise report only at load or run time: unknown or misspelt keys, deprecated keys, wrong types and values, unsupported backend and multi-transducer combinations, missing inputs and calibration entries, heating timing problems, a `path.acoustic_input` or follow-up parent that is missing or does not match, a `target_isppa_wcm2` that cannot be met or that adds a free-water simulation, and a project `config_default.yaml` that lacks keys added since it was copied. Clicking an issue opens its field. Open the panel from the error and warning counts in the header, or with Ctrl+Shift+I.

## Run Tab

- *Validate in MATLAB* only runs `load_parameters`.
- *Run locally* writes `run_<study>.m` next to the config and runs it with `matlab -batch`, with `platform: matlab` for that run only. You can run the launcher yourself with `run('run_<study>.m')`.
- *Submit to cluster* (with `sbatch` or `qsub` on this computer) writes `submit_<study>.m`, which keeps the config's `platform` and `hpc.*` settings, and runs it; PRESTUS submits the jobs. MATLAB exits once they are submitted, unless `hpc.wait_for_job` is set.
- *Plan* shows what the run will do: each stage as run, reused, or skipped, the jobs it submits, the files it writes, and the predicted load ([Run Time and Resources](doc_run-time.md)).
- *Pipeline options* (sequential follow-ups, per-stage limits of the multi-I<sub>SPPA</sub> and uncertainty pipelines) are not config keys. The Nexus keeps them per user and study and writes them into the launchers.
- **Stop** ends MATLAB together with its launcher.

The *Jobs* sub-tab reads the job IDs from the log and follows them; their logs are read from `<output>/log_hpc`, from the file system when the folder is visible, else over SSH. *Cancel selected* runs `scancel` or `qdel`. Jobs that a job submits itself (the k-Wave drive check, pipeline stages) join the list, marked ↳. The list is kept when the Nexus closes; *Clear finished* removes the jobs that have ended. When a job completes, the Results tab reads its outputs; when a run or job ends while the Nexus is in the background, it notifies you.

A progress bar above the log (local runs) and above the selected job's logs (cluster jobs, also in the *Progress* column) shows the current stage, how far its solver is, and the time left; for a whole-pipeline run also the whole run, from the predicted stage times. A job whose time left exceeds its remaining time limit is flagged in red (it will probably end in TIMEOUT). Drive checks and study runs also count their depths or runs. How exact this is depends on the solver ([backends](doc_backend.md)): the k-Wave C++ binaries and the thermal simulation report their progress, the k-Wave MATLAB solvers only a time estimate (shown with `~`), and the jwave Helmholtz solver only its GMRES residual, with live view (else elapsed time only). Cluster jobs are checked every 30 s. With [`simulation.live_view: 1`](doc_parameters_simulation.md#simulation-type-simulation), a *Live view* sub-tab appears with the first frame and shows the evolving pressure and temperature fields, one tab per field, in mm with the transducer and focus marked, and the peak temperature over time. A local run switches to it once; *Live view* on the progress card returns to it, also for the selected cluster job. `examples/example_progress_water.yaml` is a free-water run of a few minutes to try both.

??? note "Submitting over SSH"
    *Submit over SSH* submits from a computer without a scheduler, e.g. a laptop with the project drive mounted.

    - Set the login node (*Host*, e.g. `user@mentat005.dccn.nl` or an alias from `~/.ssh/config`), the folder mapping (*Local folder* `/Volumes/myproject` *is on the cluster* `/project/myproject`), and the *Setup commands* run before MATLAB (e.g. `module load slurm matlab/R2024a`). The study folder must be under the mapped folder. The settings are kept per user in `~/.config/prestus_nexus/last_session.yaml`.
    - Paths in the config under the mapped folder are translated, e.g. `/Volumes/myproject/data/sim` runs as `/project/myproject/data/sim`. They are passed to `load_parameters` as overrides; the config file stays unchanged. A local path outside the mapped folder cannot be translated, and the Nexus asks before submitting.
    - *Connect…* logs in once and keeps a shared [OpenSSH](https://www.openssh.com/) connection for 8 h of idle time; with an SSH key it is not needed. Without a key or an open connection, commands fail with "Permission denied". Connection sharing needs OpenSSH 8.4 or later on macOS or Linux; on Windows, use a key.
    - If the host's `matlab` is a site wrapper that needs a display (the log says "No DISPLAY variable set"), set *Remote MATLAB* to `<matlabroot>/bin/matlab`.

## Results Tab

The **Run** list names each run in the results folder from its output affix. For async coupling, a second list switches between the joint result and each transducer's own field (`_trans-NN`).

- *Overview* groups the run's values as the ITRUSST consensus asks TUS studies to report them, with the [non-significant risk limits](doc_outputs_metrics-table.md#non-significant-risk-limits) of the HTML report. Each row names its source (metrics CSV, free-water run, targeting CSV, placement record, equipment library, or open config). Pulse timing and transducer geometry come from the open config, which can differ from what the run used. I<sub>SPTA</sub> per tissue is that tissue's I<sub>SPPA</sub> times the run's duty cycle. For runs without the `focal_fwhm_*` columns, the Nexus measures the focal widths on the run's intensity map, which can differ from PRESTUS's by a fraction of a voxel.
- *Spatial* shows only the maps the run wrote; positions come from the run's placement record, else from the config.
- *Heating* reads `cache/sub-XXX_<medium>_HEAT<affix>.mat` (written by `modules.run_thermal_analysis`; not with the BabelBrain thermal backend). It needs [SciPy](https://scipy.org/), or h5py for `-v7.3` files.

### Follow-up Runs

**Follow-up…** creates a new study config that starts from the selected run:

| Type | What it reuses |
|---|---|
| New intensity | The run's pressure map (`path.acoustic_input`): no acoustic simulation ([Reusing outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs)). |
| New timing (thermal only) | The pressure map; only the thermal simulation runs. |
| New placement, same head | The head model, built once (`headmodel.reuse_across_placements`). |
| Pose from the time-reversal search | A pose or its refocusing phases from the run's pose table ([Pose search](doc_time_reversal.md#pose-search)). |
| Sequential | The run's end temperature and thermal dose. Added to this config's *Pipeline options* rather than run on its own ([Sequential runs](doc_sequential.md)). |

Each follow-up gets a new `io.output_affix` and is saved next to the parent as `<parent>_<type>_<affix>.yaml`. It records its parent in a `# nexus-followup:` comment line, which MATLAB ignores.

## Study Tab

The run table is the study manifest (`study.manifest`; [Study setup](doc_group_study-setup.md#study-setup)): one row per run, and settings not in the table come from the study config. An empty cell takes the study value.

- Without a *Table file*, the first save writes `<config name>_runs.tsv` next to the config and sets `study.manifest`.
- *Import from Localite…* reads the `neuronav_gui` outputs under `path.localite_post` and makes one run per target, or one per session with the targets assigned to transducers ([Importing runs from Localite](doc_group_study-setup.md#importing-runs-from-localite)).
- *Check all runs* runs `prestus_study_start(config, 'DryRun', true)`; *Run all runs locally* writes and runs `run_study_<study>.m`; *Build group report* runs `prestus_group_report_start(config)`. Their output streams to the Run tab.

## Timing Diagram and TUS Calculator

The *Sonication protocol* card draws the protocol from `timing.*` and writes back what you change in it: drag an edge, or click it and type a duration (`20 ms`, `500 µs`; the PRI also takes a frequency, `5 Hz`). Dragged values snap to two significant digits.

![Nexus Sonication tab with the sonication timing diagram](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/nexus_simulation-timing_screenshot.png)

The [TUS calculator](https://www.socsci.ru.nl/fusinitiative/tuscalculator/) (Radboud University, MIT licence) estimates derived quantities for a protocol. The Nexus does not repeat its calculations; it exchanges settings files with it. **Export for TUS calculator…** writes a JSON file for the calculator's *Upload settings*; **Import from TUS calculator…** reads a file downloaded from it (JSON or YAML) and shows the values before setting them.

Both use the same ramp convention: the ramp duration is the rise (= fall) time at each end of a pulse (`timing.pulse_ramp_dur`, `timing.pulse_ramp_shape`). PRESTUS models pulse ramps only, so the calculator's pulse-train and repetition ramps are exported as 0 and reported when an imported file sets them. A level disabled in the calculator sets the matching PRESTUS durations to 0.

The *Sonication sequence* card below sets up stances fired one after another ([Planning a Sequence in the Nexus](doc_sequential.md#planning-a-sequence-in-the-nexus)).

## MATLAB Session

The [Calibration](#calibration-page) and [Localite](#localite-page) pages bring `calibration_gui` and `neuronav_gui` into the Nexus and compute in MATLAB, with the same functions. Their first step waits for MATLAB to start.

??? note "When the session starts and ends"
    The Nexus keeps one MATLAB session open, started the first time an action needs it, e.g. *Load*, *Predict*, or *Fit* on the Calibration page or *Detect series* on the Localite page; showing a page does not start it. The MATLAB binary is the one set on the Run tab. The first call takes as long as MATLAB needs to start (15–60 s); later ones take only the computation. The pages show the session's state (*starting*, *ready*, *working*), and **Stop** ends it.

    - The session ends by itself after 30 min without use, since it holds a MATLAB licence, and when the Nexus quits or stops answering. The next call starts a new one, and the pages open their transducer or session again.
    - Long work does not use the session: batch ingestion and the k-Wave drive check are written as scripts and run on the Run tab, where the log streams and *Stop* works.

??? developer "How the session exchanges data"
    The session is `nexus_bridge`, started with `matlab -batch`. The Nexus and MATLAB exchange JSON files in a temporary folder. MATLAB runs only `nexus_calibration`, `nexus_neuronav`, and `nexus_ping` from it.

## Standalone MATLAB GUIs

The **Tools** menu opens the Localite and Calibration pages, and also opens PRESTUS's standalone MATLAB GUIs in an interactive MATLAB session (they need the MATLAB desktop):

- **neuronav_gui in MATLAB…** ([Localite ingestion](doc_placement_neuronav_workflows.md#gui-workflow)). The study config is saved first and passed to the GUI.
- **calibration_gui in MATLAB…** saves entries to `transducer.library_path` when that is set. Its *View geometry* and *New…* are on the *Equipment* sub-tab ([New transducer or driving system in the Nexus](doc_transducer_library_maintenance.md#new-transducer-or-driving-system-in-the-nexus)).
