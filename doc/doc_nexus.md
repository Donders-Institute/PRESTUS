# PRESTUS Nexus

!!! note "At a glance"
    - **What:** A desktop app (experimental) that builds a study config with live checks, shows placements and results, and runs the study locally or submits it to a cluster.
    - **When you need it:** To plan placements, check a config and submit cluster jobs from a laptop, as an alternative to the MATLAB GUI.
    - **Before you start:** The `prestus-nexus` conda environment ([Installation](#installation)). MATLAB runs the simulations.

The Nexus (Python/Qt) writes an ordinary override YAML for `load_parameters`, containing only the values that differ from `config_default.yaml`. MATLAB does all loading, validation and simulation. The Nexus is experimental; the MATLAB GUI ([doc_gui.md](doc_gui.md)) remains available.

## Starting

From MATLAB:

```matlab
prestus_nexus                                  % untitled study (File > New / Open)
prestus_nexus('config_study.yaml')             % open a study config
prestus_nexus(Python='/path/to/env/bin/python')
```

From a terminal:

```bash
cd PRESTUS/nexus
python -m prestus_nexus [config_study.yaml]
```

MATLAB's default `python3` usually lacks PyQt6, so pass the environment's Python: `prestus_nexus(Python='<conda>/envs/prestus-nexus/bin/python')`. From a terminal, `conda activate prestus-nexus` followed by `prestus-nexus` also starts it.

## Workflow

### Start tab

The Nexus opens on *Start*, which sets up a first simulation in four steps. Every other setting keeps its PRESTUS default. Each step shows *To do* until it is filled in, then *Done* or the number of problems (hover for the messages). The Start page writes the same config keys as the other tabs, so the two never disagree.

1. **Head image:** choose the T1-weighted NIfTI, or drop it onto the window. From the file name it sets `path.anat`, `path.t1_pattern` (as `sub-%1$03d_…` when the name holds `sub-<n>`) and `subject_id`. It also finds a T2 next to the T1 (`T1w`→`T2w`) and a folder with `m2m_sub-XXX` in the anatomy folder, one of its subfolders or one of its siblings (`path.seg`). If the segmentation is missing, it looks for SimNIBS (`charm` on the `PATH`, then the installers' default folders on macOS, Windows and Linux) to set `startup.simnibs_bin_path`. The results folder defaults to `prestus_sim` next to the study config, or next to the anatomy folder if the config is unsaved, never inside the anatomy folder.
2. **Transducer:** the device and driving system from the equipment library, the focal depth setting (with the device's range) and the free-water I_SPPA. The calibration entry PRESTUS will use is shown below.
3. **Target:** an MNI target, either alone (`placement.mode: heuristic`, which needs a target name) or with a transducer position (`mni`). A study using another mode is pointed to the Placement tab.
4. **Check and run:** *Run on* (`platform`) and the errors still to fix, each a link to its field. *Continue to Run* saves the study and opens the Run tab.

### Study config

The Nexus opens a given config, or an untitled one. **File ▸ New study…** creates a config and can copy `config_default.yaml` into its folder, which pins the defaults as `prestus_config_init` does. **Open…** and **Open recent** switch to another config. The subject, paths, T1/T2 patterns and SimNIBS folder are set on the Setup tab. The Issues panel checks them as you type: the T1 must resolve, and the `m2m_sub-XXX` folder must exist or SimNIBS must be set.

### Tabs

The window has one tab per area; each shows its config sections as cards:

| Tab | Sections | Also on the tab |
|---|---|---|
| Start | (the four steps above) | |
| Setup | subject, paths, environment, input/output | |
| Transducer | transducer, calibration conventions | device panel, cross-section drawing, calibration coverage; 'Calibration' sub-tab |
| Placement | placement (incl. the Localite folders), Localite ingestion | placement viewer ('Viewer' sub-tab, next to 'Settings'); 'Localite' sub-tab; Localite shortcuts |
| Simulation | medium & solvers, grid, sonication protocol, thermal, head model, pseudo-CT, tissue layers and properties, jwave, BabelBrain | timing diagram, TUS calculator exchange; 'Pseudo-CT' sub-tab |
| Study | study name, targeting & analysis | study manifest editor |
| Run | pipeline stages, where to run, HPC | run log, YAML review |

- Each field has a readable name. Hovering over it or clicking **?** shows the full help from `config_default.yaml`, the config key and the default.
- Changed values are marked ●, fields filled from the equipment library are marked *from library*, and ↺ resets a field.
- A tab's label shows a red (errors) or amber (warnings only) dot with the number of problems in its fields. Its tooltip also counts the changed values.
- Fields appear only when they apply to the current modes (e.g. the Localite settings only with `placement.mode: localite`).
- *Placement mode* marks each mode for the configured subject. **Unavailable** (greyed out; the tooltip names the missing input) means PRESTUS cannot produce an input the mode needs: no T1, no SimNIBS segmentation (PlanTUS: also the mesh `sub-XXX.msh`) and no `startup.simnibs_bin_path` to run charm, or no recorded Localite poses. **Needs setup** means charm will run first or the mode's tool paths are still empty. Tool paths never grey a mode out, since they are only shown once it is selected. The mode in use always stays selectable, and why it cannot run yet is written under the field.

### Advanced options

Each tab hides the settings most studies leave at their defaults, and says how many are hidden.

- **+ n advanced** in a card's header shows that card's advanced settings (*Hide advanced* folds them again). **Show advanced options** shows them on every card of the tab, including the fully advanced sections.
- An advanced setting that holds a value is always shown. Jumping to an issue opens the card it is on.
- **Locked cards.** The tissue properties card (`medium_properties`) is read-only. **🔒 Unlock** asks for confirmation and makes it editable for the session; **Lock** reverts. Values from a loaded config still apply.
- **Uncertainty variants.** With *Uncertainty analysis* on, the Simulation tab shows the four variant medium files and a card listing the default, Temp-H, Temp-L, Trans-H and Trans-L properties, with the brain and skull range of each pair plotted one property per column. Thermal properties show the Temp pair only. Greyed values are not in the variant file and use the default.
- **Pseudo-CT sub-tab.** Shows `pseudoCT.nii.gz` over the T1 and the in-skull distributions of the mapped skull properties against the band of the uncertainty bound variants' skull values, with the fraction of skull voxels outside it. The distributions come from the default run when it exists; otherwise **Predict skull properties** predicts them without simulating, labelled as a prediction. See [Pseudo-CT runs](doc_uncertainty.md#pseudo-ct-runs).
- **Tune for** (Speed, Accuracy, Memory, Disk space) shows the settings that affect that goal, advanced ones included, on every tab, and marks them after their label (e.g. *Grid resolution · speed, accuracy, memory*). The goal is highlighted. The tag says only what a setting affects; check its help (**?**) before changing it. Set it back to **,** to return to the normal view.

??? note "Details: what counts as advanced"
    - *Setup:* extra MATLAB paths, FSL/ANTs locations, HPC modules for pseudo-CT, the individual `io.save_*` switches, external acoustic NIfTI.
    - *Transducer:* element geometry (diameters, curvature, exit-plane distance), matrix layout, bowl-referenced depths, focal-depth reference convention. Amplitude, phases and steering stay visible.
    - *Placement:* Localite ingestion details (export subfolder, pulses per series, tracker offsets), the heuristic search criteria, PlanTUS weights and outputs, TUS_entry options, neuronav naming rules.
    - *Simulation:* precision and interactive mode, everything on the grid except its resolution, protocol time steps and cool-down, the thermal model, and the whole head-model, segmentation, tissue-layer and tissue-property sections. Debug outputs and the pseudo-CT section stay visible so the toggles that turn them on are never hidden.
    - *Study:* the target area radius.
    - *Run:* intermediate pipeline stages (grid, medium, source, analyses), report style, cluster profile, reservation, job naming and waiting.

    Locked cards are set by `locked: true` under `sections:` in `overlay.yaml`. Which keys count for which *Tune for* goal is listed under `tuning:` in `nexus/prestus_nexus/overlay.yaml`, mainly the grid, time-step, FOV, precision and solver settings and the large `io.save_*` outputs.

### Transducer tab

- Pick a device and driving system from `config/equipment`, or enter custom geometry.
- The panel shows the focal range (`min_foc`/`max_foc`) and which calibration entry PRESTUS will use.
- With several transducers, transducer 2 and later inherit transducer 1's values unless set.
- The cross-section on the right draws the elements, exit plane, geometric focus, the focal depth you set (✕) and the device's steering range.
- The *Calibration* sub-tab compares, fits and saves library entries for the config's transducers ([Calibration page](#calibration-page)).

### Placement viewer

The Placement tab has two sub-tabs: **Settings** (the `placement`/Localite cards) and **Viewer** (manual placement), which takes the full tab width.

- Axial, coronal and sagittal T1 slices with the segmentation overlay. Without a T1, the viewer shows what is missing and a button to the field.
- **Click sets** (*Focus* · *Transducer* · *Browse*) chooses what a click on a slice does: set the focus or the transducer (bowl apex), or with *Browse* move the crosshair without writing anything (scrolling always browses). The bowl is drawn from the device geometry.
- With several transducers, every bowl, beam axis and acoustic FOV box is drawn; the ones not selected are dimmed. Each FOV box is resolved from its own transducer (`'auto'` diameter and length), as PRESTUS does: coherent coupling simulates their union, async coupling one box per transducer ([details](doc_multitransducer.md#partial-acoustic-fov)).
- **+**/**−** by the transducer selector add or remove a transducer without switching to the Transducer tab.
- **Coordinates…** opens a row to type X/Y/Z directly instead of clicking. **RAS mm** or **Voxel** need no conversion. **MNI mm** is converted to this subject's space with SimNIBS when you click, which needs the subject's `m2m_sub-XXX` folder; **Set focus** uses the nonlinear warp and **Set transducer** the affine, as in [MNI placement](doc_placement_mni.md#how-a-position-is-resolved). **Go to** only moves the crosshair.
- With manual placement, the readout's red flags are also listed in the Issues panel, for every transducer: a target more than 2 mm from the focal depth setting, an exit plane (centre or rim) inside the scalp, and beam rays through an eye. They are checked once the viewer (or the Run tab) has loaded the T1 and segmentation; clicking one opens the viewer on that transducer.
- *Snap to scalp* keeps the beam axis and puts the exit plane a chosen gap outside the scalp. *Transducer at focal depth* (under **More**, with the gap) instead puts the transducer where the focal depth setting lands on the target. **View** holds the overlay opacity, the FOV size and *Centre on focus*.
- **Beam:** with a segmentation, the readout adds skull on the beam (mean ± SD over 37 rays from the aperture to the focus; the central ray), the mean skull incidence (amber above 20°), the gap between the exit-plane rim and the scalp along the beam (mean ± SD; red when the minimum is negative, i.e. the transducer cuts into the head; a large SD means uneven contact) and the share of rays that pass an eye. These are the measures of the placement comparison (`placement.beam_qc`).
- The readout keeps two depths apart: the **target depth** (exit plane → `focus_pos`, from the placement) and the **focal depth setting** (`focal_distance_ep`, the depth set on the equipment). For annular arrays the setting selects the calibrated phases, so the simulated focus lands at the setting even when the target is deeper or shallower; the readout gives the difference (red above 2 mm). *Set focal depth to target* copies the target depth into the setting.
- Trajectories can be imported from Brainsight (NIfTI coordinates) or 3D Slicer, using BabelBrain's conventions.
- Positions are written in the space set by `placement.manual.coordinate_space` (`voxel` or `ras`). See [doc_placement.md](doc_placement.md#manual-coordinate-selection).
- Other placement modes are shown read-only: the Localite pose from the `neuronav_gui` outputs, the heuristic choice from `transducer_positioning`, and for every mode the positions the last run's placement stage resolved (`<path.sim>/sub-XXX/sub-XXX_placement<affix>.json`). To get them without simulating, run the pipeline once with only `modules.run_transducer_placement` enabled.
- With `placement.mode: localite`, the Placement card offers *Localite page* and *Import into study…*; the *Localite* sub-tab detects, names and writes the recorded series of a session ([Localite page](#localite-page)).
- With `placement.mode: plantus`, the Placement card runs only the [PlanTUS placement](doc_placement_plantus.md#running-plantus): *Choose in Workbench…* opens PlanTUS's maps to pick the position, *Automatic* takes the composite-score optimum. The run shows on the Run tab, and the Viewer opens on the result.

??? note "Details: MNI coordinates in the viewer"
    **MNI mm** is converted with SimNIBS's `mni2subject_coords` CLI tool, the same one-shot external call `mni2subject_coords_LDfix.m` uses; it runs only when you click, not on every redraw. It needs the `m2m_sub-XXX` folder (under `path.seg`/`path.anat`) for both the nonlinear and affine transform. Setting `startup.simnibs_bin_path` alone is not enough, and the viewer says so (distinctly from "SimNIBS's bin folder isn't set") when the m2m folder is missing. **Set focus** uses the nonlinear warp (most accurate in-brain) and **Set transducer** the affine (`12dof`, well-defined outside the brain too). If the m2m folder exists but its nonlinear warp failed or was never produced, **Set focus** retries with the affine and says so; this does not help when SimNIBS has not been run at all. **Go to** moves the crosshair using the focus's transform, or the transducer's if *Click sets transducer* is selected. **Set focus** / **Set transducer** only work with `placement.manual.coordinate_space` set to `voxel` or `ras`. A study using `coordinate_space: mni` writes MNI values directly (PRESTUS converts them at placement-stage runtime) and is not shown in the viewer; type an MNI point and use **Go to** to look at the matching subject-space location.

### Issues panel

The panel lists what PRESTUS would otherwise report only at load or run time: unknown or misspelt keys (with suggestions), deprecated keys, wrong types and values, backend combinations and frequency mismatches, missing inputs and missing calibration entries, a `target_isppa_wcm2` that cannot be met (no library calibration and `modules.run_water_baseline: 0`) or that adds a free-water simulation (`auto`), a project `config_default.yaml` that lacks keys added to PRESTUS since it was copied, and (manual placement) the placement viewer's checks. Clicking an issue opens its tab and field, showing advanced options if needed. The panel is closed by default: the header shows the error and warning counts, and clicking them (or **View ▸ Issues…**, Ctrl+Shift+I) opens the list. It closes on Esc, a click outside, or a click on an issue.

### Run tab

- The highlighted button follows *Run on*: *Run locally* for `platform: matlab`, *Submit to cluster* for `slurm` or `qsub` (and for `auto` when a scheduler or the SSH host is available). The other one, *Write launcher* and *Save log…* are under **More**. The HPC card hides the scheduler fields while *Run on* is `matlab`.
- *Validate in MATLAB* only runs `load_parameters`.
- *Show deployment, save & resource details* lists the jobs the run submits, which outputs it saves, and a resource estimate from earlier runs of the project (`<path.sim>/run_stats.jsonl`): acoustic and thermal simulation time and peak RAM, as median and maximum. See [Resource estimates](doc_hpc.md#resource-estimates-from-earlier-runs).
- *Run locally* writes `run_<study>.m` next to the config and runs it with `matlab -batch`, streaming the log. It sets `platform: matlab` for that run only. You can also run the launcher yourself with `run('run_<study>.m')`.
- *Submit to cluster* (when `sbatch` or `qsub` is on this computer, e.g. on a cluster login node) writes `submit_<study>.m`, which keeps the config's `platform` and `hpc.*` settings, and runs it: PRESTUS submits as usual (chained pipelines submit several jobs). MATLAB exits once the jobs are submitted (unless `hpc.wait_for_job` is set). The *Run profile* presets fill in the `hpc.*` fields. The submitted jobs appear in the *Jobs* side tab (details below).
- *Submit over SSH* submits from a computer without a scheduler, e.g. a laptop with the project drive mounted. Set the login node (*Host*, e.g. `user@mentat005.dccn.nl` or an alias from `~/.ssh/config`), the folder mapping (*Local folder* `/Volumes/myproject` *is on the cluster* `/project/myproject`), the *Setup commands* run before MATLAB (e.g. `module load slurm matlab/R2024a`) and, if the cluster uses its own checkout, the *Remote PRESTUS* folder. The study folder must be under the mapped folder. The settings are kept per user (`~/.config/prestus_nexus/last_session.yaml`).
- *Review YAML* shows the exact YAML that will be written.
- *Results* shows what the runs wrote to `<path.sim>/sub-XXX/` (below).

??? note "Details: the Jobs tab"
    The Nexus reads the job IDs from the log and lists them in the *Jobs* side tab (next to *Run log* and *Review YAML*, brought forward after a submit), refreshing them every 30 s until they finish. *Pop out ↗*, or double-clicking a job, opens the same view in its own window, which can be maximised or moved to a second screen. For SLURM jobs the list shows the state, elapsed time and time limit, time left (red when under 10 % of the limit, at least 15 min and at most 1 h, i.e. a `TIMEOUT` is near), the partition with the node (or, while pending, the reason such as `(Priority)`), peak memory (`sstat` while running, `sacct` after) and the exit code; with qsub only the state. Clicking a job shows its logs from `<output>/log_hpc` in three tabs, *Output* (`sub-XXX_slurm_output_<id>.log`), *Errors* and *Report* (`sub-XXX_job<id>_report.txt`, written by the report job after the job ends: wall time, resources, CO₂e, cost), below the job table (drag the divider to share the space), refreshed every 10 s while the job runs. *Follow* keeps the end of the log in view, *Wrap lines* wraps long lines, and A−/A+ change the text size (remembered). Logs are read from the file system when the folder is visible (on the cluster, or through the SSH folder mapping), else over SSH. *Cancel selected* runs `scancel`/`qdel`.

??? note "Details: submitting over SSH"
    - *Submit to cluster* runs `load_parameters` + `prestus_pipeline_start` on the host with `matlab -batch` (in a login shell, from the mapped study folder), and the job list follows the jobs over SSH.
    - *Connect…* logs in once: ssh's prompts (password, one-time code, an unknown host key) appear in a Nexus dialog, and the login is kept as a shared OpenSSH connection that every later submit, status check and cancel reuses until it has been idle for 8 h (*Disconnect* closes it). With an SSH key, *Connect…* is not needed. The other commands never prompt (`BatchMode=yes`), so without a key or an open connection they fail with "Permission denied". Connection sharing needs OpenSSH 8.4 or later on macOS or Linux; on Windows use a key.
    - *Test connection* reports whether `sbatch`/`qsub`, MATLAB and PRESTUS are found on the host. If the host's `matlab` is a site wrapper that needs a display (at the Donders, `module load matlab/…` gives one; the log then says "No DISPLAY variable set"), set *Remote MATLAB* to the MATLAB binary itself (`<matlabroot>/bin/matlab`); `matlab -batch` runs without a display.
    - Config values under the mapped folder are translated the same way, e.g. `path.sim: /Volumes/myproject/data/sim` runs as `/project/myproject/data/sim`: the Nexus passes them to `load_parameters` as overrides (listed in the log) and leaves the config file unchanged. This covers the merged config (study, project `config_default.yaml` copy, transducers). A path on this Mac outside the mapped folder (`/Volumes/…`, `/Users/…`) cannot be translated; the Nexus lists those keys and asks before submitting.

**Results.** One folder usually holds several runs, told apart by their output affix (study conditions, uncertainty variants, async and sequential runs). The **Run** list names each one from its affix and starts on the run the open config produces (for async coupling: the joint result; a second list switches between the joint result and each transducer's own field, its `_trans-NN` run). Every sub-tab follows the selected run. Files are read when needed: only the shown sub-tab reads its files, volumes are read in the background (the viewer says *Loading …*) and kept for the next view, and opening the tab again reloads nothing unless the runs' files changed (*Refresh* rescans and reloads).

- *Acoustic* / *Thermal*: the run's maps on the T1 (only the maps it wrote are offered), centred on the map's peak (*Go to peak* re-centres). A colour bar gives the map's unit (W/cm², Pa, °C, min). Intensity, pressure, MI, temperature rise and CEM43 are coloured from 0 to the map's maximum; absolute temperatures from their 1st to 99th percentile. *Hide below* hides values under a share of the colour scale's top (10 % by default for intensity, pressure, MI and temperature rise). Clicking a slice moves the crosshair; the line under the slices gives its RAS position and the map's value there. Every transducer (square) and focus (x) is marked, as the run's placement stage recorded them (`sub-XXX_placement<affix>.json`), else from the config; points away from the slice are drawn faint.
- *Heating*: the run's thermal time courses from `cache/sub-XXX_<medium>_HEAT<affix>.mat` (written by the thermal stage). *Metric* picks the peak temperature, temperature rise, CEM43 (k-Wave) or CEM43 (ISO) of each tissue layer (max over the layer), or the temperature, rise or CEM43 on the focal axis (per tissue and at the target). Tick the tissues to show; the choice is kept across metrics and runs. The legend gives each tissue's maximum. *Log scale* helps with CEM43; *Shade sonication* marks the sonication periods (HEAT files written since 2026-10-05). The toolbar zooms, pans and saves the plot. Needs scipy (or h5py for `-v7.3` files). For an async transducer's own field the joint run's time courses are shown.
- *Report*: the run's HTML report, or any other report in the folder (sequential and multi-ISPPA summaries).
- *Data Table*: the run's metrics (last CSV row).
- *Compare runs*: one row per run with the key metrics (ISPPA, MI, focal distance and volume, temperature rise, CEM43); double-click a row to show that run.

### Study tab

Runs many subjects with per-run settings. See [Group analysis](doc_group.md#study-setup).

- The table (the study manifest, `study.manifest`) has one row per simulation run; settings not in the table come from the study config. A single subject does not need it. *How to fill the table* shows the steps and an example.
- *Add subjects…* adds one run per subject ID (`1-12, 15, sub-020`). *Add setting column…* adds a config key to vary per run (empty cell = the study value) or the run columns `include`, `config` and `notes`; hovering over a column header shows its help. *More* holds adding an empty row, removing rows or a column, and *Save table as…*.
- The table is saved to the *Table file* (relative to the config). Without one, the first save writes `<config name>_runs.tsv` next to the config and sets `study.manifest`.
- *Import from Localite…* reads the `neuronav_gui` outputs under `path.localite_post` (`neuronav_ingest_summary.csv` and the per-session `InstrumentMarker_*.json`). It makes one run per target, or one run per session with the targets assigned to transducers 1, 2, …, with a preview before adding. See [Group analysis](doc_group.md#from-localite-ingestion).
- The summary under the table checks the manifest as you type: subject IDs, `include` values, unknown keys (with suggestions), and rows that would write to the same outputs.
- *Check all runs* saves the table and runs `prestus_study_start(config, 'DryRun', true)`. *Run all runs locally* writes `run_study_<study>.m` and runs it; *Write launcher script* only writes it. These stay disabled while the table is empty. *Build group report* runs `prestus_group_report_start(config)` (group maps + report). Checks, runs and the report stream to the Run tab.

## Timing diagram and TUS calculator

The *Sonication protocol* card (Simulation tab) draws the protocol from `timing.*`: one pulse (with its ramp), one pulse train, and the train repetitions. It labels PD, PRI/PRF, PTD, PTRI and PTRD with the pulse count, duty cycles, the on-time, and the Ipa/Itp factor of the ramp (`compute_ramp_factor`). It also lists inconsistent settings (e.g. a pulse longer than its repetition interval).

The [TUS calculator](https://www.socsci.ru.nl/fusinitiative/tuscalculator/) (Radboud University, MIT licence) estimates derived quantities (derated intensities, MI, energy deposition, temperature rise) for a protocol. The Nexus does not copy its calculations. It exchanges settings files with it instead:

- **Export for TUS calculator…** writes a JSON file for the calculator's *Upload settings*: pulse duration, ramp duration and shape, PRI, PTD, PTRI, PTRD (in ms), whether the train and repetition levels are enabled, and (when set) the frequency, target Isppa, aperture and focal depth.
- **Import from TUS calculator…** reads a settings file downloaded from the calculator (JSON or YAML). It shows the values it will set (in seconds), and optionally also sets the frequency and target Isppa.
- **Open TUS calculator ↗** opens the web app.

Both use the same ramp convention: the ramp duration is the rise (= fall) time at each end of a pulse (`timing.pulse_ramp_dur`, `timing.pulse_ramp_shape` `linear` or `tukey`). PRESTUS models pulse ramps only, so the calculator's pulse-train and repetition ramps are exported as 0 and reported when an imported file sets them. A level disabled in the calculator sets the corresponding PRESTUS durations to 0.

## MATLAB session

The Calibration and Localite pages compute in MATLAB, with the same functions as `calibration_gui` and `neuronav_gui`. The Nexus keeps one MATLAB session open (`nexus_bridge`, started with `matlab -batch` the first time an action needs it, e.g. *Load*, *Predict* or *Fit* on the Calibration page or *Detect series* on the Localite page; showing a page does not start it). The MATLAB binary is the one set on the Run tab. The first call takes as long as MATLAB needs to start (15–60 s); later ones take only the computation. The pages show the session's state (*starting*, *ready*, *working*), and **Stop** ends it.

- The session ends by itself after 30 min without use, since it holds a MATLAB licence, and when the Nexus quits or stops answering. The next call starts a new one, and the pages open their transducer or session again.
- Long work does not use the session: batch ingestion and the k-Wave drive check are written as scripts and run on the Run tab, where the log streams and *Stop* works.

??? note "Details: how the session exchanges data"
    The Nexus and MATLAB exchange JSON files in a temporary folder. MATLAB runs only `nexus_calibration`, `nexus_neuronav` and `nexus_ping` from it.

## Calibration page

*Transducer ▸ Calibration* is the transducer-calibration GUI for the study's own transducers. It does what `calibration_gui` does ([Transducer library](doc_transducer_library.md)), computed by the same functions (`calibration_desired_profiles`, `calibration_evaluate`, `calibration_fit`, `calibration_save_entry`).

The controls are grouped in three steps: **1 Data** (profile CSV, measured and desired depths), **2 Fit** (where to calibrate, reference ISPPA, *Fit*, *Evaluate*) and **3 Save** (save, *Export*, k-Wave drive check). The fitter's settings (ring profile, free parameters, forward model, near-field skip, loss, Huber δ, starts) are under **Fit options**, collapsed by default, because the defaults suit most fits.

Until MATLAB has answered, the right-hand side shows the current library entry from its YAML alone:

- **Quality:** the global model's R² per calibrated depth, with the geometric and per-depth-fit medians as reference lines.
- **Drive:** the drive amplitude per depth.
- **Phases:** each ring's phase correction per depth, unwrapped; steps over 30° between neighbouring depths are marked.

The calibrated range is shaded and this study's depths are drawn as lines, coloured by whether the entry covers them. The same plot sits under the transducer drawing on the Settings sub-tab. After *Predict*, *Evaluate* or *Fit*, the *Profiles* view shows the profiles with a bar per depth coloured by R² (green ≥ 0.9, amber ≥ 0.8, red below); clicking a bar shows that depth.

- **Transducer** lists the config's transducers that use a device from `config/equipment`; the driving system and `transducer.library_path` are taken from the config. The current library entry is the one PRESTUS would use (the user folder first, then `config/calibration`).
- **Profile CSV** defaults to the combination's `char_data_path` under `local/profiles`. Loading it shows the geometric model and the library entry against the measured profiles (R² per depth).
- **Desired depths** and **Predict** model any depths through the runtime library path, with the measured, interpolated or (within *extrapolate ≤*) extrapolated profile at each depth. The desired depths are prefilled with the depths this study uses with the device: the focal depth settings (`focal_distance_ep` or `_bowl` of every matching transducer, and the manifest's `transducer[(i)].focal_distance_ep` columns) and each matching transducer's target depth in the placement viewer (exit plane → `focus_pos`, 0.1 mm). They follow the config and placement each time the page is shown, until you type your own; **Study depths** fills them in again. Without either, the depths of the entry's `meta.fit`, then the device's defaults, are used. The plot then shows focus against target per depth, and the table the phases and drive each depth resolves to.
- **Ring profile** (*keep the entry's* / *fit* / *uniform rings*) selects the within-ring velocity profile ([Transducer library](doc_calibration.md#optional-within-ring-velocity-profile)); *fit* fits it on the measured profiles before the per-depth fits. **Options…** sets the profile fit and a plane-scan file (off-axis data) fitted with the axial profiles.
- **Fit** calibrates per-depth refinements at the ticked measured depths or at the desired depths. **Reference ISPPA** defaults to `transducer.target_isppa_wcm2`. *Advanced…* holds the fitter's further options (as in `calibration_gui`; see [Calibration](doc_calibration.md#refining-a-library-entry-to-your-hydrophone-profiles)), the desired-profile tolerance and the k-Wave drive check's grid and correction, and loads a preset.
- **Save entry** writes `<serial>[_<ds>].yaml` to `transducer.library_path`. If that is not set, it asks for a folder and sets it in the config. The entry records how it was fitted (`meta.fit`), and opening it later restores those settings. The device panel on the Settings sub-tab then shows the new entry as the one PRESTUS will use.
- **k-Wave drive check…** (per-depth fits) writes the unsaved fit and `calibration/kwave_drive_<combo>.m` next to the config. The script runs `kwave_calibrate_drive` with the study config's grid and platform, then saves the drive-corrected entry. Run it on the Run tab, or yourself (e.g. on the cluster); it takes minutes per depth.
- *Export profiles…* and *Export settings…* write the same CSVs as `calibration_gui`.

Matrix transducers and new equipment entries (`calibration_gui`'s *View geometry* and *New…*) are only available in `calibration_gui`.

## Localite page

*Placement ▸ Localite* ingests Localite recordings with the study config's settings (`path.localite_raw`/`_post`, `placement.localite.*`, `io.overwrite_files`, `neuronav.*`). These are edited on the Settings sub-tab (*Edit settings*), not on the page. The config is saved before each step, because MATLAB reads it with `load_parameters` as a run does. The steps are those of `neuronav_gui` ([doc_placement_neuronav.md](doc_placement_neuronav.md)):

1. **Detect series** (`neuronav_detect_series`, `neuronav_name_series`) lists the series of the chosen subject and session with their QC (pulses, position and rotation deviation, tracking flags). Each series is named after the coil table, the planned target nearest its beam axis, or `series_<k>`. Tabs show the pulse positions in 3-D (with the beam direction and the planned targets), the deviation per pulse, the marker files found and the log.
2. **Write outputs** (`neuronav_check_target_map`, `neuronav_ingest_markers`) writes the ticked series under their (editable) names and colours to `path.localite_post/sub-XXX/ses-XX/`. When the inputs or the settings change after detection, the page asks you to detect again.

- **Use in this study** sets `placement.localite.session` and `.target` to the written session and target(s) (one per transducer).
- **Import into study…** builds manifest rows from all ingested sessions (Study tab).
- **Ingest all sessions…** scans every session (`neuronav_scan_sessions`) and lists the planned-target labels found, whose target names you can edit (*Group labels by stem*, *Use labels as found*). *Save label map and run batch* writes the map to `neuronav.target_labels` and the write mode to `neuronav.batch_write`, then runs `ingest_<study>.m` (`neuronav_ingest_batch` on the ticked sessions) on the Run tab.
- **MNI view…** maps `neuronav_ingest_summary.csv` to MNI space with each subject's SimNIBS warp (`neuronav_summary_to_mni`, needs `path.seg`). It plots transducer, focus and beam axis over the MNI152 head outline, with filters per subject, session and target.

## MATLAB GUIs

The **Tools** menu opens the Localite and Calibration pages, and also opens PRESTUS's standalone MATLAB GUIs in an interactive MATLAB session (they need the MATLAB desktop):

- **neuronav_gui in MATLAB…** ([doc_gui.md](doc_gui.md#localite-ingest-gui)). The study config is saved first and passed to the GUI.
- **calibration_gui in MATLAB…** saves entries to `transducer.library_path` when that is set. It also offers *View geometry* and *New…* (adding a transducer to `config/equipment`), which the Calibration page does not.

## Project defaults

When the study folder holds a copy of `config_default.yaml` (from *File ▸ New study…* or `prestus_config_init`), the Nexus merges it onto the toolbox default exactly as `load_parameters` does: the copy's values win, keys added to PRESTUS since the copy was made keep their toolbox defaults. The *Review YAML* page names both files. **Tools ▸ Update project defaults…** replaces the copy with the toolbox default (keeping a dated backup), optionally after writing the values the copy pinned differently into the study config so the simulation settings stay the same; `prestus_config_upgrade(folder)` does the same from MATLAB.

## Installation

One command creates a conda environment `prestus-nexus` with everything the Nexus uses (Python 3.11, PyQt6 with WebEngine, PyYAML, numpy, nibabel, matplotlib, scipy):

```bash
cd PRESTUS/nexus
conda env create -f environment.yml
```

After a PRESTUS update that changes the dependencies: `conda env update -n prestus-nexus -f environment.yml`.

??? note "Details: pip, Qt and h5py"
    Without conda, `pip install -e .` in `PRESTUS/nexus` installs the same packages into any Python ≥ 3.10 environment; the list lives in `nexus/pyproject.toml`.

    Install Qt from pip only (as `environment.yml` does). conda's `pyqt` together with pip's `PyQt6-WebEngine` fails to load (`Library not loaded: …QtPrintSupport`). In an older environment built with conda's `pyqt`, run `conda remove pyqt` and then `pip install -e .`.

    Reading `-v7.3` HEAT files on the *Heating* sub-tab also needs `h5py` (not installed by default; PRESTUS saves `-v7` unless MATLAB's preference says otherwise).

??? note "Details: local copy and cache folder"
    The Nexus runs from a local copy of its own code. On start it copies the files of `nexus/prestus_nexus` that changed since the last start into the user's cache folder (`~/Library/Caches/prestus_nexus` on macOS, `~/.cache/prestus_nexus` on Linux) and runs from there, so Python does not read its modules from a network drive (on the Donders project drive this halves the start-up time). The checkout stays the source: an edit there is copied on the next start. The parsed `config/equipment` files are cached in the same folder while they are unchanged. Set `PRESTUS_NEXUS_LOCAL_COPY=0` to run from the checkout directly, and `PRESTUS_NEXUS_CACHE` to use another cache folder.

??? note "Details: standalone app"
    `nexus/packaging/prestus_nexus.spec` builds a standalone app with PyInstaller (`pip install -e . pyinstaller`, then `pyinstaller packaging/prestus_nexus.spec` in `nexus/`). The app does not contain PRESTUS: on first start it asks for the PRESTUS folder (or use `--prestus-root` / `PRESTUS_ROOT`) and remembers it.

??? note "Details: compatibility with MATLAB's YAML loader"
    PRESTUS reads configs with `yaml.loadFile` (SnakeYAML, YAML 1.1). The Nexus:

    - quotes all strings (unquoted `yes`/`no`/`04:00:00` would change type);
    - reads floats such as `2.17e-3` the way SnakeYAML does;
    - writes every element of a transducer list with the same top-level keys, which `yaml.loadFile` requires.

    A round-trip test writes configs with the Nexus, loads them with `load_parameters`, and compares every value. It covers both a serial-based and a fully inline transducer.

??? note "Details: maintenance"
    The schema is derived from `config_default.yaml` at startup. `nexus/prestus_nexus/overlay.yaml` adds only widget kinds, allowed values, units, visibility rules, keys the code reads that the default omits, and deprecated keys. `nexus/tests/test_schema.py` fails when the overlay no longer matches the default config or the code. See `nexus/README.md` for running the tests.
