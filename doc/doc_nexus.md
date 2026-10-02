# PRESTUS Nexus

The PRESTUS Nexus is a standalone desktop app (Python/Qt) for putting a study config together before you run it. It writes an ordinary override YAML for `load_parameters`, containing only the values that differ from `config_default.yaml`. MATLAB still does all loading, validation and simulation. The Nexus is *experimental*; the MATLAB GUI ([doc_gui.md](doc_gui.md)) remains available.

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

Requirements: Python ≥ 3.10 with PyQt6 and PyYAML. The placement viewer also needs numpy, nibabel and matplotlib. `nexus/environment.yml` creates a suitable conda environment.

## Workflow

1. **Study config.** The Nexus opens a given config, or an untitled one. **File ▸ New study…** creates a config and can copy `config_default.yaml` into its folder, which pins the defaults as `prestus_config_init` does. **Open…** and **Open recent** switch to another config. The subject, paths, T1/T2 patterns and SimNIBS folder are set on the Setup tab, and the Issues panel checks them as you type: the T1 must resolve, and the `m2m_sub-XXX` folder must exist or SimNIBS must be set.
2. **Tabs.** The window has one tab per area; each shows its config sections as cards:

   | Tab | Sections | Also on the tab |
   |---|---|---|
   | Setup | subject, paths, environment, input/output | |
   | Transducer | transducer, calibration conventions | device panel, cross-section drawing |
   | Placement | placement (incl. the Localite folders), Localite ingestion | placement viewer (own 'Viewer' sub-tab, next to 'Settings'), Localite shortcuts |
   | Simulation | medium & solvers, grid, sonication protocol, thermal, head model, pseudo-CT, tissue layers and properties, jwave, BabelBrain | timing diagram, TUS calculator exchange |
   | Study | study name, targeting & analysis | study manifest editor |
   | Run | pipeline stages, where to run, HPC | run log, YAML review |

   - Each field has a readable name. Hovering over it or clicking **?** shows the full help from `config_default.yaml`, the config key and the default.
   - Changed values are marked ● (the tab label counts them), fields filled from the equipment library are marked *from library*, and ↺ resets a field.
   - Fields appear only when they apply to the current modes (e.g. the Localite settings only with `placement.mode: localite`).
3. **Advanced options.** Each tab hides the settings most studies leave at their defaults, and says how many are hidden. Tick **Show advanced options** on that tab to see them. An advanced setting that holds a value is always shown. What counts as advanced:
   - *Setup:* extra MATLAB paths, FSL/ANTs locations, HPC modules for pseudo-CT, the individual `io.save_*` switches, external acoustic NIfTI.
   - *Transducer:* element geometry (diameters, curvature, exit-plane distance), matrix layout, bowl-referenced depths, focal-depth reference convention. Amplitude, phases and steering stay visible.
   - *Placement:* Localite ingestion details (export subfolder, pulses per series, tracker offsets), the heuristic search criteria, PlanTUS weights and outputs, TUS_entry options, neuronav naming rules.
   - *Simulation:* precision and interactive mode, everything on the grid except its resolution, protocol time steps and cool-down, the thermal model, and the whole head-model, segmentation, tissue-layer and tissue-property sections. Debug outputs and the pseudo-CT section stay visible so the toggles that turn them on are never hidden.
   - *Study:* the target area radius.
   - *Run:* intermediate pipeline stages (grid, medium, source, analyses), report style, cluster profile, reservation, job naming and waiting.
4. **Transducer.**
   - Pick a device and driving system from `config/equipment`, or enter custom geometry.
   - The panel shows the focal range (`min_foc`/`max_foc`) and which calibration entry PRESTUS will use.
   - With several transducers, transducer 2 and later inherit transducer 1's values unless set.
   - The cross-section on the right draws the elements, exit plane, geometric focus, the focal depth you set (✕) and the device's steering range.
5. **Placement viewer** (manual placement). The Placement tab has two sub-tabs: **Settings** (the
   `placement`/Localite cards) and **Viewer**, which gives the viewer the full tab width instead of
   sharing it with the settings — the same nested-tab pattern the Run tab uses for Run log / Review
   YAML / Results.
   - Axial, coronal and sagittal T1 slices with the segmentation overlay.
   - Click to set the focus or the transducer (bowl apex), or switch to **Click to browse** to just
     move the crosshair without writing anything (scrolling always just browses). The bowl is drawn
     from the device geometry.
   - **+**/**−** by the transducer selector add or remove a transducer without switching to the
     Transducer tab; with several transducers, transducer 2 and later inherit transducer 1's values
     unless set (see *Transducer* below).
   - **Enter coordinates** types X/Y/Z directly instead of clicking: **RAS mm** or **Voxel** need
     no conversion; **MNI mm** is converted to this subject's space with SimNIBS's
     `mni2subject_coords` CLI tool, the same one-shot external call `mni2subject_coords_LDfix.m`
     uses — it is not run on every redraw, only when you click. This needs SimNIBS (charm) to have
     already been run for this subject: the `m2m_sub-XXX` folder (under `path.seg`/`path.anat`)
     must exist, for *both* the nonlinear and affine transform — setting `startup.simnibs_bin_path`
     alone is not enough, and the viewer says so (distinctly from "SimNIBS's bin folder isn't set")
     when the m2m folder is missing. As in
     [MNI placement](doc_placement_mni.md#how-a-position-is-resolved), **Set focus** uses the
     nonlinear warp (most accurate in-brain) and **Set transducer** the affine (`12dof`,
     well-defined outside the brain too); if the m2m folder exists but its nonlinear warp
     specifically failed or was never produced, **Set focus** automatically retries with the affine
     and says so — it does **not** help when SimNIBS hasn't been run at all. **Go to** just moves
     the crosshair (using the focus's transform, or the transducer's if 'Click sets transducer' is
     selected) without writing anything. **Set focus** / **Set transducer** only work with
     `placement.manual.coordinate_space` set to `voxel` or `ras` — a study using `coordinate_space:
     mni` still writes MNI values
     directly (PRESTUS converts them at placement-stage runtime) and is not shown in the viewer;
     type an MNI point and use **Go to** here just to look at the matching subject-space location.
   - *Snap to scalp* keeps the beam axis and puts the exit plane a chosen gap outside the scalp. *Match focal distance* then sets `focal_distance_ep` to fit.
   - Trajectories can be imported from Brainsight (NIfTI coordinates) or 3D Slicer, using BabelBrain's conventions.
   - Positions are written in the space set by `placement.manual.coordinate_space` (`voxel` or `ras`). See [doc_placement.md](doc_placement.md#manual-coordinate-selection).
   - Other placement modes are shown read-only: the Localite pose from the `neuronav_gui` outputs, the heuristic choice from `transducer_positioning`, and for every mode the positions the last run's placement stage resolved (`<path.sim>/sub-XXX/sub-XXX_placement<affix>.json`, written by `write_placement_positions`). To get them without simulating, run the pipeline once with only `modules.run_transducer_placement` enabled.
   - With `placement.mode: localite`, the Placement card offers *Open neuronav_gui* and *Import into study…*.
6. **Timing diagram and TUS calculator** (Simulation tab, *Sonication protocol*). See the next section.
7. **Issues panel.** Lists what PRESTUS would otherwise report only at load or run time: unknown or misspelt keys (with suggestions), deprecated keys, wrong types and values, backend combinations and frequency mismatches, missing inputs and missing calibration entries, a `target_isppa_wcm2` that cannot be met (no library calibration and `modules.run_water_baseline: 0`) or that adds a free-water simulation (`auto`), and a project `config_default.yaml` that lacks keys added to PRESTUS since it was copied. Clicking an issue opens its tab and field, showing advanced options if needed. The header shows the error and warning counts.
8. **Run tab.**
   - *Validate in MATLAB* only runs `load_parameters`.
   - *Run locally* writes `run_<study>.m` next to the config and runs it with `matlab -batch`, streaming the log. It sets `platform: matlab` for that run only. You can also run the launcher yourself with `run('run_<study>.m')`.
   - *Submit to cluster* (when `sbatch` or `qsub` is on this computer, e.g. on a cluster login node) writes `submit_<study>.m`, which keeps the config's `platform` and `hpc.*` settings, and runs it: PRESTUS submits as usual (`prestus_pipeline_start`; chained pipelines submit several jobs). The Nexus reads the job IDs from the log and lists them under *Cluster jobs*, refreshing their state every 30 s (`squeue`/`sacct`, or `qstat`) until they finish; *Cancel selected* runs `scancel`/`qdel`. The *Run profile* presets fill in the `hpc.*` fields.
   - *Review YAML* shows the exact YAML that will be written.
9. **Study tab.** Runs many subjects with per-run settings. See [Group analysis](doc_group.md#study-setup).
   - Set the manifest file (`study.manifest`, relative to the config) and edit the table: one row per run, with `subject_id`, `condition`, `include`, `config` and any config key as a column (*Add column…* offers the known keys).
   - *Import from Localite…* reads the `neuronav_gui` outputs under `path.localite_post` (`neuronav_ingest_summary.csv` and the per-session `InstrumentMarker_*.json`). It makes one run per target, or one run per session with the targets assigned to transducers 1, 2, …, with a preview before adding. *Open neuronav_gui* starts the MATLAB ingestion GUI with this config. See [Group analysis](doc_group.md#from-localite-ingestion).
   - The summary under the table checks the manifest as you type: subject IDs, `include` values, unknown keys (with suggestions), and rows that would write to the same outputs.
   - *Check all rows in MATLAB* runs `prestus_study_start(config, 'DryRun', true)`. *Run study locally* writes `run_study_<study>.m` and runs it. *Build group report* runs `prestus_group_report_start(config)` (group maps + report). All three stream to the Run tab.

## Timing diagram and TUS calculator

The *Sonication protocol* card draws the protocol from `timing.*`: one pulse (with its ramp), one pulse train, and the train repetitions. It labels PD, PRI/PRF, PTD, PTRI and PTRD with the pulse count, duty cycles, the on-time, and the Ipa/Itp factor of the ramp (`compute_ramp_factor`). It also lists inconsistent settings (e.g. a pulse longer than its repetition interval).

The [TUS calculator](https://www.socsci.ru.nl/fusinitiative/tuscalculator/) (Radboud University, MIT licence) estimates derived quantities (derated intensities, MI, energy deposition, temperature rise) for a protocol. The Nexus does not copy its calculations. It exchanges settings files with it instead:

- **Export for TUS calculator…** writes a JSON file for the calculator's *Upload settings*: pulse duration, ramp duration and shape, PRI, PTD, PTRI, PTRD (in ms), whether the train and repetition levels are enabled, and (when set) the frequency, target Isppa, aperture and focal depth.
- **Import from TUS calculator…** reads a settings file downloaded from the calculator (JSON or YAML). It shows the values it will set (in seconds), and optionally also sets the frequency and target Isppa.
- **Open TUS calculator ↗** opens the web app.

Both use the same ramp convention: the ramp duration is the rise (= fall) time at each end of a pulse (`timing.pulse_ramp_dur`, `timing.pulse_ramp_shape` `linear` or `tukey`). PRESTUS models pulse ramps only, so the calculator's pulse-train and repetition ramps are exported as 0 and reported when an imported file sets them. A level disabled in the calculator sets the corresponding PRESTUS durations to 0.

## MATLAB tools

The **Tools** menu opens PRESTUS's standalone MATLAB GUIs in an interactive MATLAB session. They need the MATLAB desktop, so they are not run through `-batch`:

- **Localite marker ingest** (`neuronav_gui`, see [doc_gui.md](doc_gui.md#localite-ingest-gui)). The study config is saved first and passed to the GUI, which reads `path.localite_*`, `placement.localite.*` and `neuronav.coil_map` from it. After ingestion, set `placement.mode: localite` in the Nexus, or import the ingested targets into a study manifest (Study tab).
- **Transducer calibration** (`calibration_gui`). Fitted library entries are saved to `transducer.library_path` when that is set. The Nexus's transducer panel then reports them as the calibration entry PRESTUS will use.

## Project defaults

When the study folder holds a copy of `config_default.yaml` (from *File ▸ New study…* or `prestus_config_init`), the Nexus merges it onto the toolbox default exactly as `load_parameters` does: the copy's values win, keys added to PRESTUS since the copy was made keep their toolbox defaults. The *Review YAML* page names both files. **Tools ▸ Update project defaults…** replaces the copy with the toolbox default (keeping a dated backup), optionally after writing the values the copy pinned differently into the study config so the simulation settings stay the same; `prestus_config_upgrade(folder)` does the same from MATLAB.

## Standalone app

`nexus/packaging/prestus_nexus.spec` builds a standalone app with PyInstaller (`pip install -e ".[viewer]" pyinstaller`, then `pyinstaller packaging/prestus_nexus.spec` in `nexus/`). The app does not contain PRESTUS: on first start it asks for the PRESTUS folder (or use `--prestus-root` / `PRESTUS_ROOT`) and remembers it.

## Compatibility with MATLAB's YAML loader

PRESTUS reads configs with `yaml.loadFile` (SnakeYAML, YAML 1.1). The Nexus:

- quotes all strings (unquoted `yes`/`no`/`04:00:00` would change type);
- reads floats such as `2.17e-3` the way SnakeYAML does;
- writes every element of a transducer list with the same top-level keys, which `yaml.loadFile` requires.

A round-trip test writes configs with the Nexus, loads them with `load_parameters`, and compares every value. It covers both a serial-based and a fully inline transducer.

## Maintenance

The schema is derived from `config_default.yaml` at startup. `nexus/prestus_nexus/overlay.yaml` adds only widget kinds, allowed values, units, visibility rules, keys the code reads that the default omits, and deprecated keys. `nexus/tests/test_schema.py` fails when the overlay no longer matches the default config or the code. See `nexus/README.md` for running the tests.
