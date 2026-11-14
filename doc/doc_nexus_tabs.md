# Nexus Tabs and Pages

The QuickStart tab sets up a first simulation and leaves every other setting at its PRESTUS default; the other tabs hold all of those settings. This page describes each tab and page of the Nexus, what its controls do and which config keys they write, so you can find a setting, read a check, and follow a run to its results.

On this page:

- [Tab Overview](#tab-overview): what each tab holds and how fields mark changes and problems.
- [Menus and Shortcuts](#menus-and-shortcuts): undo, help, tools and keyboard shortcuts.
- [Segmentation Step](#segmentation-step): running SimNIBS charm on its own, before anything else.
- [QuickStart Tab](#quickstart-tab): the four steps and the config keys they write.
- [Advanced Options](#advanced-options): settings hidden by default and locked cards.
- [Transducer Tab](#transducer-tab): choosing a device or entering custom geometry.
- [Placement Viewer](#placement-viewer): viewing and setting the transducer position on the T1.
- [Issues Panel](#issues-panel): the config problems found before PRESTUS loads or runs the study.
- [Run Tab](#run-tab): running locally or submitting to a cluster, and the run log.
- [Results Tab](#results-tab): browsing the runs in a results folder and starting follow-up runs.
- [Study Tab](#study-tab): the run table for many subjects or conditions.
- [Timing Diagram and TUS Calculator](#timing-diagram-and-tus-calculator): drawing the sonication protocol and exchanging it with the TUS calculator.

!!! note "At a glance"
    - **When you need it:** You change a setting beyond the QuickStart tab's four steps, run many subjects, or look up what a control writes to the config.
    - **Before you start:** The Nexus installed and started ([PRESTUS Nexus (GUI)](doc_nexus.md)); [A First Simulation](doc_nexus.md#a-first-simulation) walks through the main tabs once.

## Tab Overview

The window has one tab per area; each shows its config sections as cards:

| Tab | Sections | Also on the tab |
|---|---|---|
| QuickStart | (the four steps below) | |
| Setup | subject, paths, environment, input/output | |
| Transducer | transducer, geometry, drive and steering | device panel with the cross-section under it; 'Calibration' sub-tab (calibration coverage and the calibration conventions, e.g. the focal-depth reference) |
| Placement | placement (incl. the Localite folders), Localite ingestion | sub-tabs *Settings*, *Viewer*, *Beam axis*, *3D*, *Transducers* (with 2+ transducers), *Localite*, *PlanTUS*, *TUS entry*; Localite shortcuts |
| Medium | the medium (a row above the sub-tabs); 'Head model': segmentation, head model, tissue layers; 'Acoustic properties': tissue properties; 'Uncertainty': uncertainty analysis; 'Skull CT': CT-based skull mapping (measured CT or UTE pseudo-CT) | SimNIBS environment (installs found on the computer, and with SSH submission on the cluster: *Scan cluster*), its version, the existing segmentation with its charm report, and the re-run switch; tissue property table; uncertainty variants card with the ITRUSST bounds plot; mapping equations and the CT view |
| Sonication | sonication protocol, multi-transducer coupling (with more than one transducer), thermal model | timing diagram, TUS calculator exchange; coherent vs. async field prediction |
| Solver | acoustic and thermal solvers, their compute backend, grid, [jwave](https://github.com/ucl-bug/jwave), [BabelBrain](https://github.com/ProteusMRIgHIFU/BabelBrain), elastic skull, time reversal | 'Compute backend' card: the device setting of each solver that runs (k-Wave, also for heating and time reversal; jwave; BabelBrain), several when solvers are mixed |
| Study | study, targeting analysis | study manifest editor |
| Run | sub-tabs *Stages* (pipeline stages, debug outputs), *Where to run* (run profile, platform, HPC, pipeline options), *Plan* (what the run will do), *Review YAML*, *Launch*, *Jobs* | run buttons and log, job list |
| Results | (none) | the runs in the results folder, follow-up runs |

??? note "Fields, change marks and mode availability"
    - Each field has a readable name in words, with abbreviations spelled out or in their usual capitals (HU, CT) and intensities in ITRUSST notation (I<sub>SPPA</sub>, I<sub>SPTA</sub>; plain ISPPA in drop-downs, file names and copied text). Hovering over it or clicking **?** shows the full help from `config_default.yaml`, with the choices listed one per line, the config key, the default and a link to its documentation page. Choices that a picture explains (placement method, medium, coupling, array type, steering, focal-depth reference and others) show a small schematic before their name in the drop-down and in the help; the coupling choice also shows two larger figures under it. **↗** in a card's header opens the card's page.
    - Changed values are marked ●, fields filled from the equipment library are marked *from library*, and ↺ resets a field.
    - A tab's label shows a red (errors) or amber (warnings only) dot with the number of problems in its fields. Its tooltip also counts the changed values.
    - Fields appear only when they apply to the current modes (e.g. the Localite settings only with `placement.mode: localite`).
    - *Placement method* (the row at the top of the Placement tab) marks each mode for the configured subject. **Unavailable** (greyed out; the tooltip names the missing input) means PRESTUS cannot produce an input the mode needs: no T1, no SimNIBS segmentation (PlanTUS: also the mesh `sub-XXX.msh`) and no `startup.simnibs_bin_path` to run charm, or no recorded Localite poses. **Needs setup** means charm will run first or the mode's tool paths are still empty. Tool paths never grey a mode out, since they are only shown once it is selected. The mode in use always stays selectable, and why it cannot run yet is written next to it.

??? note "Exporting figures"
    Each figure has a download button in its top-right corner. Plots export as SVG, PDF or PNG at 300 or 600 dpi; drawn figures (timing diagram, transducer cross-section, the R² strip of the Calibration sub-tab) as SVG, PDF or PNG at 3× or 6× the screen size; 3D views as PNG only. **Copy image** puts the figure on the clipboard. Figures export in the light theme on white unless **Keep current theme** is ticked. The file name defaults to `sub-XXX_<view>` plus the output affix.

## Menus and Shortcuts

- **File:** new, open, open recent, save, save as, reload from disk, quit.
- **Edit ▸ Undo / Redo** (Ctrl+Z / Ctrl+Shift+Z) step through value changes, resets, and added or removed transducers. Typing in one field counts as one step. *Reload from disk* clears the history.
- **View:** *Dark mode*, and *Issues…* (Ctrl+Shift+I) to open the [Issues panel](#issues-panel).
- **Tools:** *Localite*, *PlanTUS*, *TUS entry* and *Calibration* open those sub-tabs; `neuronav_gui` and `calibration_gui` open in MATLAB; *Update project defaults…* ([Project defaults](doc_nexus.md#project-defaults)); *Usage statistics…* (below).
- **Help:** the Nexus overview, this page, *Docs for this tab* (F1, the page of the tab and sub-tab shown), the parameters reference and troubleshooting. *Use local documentation files* opens the Markdown files in `doc/` instead of the published pages, e.g. offline.
- **Keyboard:** Ctrl+1 to Ctrl+9 switch tabs, Ctrl+Shift+V validates in MATLAB and Ctrl+R runs locally.

??? note "Usage statistics"
    The Nexus runs MATLAB with `matlab -batch`, where PRESTUS never asks whether to share [usage statistics](doc_telemetry.md). **Tools ▸ Usage statistics…** shows the recorded decision and records a new one (*Allow*, *Decline* or *Forget my decision*) in `~/.prestus/telemetry.json`, the file MATLAB reads. The Nexus itself sends nothing.

??? note "When something goes wrong"
    - An unexpected error shows a dialog with the details instead of closing the Nexus; the traceback is also written to `nexus_errors.log` in the [cache folder](doc_nexus.md#installation-options).
    - Saving writes to a temporary file first, so a full or read-only share leaves the old config intact and shows an error.
    - The MATLAB session must start within 180 s, and a request may take up to 600 s; after a time-out the session restarts on the next request.

## Segmentation Step

The SimNIBS segmentation (`m2m_sub-XXX`) feeds the layered head model, heuristic, MNI, PlanTUS and TUS_entry placement, the Viewer's scalp and beam checks and the 3D head. charm takes an hour or more, so it can run on its own, before anything else.

- The Medium tab's *Segmentation* card shows whether it is found, incomplete (tissue labels or MNI warps missing: charm unfinished or failed) or missing, and what in the config needs it. **Run segmentation** runs only charm for this subject (`modules.segmentation_only: 1` with `simulation.medium: layered`): as a cluster job where the Run tab submits there, else in MATLAB here. The card updates when the job has written the segmentation.
- **Segment again…** asks first: charm deletes the m2m folder and segments again (`io.overwrite_simnibs: 1` for that run only), and placements, head models and MNI conversions made from the old segmentation no longer match it.
- The Placement tab shows the same row at the top while the placement method needs a segmentation that is not complete; QuickStart offers *Run segmentation now* in its first step.

## QuickStart Tab

The Nexus opens on *QuickStart*, which sets up a first simulation in four steps. Every other setting keeps its PRESTUS default. Each step shows *To do* until it is filled in, then *Done* or the number of problems (hover for the messages). The QuickStart page writes the same config keys as the other tabs, so the two never disagree.

1. **Head image:** choose the T1-weighted NIfTI, or drop it onto the window. From the file name it sets `path.anat`, `path.t1_pattern` (as `sub-%1$03d_…` when the name holds `sub-<n>`) and `subject_id`. It also finds a T2 next to the T1 (`T1w`→`T2w`) and a folder with `m2m_sub-XXX` in the anatomy folder, one of its subfolders or one of its siblings (`path.seg`). If the segmentation is missing, it looks for SimNIBS 4 (`charm` on the `PATH`, then the newest install in conda environments and the installers' default folders on macOS, Windows and Linux) to set `startup.simnibs_bin_path`, and *Run segmentation now* runs only charm ([Segmentation step](#segmentation-step)). The results folder defaults to `prestus_sim` next to the study config, or next to the anatomy folder if the config is unsaved, never inside the anatomy folder.
2. **Transducer:** the device and driving system from the equipment library, the focal depth setting (with the device's range) and the free-water I<sub>SPPA</sub>. The calibration entry PRESTUS will use is shown below.
3. **Target:** an MNI target, either alone (`placement.mode: heuristic`, which needs a target name) or with a transducer position (`mni`). A study using another mode is pointed to the Placement tab.
4. **Check and run:** *Run on* (`platform`) and the errors still to fix, each a link to its field. *Continue to Run* saves the study and opens the Run tab.

Each step's *Learn more* opens the page that explains it. *Already have results? Start a follow-up run* opens the Results tab ([Follow-up Runs](#follow-up-runs)).

## Advanced Options

Each tab hides the settings most studies leave at their defaults, and its **Advanced** toggle says how many, e.g. *Advanced (3)*. Settings the chosen setup does not read stay hidden too: with a `water` medium, the anatomy, SimNIBS, MNI-saving, defacing, targeting and free-water-comparison settings; settings the pipelines fill in themselves (such as the `io.adopted_*` start maps of a sequential follow-up) until they hold a value. Every card has an icon and a line saying what it decides; every sub-tab has an icon.

??? note "Advanced controls and locked cards"
    - The advanced icon in a card's header, with the number of hidden settings next to it in grey, shows that card's advanced settings; it is highlighted while they are shown, and a second click folds them again. **Advanced** at the top of the tab shows them on every card of the tab, including the fully advanced sections.
    - An advanced setting that holds a value is always shown. Jumping to an issue opens the card it is on.
    - **Locked cards.** *Tissue layers* (Medium tab, *Head model*) is greyed and read-only until unlocked, as below. The *Tissue properties* card (Medium tab, *Acoustic properties*) shows `medium_properties` as one table, tissues as columns (each headed by a glyph of its head layer) and properties as rows, greyed and read-only. With advanced options shown, **🔒 Unlock** asks for confirmation and makes the cells editable for the session; **Lock** reverts. Changed values are bold, and emptying a cell resets it. Hovering a cell shows its default and the uncertainty-variant values. Values from a loaded config still apply.
    - **Which tissues.** A `water` medium shows the water column only and hides the *Head model* and *Skull CT* sub-tabs. A `phantom` shows every tissue, with a note that only the tissues its layers assign are modelled, and hides *Skull CT*.
    - **Thermal rows and absorption.** Absorption follows attenuation; the absorption rows and those from thermal conductivity on are shaded: they matter for the thermal simulation only. A tissue sets either the absorption fraction or the absorption coefficient; the table shows the other in (brackets) as implied. Both follow the attenuation power law (y = `alpha_power`). Typing a value into an implied cell switches the tissue to that field.
    - **ITRUSST reference.** **Open ITRUSST tissue properties ↗** in the card header opens the [ITRUSST tissue-property overview](https://itrusst.github.io/itrusst_properties/) with the reference values and bounds.
    - **Uncertainty sub-tab.** The *Uncertainty analysis* card turns on `simulation.uncertainty`; the property card's note links to it while it is off. With it on, the card shows the four variant medium files, followed by the *Uncertainty variants* card: the brain and skull values this config uses (line) plotted against the bound values of the variants, one property per column (Temp-H/L bound brain heating, Trans-H/L brain intensity, acoustic properties only), then the default, Temp-H, Temp-L, Trans-H and Trans-L values. Greyed values are not in the variant file and use the default. The sub-tab is hidden for a water medium.
    - **Skull CT sub-tab.** Holds the CT mapping switch, the measured CT file pattern (`path.ct_pattern`; empty = pseudo-CT from the UTE image) and the mapping settings. With the BabelBrain backend it also holds BabelBrain's own HU mapping and internal-air settings (`babelbrain.ct_mapping`, `babelbrain.air_*`); the *BabelBrain options* card links here. *From image to skull properties* lists the steps from the input image to the skull properties, each marked as run, reused or skipped: for a pseudo-CT, the UTE image, bias correction and normalisation, the UTE→pseudo-HU line (plotted, read from `pct_create_pseudoCT.sh` or the subject's fit) and the pseudo-CT; for a measured CT, the CT and its resampling. It then plots the HU→density, sound speed and attenuation curves of the selected mappings, evaluated by MATLAB (`pct_mapping_curves`) when MATLAB runs or on *Draw the mapping curves*. The equations of the selected mappings (from [Skull Property Mappings](doc_pseudoCT_property-mappings.md)) fold out below. It shows the measured CT or `pseudoCT.nii.gz` over the T1 and the in-skull distributions of the mapped skull properties against the band of the uncertainty bound variants' skull values, with the fraction of skull voxels outside it. The distributions come from the default run when it exists; otherwise **Predict skull properties** predicts them without simulating, labelled as a prediction. See [Pseudo-CT runs](doc_uncertainty.md#pseudo-ct-runs).

??? note "What counts as advanced"
    - *Setup:* extra MATLAB paths, FSL/ANTs locations, HPC modules for pseudo-CT, the individual `io.save_*` switches, external acoustic NIfTI.
    - *Transducer:* element geometry (diameters, curvature, exit-plane distance), matrix layout, bowl-referenced depths, focal-depth reference convention. Amplitude, phases and steering stay visible.
    - *Placement:* Localite ingestion details (export subfolder, pulses per series, tracker offsets), the heuristic search criteria, PlanTUS weights and outputs, [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) options, neuronav naming rules.
    - *Medium:* the *Tissue layers* card (the SimNIBS charm labels each tissue takes; locked like the property table), unlocking the tissue property table, charm debug output, `--noneck` and custom charm settings, and the tuning parameters of *Segmentation refinement*, *Smoothing* (thresholds and width, smoothing of property maps) and *Grid refinement* (padding), and the rubber-wrap display.
    - *Sonication:* protocol time steps and cool-down, the thermal model.
    - *Solver:* precision and interactive mode, everything on the grid except its resolution. Debug outputs stay visible.
    - *Study:* the target area radius.
    - *Run:* intermediate pipeline stages (grid, medium, source, analyses), report style, cluster profile, reservation, job naming and waiting.


??? developer "Where locked cards are defined"
    Locked cards are set by `locked: true` under `sections:` in `nexus/prestus_nexus/overlay.yaml`. The `tuning:` list there (which keys affect speed, accuracy, memory or disk) is kept for a later *Tune for* control and is not shown.

## Transducer Tab

Pick a device and driving system from `config/equipment`, or enter custom geometry.

??? note "Device panel, drawing and Calibration sub-tab"
    - The panel shows the focal range (`min_foc`/`max_foc`) and which calibration entry PRESTUS will use.
    - With several transducers, transducer 2 and later inherit transducer 1's values unless set.
    - Once a transducer is chosen, the cross-section under the device panel draws the elements, exit plane, geometric focus, the focal depth you set (✕) and the device's steering range.
    - *Geometry* (array type, frequency, element geometry) starts collapsed when the equipment library supplies it; open it to check or override a value. *Drive and steering* holds the values set per run: the target free-water I<sub>SPPA</sub>, element amplitudes and phases, and matrix steering. *Sweep…* next to the target I<sub>SPPA</sub> edits a list of targets, which runs the [multi-ISPPA pipeline](doc_multiintensity.md). Both cards follow the transducer chosen in the *Transducer* card.
    - The *Calibration* sub-tab shows what the library entry covers (fit quality, drive and phase corrections per depth) and compares, fits and saves library entries for the config's transducers ([Calibration page](doc_nexus_calibration-localite.md#calibration-page)).

## Placement Viewer

The **Placement method** row at the top of the Placement tab chooses `placement.mode`, as the Medium row does on the Medium tab. The tab opens on **Settings** (the method's paths and options, and per transducer the transducer and target positions and the focal depth setting); its other sub-tabs are **Viewer** (manual placement), **Beam axis** (the medium and free-water profile along the beam), **3D** (head and transducers in 3D), **Transducers** (with two or more transducers), **Localite**, **PlanTUS** and **TUS entry**. The sub-tab of the method in use is labelled in bold (the Viewer for manual placement). **Localite**, **PlanTUS** and **TUS entry** carry the icon of their method, as in the *Placement method* drop-down.

??? note "Using the viewer"
    - Axial, coronal and sagittal T1 slices with the segmentation overlay; a line under them names the markers. Without a T1, the viewer shows what is missing and a button to the field.
    - **Click sets** (*Focus* · *Transducer* · *Browse*) chooses what a click on a slice does: set the focus or the transducer (bowl apex), or with *Browse* move the crosshair without writing anything (scrolling always browses). The bowl is drawn from the device geometry.
    - **Undo** / **Redo** (the arrows next to *Click sets*, or Ctrl+Z / Ctrl+Shift+Z while the viewer has focus) step back and forward through placement changes: positions, the focal depth setting, the placement mode, the acoustic FOV, and added or removed transducers. Edits made on other tabs in between are kept.
    - With several transducers, every bowl, beam axis and acoustic FOV box is drawn; the ones not selected are dimmed. Each FOV box is resolved from its own transducer (`'auto'` diameter and length), as PRESTUS does: coherent coupling simulates their union, async coupling one box per transducer ([details](doc_multitransducer.md#partial-acoustic-fov)).
    - **+**/**−** by the transducer selector add or remove a transducer without switching to the Transducer tab. Hovering over the selector or one of its entries shows the transducer's device (`transducer.serial`, with its equipment-library name).
    - **Coordinates…** opens a row to type X/Y/Z directly instead of clicking. **RAS mm** or **Voxel** need no conversion. **MNI mm** is converted to this subject's space with SimNIBS when you click, which needs the subject's `m2m_sub-XXX` folder; **Set focus** uses the nonlinear warp and **Set transducer** the affine, as in [MNI placement](doc_placement_mni.md#how-a-position-is-resolved). **Go to** only moves the crosshair.
    - With manual placement, the readout's red flags are also listed in the Issues panel, for every transducer: a target more than 2 mm from the focal depth setting, an exit plane (centre or rim) inside the scalp, and beam rays through an eye. They are checked once the viewer (or the Run tab) has loaded the T1 and segmentation; clicking one opens the viewer on that transducer.
    - *Snap to scalp* keeps the beam axis and puts the exit plane a chosen gap outside the scalp. *Transducer at focal depth* (under **More**, with the gap) instead puts the transducer where the focal depth setting lands on the target. **View** holds the overlay opacity, the FOV size and *Centre on focus*.
    - **Beam:** with a segmentation, the readout adds skull on the beam (mean ± SD over 37 rays from the aperture to the focus; the central ray), the mean skull incidence (amber above 20°), the gap between the exit-plane rim and the scalp along the beam (mean ± SD; red when the minimum is negative, i.e. the transducer cuts into the head; a large SD means uneven contact) and the share of rays that pass an eye. These are the measures of the placement comparison (`placement.beam_qc`).
    - The readout keeps two depths apart: the **target depth** (exit plane → `focus_pos`, from the placement) and the **focal depth setting** (`focal_distance_ep`, the depth set on the equipment). For annular arrays the setting selects the calibrated phases, so the simulated focus lands at the setting even when the target is deeper or shallower; the readout gives the difference (red above 2 mm). *Set focal depth to target* copies the target depth into the setting.
    - Trajectories can be imported from [Brainsight](https://www.rogue-research.com/tms/brainsight-tms/) (NIfTI coordinates) or [3D Slicer](https://www.slicer.org/), using BabelBrain's conventions.
    - Positions are written in the space set by `placement.manual.coordinate_space` (`voxel` or `ras`). See [doc_placement.md](doc_placement.md#manual-coordinate-selection).

??? note "Beam axis"
    - A section through the beam axis of the transducer selected on the Viewer, drawn like the transducer schematic: bowl on the left, depth from the bowl apex to the right, 1:1. It shows the simulation medium: the last layered run's `medium_masks` (`nii/sub-XXX_layered_T1w<affix>_medium_masks.nii.gz`), or before a run the SimNIBS labels mapped through `layers` (without PRESTUS's smoothing). On it are the bowl, the exit plane, the scalp on the axis with its distance to the exit plane, the target, the focal depth setting, the device's steering range and the acoustic FOV; for heuristic, MNI, PlanTUS and TUS-entry placement, a dashed line marks where `placement.skin_gap_mm` puts the exit plane. The key beside it names every marker and line and the medium layers.
    - **⟲ ⟳** or scrolling turn the section about the axis in 15° steps (0°: lateral axis towards superior); a click moves the Viewer's crosshair.
    - Under the section, on the same depth axis, is the **free-water profile**. After a run it is read from the run's axial profile (`img/sub-XXX_<medium>_axial_profile<affix>.csv`): the O'Neil free-water intensity of the settings the run used, the k-Wave free-water baseline if it ran, and the simulated intensity in the head. Before a run it is the device's library model (or the geometric model without a calibration entry) at the focal depth setting, predicted in MATLAB as the Calibration page does: with MATLAB running it appears by itself, otherwise use *Compute free-water profile* (this starts MATLAB; annular library devices only). A ▼ marks the free-water focus (the peak for a run, the FLHM centre for the prediction).
    - **Profile: Scaled · Absolute** scales each profile to its own maximum, to compare where they peak, or shows W/cm², to compare how much reaches the target. *Absolute* needs a run: the prediction has no absolute scale.
    - The readout gives the distance from the exit plane to the scalp and to the first tissue of the medium (on the axis, and the smallest across the aperture in this section; ✗ when negative: tissue sits behind the exit plane), the free-water focus, and after a run the simulated peak in brain. It says when the run used another placement or focal depth than the config now, and what is missing (segmentation, medium, profile). Readout values here and on the Viewer carry ✓, ▲ or ✗ as well as their colour. Without a T1, a placement or a transducer geometry, the tab says which and offers the field.

??? note "3D"
    - The head (scalp and brain, from the SimNIBS segmentation) with every transducer as it is simulated: one curved ring per element on the bowl, the exit plane, the beam axis, the target (focus position, a sphere with a crosshair), the geometric focus and where the focal depth setting puts the focus. The transducer selected on the Viewer is in colour, others grey.
    - With PlanTUS or TUS_entry placement, the run's target mask is drawn as well (`placement.tus_entry.target_mask`, or the mask the run wrote to its placement run folder); **Mask** shows or hides it.
    - Each transducer is drawn in a cylindrical housing, as [PlanTUS](https://github.com/mlueckel/PlanTUS) and BabelBrain draw it. Its front face lies in the exit plane; in practice a gel pad or membrane can leave a gap to the scalp. The size comes from the device's equipment entry (`housing:` in `config/equipment/<SERIAL>.yaml`, or its model's). Without a size, the cylinder has the aperture's diameter and reaches 15 mm behind the bowl apex, and the legend says so.
    - A slider next to **Scalp**, **Brain**, **Mask** and **Housing** sets each surface's opacity. The brain is see-through by default, so the target inside it shows; the housing is opaque and hides the element rings until you lower its slider.
    - The face is removed by default: the scalp in front of and below the eyes (SimNIBS eye label; the brain's front without one) is not drawn, so the view does not identify the participant. **Face** draws it.
    - Drag to turn, scroll to zoom; **Reset view** looks at the target from the transducer's side. **Scalp**, **Face**, **Brain**, **Mask** and **Housing** show or hide the surfaces.
    - It needs a display with graphics support (Qt Quick 3D); without one, the tab says so.

??? note "Transducers"
    - Shown with two or more transducers. A read-only table lists each transducer's device, positions, focal depth setting, target ISPPA and tilt from the grid z axis; *Edit in placement viewer* opens the Viewer on the selected one.
    - *Grid orientation* sets `grid.mode` (`focal_axis`, `transducer_axis` or `ras_plus`; [Grid Alignment](doc_multitransducer.md#grid-alignment)) and gives the resulting grid z axis and each beam's tilt from it. It warns when `focal_axis` falls back to transducer 1's axis because the beams are too spread (`grid.focal_axis_min_norm`, *Mean-axis spread limit*), and when a beam tilts 20° or more on a grid with fewer than 6 points per wavelength in water.
    - *Coupling* shows `simulation.transducer_coupling` and, for async coupling, `simulation.async_grid_mode`. *Timing* explains what the single `timing` block means for several transducers.
    - The Viewer's **View** menu draws the grid z axis (*Grid z axis (grid.mode)*), on by default with two or more transducers.

??? note "Other placement modes"
    - Other placement modes are shown read-only: the Localite pose from the `neuronav_gui` outputs, the heuristic choice from `transducer_positioning`, and for every mode the positions the last run's placement stage resolved (`<path.sim>/sub-XXX/sub-XXX_placement<affix>.json`). To get them without simulating, run the pipeline once with only `modules.run_transducer_placement` enabled.
    - With `placement.mode: localite`, the Placement card offers *Localite page* and *Import into study…*; the *Localite* sub-tab detects, names and writes the recorded series of a session ([Localite page](doc_nexus_calibration-localite.md#localite-page)).
    - With `placement.mode: plantus`, the Placement card runs only the [PlanTUS placement](doc_placement_plantus.md#running-plantus): *Compute maps* scores every skin position and takes the composite-score optimum, *Choose in Workbench…* picks the position on PlanTUS's maps in [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench), and *PlanTUS page* opens the sub-tab below. The run shows on the Run tab, and the Viewer opens on the result.

??? note "PlanTUS"
    - The skin of the last PlanTUS run for this subject and target (`plantus_metrics.mat` in the placement run folder), coloured by the composite score or, under **Map**, one criterion: bright is better, grey is excluded (avoidance area, or no value). Markers show the target, the composite optimum and the vertex placed.
    - A click (left or right) on the skin chooses a position. The panel lists its criteria, composite score and rank, or why it is not eligible; the view draws the beam axis to the target and the transducer where PlanTUS puts it. **Go to the optimum** chooses the best composite score.
    - **Place here** sets `placement.plantus.vertex` and places the transducer there on the same maps (Run tab, no recomputation); a position that is not eligible asks first. **Compute maps** runs PlanTUS again and clears the choice.
    - The face is scrambled as on the 3D tab; **Face** draws it. Without maps, the tab offers *Compute maps*; without a display, it says so.
    - *Compute maps* and *Place here* run in MATLAB on this computer, or, where the Run tab submits to a cluster (scheduler or *Submit over SSH*), as a placement-only job; the page reloads the maps when the job has written them. *Choose in Workbench…* needs a display and always runs here.

??? note "Heuristic"
    - With `placement.mode: heuristic`: the candidates of the last heuristic run (`tpars_sub-XXX_<target>.csv` in the placement run folder), plotted as distance to the target against the share of the beam path through excluded tissue, and listed in a table (sortable). Candidates above `criterion_intersection` are grey; ★ marks the heuristic's choice, ◆ `placement.heuristic.candidate`, ○ the candidate selected here.
    - A click in the plot or the table shows the candidate's position and criteria. **Place here** sets `placement.heuristic.candidate` and places it (placement only, reusing the candidate table); **Use the heuristic's choice** clears it; **Compute candidates** runs the heuristic again. Runs go to the cluster where the Run tab submits there, as on the PlanTUS page. Distances are in voxels, as the heuristic computes them.

??? note "TUS entry"
    - With `placement.mode: tus_entry`: the TUS_entry settings, *Run TUS_entry* (placement only; here or as a cluster job, as on the PlanTUS page) and the result from `tus_entry_report.csv` in the placement run folder: scalp entry and target centre (voxels, 1-based), entry–target distance and angle score. *Open in Viewer* shows the transducer it placed ([TUS_entry placement](doc_placement_tus_entry.md)).

??? note "MNI coordinates in the viewer"
    **MNI mm** is converted with SimNIBS when you click. It needs the `m2m_sub-XXX` folder (under `path.seg`/`path.anat`) for both the nonlinear and affine transform. Setting `startup.simnibs_bin_path` alone is not enough, and the viewer says so (distinctly from "SimNIBS's bin folder isn't set") when the m2m folder is missing. **Set focus** uses the nonlinear warp (most accurate in-brain) and **Set transducer** the affine (`12dof`, well-defined outside the brain too). If the m2m folder exists but its nonlinear warp failed or was never produced, **Set focus** retries with the affine and says so; this does not help when SimNIBS has not been run at all. **Go to** moves the crosshair using the focus's transform, or the transducer's if *Click sets transducer* is selected. **Set focus** / **Set transducer** only work with `placement.manual.coordinate_space` set to `voxel` or `ras`. A study using `coordinate_space: mni` writes MNI values directly (PRESTUS converts them at placement-stage runtime) and is not shown in the viewer; type an MNI point and use **Go to** to look at the matching subject-space location.

??? developer "MNI conversion call"
    **MNI mm** is converted with SimNIBS's `mni2subject_coords` CLI tool, the same one-shot external call `mni2subject_coords_LDfix.m` uses; it runs only when you click, not on every redraw.

## Issues Panel

The panel lists config problems before PRESTUS loads or runs the study; open it from the error and warning counts in the header.

??? note "What the panel checks"
    The panel lists what PRESTUS would otherwise report only at load or run time: unknown or misspelt keys (with suggestions), deprecated keys, wrong types and values, backend combinations and frequency mismatches, missing inputs and missing calibration entries, multi-transducer combinations PRESTUS refuses or reduces, heating timing and thermal FOV problems, too few points per shear wavelength for the elastic solver, a `path.acoustic_input`, `io.adopted_*` file or follow-up parent config that is missing or does not match, a `target_isppa_wcm2` that cannot be met (no library calibration and `modules.run_water_baseline: 0`) or that adds a free-water simulation (`auto`), a project `config_default.yaml` that lacks keys added to PRESTUS since it was copied, and (manual placement) the placement viewer's checks. Clicking an issue opens its tab and field, showing advanced options if needed. The panel is closed by default: the header shows the error and warning counts, and clicking them (or **View ▸ Issues…**, Ctrl+Shift+I) opens the list. It closes on Esc, a click outside, or a click on an issue.

## Run Tab

The sub-tabs follow the steps of a run: *Stages* (pipeline stages and debug outputs), *Where to run* (run profile, *Run on* and the scheduler settings), *Plan* (what the run will do), *Review YAML*, *Launch* (the run buttons and log) and *Jobs*.

- The highlighted button follows *Run on*: *Run locally* for `platform: matlab`, *Submit to cluster* for `slurm` or `qsub` (and for `auto` when a scheduler or the SSH host is available). The other one, *Write launcher* and *Save log…* are under **More**. The *Where to run* card shows the scheduler fields for `slurm` or `qsub`, and for `auto` only on a computer with `sbatch`/`qsub` or while *Submit over SSH* is on.
- *Validate in MATLAB* only runs `load_parameters`.
- *Plan* follows the config as you edit it: the pipeline mode (single, async, multi-ISPPA or uncertainty); each stage marked as run, reused (e.g. an existing segmentation, pseudo-CT or measured CT), skipped or stopping the run, with one line on what it does; the jobs the run submits; the files and folders it always writes; the optional outputs (`io.save_*`, `simulation.debug`) split into saved and not saved, each with its location and the flag or default that decides it, and which of them already have files on disk for this subject, medium and output affix; and a resource estimate from earlier runs of the project (`<path.sim>/run_stats.jsonl`): acoustic and thermal simulation time and peak RAM, as median and maximum. See [Resource estimates](doc_hpc.md#resource-estimates-from-earlier-runs).
- *Run locally* writes `run_<study>.m` next to the config and runs it with `matlab -batch`, streaming the log. It sets `platform: matlab` for that run only. You can also run the launcher yourself with `run('run_<study>.m')`.
- *Submit to cluster* (when `sbatch` or `qsub` is on this computer, e.g. on a cluster login node) writes `submit_<study>.m`, which keeps the config's `platform` and `hpc.*` settings, and runs it: PRESTUS submits as usual (chained pipelines submit several jobs). MATLAB exits once the jobs are submitted (unless `hpc.wait_for_job` is set). *Run profile* shows which preset the `platform` and `hpc.*` fields match (or that none does); *Apply* fills them in from another (*Local* sets `platform: matlab`). The submitted jobs appear in the *Jobs* sub-tab (details below).
- *Submit over SSH* submits from a computer without a scheduler, e.g. a laptop with the project drive mounted. Set the login node (*Host*, e.g. `user@mentat005.dccn.nl` or an alias from `~/.ssh/config`), the folder mapping (*Local folder* `/Volumes/myproject` *is on the cluster* `/project/myproject`), the *Setup commands* run before MATLAB (e.g. `module load slurm matlab/R2024a`) and, if the cluster uses its own checkout, the *Remote PRESTUS* folder. The study folder must be under the mapped folder. The settings are kept per user (`~/.config/prestus_nexus/last_session.yaml`).
- *Review YAML* shows the exact YAML that will be written.
- *Pipeline options* (on *Where to run*, folded) holds the sequential follow-ups run after this config with heat carry-over ([Sequential runs](doc_sequential.md)), per-stage time and memory limits for the multi-ISPPA and uncertainty pipelines, and the multi-ISPPA summary report switch. They are not config keys: the Nexus keeps them per user and study and writes them into `run_<study>.m`, `submit_<study>.m` and SSH submissions.
- **Stop** ends MATLAB together with its launcher.

??? note "The Jobs tab"
    The Nexus reads the job IDs from the log and lists them in the *Jobs* sub-tab (brought forward after a submit), refreshing them every 30 s until they finish. *Pop out ↗*, or double-clicking a job, opens the same view in its own window, which can be maximised or moved to a second screen. For [SLURM](https://slurm.schedmd.com/) jobs the list shows the state, elapsed time and time limit, time left (red when under 10 % of the limit, at least 15 min and at most 1 h, i.e. a `TIMEOUT` is near), the partition with the node (or, while pending, the reason such as `(Priority)`), peak memory (`sstat` while running, `sacct` after) and the exit code; with qsub only the state. Clicking a job shows its logs from `<output>/log_hpc` in three tabs, *Output* (`sub-XXX_slurm_output_<id>.log`), *Errors* and *Report* (`sub-XXX_job<id>_report.txt`, written by the report job after the job ends: wall time, resources, CO₂e, cost), below the job table (drag the divider to share the space), refreshed every 10 s while the job runs. *Follow* keeps the end of the log in view, *Wrap lines* wraps long lines, and A−/A+ change the text size (remembered). Logs are read from the file system when the folder is visible (on the cluster, or through the SSH folder mapping), else over SSH. *Cancel selected* runs `scancel`/`qdel`.

??? note "Submitting over SSH"
    - *Submit to cluster* runs `load_parameters` + `prestus_pipeline_start` on the host with `matlab -batch` (in a login shell, from the mapped study folder), and the job list follows the jobs over SSH.
    - *Connect…* logs in once: ssh's prompts (password, one-time code, an unknown host key) appear in a Nexus dialog, and the login is kept as a shared [OpenSSH](https://www.openssh.com/) connection that every later submit, status check and cancel reuses until it has been idle for 8 h (*Disconnect* closes it). With an SSH key, *Connect…* is not needed. The other commands never prompt (`BatchMode=yes`), so without a key or an open connection they fail with "Permission denied". Connection sharing needs OpenSSH 8.4 or later on macOS or Linux; on Windows use a key.
    - *Test connection* reports whether `sbatch`/`qsub`, MATLAB and PRESTUS are found on the host. If the host's `matlab` is a site wrapper that needs a display (at the Donders, `module load matlab/…` gives one; the log then says "No DISPLAY variable set"), set *Remote MATLAB* to the MATLAB binary itself (`<matlabroot>/bin/matlab`); `matlab -batch` runs without a display.
    - Config values under the mapped folder are translated the same way, e.g. `path.sim: /Volumes/myproject/data/sim` runs as `/project/myproject/data/sim`: the Nexus passes them to `load_parameters` as overrides (listed in the log) and leaves the config file unchanged. This covers the merged config (study, project `config_default.yaml` copy, transducers). A path on this Mac outside the mapped folder (`/Volumes/…`, `/Users/…`) cannot be translated; the Nexus lists those keys and asks before submitting.

## Results Tab

One folder usually holds several runs, told apart by their output affix (study conditions, uncertainty variants, async and sequential runs). The **Run** list names each one from its affix and starts on the run the open config produces (for async coupling: the joint result; a second list switches between the joint result and each transducer's own field, its `_trans-NN` run). Every sub-tab follows the selected run. Files are read when needed: only the shown sub-tab reads its files, volumes are read in the background (the viewer says *Loading …*) and kept for the next view, and opening the tab again reloads nothing unless the runs' files changed (*Refresh* rescans and reloads).

- *Overview* (opens first): the selected run grouped as the ITRUSST consensus asks TUS studies to report them: the non-significant-risk check, the transducer, placement and targeting, pulse timing, free-water acoustic parameters, the estimated in situ acoustic parameters and the thermal results. The risk check uses the limits and colours of the HTML report ([non-significant risk limits](doc_outputs_metrics-table.md#non-significant-risk-limits)): ✓ below half the limit, ▲ from half the limit up to it, ✗ above it; above a limit means inspecting adherence to the NSR criteria that govern the protocol. Each row names its source: the run's metrics CSV, the free-water run with the same output affix (`modules.run_water_baseline`), the targeting CSV (`modules.run_targeting_analysis`), the placement record, the equipment library or the open config. Pulse timing and transducer geometry come from the open config, which can differ from what the run used; PRF and duty cycles are derived from the timing settings. ISPTA is averaged over the pulse repetition within a pulse train (PD/PRI), as PRESTUS computes it; ISPTA per tissue is that tissue's ISPPA times the run's duty cycle (ISPTA/ISPPA of its CSV, else PD/PRI of the open config). The axial and lateral −3 dB focal widths come from the CSV ([`focal_fwhm_*`](doc_outputs_metrics-table.md#acoustic-metrics)); for runs written before these columns existed, the Nexus measures them in the background on the run's T1w intensity map, along the transducer axis of its placement record and through the intensity peak inside the −6 dB focal mask, and the row names that map as its source. These can differ from the grid's by a fraction of a voxel. What no output holds reads *not available*. **Copy as Markdown** copies the overview as tables.
- *Spatial* ▸ *Acoustic* / *Thermal*: the run's maps on the T1 (only the maps it wrote are offered), centred on the map's peak (*Go to peak* re-centres). A colour bar gives the map's unit (W/cm², Pa, °C, min). Intensity, pressure, MI, temperature rise and CEM43 are coloured from 0 to the map's maximum; absolute temperatures from their 1st to 99th percentile. *Overlay colormap* offers magma (the default), inferno, viridis, plasma and hot. *Hide below* hides values under a share of the colour scale's top (10 % by default for intensity, pressure, MI and temperature rise). Clicking a slice moves the crosshair; the line under the slices gives its RAS position and the map's value there. Every transducer (square) and focus (x) is marked, as the run's placement stage recorded them (`sub-XXX_placement<affix>.json`), else from the config; points away from the slice are drawn faint.
- *Heating*: the run's thermal time courses from `cache/sub-XXX_<medium>_HEAT<affix>.mat` (written by the thermal analysis, `modules.run_thermal_analysis`; not with the BabelBrain thermal backend). *Metric* picks the peak temperature, temperature rise, CEM43 ([k-Wave](https://github.com/ucl-bug/k-wave)) or CEM43 (ISO) of each tissue layer (max over the layer), or the temperature, rise or CEM43 on the focal axis (per tissue and at the target). Tick the tissues to show; the choice is kept across metrics and runs. The legend gives each tissue's maximum. *Log scale* helps with CEM43; *Shade sonication* marks the sonication periods. Needs [SciPy](https://scipy.org/) (or h5py for `-v7.3` files). For an async transducer's own field the joint run's time courses are shown.
- *Report*: the run's HTML report, or any other report in the folder (sequential and multi-ISPPA summaries). *Export report as PDF…* saves the shown report.
- *Data table*: the run's metrics (last CSV row).
- *Poses*: shown when the run wrote a time-reversal pose table (`_timereversal_poses_T<i>.csv`); *Use this pose…* starts a follow-up at the selected pose.
- *Compare runs*: one row per run with the key metrics (ISPPA, MI, focal distance and volume, temperature rise, CEM43); double-click a row to show that run.

### Follow-up Runs

**Follow-up…** in the Results header creates a new study config that starts from the selected run, and opens it. Pick one of five types:

- *New intensity*: another target ISPPA on the same beam. It reuses the run's pressure map (`path.acoustic_input`), so no acoustic simulation runs ([Reusing outputs](doc_outputs_reuse.md#starting-from-earlier-stage-outputs)).
- *New timing (thermal only)*: other pulse timing on the same pressure map; only the thermal simulation runs.
- *New placement, same head*: another transducer position; the head model is built once and reused (`headmodel.reuse_across_placements`).
- *Pose from the time-reversal search*: a pose from the run's pose table as a manual placement, or its refocusing phases ([Pose search](doc_time_reversal.md#pose-search)).
- *Sequential*: a run that continues from this run's end temperature and thermal dose. It is added to this config's *Pipeline options* ([Run tab](#run-tab)) rather than run on its own ([Sequential runs](doc_sequential.md)).

The dialog lists problems as you edit, and errors block *Create*. Each follow-up gets a new `io.output_affix` and is saved next to the parent as `<parent>_<type>_<affix>.yaml`. The config records its parent in a `# nexus-followup:` comment line, which MATLAB ignores and the Nexus keeps on every save. The Results header names a run's parent and its follow-ups.

## Study Tab

Runs many subjects with per-run settings. See [Group analysis](doc_group.md#study-setup).

??? note "The run table"
    - The table (the study manifest, `study.manifest`) has one row per simulation run; settings not in the table come from the study config. A single subject does not need it. *How to fill the table* shows the steps and an example.
    - *Add subjects…* adds one run per subject ID (`1-12, 15, sub-020`). *Add setting column…* adds a config key to vary per run (empty cell = the study value) or the run columns `include`, `config` and `notes`; hovering over a column header shows its help. *More* holds adding an empty row, removing rows or a column, and *Save table as…*.
    - The table is saved to the *Table file* (relative to the config). Without one, the first save writes `<config name>_runs.tsv` next to the config and sets `study.manifest`.
    - *Import from Localite…* reads the `neuronav_gui` outputs under `path.localite_post` (`neuronav_ingest_summary.csv` and the per-session `InstrumentMarker_*.json`). It makes one run per target, or one run per session with the targets assigned to transducers 1, 2, …, with a preview before adding. See [Group analysis](doc_group.md#importing-runs-from-localite).
    - The summary under the table checks the manifest as you type: subject IDs, `include` values, unknown keys (with suggestions), and rows that would write to the same outputs.
    - *Check all runs* saves the table and runs `prestus_study_start(config, 'DryRun', true)`. *Run all runs locally* writes `run_study_<study>.m` and runs it; *Write launcher script* only writes it. These stay disabled while the table is empty. *Build group report* runs `prestus_group_report_start(config)` (group maps + report). Checks, runs and the report stream to the Run tab.

## Timing Diagram and TUS Calculator

The *Sonication protocol* card (Sonication tab) draws the protocol from `timing.*`: one pulse (with its ramp), one pulse train, and the train repetitions. It labels PD, PRI/PRF, PTD, PTRI and PTRD with the pulse count, duty cycles, the on-time, and the Ipa/Itp factor of the ramp (`compute_ramp_factor`). It also lists inconsistent settings (e.g. a pulse longer than its repetition interval).

You can edit the protocol in the diagram. Drag an edge to change its value: the pulse end (PD), the end of the rise (ramp), the dashed start of the next pulse (PRI), the train end (PTD), the dashed start of the next train (PTRI) or the end of the repetitions (PTRD). Dragged values snap to two significant digits and stay within the other settings (a pulse cannot outgrow its PRI). To type an exact value, click an edge or its label and enter a duration with its unit (`20 ms`, `500 µs`); the PRI also accepts a frequency (`5 Hz`). Enter applies the value and Esc cancels. The fields for these six values are folded under *Show the values as fields*; a field whose value is invalid stays in view.

![Nexus Sonication tab with the sonication timing diagram](https://github.com/jkosciessa/PRESTUS_bin/raw/main/img/nexus_simulation-timing_screenshot.png)

The [TUS calculator](https://www.socsci.ru.nl/fusinitiative/tuscalculator/) (Radboud University, MIT licence) estimates derived quantities (derated intensities, MI, energy deposition, temperature rise) for a protocol. The Nexus does not copy its calculations. It exchanges settings files with it instead:

- **Export for TUS calculator…** writes a JSON file for the calculator's *Upload settings*: pulse duration, ramp duration and shape, PRI, PTD, PTRI, PTRD (in ms), whether the train and repetition levels are enabled, and (when set) the frequency, target I<sub>SPPA</sub>, aperture and focal depth.
- **Import from TUS calculator…** reads a settings file downloaded from the calculator (JSON or YAML). It shows the values it will set (in seconds), and optionally also sets the frequency and target I<sub>SPPA</sub>.
- **Open TUS calculator ↗** opens the web app.

??? note "Ramp convention"
    Both use the same ramp convention: the ramp duration is the rise (= fall) time at each end of a pulse (`timing.pulse_ramp_dur`, `timing.pulse_ramp_shape` `linear` or `tukey`). PRESTUS models pulse ramps only, so the calculator's pulse-train and repetition ramps are exported as 0 and reported when an imported file sets them. A level disabled in the calculator sets the corresponding PRESTUS durations to 0.
