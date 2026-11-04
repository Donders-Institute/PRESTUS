# Concepts and glossary

This page explains what PRESTUS does in one picture and defines the terms the rest of the documentation uses. Read it once before your first head simulation, and come back to the glossary whenever a term is unclear.

## What PRESTUS does

PRESTUS estimates where an ultrasound transducer on the head focuses its energy, how strong the field is inside the brain, and how much it heats the tissue. It does so in a sequence of steps, each building on the previous one:

```
  T1 (+ T2/UTE) MRI of the participant
        │
  1. Segmentation        SimNIBS charm labels skin, skull, brain, …
        │
  2. Placement           where the transducer sits and where it aims
        │
  3. Grid & medium       a 3D grid of tissue properties (sound speed, density, attenuation)
        │
  4. Acoustic simulation the pressure field of the transducer through the head (k-Wave),
        │                scaled to the intensity you asked for
  5. Thermal simulation  heating over the sonication protocol (optional)
        │
  6. Outputs             metrics table, NIfTI maps, HTML report with safety checks
```

Every step reads its settings from a config file ([Configuring a study](doc_configuration.md)), and each can be switched on or off ([Modules](doc_modules.md)). A simulation in **free water** skips steps 1 to 3: no MRI is needed, which makes it the quickest way to try PRESTUS ([Your first simulation](doc_getting-started.md)).

Simulations are estimates. Tissue properties, especially of the skull, are uncertain, and results depend on the segmentation, the grid and the placement. PRESTUS is a research tool; it is not approved for clinical use.

## Glossary

### Transducer and equipment

Transducer
:   The device that emits the ultrasound. PRESTUS models focused, bowl-shaped transducers.

Annular array
:   A transducer made of concentric rings. Changing the phase of each ring moves the focus along the beam axis (deeper or shallower), but not sideways. Most TUS transducers are annular. See [Transducer Definition](doc_transducer.md).

Matrix array
:   A transducer made of many small elements that can steer the focus in 3D.

Driving system
:   The electronics that drive the transducer (e.g. Sonic Concepts TPO, IGT). The same transducer behaves slightly differently with different driving systems, so calibrations are per transducer and driving-system pair.

Serial, combo
:   The identifiers of a transducer (`transducer.serial`) and its driving system (`transducer.combo.ds_serial`) in the shipped [transducer library](doc_transducer_library.md). Setting them fills in geometry, phases and drive.

Bowl apex
:   The centre of the back of the transducer bowl. PRESTUS positions the transducer by this point (`transducer.trans_pos`).

Exit plane
:   The plane of the transducer's front face, where the ultrasound leaves the device. Focal depths are usually given from here.

Geometric focus
:   The natural focus of the bowl, at its centre of curvature, without any phase steering.

Focal depth setting (`transducer.focal_distance_ep`)
:   The depth, from the exit plane, that the device is set to focus at. For annular arrays it selects the ring phases.

Target depth
:   The distance from the exit plane to the target you planned. Ideally it matches the focal depth setting; if not, the simulated focus lands at the setting, not at the target.

FLHM
:   Focal Length at Half Maximum: the stretch of the beam axis where the intensity is above half its peak. Manufacturers quote focal depths as the centre of this stretch. See [Calibration](doc_calibration.md#focal-depth-convention-flhm).

Calibration
:   Fitting the phases and drive of the simulated transducer so that its field in water matches measurements (hydrophone profiles) of the real device. See [Calibration](doc_calibration.md).

### Intensity and safety metrics

Free water
:   A simulation in water only, without a head. Measurements and manufacturer specifications refer to free water.

ISPPA
:   Spatial-peak pulse-average intensity (W/cm²): the highest intensity in the field, averaged over the duration of a pulse. PRESTUS reports it in free water and in each tissue (e.g. `Isppa_brain`).

Target ISPPA (`transducer.target_isppa_wcm2`, default 30 W/cm²)
:   The free-water ISPPA you want to deliver. PRESTUS scales the simulated head field so that the same drive would give this ISPPA in free water ([ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)).

Free-water baseline
:   A short free-water simulation of the same transducer that PRESTUS runs to know what ISPPA the simulated drive gives. It is the reference for that scaling.

In situ
:   In the tissue, after the skull. The in-situ intensity at the target is usually much lower than the free-water value, mainly because of the skull.

Ipa, Ita
:   Pulse-average and time-average intensity at a location. Ita includes the off-time of the protocol (the duty cycle). See [Outputs](doc_outputs.md#intensity-nomenclature).

MI
:   Mechanical index, a measure of the risk of cavitation (peak negative pressure in MPa divided by the square root of the frequency in MHz). `MI_tc` is the transcranial MI inside the skull.

Temperature rise
:   The increase over the starting temperature (37 °C by default) during the protocol.

CEM43
:   Cumulative equivalent minutes at 43 °C: a thermal dose that combines how hot tissue got and for how long. See [Thermal](doc_simulations-thermal.md#thermal-dose-cem43-calculation).

ITRUSST
:   The International Transcranial Ultrasonic Stimulation Safety and Standards consortium. Its consensus documents define the non-significant-risk (NSR) levels the PRESTUS report compares against (e.g. MI ≤ 1.9, temperature rise ≤ 2 °C). See [Reading your results](doc_results.md).

### Head model

Segmentation
:   Dividing the MRI into tissue types. PRESTUS uses SimNIBS `charm`.

m2m folder
:   The SimNIBS output folder of a subject (`m2m_sub-XXX`, under `path.seg`). PRESTUS reads it but never writes into it.

Medium
:   What the ultrasound travels through: `water` (free water), `layered` (a head from the segmentation, with layers such as skin, skull and brain) or `phantom` (a designed test object). See [Medium Setup](doc_medium.md).

Medium properties
:   Sound speed, density, attenuation and thermal properties of each tissue (`medium_properties`).

Pseudo-CT
:   An image of skull density estimated from a special MRI (UTE/PETRA). It gives each skull voxel its own properties instead of one value per layer. See [pseudo-CT](doc_pseudoCT.md).

Placement mode (`placement.mode`)
:   How the transducer position and target are obtained: typed in (`manual`), given in MNI coordinates (`mni`), searched automatically (`heuristic`, `plantus`, `tus_entry`) or read from neuronavigation (`localite`). See [Placement](doc_placement.md).

Voxel, RAS, MNI
:   Three ways to give a position: as an index into the T1 image (voxel), in millimetres in the participant's scanner space (RAS), or in millimetres in a standard brain template (MNI). See [Coordinate systems](doc_coordinate_systems.md).

### Simulation

k-Wave
:   The open-source acoustic simulation toolbox that PRESTUS uses by default. BabelBrain and jwave are alternatives ([Backends](doc_backend.md)).

Grid
:   The 3D box of points on which the simulation is computed. `grid.resolution_mm` is the spacing between points: finer grids are more accurate but slower and need more memory.

PPW
:   Points per wavelength: how finely the grid samples the ultrasound wave. PRESTUS warns below 6.

PML
:   Perfectly matched layer: an absorbing border around the grid that stops waves from reflecting off its edges.

Axisymmetric
:   A 2D simulation that assumes the setup is rotationally symmetric around the beam axis. Much faster, for free water and simple tests.

Acoustic FOV
:   An optional smaller sub-grid around the beam for the acoustic simulation, to save time and memory ([Limited acoustic FOV](doc_limited_fov.md)).

Sonication protocol (`timing.*`)
:   How the ultrasound is switched on and off: pulse duration (PD), pulse repetition interval (PRI), pulse train duration (PTD) and its repetitions (PTRI, PTRD). Needed for thermal simulations.

`code_type`
:   Which k-Wave implementation runs the simulation: MATLAB or C++, on CPU or GPU. A GPU makes head simulations much faster.

Platform (`platform`)
:   Where the run happens: in your MATLAB session (`matlab`) or as a job on a cluster (`slurm`, `qsub`). See [Running on a cluster](doc_hpc.md).

### Files and runs

Study config
:   Your YAML file with the settings that differ from the defaults ([Configuring a study](doc_configuration.md)).

Output affix (`io.output_affix`)
:   A label added to output file names, to keep several runs of the same subject apart (it appears as `_desc-<label>`).

Cache
:   Intermediate results (`cache/` in the output folder) that PRESTUS reuses to skip work on re-runs.

Study manifest
:   A table with one row per run (subject, condition and per-run settings) for running many subjects at once ([Group analysis](doc_group.md#study-setup)).
