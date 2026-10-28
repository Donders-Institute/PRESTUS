# Troubleshooting

Known problems and their fixes, grouped by where they show up. Most entries link to the page that explains the behaviour. If you solve an issue that is not listed here, please add it. If you find an error in the pipeline or have a suggestion, please [open an issue on GitHub](https://github.com/Donders-Institute/PRESTUS/issues).

## Setup and starting a run

| Symptom | Cause and fix |
|---|---|
| `makeBowl` (or another k-Wave function) is not found | k-Wave is not on the MATLAB path. Clone with submodules and add `external/` ([Installation](doc_installation.md)). |
| `No HPC system detected (sbatch/qsub)` | `platform: auto` (the default) only detects SLURM or qsub. On a workstation set `platform: matlab` ([HPC](doc_hpc.md#where-a-run-executes)). |
| Warning `prestus:unknownConfigKey` | A key in your YAML is not a PRESTUS setting, often a typo; the warning suggests the closest known key ([Configuring a study](doc_configuration.md#how-a-config-is-assembled)). |
| Warning `prestus:staleProjectDefault` | Your project copy of `config_default.yaml` predates settings added since; they take their toolbox default. Refresh it with `prestus_config_upgrade(folder)` ([Configuring a study](doc_configuration.md#a-project-config-folder)). |
| macOS: files cannot be read when MATLAB is started from the Finder or Dock | Give `MATLAB.app` Full Disk Access ([GUI](doc_gui.md#macos-full-disk-access-required)). |
| Nexus: `Library not loaded: …QtPrintSupport` | conda's `pyqt` mixed with pip's `PyQt6-WebEngine`; install Qt from pip only ([Nexus](doc_nexus.md#installation)). |

## Segmentation (SimNIBS `charm`)

| Symptom | Cause and fix |
|---|---|
| `ValueError: The qform and sform of do not match. Please run charm with the --forceqform option` | Set `segmentation.use_qform: 1`. |
| `undefined symbol` errors from `charm` | Set `hpc.ld_library_path` to the shared libraries SimNIBS needs (e.g. `/opt/gcc/7.2.0/lib64`; [Installation](doc_installation.md#specifying-simnibs-paths-in-prestus)). |
| Two runs of the same subject (e.g. left and right transducer) both start a segmentation | Neither finds a finished `m2m_sub-XXX` folder. Segment first (`modules.segmentation_only: 1`), then start the simulations. |
| A segmentation is not redone after changing the input images | Existing segmentations are reused unless `io.overwrite_simnibs: 1` (independent of `io.overwrite_files`). |

## Placement

| Symptom | Cause and fix |
|---|---|
| PlanTUS interactive mode runs headless | No display was found (e.g. a cluster job or an SSH session without X forwarding); run it on a desktop or in a VNC session. |
| PlanTUS places the last clicked vertex without a review dialog | PyQt5 is missing in the SimNIBS Python; it ships with SimNIBS, so check `placement.plantus.env_path`. |
| PlanTUS stops because intermediate files are missing | `wb_command` (Connectome Workbench) is not found; install it in the SimNIBS environment or set `placement.plantus.connectome_wb_path` ([Installation](doc_installation.md#optional-plantus-dependencies)). |
| TUS_entry: `replacement has 1 row, data has 0` | TUS_entry needs a solid head mask; check `placement.tus_entry.scalp_mask` ([TUS_entry](doc_placement_tus_entry.md#required-configuration)). |
| MNI or RAS positions end up in the wrong place | Manual `ras`/`mni` coordinates are converted in the placement stage, so `modules.run_transducer_placement` must be `1` ([Placement](doc_placement.md#manual-coordinate-selection)). For MNI and Localite conventions see [Coordinate systems](doc_coordinate_systems.md). |

## Simulation

| Symptom | Cause and fix |
|---|---|
| `prestus:isppa:noBaseline` or `prestus:isppa:unsimulatedBaseline` | `transducer.target_isppa_wcm2` needs a free-water baseline that is not available, e.g. `modules.run_water_baseline: 0` with an uncalibrated transducer. Use `auto`, or set the target to `[]` to report the field at `elem_p_amp` ([ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)). |
| Warning `prestus:isppa:largeScale` | The target is more than 4× (or less than 0.25×) the simulated drive; check that this is still in the linear regime. |
| Warning `load_transducer_from_library:uncalibratedAmplitude` | The library drive at this depth was fitted analytically and not checked in k-Wave; see [k-Wave drive calibration](doc_transducer_library.md#optional-free-water-k-wave-drive-calibration). |
| `prestus:cache:signatureMismatch` | A reused grid, medium or acoustic cache was built with other settings. Recompute (`io.overwrite_files: always` or `force`) or delete the cache files ([Cache reuse checks](doc_outputs.md#cache-reuse-checks)). |
| Changed source or medium settings have no effect on a re-run | The acoustic cache only checks the grid size. Use `io.overwrite_files: force` or delete the relevant `cache/` files ([Overwrite behaviour](doc_outputs.md#overwrite-behaviour)). |
| Points-per-wavelength (PPW) warning | The grid is too coarse for the frequency; reduce `grid.resolution_mm` ([Spatial and temporal sampling](doc_simulations-acoustic.md#spatial-and-temporal-sampling)). |
| Out of memory in the thermal simulation | Keep `thermal.record_t_at_every_step: 0`, reduce `thermal.sensor_xy_halfsize`, coarsen `grid.thermal_resolution_mm`, or use a GPU with more memory ([Thermal](doc_simulations-thermal.md#additional-parameters), [HPC](doc_hpc.md#gpu-support)). |
| Out of memory in the acoustic simulation | Use a coarser grid, a [limited acoustic FOV](doc_limited_fov.md), or the C++ binaries ([Memory](doc_simulations-acoustic.md#memory-and-source-signal-optimisation)). |
