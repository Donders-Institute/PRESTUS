# Troubleshooting

This page lists known problems by the message or symptom you see, each with its cause and fix.

On this page:

- [Installation and Setup](#installation-and-setup): missing [k-Wave](https://github.com/ucl-bug/k-wave), platform detection, config warnings, macOS file access and Nexus installation.
- [Segmentation (SimNIBS charm)](#segmentation-simnibs-charm): [`charm`](https://simnibs.github.io/simnibs/build/html/documentation/command_line/charm.html) errors and segmentations that are started twice or not redone.
- [Placement](#placement): [PlanTUS](https://github.com/mlueckel/PlanTUS), [TUS_entry](https://github.com/CyrilAtkinson/TUS_entry) and positions that land in the wrong place.
- [Simulation](#simulation): ISPPA targeting errors, cache mismatches, sampling warnings and running out of memory.

!!! note "At a glance"
    - **When you need it:** When a run stops with one of the messages below, or a result lands in the wrong place.
    - **Before you start:** The error identifier or warning text from the [MATLAB](https://www.mathworks.com/products/matlab.html) console or the `log/` folder.

Most entries link to the page that explains the behaviour. If you solve an issue that is not listed here, add it. If you find an error in the pipeline or have a suggestion, [open an issue on GitHub](https://github.com/Donders-Institute/PRESTUS/issues).

## Installation and Setup

| Symptom | Cause and fix |
|---|---|
| `makeBowl` (or another k-Wave function) is not found | k-Wave is not on the MATLAB path. Clone with submodules and add `external/` ([Installation](doc_installation.md)). |
| `No HPC system detected (sbatch/qsub)` | `platform: auto` (the default) only detects [SLURM](https://slurm.schedmd.com/) or qsub. On a workstation set `platform: matlab` ([HPC](doc_hpc.md#where-a-run-executes)). |
| Warning `prestus:unknownConfigKey` | A key in your YAML is not a PRESTUS setting, often a typo; the warning suggests the closest known key ([Configuring a study](doc_configuration.md#how-a-config-is-assembled)). |
| Warning `prestus:staleProjectDefault` | Your project copy of `config_default.yaml` predates settings added since; they take their toolbox default. Refresh it with `prestus_config_upgrade(folder)` ([Configuring a study](doc_configuration.md#project-config-folder)). |
| macOS: files cannot be read when MATLAB is started from the Finder or Dock | Add `MATLAB.app` under **System Settings → Privacy & Security → Full Disk Access** and restart MATLAB. MATLAB started from a terminal is not affected. |
| Nexus: `Library not loaded: …QtPrintSupport` | [conda](https://docs.conda.io/)'s `pyqt` mixed with [pip](https://pip.pypa.io/)'s `PyQt6-WebEngine`; install Qt from pip only ([Nexus](doc_nexus.md#installation)). |

## Segmentation (SimNIBS charm)

| Symptom | Cause and fix |
|---|---|
| `ValueError: The qform and sform of do not match. Please run charm with the --forceqform option` | Set `segmentation.use_qform: 1`. |
| `undefined symbol` errors from `charm` | Empty, `hpc.ld_library_path` already gives SimNIBS tools the install's own libraries. Otherwise set it to the shared libraries [SimNIBS](https://simnibs.github.io/simnibs/build/html/index.html) needs (e.g. `/opt/gcc/7.2.0/lib64`; [Installation](doc_installation.md#point-prestus-to-simnibs)). |
| Two runs of the same subject (e.g. left and right transducer) both start a segmentation | Neither finds a finished `m2m_sub-XXX` folder. Segment first (`modules.segmentation_only: 1`), then start the simulations. |
| A segmentation is not redone after changing the input images | Existing segmentations are reused unless `io.overwrite_simnibs: 1` (independent of `io.overwrite_files`). |

## Placement

| Symptom | Cause and fix |
|---|---|
| PlanTUS interactive mode runs headless | No display was found (e.g. a cluster job or an SSH session without X forwarding); run it on a desktop or in a VNC session. |
| PlanTUS places the last clicked vertex without a review dialog | [PyQt](https://www.riverbankcomputing.com/software/pyqt/)5 is missing in the SimNIBS Python; it ships with SimNIBS, so check `placement.plantus.env_path`. |
| PlanTUS stops because intermediate files are missing | `wb_command` ([Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench)) is not found; install it in the SimNIBS environment or set `placement.plantus.connectome_wb_path` ([Installation](doc_installation.md#plantus-dependencies)). |
| TUS_entry: `replacement has 1 row, data has 0` | TUS_entry needs a solid head mask; check `placement.tus_entry.scalp_mask` ([TUS_entry](doc_placement_tus_entry.md#required-configuration)). |
| MNI or RAS positions end up in the wrong place | Manual `ras`/`mni` coordinates are converted in the placement stage, so `modules.run_transducer_placement` must be `1` ([Placement](doc_placement.md#manual-coordinate-selection)). For MNI and [Localite](https://www.localite.de/) conventions see [Coordinate systems](doc_coordinate_systems.md). |

## Simulation

| Symptom | Cause and fix |
|---|---|
| `prestus:isppa:noBaseline` or `prestus:isppa:unsimulatedBaseline` | `transducer.target_isppa_wcm2` needs a free-water baseline that is not available, e.g. `modules.run_water_baseline: 0` with an uncalibrated transducer. Use `auto`, or set the target to `[]` to report the field at `elem_p_amp` ([ISPPA targeting](doc_simulations-acoustic.md#isppa-targeting)). |
| Warning `prestus:isppa:largeScale` | The target is more than 4× (or less than 0.25×) the simulated drive; check that this is still in the linear regime. |
| Warning `load_transducer_from_library:uncalibratedAmplitude` | The library drive at this depth was fitted analytically and not checked in k-Wave; see [k-Wave drive calibration](doc_calibration_refinement.md#free-water-k-wave-drive-calibration). |
| `prestus:cache:signatureMismatch` | A reused grid, medium or acoustic cache was built with other settings. Recompute (`io.overwrite_files: always` or `force`) or delete the cache files ([Cache reuse checks](doc_outputs_reuse.md#cache-reuse-checks)). |
| Changed source or medium settings have no effect on a re-run | A grid, medium or acoustic cache was reused although its settings differ: `io.cache_signature_check` is `warn` (reuse with a warning) or `off`, or a chained-pipeline sub-job loaded a cache another stage wrote, which is not checked ([Cache reuse checks](doc_outputs_reuse.md#cache-reuse-checks)). Keep `io.cache_signature_check: error`, use `io.overwrite_files: force` or delete the relevant `cache/` files ([Overwrite behaviour](doc_outputs_reuse.md#overwrite-behaviour)). |
| Points-per-wavelength (PPW) warning | The grid is too coarse for the frequency; reduce `grid.resolution_mm` ([Spatial and temporal sampling](doc_simulations-acoustic.md#spatial-and-temporal-sampling)). |
| Out of memory in the thermal simulation | Coarsen `grid.thermal_resolution_mm`, or use a GPU with more memory ([Thermal](doc_simulations-thermal.md#step-temperature-and-sensor-settings), [HPC](doc_hpc.md#gpu-support)). |
| Out of memory in the acoustic simulation | Use a coarser grid, a [limited acoustic FOV](doc_limited_fov.md), or the C++ binaries ([Memory](doc_simulations-acoustic.md#memory-use)). |
