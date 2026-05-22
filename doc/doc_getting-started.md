# Quick Start Guide

### Install PRESTUS and depencencies 

See [Installation](doc_installation.md).

### [Optional] Explore demos

A simplified 2D benchmarking example can be found [here](https://github.com/jkosciessa/PRESTUS_2D_demo). It is a lightweight example that is designed to run on local CPU ressources.   

A 3D demo using SimNIBS' Ernie template is provided as a [DataLad dataset](https://gin.g-node.org/PRESTUS/sim_ernie/). [Note: This example is based on a prior PRESTUS version. It still has to be updated, but provides a practical example of a possible data management scheme.]     

A (currently outdated) demo that focuses on transducer calibration can be found [here](https://github.com/jkosciessa/PRESTUS_bin/blob/main/tutorial/PRESTUS_intro_tutorial.md).

### Create virtual transducer settings

See [Transducer] and [Calibration](doc_calibration.md).

### Create a study-specific config file

Configuration metadata must be fully specified in (a set of) config files (in the .yaml format). See [the configuration overview](doc_parameters.md) for all available parameters.

**Recommended:** Keep a project-specific config directory outside the PRESTUS toolbox. Use `prestus_config_init` to copy the current toolbox defaults:

```matlab
% Copy config_default.yaml to your project config dir and create a project override file
prestus_config_init('/path/to/project/config', project_name='projectX');
```

This copies `config_default.yaml` from the toolbox into your project directory and creates an empty `config_projectX.yaml` for your overrides. Edit `config_projectX.yaml` to specify only the parameters that differ from the defaults (see [the configuration overview](doc_parameters.md)). To review and adjust parameters interactively, pass `open_gui=true`:

```matlab
prestus_config_init('/path/to/project/config', project_name='projectX', open_gui=true);
```

To load parameters in your scripts:

```matlab
parameters = load_parameters('config_projectX.yaml', '/path/to/project/config');
```

`load_parameters` always uses `config_default.yaml` as the base (resolved from your project directory if present, otherwise from the toolbox), then deep-merges the project-specific config on top. This means your project config only needs to contain the parameters you want to override.

> **Note:** `config_default.yaml` is not intended to be edited directly in the toolbox. Keeping a project-specific copy also pins the baseline parameter set for reproducibility across toolbox updates.

### Choose the Simulation Medium

See [Medium Setup](doc_medium.md).

### Choose the Simulation Backend

See [Backend](doc_backend.md).

### Specify I/O and Planning Images

Specify folder management for inputs, segmentations, and simulations:

```
parameters.ld_library_path
parameters.data_path
parameters.seg_path
parameters.sim_path
parameters.paths_to_add = {pn.kwave, pn.minimize};
parameters.simnibs_bin_path = fullfile('.conda', 'envs', 'simnibs_env', 'bin');
```

For `layered` simulations, specify at least one T1w planning image (and optionally a T2w/UTE image) for the segmentation.

```
parameters.path.t1_pattern = 'sub-%03d_T1w.nii.gz';
parameters.path.t2_pattern = 'sub-%03d_UTE.nii.gz';  % optional: T2w or PETRA UTE image
```

### Specify a transducer and target location

For `layered` simulations, specify locations via `parameters.transducer.trans_pos` and `parameters.transducer.focus_pos`. 

See [Placement](doc_placement.md).

### [Optional] One-shot SimNIBS segmentation

It may be desirable to run an intial SimNIBS call prior to running the full pipeline incl. segmentation postprocessing, source setup, acoustic and thermal simulations. By default, existing segmentations will be reused and not overwritten unless explicitly requested with `overwrite_simnibs` (regardless of `overwrite_files`).

Set `modules.segmentation_only = 1` to stop the pipeline immediately after segmentation, skipping all subsequent steps.

### [Optional] Use a pseudoCT for subject-specific skull properties

Point `path.t2_pattern` to a PETRA UTE image and set `pct.enabled = 1`. PRESTUS will pass the UTE image to SimNIBS during segmentation and then automatically generate a pseudoCT (`pseudoCT.nii.gz`) in the SimNIBS output folder. If a pseudoCT already exists in the SimNIBS `m2m` folder, generation is skipped.

See [pseudoCT](doc_pseudoCT.md) for algorithm choices.

### Run the prestus_pipeline | [Optional] iterate across parameters

When the configuration is fully specified, the full pipeline only needs a subject id (a number consisting of no more than 3 digits). To run multiple subjects or setups in parallel, the `prestus_pipeline` can be submitted using high performance computing jobs as a part of `prestus_pipeline_start` (see [HPC documentation](doc_hpc.md)).

PRESTUS can be parallelized across-subject and within-subject across treatment setups, e.g., by iterating across different transducers, transducer placement and target coordinates, or amplitudes and/or focal depth settings (requiring either dynamic [recalibration](doc_calibration.md) using `transducer_calibration` or a precomputed lookup table) via custom MATLAB scripts.
