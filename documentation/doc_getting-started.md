# Getting started

### Install packages 

See [Installation](doc_installation.md).

### [DEPRECATED] Run the tutorial (documentation/PRESTUS_intro_tutorial.md)
- The acoustic profile is derived from the manufacturers measurements and adjusted to more closely match the actual characteristics of the transducer.
- The structural data is segemented with SimNIBS 4's charm.
- The segmented data is translated to a matrix of acoustic properties.
- k-wave performs acoustic and thermal simulations.
- Select outputs are plotted.

### [Optional] Explore demos

A simplified 2D benchmarking example can be found [here](https://github.com/jkosciessa/PRESTUS_2D_demo). It is a lightweight example that is designed to run on local CPU ressources.

A 3D demo using SimNIBS' Ernie template is provided as a [DataLad dataset](https://gin.g-node.org/PRESTUS/sim_ernie/). [Note: This example is based on a prior PRESTUS version. It still has to be updated, but provides a practical example of a possible data management scheme.]

### Create a study-specific config file

Configuration metadata must be fully specified in (a set of) config files (in the .yaml format) placed in the 'configs' folder (see [the configuration overview](doc_parameters.md). Some parameters found in the 'default_config' are mandatory for the pipeline to run while others can be left out or are expected to be changed based on the requirements of the analysis. The different parameters that can be used are found in the `default_config.yaml` and in [the configuration overview](doc_parameters.md). The default config is not intended to be altered. Instead, parameter changes should be specified via an application-specific config file (or after loading the default config and before performing the subject call in MATLAB).

The `load_parameters.m` function imports the configuration. Multiple configurations can be read in sequentially to overwrite specific portions of the default configuration.

### Specify a transducer and target location

Specify locations via `parameters.transducer.trans_pos` and `parameters.transducer.focus_pos`. 

See [Placement](doc_placement.md).

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
parameters.t1_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T1.nii.gz");
parameters.t2_path_template = fullfile(sprintf('m2m_sub-%03d', subject_id), "T2_reg.nii.gz");
```

### [Optional] Advance SimNIBS segmentation

It may be desirable to run an advance SimNIBS call prior to running the full pipeline incl. segmentation postprocessing, source setup, acoustic and thermal simulations. By default, existing segmentations will be reused and not overwritten unless explicitly requested with `overwrite_simnibs` (regardless of `overwrite_files`).

To this end, `run_source_setup`, `run_acoustic_sims`, `run_heating_sims`, and `run_posthoc_water_sims` can be deactivated. This separate step is required to inform the skull layer [using pseudoCTs](doc_pseudoCT.md).

### Run the single_subject_pipeline

When the configuration is fully specified, the full pipeline only needs a subject id (a number consisting of no more than 3 digits). To run multiple subjects or setups in parallel, the `single_subject_pipeline` can be submitted using high performance computing jobs (see [HPC documentation](doc_hpc.md)).

PRESTUS can be parallelized across-subject and within-subject across treatment setups, e.g., by iterating across different transducers, transducer placement and target coordinates, or amplitudes and/or focal depth settings (requiring either dynamic [recalibration](doc_calibration.md) using `transducer_calibration` or a precomputed lookup table) via custom MATLAB scripts.
