# Getting started

### Install packages (ensure the paths and subfolders are added in Matlab)
- SimNIBS 4
- k-Wave 1.4.1 (place in PRESTUS' toolbox folder)
    - Note that k-Wave 1.4 is supported, but version 1.4.1 (currently [GitHub exclusive](https://github.com/ucl-bug/k-wave/releases/tag/v1.4.1)) introduced GPU support for thermal simulations with kWaveDiffusion.

### [Optional] Start with the tutorial (documentation/PRESTUS_intro_tutorial.md)
- The acoustic profile is derived from the manufacturers measurements and adjusted to more closely match the actual characteristics of the transducer.
- The structural data is segemented with SimNIBS 4's charm.
- The segmented data is translated to a matrix of acoustic properties.
- k-wave performs acoustic and thermal simulations.
- Select outputs are plotted.

### Create a study-specific config file
- Configuration metadata must be fully specified in (a set of) config files (in the .yaml format) placed in the 'configs' folder (see ```doc_config.md```). 
- The ```load_parameters.m``` function imports the configuration. 
- Some of the parameters found in the 'default_config' are mandatory for the pipeline to run while others can be left out based on the requirements of the analysis.
- The different parameters that can be used are found in the 'default_config.yaml'. The default config should generally not be altered. Instead, parameter changes should be specified via an application-specific config file or before performing the subject call in MATLAB).

### Specify a transducer and target location
Specify locations via `parameters.transducer.trans_pos` and `parameters.transducer.focus_pos`. 

See [placement documentation](doc_placement.md).

### Choose simulation medium

See [Medium Setup](doc_medium.md).

### Choose simulation type

See [Backend](doc_backend.md).

### Run the single_subject_pipeline

When the configuration is completely specified, the full pipeline only needs a subject id (a number consisting of no more than 3 digits). To run multiple subjects or setups in parallel, the `single_subject_pipeline` can be submitted using high performance computing jobs (see [HPC documentation](doc_hpc.md)).
