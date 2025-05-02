# How to design your own processing pipeline

### install packages (ensure the paths and subfolders are added in Matlab)
- SimNIBS 4
- k-wave 1.4 (place in PRESTUS' toolbox folder)

### start with the tutorial (documentation/PRESTUS_intro_tutorial.md)
- The acoustic profile is derived from the manufacturers measurements and adjusted to more closely match the actual characteristics of the transducer.
- The structural data is segemented with SimNIBS 4's charm.
- The segmented data is translated to a matrix of acoustic properties.
- k-wave performs acoustic and thermal simulations.
- Select outputs are plotted.

### create a study-specific config file
- Configuration metadata must be fully specified in (a set of) config files (in the .yaml format) placed in the 'configs' folder (see ```doc_config.md```). 
- The ```load_parameters.m``` function imports the configuration. 
- Some of the parameters found in the 'default_config' are mandatory for the pipeline to run while others can be left out based on the requirements of the analysis.
- The different parameters that can be used are found in the 'default_config.yaml'. The default config should generally not be altered. Instead, parameter changes should be specified via an application-specific config file or before performing the subject call in MATLAB).

### specify or import a transducer location
- These can be specified in your config [T1_grid_space].
- Another option is to import transducer locations from Localite and translate them to the correct space using the `get_trans_pos_from_trigger_markers` and `ras_to_grid` functions.

### choose a simulation medium

- see doc_medium.md

### run the single_subject_pipeline
- Ideally, the full pipeline only needs a subject id (a number consisting of no more than 3 digits).
- To run multiple subjects or setups in parallel, the single_subject_pipeline can be submitted using high performance computing jobs (see doc_hpc.md).
