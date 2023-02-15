# How to design your own processing pipeline
## Packages to install (ensure the paths and subfolders are added in Matlab)
- SimNIBS 4 or higher.
- k-wave 1.4 or higher.

## Follow the tutorial (examples/tuSIM_intro_tutorial.mlx) which explains how:
- The steps being taken by the wrapper.
- The acoustic profile is derived from the manufacturer specs and adjusted to more closely match the actual characteristics of the transducer.
- Preprocessing of the structural data is done using SimNIBS.
- All parameters are loaded into k-wave.

## Run the single_subject_pipeline
- The only things the pipeline needs is a subject id (a number consisting of no more than 3 digits) and a set of parameters presented in a config file (placed in the 'configs' folder) and imported using 'load_parameters.m'.
- The different parameters that can be used are found in the 'default_config.yaml'.

### Making a config file
- Config files are created in the .yaml format.
- Some of the parameters found in the 'default_config' are mandatory for the pipeline to run while others can be left out based on the requirements of the analysis.
