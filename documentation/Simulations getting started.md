# How to design your own processing pipeline
## Packages to install (ensure the paths and subfolders are added in Matlab)
- SimNIBS 4.
- k-wave 1.4 (place in PRESTUS' toolbox folder).

## Start with the tutorial (examples/tuSIM_intro_tutorial.mlx) which visualises the steps you and the wrapper completes to do the simulations:
- The acoustic profile is derived from the manufacturers measurements and adjusted to more closely match the actual characteristics of the transducer.
- Segmentation of the structural data is then conducted by SimNIBS.
- Finally, the segmented data is translated to a matrix of acoustic properties by k-wave, after which the same program conducts the simulations.

## Run the single_subject_pipeline
- The only things the pipeline needs is a subject id (a number consisting of no more than 3 digits) and a set of parameters presented in a config file (placed in the 'configs' folder) and imported using 'load_parameters.m'.
- The different parameters that can be used are found in the 'default_config.yaml'.
    - It should be noted that the default config should not be altered. Instead, when you want to change a parameter copy that parameter to your own config file and change it there.

### Making a config file
- Config files are created in the .yaml format.
- Some of the parameters found in the 'default_config' are mandatory for the pipeline to run while others can be left out based on the requirements of the analysis.
    - The default config shows which ones should and don't have to be changed.

## Importing a transducer location
- These can be specified in your config [T1_grid_space].
- Another option is to import transducer locations from Localite and translate them to the correct space using the `get_trans_pos_from_trigger_markers` and `ras_to_grid` functions.