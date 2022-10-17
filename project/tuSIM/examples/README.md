# How to design your own processing pipeline
## Packages to install[^Ensure the paths are created in Matlab]
- SimNIBS.
- k-wave.

# Making a config file
- Config files are created in the .yaml format.
- The 'default_config' contains all possible parameters that are accepted by the single_subject_pipeline.

## Some fixes for possible issues that may arise:
- Issues that arise in the preprocessing_brain function can arise when the SimNIBS file has not been manually adjusted.
	- To fix this, update SimNIBS to include our changes to the 'headreco' file.