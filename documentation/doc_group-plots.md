# How to use the 'create_group_MNI_plots.m' script
## Reasons for using this function
- This will ensure that all figures will use the same scale
- Allows one to add a ROI MNI mask to visualise targeting accuracy
- Adds additional columns to each csv with information about intensity in the ROI mask, percentage of fwhm voxels within ROI etc.
- When using structural MRI's, it also allows one to normalise the brightness so the discrepancy in average brightness between individual T1's is reduced.

## Necessary input:
- subject_list: list of all subjects that you want to include in the figures.
- parameters: load the '.yaml' file using 'load_parameters.m' just as you would do in a pipeline.
- options: a structure similar to 'parameters'. Most of these are (as the name would suggest) optional. only one of the two following parameters needs to be selected:
	- options.slice_to_plot = 0 (give the number of the slice)
	- options.plot_max_intensity = 0 (turn on by changing to 1)

## Plotting a ROI mask
- Use a mask saved in MNI space
- Load it into the function under the option 'ROI_MNI_mask'

## Optional parameters:
    options.ROI_MNI_mask (:,:,:) # Loads in a 3d matrix of your ROI mask (in MNI space)
    options.slice_label = 'y' # Selects the axis along which your slice is made
    options.rotation = 90;  # Rotates your structural background figure
    options.plot_heating = 1 # Binary option to enable or disable heating plots
    options.outputs_suffix = '' # Allows you to add a string at the end of the name of your output to say, differentiate between plots of maximum intensity and a given slice
    options.isppa_thresholds = [] # Manually set the Ipa range that is to be plotted
    options.add_FWHM_boundary = 0 # Binary option to add a dotted boundary around the FWHM of intensity
    options.add_ROI_boundary = 1 # Binary option to add a red boundary around your ROI mask
    options.skip_missing = 0 # Binary option to skip subjects with missing files
    options.brightness_correction = 0 # Binary option to normalise the brightness between structural T1's
    options.average_target_brightness = 100 # Average targetted brightness, only used when brightness correction is enabled
