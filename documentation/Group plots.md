# How to use the 'create_group_MNI_plots.m' script
## Reasons for using this function
- This will ensure that all figures will use the same scale
- Allows one to add a ROI MNI mask to visualise targeting accuracy
- Adds additional columns to each csv with information about intensity in the ROI mask, percentage of fwhm voxels within ROI etc.
- When using structural MRI's, it also allows to normalise the brightness so the average brightness of the individual T1's is reduced.

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
    options.ROI_MNI_mask (:,:,:)
    options.slice_label = 'y'
    options.rotation = 90;
    options.plot_heating = 1
    options.outputs_suffix = ''
    options.isppa_thresholds = []
    options.add_FWHM_boundary = 0
    options.add_ROI_boundary = 1
    options.skip_missing = 0
