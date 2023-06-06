# How to plot your entire group together
- This will ensure that all figures will be adjusted to the legend
- All figures are transformed into MNI space

## Necessary input:
- outputs_path: full path to the location of your output.
- subject_list: list of all subjects that you want to include in the figures.
- parameters: load the '.yaml' file using 'load_parameters.m' just as you would do in a pipeline.
- options: a structure similar to 'parameters'. Most of these are (as the name would suggest) optional. only one of the two following parameters needs to be selected:
	- options.slice_to_plot = 0 (give the number of the slice)
	- options.plot_max_intensity = 0 (turn on by changing to 1)

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
