function prestus_group_start(subject_list, config, options) %#ok<INUSD>
% PRESTUS_GROUP_START  DEPRECATED entry point for group-level MNI-space analysis
%
% Deprecated: its back-end, CREATE_GROUP_MNI_PLOTS, predates the current
% output layout (NIfTI names in nii/, sub-XXX_<medium><affix>.csv) and no
% longer finds the per-subject files. This function now warns and builds
% the group MNI maps (CREATE_GROUP_MAPS: focal-volume overlap, mean/SD
% intensity, mean heating) and the group HTML report with targeting
% accuracy via PRESTUS_GROUP_REPORT_START. The old plot options are
% accepted but ignored; for an ROI, set analysis.target_roi_mni to an MNI
% NIfTI mask in the config.
%
% Use as:
%   prestus_group_start(subject_list, config)
%
% Input:
%   subject_list  - numeric array of subject IDs, e.g. [1 2 3 5 7]
%   config        - YAML config file path (char) OR a pre-loaded parameters struct
%   options       - legacy name-value options (ignored)
%
% See also: PRESTUS_GROUP_REPORT_START, CREATE_GROUP_MAPS, PRESTUS_STUDY_START

arguments
    subject_list (1,:) {mustBeNumeric}
    config
    options.ROI_MNI_mask = []
    options.slice_to_plot = 0
    options.plot_max_intensity = 0
    options.slice_label = 'y'
    options.rotation = 90
    options.plot_heating = 1
    options.outputs_suffix = ''
    options.intensity_thresholds = []
    options.add_FWHM_boundary = 0
    options.add_ROI_boundary = 0
    options.skip_missing = 1
    options.brightness_correction = 0
    options.average_target_brightness = []
end

    warning('prestus_group_start:deprecated', ...
        ['prestus_group_start / create_group_MNI_plots are deprecated; running ' ...
         'prestus_group_report_start (group MNI maps + group report) instead. ' ...
         'Plot options are ignored; set analysis.target_roi_mni for an ROI.']);
    prestus_group_report_start(config, subject_list);
end
