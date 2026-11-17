function prestus_group_start(subject_list, config, options)
% PRESTUS_GROUP_START  Entry point for group-level analyses in MNI space
%
% Builds the group MNI maps (CREATE_GROUP_MAPS: focal-volume overlap,
% mean/SD intensity, mean heating) and the group HTML report with
% targeting accuracy via PRESTUS_GROUP_REPORT_START. The plot options of
% earlier versions are accepted but ignored; for an ROI, set
% analysis.target_roi_mni to an MNI NIfTI mask in the config.
%
% Use as:
%   prestus_group_start(subject_list, config)
%
% Input:
%   subject_list  - numeric array of subject IDs, e.g. [1 2 3 5 7]
%   config        - YAML config file path (char) OR a pre-loaded parameters struct
%   options       - plot options of earlier versions (ignored)
%
% See also: PRESTUS_GROUP_REPORT_START, CREATE_GROUP_MAPS, PRESTUS_STUDY_START

arguments
    subject_list (1,:) {mustBeNumeric}
    config
    options.ROI_MNI_mask
    options.slice_to_plot
    options.plot_max_intensity
    options.slice_label
    options.rotation
    options.plot_heating
    options.outputs_suffix
    options.intensity_thresholds
    options.add_FWHM_boundary
    options.add_ROI_boundary
    options.skip_missing
    options.brightness_correction
    options.average_target_brightness
end

    if ~isempty(fieldnames(options))
        warning('prestus_group_start:ignoredOptions', ...
            'Plot options are ignored; set analysis.target_roi_mni for an ROI.');
    end
    prestus_group_report_start(config, subject_list);
end
