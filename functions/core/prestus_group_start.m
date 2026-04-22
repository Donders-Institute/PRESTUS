function prestus_group_start(subject_list, config, options)
% PRESTUS_GROUP_START  Entry point for group-level MNI-space analysis
%
% Loads a PRESTUS configuration (YAML file or already-loaded struct),
% applies any subject-list-level overrides, and calls create_group_MNI_plots
% to generate MNI-space summary plots and group CSV statistics across all
% specified subjects.
%
% This function is the group-level analogue of prestus_pipeline_start.
% It always runs locally (no HPC dispatch) because it only reads existing
% per-subject output files and does not perform heavy computation.
%
% Use as:
%   prestus_group_start(subject_list, config)
%   prestus_group_start(subject_list, config, Name=Value, ...)
%
% Input:
%   subject_list  - numeric array of subject IDs, e.g. [1 2 3 5 7]
%   config        - YAML config file path (char) OR a pre-loaded parameters
%                   struct. The loaded config must define at minimum:
%                     path.sim            — simulation output folder
%                     simulation.medium   — medium type ('layered'|'water')
%                     layers.brain        — SimNIBS label indices for brain
%                     layers.water        — SimNIBS label indices for water
%
% Name-value options (passed through to create_group_MNI_plots):
%   ROI_MNI_mask          - (3D logical) binary ROI in MNI space
%   slice_to_plot         - (int) fixed slice index to display
%   plot_max_intensity    - (logical) use slice with peak ISPPA
%   slice_label           - ('x'|'y'|'z') axis to slice along (default 'y')
%   rotation              - (deg) image rotation for display (default 90)
%   plot_heating          - (logical) include heating maps (default true)
%   outputs_suffix        - (char) suffix appended to output filenames
%   intensity_thresholds  - ([lo hi]) manual ISPPA colorbar range
%   add_FWHM_boundary     - (logical) outline ISPPA > half-max region
%   add_ROI_boundary      - (logical) overlay ROI boundary and compute stats
%   skip_missing          - (logical) skip subjects with missing files (default true)
%   brightness_correction - (logical) normalise brightness across subjects
%   average_target_brightness - (float) target mean brightness value
%
% Example:
%   % Using a YAML file
%   prestus_group_start(1:10, 'configs/config_study.yaml', ...
%       slice_label='y', plot_max_intensity=true, skip_missing=true);
%
%   % Using an existing parameters struct
%   prestus_group_start([1 3 5], parameters, add_FWHM_boundary=true);
%
% See also: CREATE_GROUP_MNI_PLOTS, PRESTUS_PIPELINE_START, LOAD_PARAMETERS

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

    %% Load parameters
    if ischar(config) || isstring(config)
        fprintf('Loading config: %s\n', config);
        parameters = load_parameters(config);
    elseif isstruct(config)
        parameters = config;
    else
        error('config must be a YAML file path (char) or a parameters struct.');
    end

    %% Print summary
    fprintf('========================================\n');
    fprintf('PRESTUS GROUP ANALYSIS\n');
    fprintf('Subjects: %s\n', num2str(subject_list));
    fprintf('Medium:   %s\n', parameters.simulation.medium);
    fprintf('Output:   %s\n', parameters.path.sim);
    fprintf('========================================\n\n');

    %% Build pass-through options struct for create_group_MNI_plots
    % Only pass options that were explicitly set (non-empty / non-default)
    % by forwarding the full options struct as name-value pairs.
    fnames = fieldnames(options);
    passthrough = {};
    for i = 1:numel(fnames)
        val = options.(fnames{i});
        if ~isempty(val)
            passthrough{end+1} = fnames{i}; %#ok<AGROW>
            passthrough{end+1} = val;       %#ok<AGROW>
        end
    end

    %% Run group analysis
    if isempty(passthrough)
        create_group_MNI_plots(subject_list, parameters);
    else
        create_group_MNI_plots(subject_list, parameters, passthrough{:});
    end

    fprintf('\nGroup analysis complete.\n');
end
