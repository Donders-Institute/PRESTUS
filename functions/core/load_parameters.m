function parameters = load_parameters(varargin)

% LOAD_PARAMETERS Loads and merges configuration files for simulation parameters.
%
% This function loads a default configuration file and optionally merges it with
% additional configuration files or structures provided as input arguments. It
% performs checks on the loaded parameters, calculates derived values, and sanitizes
% the output file affix. The function also verifies paths and transducer settings
% required for simulations.
%
% Input:
%   varargin - Optional input arguments:
%              * Single argument: Path to an extra configuration file or a struct.
%              * Two arguments: Path to an extra configuration file and its location.
%
% Output:
%   parameters - Struct containing merged simulation parameters.

    %% Load default configuration file

    parameters = yaml.loadFile('default_config.yaml', "ConvertToArray", true);

    %% Merge with additional configuration files or structures

    if nargin == 1
        extra_config_file = varargin{1};
        if isstruct(extra_config_file)
            extra_parameters = extra_config_file; % Use provided struct directly
        else
            extra_parameters = yaml.loadFile(fullfile(extra_config_file), "ConvertToArray", true);
        end
        parameters = MergeStruct(parameters, extra_parameters);
    elseif nargin == 2
        extra_config_file = varargin{1};
        extra_config_location = varargin{2};
        extra_parameters = yaml.loadFile(fullfile(extra_config_location, extra_config_file), "ConvertToArray", true);
        parameters = MergeStruct(parameters, extra_parameters);
    elseif nargin > 2
        % Check that extra inputs come in pairs (filename, location)
        if mod(numel(varargin), 2) ~= 0
            error('Extra parameters must be provided as name/value pairs.');
        end
        % Loop through each pair and assign to the structure
        for i = 1:2:numel(varargin)
            extra_config_file = varargin{i};
            extra_config_location = varargin{i+1};
            extra_parameters = yaml.loadFile(fullfile(extra_config_location, extra_config_file), "ConvertToArray", true);
            parameters = MergeStruct(parameters, extra_parameters);
        end
    end

    %% Check interactive mode requirements

    assert(parameters.simulation.interactive == 0 || usejava('desktop'), ...
           'MATLAB should run in desktop mode if parameters.simulation.interactive is enabled in PRESTUS config');

    %% Transducer settings validation and derived calculations

    for t_i = 1:numel(parameters.transducer)

        % Ensure source phase is set in radians or degrees
        if ~isfield(parameters.transducer(t_i), 'source_phase_rad')
            assert(isfield(parameters.transducer(t_i), 'source_phase_deg'), ...
                   'Source phase should be set in transducer parameters as source_phase_rad or source_phase_deg');
            parameters.transducer(t_i).source_phase_rad = parameters.transducer(t_i).source_phase_deg / 180 * pi;
        end

        % Calculate distance to transducer plane if not provided
        if ~isfield(parameters.transducer(t_i), 'dist_to_plane_mm')
            parameters.transducer(t_i).dist_to_plane_mm = sqrt(parameters.transducer(t_i).curv_radius_mm^2 - ...
                                                          (max(parameters.transducer(t_i).Elements_OD_mm) / 2)^2);
            fprintf('Distance to transducer plane is not provided, calculated as %.2f mm\n', ...
                    parameters.transducer(t_i).dist_to_plane_mm);
        end

        % Calculate distance between target and ep/bowl is not provided
        if ~isfield(parameters, 'expected_focal_distance_bowl') || ~isfield(parameters, 'expected_focal_distance_ep')
            parameters = focal_distance_calculation(parameters);
        end

        % Ensure source amplitude matches number of transducer elements
        if length(parameters.transducer(t_i).source_amp) == 1 && parameters.transducer(t_i).n_elements > 1
            parameters.transducer(t_i).source_amp = repmat(parameters.transducer(t_i).source_amp, [1, parameters.transducer(t_i).n_elements]);
        end

        % Evaluate source phase expressions if stored as cell arrays
        if iscell(parameters.transducer(t_i).source_phase_rad)
            for i = 1:length(parameters.transducer(t_i).source_phase_rad)
                if ~isnumeric(parameters.transducer(t_i).source_phase_rad{i})
                    parameters.transducer(t_i).source_phase_rad{i} = eval(parameters.transducer(t_i).source_phase_rad{i});
                end
            end
            parameters.transducer(t_i).source_phase_rad = cell2mat(parameters.transducer(t_i).source_phase_rad);
        end

        % Ensure source phase degrees are calculated from radians if not provided
        if ~isfield(parameters.transducer(t_i), 'source_phase_deg')
            parameters.transducer(t_i).source_phase_deg = parameters.transducer(t_i).source_phase_rad / pi * 180;
        end

    end

    %% Validate thermal simulation settings
    % Request no timing overview here

    if parameters.modules.run_heating_sims
        thermal_parameters(parameters, true);
    end

    %% Output file settings validation and sanitization

    % Sanitize output file affix to ensure valid characters only
    sanitized_affix = regexprep(parameters.io.output_affix, '[^a-zA-Z0-9_]', '_');
    if ~strcmp(sanitized_affix, parameters.io.output_affix)
        fprintf('The original `results_filename_affix` was sanitized. "%s" will be used instead of "%s"\n', ...
                sanitized_affix, parameters.io.output_affix);
        parameters.io.output_affix = sanitized_affix;
    end

    % Set output directory based on absolute or relative path
    if isfield(parameters.path, 'sim') && ~strcmp(parameters.path.sim, '')
        javaFileObj = java.io.File(parameters.path.sim); % Check path type (absolute/relative)
        if javaFileObj.isAbsolute()
            parameters.path.sim = fullfile(parameters.path.sim);
        else
            parameters.path.sim = fullfile(parameters.path.anat, parameters.path.sim);
        end
    else
        % Default output directory within data path
        parameters.path.sim = fullfile(parameters.path.anat, 'tussim');
    end

    %% Validate paths for required libraries and binaries

    % Check LD_LIBRARY_PATH existence and warn user if missing
    if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'ld_library_path') && ...
       ~strcmp(parameters.hpc.ld_library_path, "") && ~exist(parameters.hpc.ld_library_path, 'dir')
        assert(all(confirmation_dlg('The path `hpc.ld_library_path` has been specified but does not exist. Do you want to continue?', ...
                                    'Yes', 'No')), 'Exiting');
    end

    % Check segmentation software path existence and warn user if missing
    if isfield(parameters.startup, 'simnibs_bin_path') && ~strcmp(parameters.startup.simnibs_bin_path, "") && ...
       ~exist(fullfile(parameters.startup.simnibs_bin_path, 'charm'), 'file')
        assert(all(confirmation_dlg(sprintf('charm does not exist at %s. Do you want to continue?', ...
                                            parameters.startup.simnibs_bin_path), ...
                                    'Yes', 'No')), 'Exiting');
    elseif (~isfield(parameters.startup, 'simnibs_bin_path') || strcmp(parameters.startup.simnibs_bin_path, "")) && contains(parameters.simulation.medium, {'layered'})
        warning('No path to SimNIBS binaries provided. Segmentation and MNI-conversion may fail...');
    end

    %% Default segmentation path fallback

    if ~isfield(parameters.path, 'seg') || isempty(parameters.path.seg) || strcmp(parameters.path.seg, '')
        parameters.path.seg = parameters.path.anat;
    end

    %% Convert additional paths into cell arrays for processing

    % Convert `paths_to_add` into cell array format (split by semicolon)
    if isfield(parameters.startup, 'paths_to_add') && ~isempty(parameters.startup.paths_to_add)
        parameters.startup.paths_to_add = cellstr(strsplit(parameters.startup.paths_to_add, ';'));
    end

    % Convert `subpaths_to_add` into cell array format (split by semicolon)
    if isfield(parameters.startup, 'subpaths_to_add') && ~isempty(parameters.startup.subpaths_to_add)
        parameters.startup.subpaths_to_add = cellstr(strsplit(parameters.startup.subpaths_to_add, ';'));
    end

    %% MATLAB version check for compatibility

    % Warn user about outdated MATLAB versions (< R2022b)
    if verLessThan('matlab', '9.13')
       assert(all(confirmation_dlg('MATLAB appears to be outdated. Please update before continuing. Do you want to continue?', ...
                                   'Yes', 'No')), 'Exiting');
    end

    %% [DEBUG] Summary of requested parameters
    % With debugging off, the parameters are saved to the log at the
    % start of pipeline execution (see path_log_setup).

    if parameters.simulation.debug == 1
        fprintf('PRELIMINARY specified parameters. Note that this may not be the final specification...\n');
        print_parameter_summary(parameters)
    end

end
