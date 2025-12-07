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
    parameters = yaml.loadFile('configs/default_config.yaml', "ConvertToArray", true);

    %% Merge with additional configuration files or structures
    if nargin == 1
        extra_config_file = varargin{1};
        if isstruct(extra_config_file)
            extra_parameters = extra_config_file; % Use provided struct directly
        else
            extra_parameters = yaml.loadFile(fullfile('configs', extra_config_file), "ConvertToArray", true); 
        end        
        parameters = MergeStruct(parameters, extra_parameters);
    elseif nargin == 2
        extra_config_file = varargin{1};
        extra_config_location = varargin{2};
        extra_parameters = yaml.loadFile(fullfile(extra_config_location, extra_config_file), "ConvertToArray", true);
        parameters = MergeStruct(parameters, extra_parameters);
    elseif nargin > 2
        % Replace files in the config from the extra config
        extra_config_file = varargin{1};
        extra_config_location = varargin{2};
        extra_parameters = yaml.loadFile(fullfile(extra_config_location, extra_config_file), "ConvertToArray", true);
        parameters = MergeStruct(parameters, extra_parameters);

        % Process any additional name/value pairs to overwrite parameters
        % Get the extra arguments (starting at the 3rd argument)
        extraArgs = varargin(3:end);
        % Check that they come in pairs
        if mod(numel(extraArgs), 2) ~= 0
            error('Extra parameters must be provided as name/value pairs.');
        end
        % Loop through each pair and assign to the structure
        for i = 1:2:numel(extraArgs)
            fieldName = extraArgs{i};
            fieldValue = extraArgs{i+1};
            if contains(fieldName, '.')
                % Split the field name into its components
                parts = strsplit(fieldName, '.');
                % Update the nested field
                parameters = setfield(parameters, parts{:}, fieldValue);
            else
                % Update or add the top-level field
                parameters.(fieldName) = fieldValue;
            end
        end
    end

    %% Check interactive mode requirements
    assert(parameters.interactive == 0 || usejava('desktop'), ...
           'MATLAB should run in desktop mode if parameters.interactive is enabled in PRESTUS config');

    %% Transducer settings validation and derived calculations
    
    % To enable n transducers while maintaining backward comp of configs,
    % a single "transducer" field in parameters becomes "transducers(1)"

    if isfield(parameters, 'transducer')

        if isfield(parameters, 'transducers')
            error('the parameter file(s) include both fields transducer as well as transducers - only one of those fields is expected!');
        end
        parameters.transducers(1) = parameters.transducer;
        parameters.transducer = [];

    elseif nargin == 1
        % Warn user about missing transducer information
        assert(all(confirmation_dlg('The transducer info is missing in the configuration file. Do you want to continue?', 'Yes', 'No')), ...
               'Exiting');
        parameters.transducers = [];

    end

    for t_i = 1:numel(parameters.transducers)
        
        % Ensure source phase is set in radians or degrees
        if ~isfield(parameters.transducers(t_i), 'source_phase_rad')
            assert(isfield(parameters.transducers(t_i), 'source_phase_deg'), ...
                   'Source phase should be set in transducer parameters as source_phase_rad or source_phase_deg');
            parameters.transducers(t_i).source_phase_rad = parameters.transducers(t_i).source_phase_deg / 180 * pi;
        end

        % Calculate distance to transducer plane if not provided
        if ~isfield(parameters.transducers(t_i), 'dist_to_plane_mm')
            parameters.transducers(t_i).dist_to_plane_mm = sqrt(parameters.transducers(t_i).curv_radius_mm^2 - ...
                                                          (max(parameters.transducers(t_i).Elements_OD_mm) / 2)^2);
            fprintf('Distance to transducer plane is not provided, calculated as %.2f mm\n', ...
                    parameters.transducers(t_i).dist_to_plane_mm);
        end

        % Ensure source amplitude matches number of transducer elements
        if length(parameters.transducers(t_i).source_amp) == 1 && parameters.transducers(t_i).n_elements > 1
            parameters.transducers(t_i).source_amp = repmat(parameters.transducers(t_i).source_amp, [1, parameters.transducers(t_i).n_elements]);
        end

        % Evaluate source phase expressions if stored as cell arrays
        if iscell(parameters.transducers(t_i).source_phase_rad)
            for i = 1:length(parameters.transducers(t_i).source_phase_rad)
                if ~isnumeric(parameters.transducers(t_i).source_phase_rad{i})
                    parameters.transducers(t_i).source_phase_rad{i} = eval(parameters.transducers(t_i).source_phase_rad{i});
                end
            end
            parameters.transducers(t_i).source_phase_rad = cell2mat(parameters.transducers(t_i).source_phase_rad);
        end

        % Ensure source phase degrees are calculated from radians if not provided
        if ~isfield(parameters.transducers(t_i), 'source_phase_deg')
            parameters.transducers(t_i).source_phase_deg = parameters.transducers(t_i).source_phase_rad / pi * 180;
        end

    end

    %% Derived grid settings
    parameters.grid_step_m = parameters.grid_step_mm / 1e3; % Convert grid step from mm to meters

    % Set simulation dimensions based on default grid dimensions or fallback to 3D
    if ~isfield(parameters, 'n_sim_dims')
        if isfield(parameters, 'default_grid_dims')
            parameters.n_sim_dims = length(parameters.default_grid_dims);
        else
            parameters.n_sim_dims = 3;
        end            
    end

    % Set default grid dimensions based on grid size if not explicitly defined
    if ~isfield(parameters, 'default_grid_dims')
        parameters.default_grid_dims = repmat(parameters.default_grid_size, [1, parameters.n_sim_dims]);
    end

    %% Validate thermal simulation settings
    if parameters.run_heating_sims
        check_thermal_parameters(parameters);
    end

    %% Output file settings validation and sanitization
    
    % Sanitize output file affix to ensure valid characters only
    sanitized_affix = regexprep(parameters.results_filename_affix, '[^a-zA-Z0-9_]', '_');
    if ~strcmp(sanitized_affix, parameters.results_filename_affix)
        fprintf('The original `results_filename_affix` was sanitized. "%s" will be used instead of "%s"\n', ...
                sanitized_affix, parameters.results_filename_affix);
        parameters.results_filename_affix = sanitized_affix;
    end

    % Set output directory based on absolute or relative path
    if isfield(parameters, 'output_location')
        javaFileObj = java.io.File(parameters.output_location); % Check path type (absolute/relative)
        if javaFileObj.isAbsolute()
            parameters.sim_path = fullfile(parameters.output_location);
        else
            parameters.sim_path = fullfile(parameters.data_path, parameters.output_location);
        end
    else
        % Default output directory within data path
        parameters.sim_path = fullfile(parameters.data_path, 'sim_outputs/');
    end
 
    %% Validate paths for required libraries and binaries
    
    % Check LD_LIBRARY_PATH existence and warn user if missing
    if isfield(parameters, 'ld_library_path') && ~exist(parameters.ld_library_path, 'dir')
        assert(all(confirmation_dlg('The path in `parameters.ld_library_path` does not exist. Do you want to continue?', ...
                                    'Yes', 'No')), 'Exiting');
    end

    % Check segmentation software path existence and warn user if missing
    if isfield(parameters, 'simnibs_bin_path') && ~exist(fullfile(parameters.simnibs_bin_path, parameters.segmentation_software), 'file')
        assert(all(confirmation_dlg(sprintf('The segmentation software (%s) does not exist at %s. Do you want to continue?', ...
                                            parameters.segmentation_software, parameters.simnibs_bin_path), ...
                                    'Yes', 'No')), 'Exiting');
    end

    %% Default segmentation path fallback
    if ~isfield(parameters, 'seg_path') || isempty(parameters.seg_path)
        parameters.seg_path = parameters.data_path;
    end

    %% Default: deactivate pseudoCT unless specified
    if ~isfield(parameters, 'usepseudoCT')
        parameters.usepseudoCT = 0;
    end
    
    %% Convert additional paths into cell arrays for processing
    
    % Convert `paths_to_add` into cell array format (split by semicolon)
    if isfield(parameters, 'paths_to_add') && ~isempty(parameters.paths_to_add)
        parameters.paths_to_add = cellstr(strsplit(parameters.paths_to_add, ';'));
    end
    
    % Convert `subpaths_to_add` into cell array format (split by semicolon)
    if isfield(parameters, 'subpaths_to_add') && ~isempty(parameters.subpaths_to_add)
        parameters.subpaths_to_add = cellstr(strsplit(parameters.subpaths_to_add, ';'));
    end
    
    %% MATLAB version check for compatibility
    
    % Warn user about outdated MATLAB versions (< R2022b)
    if verLessThan('matlab', '9.13')
       assert(all(confirmation_dlg('MATLAB appears to be outdated. Please update before continuing. Do you want to continue?', ...
                                   'Yes', 'No')), 'Exiting');
    end

end
