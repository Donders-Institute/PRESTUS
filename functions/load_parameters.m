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
    if isfield(parameters, 'transducer')
        
        % This section validates and initializes transducer element parameters
        
        % Check if older config that only supports annular arrays is used
        if isfield(parameters.transducer, 'array_shape')
            % It handles two primary types of transducers: matrix and annular arrays

            % Ensure the array_shape type field exists in the parameters structure
            assert(isfield(parameters.transducer.array_shape, 'type'), ...
                'ERROR: Missing type field. Please specify either "matrix" or "annular" in transducer.array_shape.');

            if strcmp(parameters.transducer.array_shape.type, 'matrix')
                % Handle matrix transducer configuration

                fprintf('Matrix transducer detected. Using kWaveArray (set to 1).\n');
                parameters.use_kWaveArray = 1;

                assert(isfield(parameters.transducer.array_shape.matrix, 'steering'), ...
                        'ERROR: Missing steering definition. Choose "1D" or "3D".');

                if strcmp(parameters.transducer.array_shape.matrix.steering, '1D')
                    parameters.transducer.perform_focus_rotation = true;
                elseif strcmp(parameters.transducer.array_shape.matrix.steering, '3D')
                    parameters.transducer.perform_focus_rotation = false;
                else
                   error('Steering option "%s" is not implemented.', ...
                    parameters.transducer.array_shape.matrix.steering);
                end
                
                % Initialize phase to zero degrees by default
                % Later phases are calculated based on set focus
                parameters.transducer.array_shape.matrix.source_phase_deg = 0;
                parameters.transducer.source_phase_deg = parameters.transducer.array_shape.matrix.source_phase_deg;

                % Verify that matrix_shape parameter exists
                assert(isfield(parameters.transducer.array_shape.matrix, 'matrix_shape'), ...
                    'ERROR: Missing matrix_shape field for matrix transducer. Please specify.');

                assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape, 'type'), ...
                    'ERROR: Missing matrix_shape type field for matrix transducer. Please specify "define_here"');

                assert(isfield(parameters.transducer.array_shape.matrix, 'elem_height'), ...
                    'ERROR: Missing elem_height parameter for grid. This defines height of each element.');

                assert(isfield(parameters.transducer.array_shape.matrix, 'elem_width'), ...
                    'ERROR: Missing elem_width parameter for grid. This defines width of each element.');

                parameters.transducer.elem_height = parameters.transducer.array_shape.matrix.elem_height;
                parameters.transducer.elem_width = parameters.transducer.array_shape.matrix.elem_width;

                if strcmp(parameters.transducer.array_shape.matrix.matrix_shape.type, 'define_here')
                    % Custom matrix definition path

                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape, 'define_here'), ...
                        'ERROR: Missing additional parameters for define_here.');

                    % Check for required outer diameter parameter
                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here, 'outer_diameter_mm'), ...
                        'ERROR: Missing outer_diameter_mm for custom matrix definition. This parameter is required.');

                    end_elem = parameters.transducer.array_shape.matrix.matrix_shape.define_here.outer_diameter_mm;

                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here, 'element'), ...
                        'ERROR: Missing array_shape field for custom matrix. Please specify "rect", "disc" or "bowl".');

                    parameters.transducer.element = parameters.transducer.array_shape.matrix.matrix_shape.define_here.element;

                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here, 'is_curved'), ...
                        'ERROR: Missing curved field for custom matrix. Please specify.');

                    if parameters.transducer.array_shape.matrix.matrix_shape.define_here.is_curved
                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here, 'curved'), ...
                        'ERROR: Missing additional information of how it is curved.');

                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.curved, 'curv_radius_mm'), ...
                            'ERROR: Missing curv_radius_mm field for custom matrix. Please specify radius of curvature.');
                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.curved, 'dist_to_plane_mm'), ...
                            'ERROR: Missing dist_to_plane_mm field for custom matrix. Please specify distance from exit plane to geometric focus.');

                        parameters.transducer.curv_radius_mm = parameters.transducer.array_shape.matrix.matrix_shape.define_here.curved.curv_radius_mm;
                        parameters.transducer.dist_to_plane_mm = parameters.transducer.array_shape.matrix.matrix_shape.define_here.curved.dist_to_plane_mm;
                    else
                        parameters.transducer.curv_radius_mm = 0;

                        % Officially distance to plane when the transducer is
                        % flat approaches infinity, so set a fixed value for
                        % visualization purposes.
                        parameters.transducer.dist_to_plane_mm = 70;
                    end

                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here, 'grid_shape'), ...
                        'ERROR: Missing grid_shape for grid configuration.');

                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape, 'type'), ...
                        'ERROR: Missing additional grid_shape parameters like type for grid configuration.');

                    if strcmp(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.type, 'rect')
                        % Validate required parameters for grid
                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect, 'n_elem_row'), ...
                            'ERROR: Missing n_elem_row parameter for grid. This defines number of rows.');

                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect, 'n_elem_col'), ...
                            'ERROR: Missing n_elem_col parameter for grid. This defines number of columns.');

                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect, 'elem_spacing_height'), ...
                            'ERROR: Missing elem_spacing_height parameter for grid. This defines height spacing between elements.');

                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect, 'elem_spacing_width'), ...
                            'ERROR: Missing elem_spacing_width parameter for grid. This defines width spacing between elements.');

                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect, 'sparsity_factor'), ...
                            'ERROR: Missing sparsity_factor parameter for grid. This defines % of used elements.');

                        % Get element count directly from configuration
                        parameters.transducer.sparsity_factor = parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect.sparsity_factor;

                        % Extract grid dimensions
                        n_elem_row = parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect.n_elem_row;
                        n_elem_col = parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect.n_elem_col;

                        % Store element dimensions in standard location
                        parameters.transducer.n_elem_row = n_elem_row;
                        parameters.transducer.n_elem_col = n_elem_col;

                        % Calculate initial element count (will be adjusted later for circular cutout)
                        n_elem = n_elem_col * n_elem_row;

                        parameters.transducer.elem_spacing_height = parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect.elem_spacing_height;
                        parameters.transducer.elem_spacing_width = parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.rect.elem_spacing_width;

                    elseif strcmp(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.type, 'fibonacci')
                        % Sparse spiral grid configuration

                        % Validate required parameters for sparser grid
                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape, 'fibonacci'), ...
                            'ERROR: Missing fibonacci field in grid_shape for sparser grid configuration.');

                        assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.fibonacci, 'n_elements'), ...
                            'ERROR: Missing sparsity_factor parameter for grid. This defines % of used elements.');

                        n_elem = parameters.transducer.array_shape.matrix.matrix_shape.define_here.grid_shape.fibonacci.n_elements;
                    end
                    
                    if isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here, 'is_clover')
                        parameters.transducer.is_clover = parameters.transducer.array_shape.matrix.matrix_shape.define_here.is_clover;
                        
                        if parameters.transducer.is_clover
                            assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here, 'clover'), ...
                                'ERROR: Missing additional information of Clover setup.');

                            assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.clover, 'n_leaves'), ...
                                'ERROR: Missing additional information of Clover setup. Define number of leaves.');

                            assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.define_here.clover, 'ROC_parent'), ...
                                'ERROR: Missing additional information of Clover setup. Define ROC of parent.');

                            parameters.transducer.n_leaves = parameters.transducer.array_shape.matrix.matrix_shape.define_here.clover.n_leaves;
                            parameters.transducer.ROC_parent = parameters.transducer.array_shape.matrix.matrix_shape.define_here.clover.ROC_parent;
                        end
                    else
                        parameters.transducer.is_clover = false;
                    end

                    % Store element count in standard location for compatibility with existing code
                    parameters.transducer.n_elements = n_elem;

                    % Calculate inner and outer diameters for elements based on total diameter and count
                    [id, od] = calc_elements_id_od_mm(end_elem, parameters.transducer.n_elements);

                    % Store element dimensions in standard location for visualization compatibility
                    parameters.transducer.Elements_ID_mm = id;
                    parameters.transducer.Elements_OD_mm = od;


                elseif strcmp(parameters.transducer.array_shape.matrix.matrix_shape.type, 'stanford')
                    % Stanford transducer configuration (predefined matrix layout)

                    % Validate required parameters for Stanford configuration
                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape, 'stanford'), ...
                        'ERROR: Missing stanford field for Stanford transducer configuration.');

                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.stanford, 'n_elements'), ...
                        'ERROR: Missing n_elements parameter for Stanford transducer. This defines element count.');

                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.stanford, 'outer_diameter_mm'), ...
                        'ERROR: Missing outer_diameter_mm parameter for Stanford transducer. This defines overall size.');

                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.stanford, 'curv_radius_mm'), ...
                        'ERROR: Missing curv_radius_mm field for Stanford transducer. Please specify radius of curvature.');

                    assert(isfield(parameters.transducer.array_shape.matrix.matrix_shape.stanford, 'dist_to_plane_mm'), ...
                        'ERROR: Missing dist_to_plane_mm field for Stanford transducer. Please specify distance from exit plane to geometric focus.');

                    % Store element count in standard location for compatibility with existing code
                    parameters.transducer.n_elements = parameters.transducer.array_shape.matrix.matrix_shape.stanford.n_elements;

                    end_elem = parameters.transducer.array_shape.matrix.matrix_shape.stanford.outer_diameter_mm;

                    % Calculate inner and outer diameters for elements
                    [id, od] = calc_elements_id_od_mm(end_elem, parameters.transducer.n_elements);

                    % Store element dimensions in standard location for visualization compatibility
                    parameters.transducer.Elements_ID_mm = id;
                    parameters.transducer.Elements_OD_mm = od;

                    parameters.transducer.curv_radius_mm = parameters.transducer.array_shape.matrix.matrix_shape.stanford.curv_radius_mm;
                    parameters.transducer.dist_to_plane_mm = parameters.transducer.array_shape.matrix.matrix_shape.stanford.dist_to_plane_mm;
                else
                    error('Matrix shape option "%s" is not implemented.', ...
                          parameters.transducer.array_shape.matrix.matrix_shape.type);
                end
            elseif strcmp(parameters.transducer.array_shape.type, 'annular')
                % Handle annular transducer configuration
                assert(isfield(parameters.transducer.array_shape, 'annular'), ...
                    'ERROR: Appropriate configuration for annular transducer is missing.');

                assert(isfield(parameters.transducer.array_shape.annular, 'Elements_ID_mm'), ...
                    'ERROR: Missing Elements_ID_mm for annular transducer. This defines inner diameters of elements.');

                assert(isfield(parameters.transducer.array_shape.annular, 'Elements_OD_mm'), ...
                    'ERROR: Missing Elements_OD_mm for annular transducer. This defines outer diameters of elements.');

                assert(isfield(parameters.transducer.array_shape.annular, 'n_elements'), ...
                    'ERROR: Missing n_elements parameter for annular transducer. This defines element count.');

                assert(isfield(parameters.transducer.array_shape.annular, 'curv_radius_mm'), ...
                    'ERROR: Missing curv_radius_mm field for annular transducer. Please specify radius of curvature.');

                assert(isfield(parameters.transducer.array_shape.annular, 'dist_to_plane_mm'), ...
                    'ERROR: Missing dist_to_plane_mm field for annular transducer. Please specify distance from exit plane to geometric focus.');

                assert(isfield(parameters.transducer.array_shape.annular, 'source_phase_deg'), ...
                    'ERROR: Missing source_phase_deg field for annular transducer. Please specify phases.');

                % Copy element dimensions to standard location for compatibility
                parameters.transducer.Elements_ID_mm = parameters.transducer.array_shape.annular.Elements_ID_mm;
                parameters.transducer.Elements_OD_mm = parameters.transducer.array_shape.annular.Elements_OD_mm;

                % Copy element count to standard location for compatibility
                parameters.transducer.n_elements = parameters.transducer.array_shape.annular.n_elements;

                parameters.transducer.curv_radius_mm = parameters.transducer.array_shape.annular.curv_radius_mm;
                parameters.transducer.dist_to_plane_mm = parameters.transducer.array_shape.annular.dist_to_plane_mm;

                parameters.transducer.source_phase_deg = parameters.transducer.array_shape.annular.source_phase_deg;
            else
                error('Element shape option "%s" is not implemented.', ...
                      parameters.transducer.array_shape.type);
            end
        else
            parameters.transducer.array_shape = 'annular';
        end

        if strcmp(parameters.transducer.array_shape, 'annular')
            % Ensure source phase is set in radians or degrees
            if ~isfield(parameters.transducer, 'source_phase_rad')
                assert(isfield(parameters.transducer, 'source_phase_deg'), ...
                       'Source phase should be set in transducer parameters as source_phase_rad or source_phase_deg');
                parameters.transducer.source_phase_rad = parameters.transducer.source_phase_deg / 180 * pi;
            end
    
            % Calculate distance to transducer plane if not provided
            if ~isfield(parameters.transducer, 'dist_to_plane_mm')
                parameters.transducer.dist_to_plane_mm = sqrt(parameters.transducer.curv_radius_mm^2 - ...
                                                              (max(parameters.transducer.Elements_OD_mm) / 2)^2);
                fprintf('Distance to transducer plane is not provided, calculated as %.2f mm\n', ...
                        parameters.transducer.dist_to_plane_mm);
            end
    
            % Ensure source amplitude matches number of transducer elements
            if length(parameters.transducer.source_amp) == 1 && parameters.transducer.n_elements > 1
                parameters.transducer.source_amp = repmat(parameters.transducer.source_amp, [1, parameters.transducer.n_elements]);
            end
    
            % Evaluate source phase expressions if stored as cell arrays
            if iscell(parameters.transducer.source_phase_rad)
                for i = 1:length(parameters.transducer.source_phase_rad)
                    if ~isnumeric(parameters.transducer.source_phase_rad{i})
                        parameters.transducer.source_phase_rad{i} = eval(parameters.transducer.source_phase_rad{i});
                    end
                end
                parameters.transducer.source_phase_rad = cell2mat(parameters.transducer.source_phase_rad);
            end
    
            % Ensure source phase degrees are calculated from radians if not provided
            if ~isfield(parameters.transducer, 'source_phase_deg')
                parameters.transducer.source_phase_deg = parameters.transducer.source_phase_rad / pi * 180;
            end
        end
    else
        % Warn user about missing transducer information
        assert(all(confirmation_dlg('The transducer info is missing in the configuration file. Do you want to continue?', 'Yes', 'No')), ...
               'Exiting');
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

    if ~isfield(parameters, 'use_kWaveArray')
        parameters.use_kWaveArray = 0;
    end 

end
