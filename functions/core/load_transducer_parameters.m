function parameters = load_transducer_parameters(parameters)

% LOAD_TRANSDUCER_PARAMETERS Validates and normalizes transducer parameters.
%
% This function processes the 'transducer' field of the input parameters struct.
% It validates required fields, ensures backward compatibility with legacy 
% annular-only configurations, and normalizes parameters into a consistent format 
% for simulation use. For matrix transducers, it validates element geometry, 
% curvature, steering, grid definitions, and optional Clover or sparse grid setups. 
% For annular transducers, it validates element diameters, curvature, and phases. 
% The function also calculates derived values such as source phases in radians, 
% element counts, and distances to the transducer plane.
%
% Input:
%   parameters - Struct containing simulation parameters, expected to include
%                a 'transducer' field with one or more transducer definitions.
%
% Output:
%   parameters - Struct with validated and normalized transducer parameters,
%                including additional derived fields such as:
%                * n_elements, Elements_ID_mm, Elements_OD_mm
%                * source_phase_rad, source_phase_deg
%                * curv_radius_mm, dist_to_plane_mm
%                * n_elem_row, n_elem_col (for matrix grids)
%                * sparsity_factor, clover setup fields, etc.

    if isfield(parameters, 'transducer')
    
        for t_i = 1:numel(parameters.transducer)
    
            tr = parameters.transducer(t_i);
    
            % This section validates and initializes transducer element parameters
    
            % Check if older config that only supports annular arrays is used
            if isfield(tr, 'array_shape')
                % It handles two primary types of transducers: matrix and annular arrays
    
                % Ensure the array_shape type field exists in the parameters structure
                if ~isfield(tr.array_shape, 'type')
                    error('Transducer %i; Missing type field. Please specify either "matrix" or "annular".', t_i);
                end
    
    
                switch tr.array_shape.type
    
                    case 'matrix'
                        % Handle matrix transducer configuration
    
                        matrix_tr = tr.array_shape.matrix;
    
                        fprintf('Matrix transducer detected. Using kWaveArray (set to 1).\n');
                        parameters.use_kWaveArray = 1;
    
                        % Steering
                        if ~isfield(matrix_tr, 'steering')
                            error('Transducer %i; Missing steering definition. Choose "1D" or "3D".', t_i);
                        end
    
                        if strcmp(matrix_tr.steering, '1D')
                            tr.align_transducer_with_focus = true;
                        elseif strcmp(matrix_tr.steering, '3D')
                            tr.align_transducer_with_focus = false;
                        else
                            error('Transducer %i; Steering option "%s" is not implemented.', ...
                                t_i, matrix_tr.steering);
                        end
    
                        % Initialize phase to zero degrees by default
                        % Later phases are calculated based on set focus
                        matrix_tr.source_phase_deg = 0;
                        tr.source_phase_deg = matrix_tr.source_phase_deg;
    
                        % Element shape & other properties
                        if ~isfield(matrix_tr, 'element_shape')
                            error('Transducer %i; Missing element_shape field for matrix array. Please specify "rect", "disc" or "bowl".', t_i);
                        end
    
                        if ~isfield(matrix_tr, 'elem_height_mm')
                            error('Transducer %i; Missing elem_height parameter. Please specify to define height of each element.', t_i);
                        end
    
                        if ~isfield(matrix_tr, 'elem_width_mm')
                            error('Transducer %i; Missing elem_width parameter. Please specify to define width of each element.', t_i);
                        end
                        tr.elem_height_mm = matrix_tr.elem_height_mm;
                        tr.elem_width_mm = matrix_tr.elem_width_mm;
    
                        % Outer diameter
                        if ~isfield(matrix_tr, 'outer_diameter_mm')
                            error('Transducer %i; Missing outer_diameter_mm parameter for matrix array. Please specify to define overall size.', t_i);
                        end
    
                        % Curvature
                        if ~isfield(matrix_tr, 'is_curved')
                            error('Transducer %i; Missing is_curved field for matrix transducer. Please specify.', t_i);
                        end
    
                        if matrix_tr.is_curved
                            if ~isfield(matrix_tr.is_curved, 'curv_radius_mm')
                                error('Transducer %i; Missing curv_radius_mm field for matrix transducer. Please specify radius of curvature.', t_i);
                            end
    
                            if ~isfield(matrix_tr.is_curved, 'dist_to_plane_mm')
                                error('Transducer %i; Missing dist_to_plane_mm field for matrix transducer. Please specify distance from exit plane to geometric focus.', t_i);
                            end
    
                            tr.curv_radius_mm = matrix_tr.is_curved.curv_radius_mm;
                            tr.dist_to_plane_mm = matrix_tr.is_curved.dist_to_plane_mm;
    
                        else
                            tr.curv_radius_mm = 0;
    
                            % Distance to plane when transducer is flat
                            % approaches infinity, so set a fixed value for
                            % visualization purposes.
                            tr.dist_to_plane_mm = 70;
                        end
    
                        % Clover setup
                        if isfield(matrix_tr, 'is_clover_setup')
                            tr.is_clover_setup = matrix_tr.is_clover_setup;
    
                            if tr.is_clover_setup
                                if ~isfield(matrix_tr, 'clover')
                                    error('Transducer %i; Missing additional information of Clover setup.', t_i);
                                end
    
                                if ~isfield(matrix_tr.clover, 'n_leaves')
                                    error('Transducer %i; Missing additional information of Clover setup. Define number of leaves.', t_i);
                                end
    
                                if ~isfield(matrix_tr.clover, 'ROC_parent')
                                    error('Transducer %i; Missing additional information of Clover setup. Define ROC of parent.', t_i);
                                end
    
                                tr.n_leaves = matrix_tr.clover.n_leaves;
                                tr.ROC_parent = matrix_tr.clover.ROC_parent;
                            end
                            % Parameter isn't found, so don't create Clover setup
                        else
                            tr.is_clover_setup = false;
                        end
    
                        % Matrix grid definition
                        if ~isfield(matrix_tr.matrix_shape, 'type')
                            error('Transducer %i; Missing matrix_shape type field for matrix transducer. Please specify "define_here" or "extract_from_file"', t_i);
                        end
    
                        switch matrix_tr.matrix_shape.type
                            case 'define_here'
                                if ~isfield(matrix_tr.matrix_shape, 'define_here')
                                    error('Transducer %i; Missing define_here field and properties for matrix transducer. Please specify.', t_i);
                                end

                                define_shape_here = matrix_tr.matrix_shape.define_here;

                                if ~isfield(define_shape_here, 'grid_shape')
                                    error('Transducer %i; Missing grid_shape field and properties for matrix transducer. Please specify.', t_i);
                                end
    
                                if ~isfield(define_shape_here.grid_shape, 'type')
                                    error('Transducer %i; Missing grid_shape type field and properties for matrix transducer. Please specify.', t_i);
                                end
    
                                switch define_shape_here.grid_shape.type
                                    case 'rect'
                                        if ~isfield(define_shape_here.grid_shape, 'rect')
                                            error('Transducer %i; Missing rect field in grid_shape for rectangular grid configuration.', t_i);
                                        end

                                        rect_grid = define_shape_here.grid_shape.rect;

                                        % Validate required parameters for grid
                                        if ~isfield(rect_grid, 'n_elem_row'), ...
                                                error('Transducer %i;: Missing n_elem_row parameter for grid. Please specify to define number of rows.', t_i);
                                        end
    
                                        if ~isfield(rect_grid, 'n_elem_col'), ...
                                                error('Transducer %i; Missing n_elem_col parameter for grid. Please specify to define number of columns.', t_i);
                                        end
    
                                        if ~isfield(rect_grid, 'elem_spacing_height_mm'), ...
                                                error('Transducer %i; Missing elem_spacing_height parameter for grid. Please specify to define height spacing between elements.', t_i);
                                        end
    
                                        if ~isfield(rect_grid, 'elem_spacing_width_mm'), ...
                                                error('Transducer %i; Missing elem_spacing_width parameter for grid. Please specify to define width spacing between elements.', t_i);
                                        end
    
                                        if ~isfield(rect_grid, 'sparsity_factor'), ...
                                                error('Transducer %i; Missing sparsity_factor parameter for grid. Please specify to define % of used elements.', t_i);
                                        end
    
                                        % Extract grid dimensions
                                        n_elem_row = rect_grid.n_elem_row;
                                        n_elem_col = rect_grid.n_elem_col;
    
                                        % Store element dimensions in standard location
                                        tr.n_elem_row = n_elem_row;
                                        tr.n_elem_col = n_elem_col;
    
                                        % Calculate initial element count (will be adjusted later for circular cutout)
                                        tr.n_elements = n_elem_col * n_elem_row;
    
                                        tr.elem_spacing_height_mm = rect_grid.elem_spacing_height_mm;
                                        tr.elem_spacing_width_mm = rect_grid.elem_spacing_width_mm;
    
                                        % Get element count directly from configuration
                                        tr.sparsity_factor = rect_grid.sparsity_factor;
    
                                    case 'fibonacci'
                                        % Sparse spiral grid configuration
    
                                        % Validate required parameters for sparser grid
                                        if ~isfield(define_shape_here.grid_shape, 'fibonacci')
                                            error('Transducer %i; Missing fibonacci field in grid_shape for sparser grid configuration.', t_i);
                                        end
                                        
                                        fibo_grid = define_shape_here.grid_shape.fibonacci;
    
                                        if ~isfield(fibo_grid, 'n_elements')
                                            error('Transducer %i; Missing sparsity_factor parameter for grid. This defines % of used elements.', t_i);
                                        end
    
                                        tr.n_elements = fibo_grid.n_elements;
                                    otherwise
                                        error('Transducer %i; Grid shape type "%s" is not implemented.', ...
                                            t_i, define_shape_here.grid_shape.type);
                                end
    
                            case 'extract_from_file'
                                if ~isfield(matrix_tr.matrix_shape, 'extract_from_file')
                                    error('Transducer %i; Missing extract_from_file field and properties for matrix transducer. Please specify.', t_i);
                                end

                                extract_shape_from_file = matrix_tr.matrix_shape.extract_from_file;

                                if ~isfield(extract_shape_from_file, 'file_path')
                                    error('Transducer %i; Missing file_path to extract element positions from. Please specify.', t_i);
                                end
    
                                file_path = extract_shape_from_file.file_path;
    
                                if ~isfile(file_path)
                                    error('Transducer %i; File "%s" does not exist or is not readable.', t_i, file_path);
                                end
    
                                tr.file_path = file_path;
    
                                if ~isfield(extract_shape_from_file, 'start_row')
                                    error('Transducer %i; Missing start_row parameter. Please specify start row of element positions.', t_i);
                                end
    
                                if ~isfield(extract_shape_from_file, 'start_col')
                                    error('Transducer %i; Missing start_col parameter. Please specify start column of element positions.', t_i);
                                end
    
                                tr.start_row = extract_shape_from_file.start_row;
                                tr.start_col = extract_shape_from_file.start_col;
    
                                if ~isfield(extract_shape_from_file, 'n_elements')
                                    error('Transducer %i; Missing n_elements. Please specify number of elements.', t_i);
                                end
    
                                tr.n_elements = extract_shape_from_file.n_elements;
    
                                if isfield(extract_shape_from_file, 'select_random_subset')
    
                                    tr.select_random_subset = extract_shape_from_file.select_random_subset;
                                    if tr.select_random_subset
                                        if ~isfield(extract_shape_from_file, 'subset')
                                            error('Transducer %i; Missing subset parameters. Please specify subset properties.', t_i);
                                        end
    
                                        if ~isfield(extract_shape_from_file.subset, 'random_seed')
                                            error('Transducer %i; Missing random_seed parameter. Please specify if a random subseed for subset selection is desired.', t_i);
                                        end
    
                                        if ~isfield(extract_shape_from_file.subset, 'subset_n_elements')
                                            error('Transducer %i; Missing subset_n_element parameter. Please specify subset size.', t_i);
                                        end
    
                                        tr.random_seed = extract_shape_from_file.subset.random_seed;
                                        tr.subset_n_elements = extract_shape_from_file.subset.subset_n_elements;
                                    end
    
                                    % Parameter isn't found, so don't create random subset
                                else
                                    tr.select_random_subset = false;
                                end
    
                                if isfield(extract_shape_from_file, 'project_on_new_ROC')
                                    tr.project_on_new_ROC = extract_shape_from_file.project_on_new_ROC;
    
                                    if tr.project_on_new_ROC
                                        if ~isfield(extract_shape_from_file.project_on_new_ROC, 'new_ROC_mm')
                                            error('Transducer %i; Missing new_ROC_mm parameter. Please specify new desired ROC.', t_i);
                                        end
    
                                        tr.new_ROC_mm = extract_shape_from_file.new_ROC_mm;
    
                                    end
    
                                    % Parameter isn't found, so don't create projection
                                else
                                    tr.project_on_new_ROC = false;
                                end
    
                            otherwise
                                error('Transducer %i; Matrix shape option "%s" is not implemented.', ...
                                    t_i, matrix_tr.matrix_shape.type);
                        end
    
                        % Define annular element dimensions for visualization
                        % compatibility
                        end_elem = matrix_tr.outer_diameter_mm;
    
                        % Calculate inner and outer diameters for elements
                        [id, od] = calc_elements_id_od_mm(end_elem, tr.n_elements);
    
                        % Store element dimensions in standard location for visualization compatibility
                        tr.Elements_ID_mm = id;
                        tr.Elements_OD_mm = od;
    
                    case 'annular'
                        % Handle annular transducer configuration
                        if ~isfield(tr.array_shape, 'annular')
                            error('Transducer %i; Appropriate configuration for annular transducer is missing.', t_i);
                        end

                        annular_tr = tr.array_shape.annular;
    
                        if ~isfield(annular_tr, 'Elements_ID_mm')
                            error('Transducer %i; Missing Elements_ID_mm for annular transducer. This defines inner diameters of elements.', t_i);
                        end
    
                        if ~isfield(annular_tr, 'Elements_OD_mm')
                            error('Transducer %i; Missing Elements_OD_mm for annular transducer. This defines outer diameters of elements.', t_i);
                        end
    
                        if ~isfield(annular_tr, 'n_elements')
                            error('Transducer %i; Missing n_elements parameter for annular transducer. This defines element count.', t_i);
                        end
    
                        if ~isfield(annular_tr, 'curv_radius_mm')
                            error('Transducer %i; Missing curv_radius_mm field for annular transducer. Please specify radius of curvature.', t_i);
                        end
    
                        if isfield(annular_tr, 'dist_to_plane_mm')
                            tr.dist_to_plane_mm = annular_tr.dist_to_plane_mm;                    
                        end
    
                        if ~isfield(annular_tr, 'source_phase_deg')
                            error('Transducer %i;: Missing source_phase_deg field for annular transducer. Please specify phases.', t_i);
                        end
    
                        % Copy element count to standard location for compatibility
                        tr.n_elements = annular_tr.n_elements;
    
                        tr.Elements_ID_mm = annular_tr.Elements_ID_mm;
                        tr.Elements_OD_mm = annular_tr.Elements_OD_mm;
    
                        tr.curv_radius_mm = annular_tr.curv_radius_mm;
    
                        tr.source_phase_deg = annular_tr.source_phase_deg;
    
                    otherwise
                        error('Transducer %i; Element shape option "%s" is not implemented.', ...
                            tr.array_shape.type);
                end
            else
                tr.array_shape = 'annular';
            end
    
            if strcmp(tr.array_shape, 'annular')
    
                % Ensure source phase is set in radians or degrees
                if ~isfield(tr, 'source_phase_rad')
                    if ~isfield(tr, 'source_phase_deg'), ...
                            error('Transducer %i; Source phase should be set in transducer parameters as source_phase_rad or source_phase_deg', t_i);
                    end
                    tr.source_phase_rad = deg2rad(tr.source_phase_deg);
                end
    
                % Calculate distance to transducer plane if not provided
                if ~isfield(tr, 'dist_to_plane_mm')
                    tr.dist_to_plane_mm = sqrt(tr.curv_radius_mm^2 - ...
                        (max(tr.Elements_OD_mm) / 2)^2);
                    fprintf('Transducer %i; Distance to transducer plane is not provided, calculated as %.2f mm\n', ...
                        t_i, tr.dist_to_plane_mm);
                end
    
                % Ensure source amplitude matches number of transducer elements
                if length(tr.source_amp) == 1 && tr.n_elements > 1
                    tr.source_amp = repmat(tr.source_amp, [1, tr.n_elements]);
                end
    
                % Evaluate source phase expressions if stored as cell arrays
                if iscell(tr.source_phase_rad)
                    for i = 1:length(tr.source_phase_rad)
                        if ~isnumeric(tr.source_phase_rad{i})
                            tr.source_phase_rad{i} = str2double(tr.source_phase_rad{i});
                        end
                    end
                    tr.source_phase_rad = cell2mat(tr.source_phase_rad);
                end
    
                % Ensure source phase degrees are calculated from radians if not provided
                if ~isfield(tr, 'source_phase_deg')
                    tr.source_phase_deg = rad2deg(tr.source_phase_rad);
                end
            end
        end
    else
        % Warn user about missing transducer information
        assert(all(confirmation_dlg('The transducer info is missing in the configuration file. Do you want to continue?', 'Yes', 'No')), ...
            'Exiting');
    end

    parameters.transducer(t_i) = tr;
end