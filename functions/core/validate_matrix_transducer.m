function [parameters, tr] = validate_matrix_transducer(parameters, tr, t_i)
% VALIDATE_MATRIX_TRANSDUCER Validates configuration of a matrix transducer.
%
% This function checks that all required fields for a matrix transducer
% definition are present in the input transducer structure. 

% The function also supports optional configurations such as Clover array
% layouts, sparse Fibonacci grids, random element subsets, and projection
% of element positions onto a new radius of curvature.
%
% After validation, the relevant parameters are copied to standardized
% fields used throughout the simulation pipeline. Additional derived values
% such as the number of elements and equivalent annular element diameters
% (for visualization compatibility) are calculated.
%
% Input:
%   parameters - Struct containing global simulation parameters.
%   tr         - Struct containing a single transducer definition.
%   t_i        - Index of the transducer in the configuration (used for
%                informative error messages).
%
% Output:
%   parameters - Updated parameters struct (e.g., enabling kWaveArray).
%   tr         - Transducer struct with validated and normalized fields,
%                including geometry, grid configuration, and derived values.

    matrix_tr = tr.array_shape.matrix;

    fprintf('Matrix transducer detected. Using kWaveArray (set to 1).\n');
    parameters.use_kWaveArray = 1;

    % ---------------------------------------------------------------------
    % Steering configuration
    % Determines whether the array orientation follows the acoustic focus.
    % ---------------------------------------------------------------------
    assert(isfield(matrix_tr, 'steering'), ...
        'Transducer %i; Missing steering definition. Choose "1D" or "3D".', t_i);

    switch matrix_tr.steering
        case '1D'
            matrix_tr.align_transducer_with_focus = true;
        case '3D'
            matrix_tr.align_transducer_with_focus = false;
        otherwise
            error('Transducer %i; Steering option "%s" is not implemented.', ...
                t_i, matrix_tr.steering);
    end

    % Initialize element phases to zero degrees by default.
    % Phase delays are later adjusted based on the defined focus.
    matrix_tr.source_phase_deg = 0;

    % ---------------------------------------------------------------------
    % Element geometry
    % Defines the physical dimensions of individual matrix elements.
    % ---------------------------------------------------------------------
    assert(isfield(matrix_tr, 'element_shape'), ...
        'Transducer %i; Missing element_shape field for matrix array. Please specify "rect", "disc" or "bowl".', t_i);

    assert(isfield(matrix_tr, 'elem_height_mm'), ...
        'Transducer %i; Missing elem_height parameter. Please specify to define height of each element.', t_i);

    assert(isfield(matrix_tr, 'elem_width_mm'), ...
        'Transducer %i; Missing elem_width parameter. Please specify to define width of each element.', t_i);

    % Validate outer diameter
    assert(isfield(matrix_tr, 'outer_diameter_mm'), ...
        'Transducer %i; Missing outer_diameter_mm parameter for matrix array. Please specify to define overall size.', t_i);

    assert(matrix_tr.outer_diameter_mm > 0)

    % ---------------------------------------------------------------------
    % Curvature definition
    % Determines whether the array lies on a spherical surface.
    % ---------------------------------------------------------------------
    assert(isfield(matrix_tr, 'is_curved'), ...
        'Transducer %i; Missing is_curved field for matrix transducer. Please specify.', t_i);

    if matrix_tr.is_curved
        assert(isfield(matrix_tr.curved, 'curv_radius_mm'), ...
            'Transducer %i; Missing curv_radius_mm field for matrix transducer. Please specify radius of curvature.', t_i);

        % Calculate distance to transducer plane if not provided
        if ~isfield(matrix_tr.curved, 'dist_to_plane_mm')
            assert(matrix_tr.curved.curv_radius_mm > matrix_tr.outer_diameter_mm/2, ...
                'Transducer %i; curv_radius_mm must exceed aperture radius.', t_i);

            matrix_tr.curved.dist_to_plane_mm = sqrt(matrix_tr.curved.curv_radius_mm^2 - ...
                (matrix_tr.outer_diameter_mm / 2)^2);

            fprintf('Transducer %i; Distance to transducer plane is not provided, calculated as %.2f mm\n', ...
                t_i, matrix_tr.curved.dist_to_plane_mm);
        end

    else
        matrix_tr.curved.curv_radius_mm = inf;

        % For a flat transducer the distance to the focal plane approaches
        % infinity. A finite value is assigned here for visualization purposes.
        matrix_tr.curved.dist_to_plane_mm = 70;
    end

    % ---------------------------------------------------------------------
    % Optional Clover multi-aperture configuration
    % ---------------------------------------------------------------------
    if isfield(matrix_tr, 'is_clover_setup')

        if matrix_tr.is_clover_setup
            assert(isfield(matrix_tr, 'clover'), ...
                'Transducer %i; Missing additional information of Clover setup.', t_i);

            assert(isfield(matrix_tr.clover, 'n_leaves'), ...
                'Transducer %i; Missing additional information of Clover setup. Define number of leaves.', t_i);

            assert(isfield(matrix_tr.clover, 'ROC_parent'), ...
                'Transducer %i; Missing additional information of Clover setup. Define ROC of parent.', t_i);
        end
    else
        % Clover array configuration not provided; skip Clover setup
        matrix_tr.is_clover_setup = false;
    end

    % ---------------------------------------------------------------------
    % Matrix grid definition
    % Elements can be defined directly or extracted from an external file.
    % ---------------------------------------------------------------------
    assert(isfield(matrix_tr.matrix_shape, 'type'), ...
        'Transducer %i; Missing matrix_shape type field for matrix transducer. Please specify "define_here" or "extract_from_file"', t_i);

    switch matrix_tr.matrix_shape.type
        case 'define_here'
            assert(isfield(matrix_tr.matrix_shape, 'define_here'), ...
                'Transducer %i; Missing define_here field and properties for matrix transducer. Please specify.', t_i);

            define_shape_here = matrix_tr.matrix_shape.define_here;

            assert(isfield(define_shape_here, 'grid_shape'), ...
                'Transducer %i; Missing grid_shape field and properties for matrix transducer. Please specify.', t_i);

            grid_shape = define_shape_here.grid_shape;

            assert(isfield(grid_shape, 'type'), ...
                'Transducer %i; Missing grid_shape type field and properties for matrix transducer. Please specify.', t_i);

            switch grid_shape.type
                case 'rect'
                    assert(isfield(grid_shape, 'rect'), ...
                        'Transducer %i; Missing rect field in grid_shape for rectangular grid configuration.', t_i);

                    rect_grid = grid_shape.rect;

                    % Validate required parameters for grid
                    assert(isfield(rect_grid, 'n_elem_row'), ...
                            'Transducer %i;: Missing n_elem_row parameter for grid. Please specify to define number of rows.', t_i);

                    assert(isfield(rect_grid, 'n_elem_col'), ...
                            'Transducer %i; Missing n_elem_col parameter for grid. Please specify to define number of columns.', t_i);

                    assert(isfield(rect_grid, 'elem_spacing_height_mm'), ...
                            'Transducer %i; Missing elem_spacing_height parameter for grid. Please specify to define height spacing between elements.', t_i);

                    assert(isfield(rect_grid, 'elem_spacing_width_mm'), ...
                            'Transducer %i; Missing elem_spacing_width parameter for grid. Please specify to define width spacing between elements.', t_i);

                    % Validate spacing positivity
                    assert(rect_grid.elem_spacing_height_mm > 0)
                    assert(rect_grid.elem_spacing_width_mm > 0)

                    % Total physical span of the array
                    tran_width = (rect_grid.n_elem_width * matrix_tr.elem_width_mm) + ...
                        ((rect_grid.n_elem_width - 1) * rect_grid.elem_spacing_width_mm);

                    tran_height = (rect_grid.n_elem_height * matrix_tr.elem_height_mm) + ...
                        ((rect_grid.n_elem_height - 1) * rect_grid.elem_spacing_height_mm);

                    % Distance of the array side from the center
                    grid_radius = max(tran_width/2, tran_height/2);

                    aperture_radius = matrix_tr.outer_diameter_mm / 2;

                    assert(grid_radius <= aperture_radius, ...
                        ['Transducer %i; Rectangular grid exceeds outer_diameter_mm. ' ...
                        'Grid span: %.2f x %.2f mm, aperture diameter: %.2f mm.'], ...
                        t_i, tran_width, tran_height, matrix_tr.outer_diameter_mm);

                    assert(isfield(rect_grid, 'sparsity_factor'), ...
                            'Transducer %i; Missing sparsity_factor parameter for grid. Please specify to define % of used elements.', t_i);

                    assert(rect_grid.sparsity_factor > 0 && rect_grid.sparsity_factor <= 1)

                    % Extract grid dimensions
                    n_elem_row = rect_grid.n_elem_row;
                    n_elem_col = rect_grid.n_elem_col;

                    % Calculate initial element count (will be adjusted later for circular cutout)
                    matrix_tr.n_elements = n_elem_col * n_elem_row;

                case 'fibonacci'
                    % Sparse spiral grid configuration

                    % Validate required parameters for sparser grid
                    assert(isfield(grid_shape, 'fibonacci'), ...
                        'Transducer %i; Missing fibonacci field in grid_shape for sparser grid configuration.', t_i);
                 
                    assert(isfield(grid_shape.fibonacci, 'n_elements'), ...
                       'Transducer %i; Missing n_elements parameter for grid. Please specify to define number of elements.', t_i);

                    assert(isfield(grid_shape.fibonacci, 'kerf_mm'), ...
                       'Transducer %i; Missing kerf_mm parameter for grid. Please specify.', t_i);

                    matrix_tr.n_elements = grid_shape.fibonacci.n_elements;

                case 'fermat'
                    % Sparse spiral grid configuration

                    % Validate required parameters for sparser grid
                    assert(isfield(grid_shape, 'fermat'), ...
                        'Transducer %i; Missing fermat field in grid_shape for sparser grid configuration.', t_i);

                    assert(isfield(grid_shape.fermat, 'n_elements'), ...
                        'Transducer %i; Missing n_elements parameter for grid. Please specify to define number of elements.', t_i);
                    
                    matrix_tr.n_elements = grid_shape.fermat.n_elements;

                otherwise
                    error('Transducer %i; Grid shape type "%s" is not implemented.', ...
                        t_i, grid_shape.type);
            end

        case 'extract_from_file'
            assert(isfield(matrix_tr.matrix_shape, 'extract_from_file'), ...
                'Transducer %i; Missing extract_from_file field and properties for matrix transducer. Please specify.', t_i);

            extract_shape_from_file = matrix_tr.matrix_shape.extract_from_file;

            assert(isfield(extract_shape_from_file, 'file_path'), ...
                'Transducer %i; Missing file_path to extract element positions from. Please specify.', t_i);

            file_path = extract_shape_from_file.file_path;

            resolved_path = which(file_path);

            assert(~isempty(resolved_path), ...
                'Transducer %i; File "%s" does not exist or is not on the MATLAB path.', t_i, file_path);

            assert(isfield(extract_shape_from_file, 'start_row'), ...
                'Transducer %i; Missing start_row parameter. Please specify start row of element positions.', t_i);

            assert(isfield(extract_shape_from_file, 'start_col'), ...
                'Transducer %i; Missing start_col parameter. Please specify start column of element positions.', t_i);

            assert(isfield(extract_shape_from_file, 'n_elements'), ...
                'Transducer %i; Missing n_elements. Please specify number of elements.', t_i);

            matrix_tr.n_elements = extract_shape_from_file.n_elements;

            if isfield(extract_shape_from_file, 'select_random_subset')

                if extract_shape_from_file.select_random_subset
                    assert(isfield(extract_shape_from_file, 'subset'), ...
                        'Transducer %i; Missing subset parameters. Please specify subset properties.', t_i);

                    assert(isfield(extract_shape_from_file.subset, 'random_seed'), ...
                        'Transducer %i; Missing random_seed parameter. Please specify if a random subseed for subset selection is desired.', t_i);

                    assert(isfield(extract_shape_from_file.subset, 'subset_n_elements'), ...
                        'Transducer %i; Missing subset_n_element parameter. Please specify subset size.', t_i);

                end
            else
                % Random element subset not specified; use full element set
                extract_shape_from_file.select_random_subset = false;
            end

            if isfield(extract_shape_from_file, 'project_on_new_ROC')

                if extract_shape_from_file.project_on_new_ROC
                    assert(isfield(extract_shape_from_file, 'ROC_projection'), ...
                        'Transducer %i; Missing ROC_projection parameters. Please specify to define ROC to project on.', t_i);

                    assert(isfield(extract_shape_from_file.ROC_projection, 'new_ROC_mm'), ...
                        'Transducer %i; Missing new_ROC_mm parameter. Please specify new desired ROC.', t_i);

                end
            else
                % No projection configuration provided; element positions remain unchanged
                extract_shape_from_file.project_on_new_ROC = false;
            end

        otherwise
            error('Transducer %i; Matrix shape option "%s" is not implemented.', ...
                t_i, matrix_tr.matrix_shape.type);
    end
    matrix_tr.matrix_shape.extract_from_file = extract_shape_from_file;

    % ---------------------------------------------------------------------
    % Derived visualization parameters
    % Convert matrix array into equivalent annular representation.
    % ---------------------------------------------------------------------
    transducer_diameter_mm = matrix_tr.outer_diameter_mm;

    % Calculate inner and outer diameters for elements
    [id, od] = calc_elements_id_od_mm(transducer_diameter_mm, matrix_tr.n_elements);

    % Store element dimensions in standard location for visualization compatibility
    matrix_tr.Elements_ID_mm = id;
    matrix_tr.Elements_OD_mm = od;

    tr.array_shape.matrix = matrix_tr;
end