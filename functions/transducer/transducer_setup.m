function [transducer_mask, source_label, transducer_pars] = ...
    transducer_setup(transducer_pars, trans_pos, focus_pos, grid_dims, grid_res_mm)

% TRANSDUCER_SETUP Creates a transducer mask and label matrix for a computational grid.
%
% This function generates a binary mask (`transducer_mask`) and a labeled matrix 
% (`source_label`) for a multi-element ultrasound transducer. It computes these 
% based on the transducer's geometric parameters and its position relative to 
% a computational grid. Additionally, it updates the `transducer_pars` structure 
% with converted values (e.g., diameters in grid points).
%
% Use this function to define the geometry of a transducer and map it to a 
% computational grid for numerical simulations.
%
% Input:
%   transducer_pars   - Struct containing geometric parameters of the transducer:
%                       * n_elements: Number of elements in the transducer.
%                       * Elements_OD_mm: Outer diameter of each element (in mm).
%                       * Elements_ID_mm: Inner diameter of each element (in mm).
%                       * curv_radius_mm: Curvature radius of the elements (in mm).
%
%   trans_pos         - [1x2/3] or [2/3x1] array specifying the Cartesian coordinates 
%                       of the transducer position.
%
%   focus_pos         - [1x2/3] or [2/3x1] array specifying the Cartesian coordinates 
%                       of the geometric focus position.
%
%   grid_dims         - [Nx, Ny, Nz] array defining the dimensions of the computational grid.
%
%   grid_res_mm       - Scalar specifying the step size of the computational grid (in mm).
%
% Output:
%   transducer_mask   - Binary matrix of size `grid_dims`. Non-zero values represent 
%                       active regions of specific elements in the transducer.
%
%   source_label      - Labeled matrix of size `grid_dims`. Each non-zero value uniquely 
%                       identifies an individual element in the transducer.
%
%   transducer_pars   - Updated struct with additional fields:
%                       * Elements_OD: Outer diameter in grid points.
%                       * Elements_ID: Inner diameter in grid points.
%                       * radius_grid: Curvature radius in grid points.

    % Check if the transducer positions and focus positions have matching dimensions
    if ~isequal(size(focus_pos), size(trans_pos))
        error('Transducer and focus positions should be arrays of equal size')
    end

    % Ensure the positions are formatted as row vectors with size [1, 3]
    if isequal(size(focus_pos), [3 1]) | isequal(size(focus_pos), [2 1])
        focus_pos = focus_pos'; % Convert column vector to row vector
        trans_pos = trans_pos'; % Convert column vector to row vector
    elseif ~isequal(size(focus_pos), [1 3]) & ~isequal(size(focus_pos), [1 2])
        error('Transducer and focus positions should have the size [1 2] or [1 3]')
    end

    % Initialize computational grids for the transducer mask and source label matrix
    transducer_mask = zeros(grid_dims); % Binary mask representing active transducer regions
    source_label = zeros(grid_dims);   % Label matrix for identifying individual elements
    switch transducer_pars.array_shape.type
        case 'annular'
            % Convert element diameters from millimeters to grid points and ensure they are odd integers
			transducer_pars.array_shape.annular.Elements_OD = ...
                2*floor(transducer_pars.array_shape.annular.Elements_OD_mm / grid_res_mm / 2) + 1; % Outer diameter in grid points
			transducer_pars.array_shape.annular.Elements_ID = ...
                2*floor(transducer_pars.array_shape.annular.Elements_ID_mm / grid_res_mm / 2) + 1; % Inner diameter in grid points

			% Handle cases where inner diameter is zero (e.g., for flat elements)
			transducer_pars.array_shape.annular.Elements_ID(transducer_pars.array_shape.annular.Elements_ID_mm == 0) = 0;

			% Convert the curvature radius from millimeters to grid points
			transducer_pars.array_shape.annular.radius_grid = round(transducer_pars.array_shape.annular.curv_radius_mm / grid_res_mm); % Radius in grid points

            % Loop through each transducer element to create its geometry
            for el_i = 1:transducer_pars.array_shape.annular.n_elements
                % Create the outer bowl geometry for the current element
                bowl = makeBowl(...
                    grid_dims, ...
                    trans_pos, ...
                    transducer_pars.array_shape.annular.radius_grid,...
                    transducer_pars.array_shape.annular.Elements_OD(el_i), ...
                    focus_pos);

                % If the inner diameter is greater than zero, subtract the inner bowl geometry
                if transducer_pars.array_shape.annular.Elements_ID(el_i) > 0
                    bowl = bowl - makeBowl(...
                        grid_dims, ...
                        trans_pos, ...
                        transducer_pars.array_shape.annular.radius_grid, ...
                        transducer_pars.array_shape.annular.Elements_ID(el_i), ...
                        focus_pos);
                end

                % Add the current element's bowl geometry to the binary mask
                transducer_mask = transducer_mask + bowl;

                % Assign a unique label to this element in the source label matrix
                source_label = source_label + el_i * bowl;
            end
        case 'matrix'
            matrix_tp = transducer_pars.array_shape.matrix;

            trans_pos_m = trans_pos * grid_step_mm / 1e3;
            focus_pos_m = focus_pos * grid_step_mm / 1e3;

            switch matrix_tp.matrix_shape.type
                case 'define_here'
                    [elem_pos_m, transducer_pars] = convert_to_element_pos(parameters, transducer_pars,  trans_pos_m, focus_pos_m);
                case 'extract_from_file'
                    elem_pos_m = extract_element_pos(parameters, transducer_pars,  trans_pos_m);
                otherwise
                    error('Matrix shape %s is unknown or not implemented.', matrix_tp.matrix_shape.type)
            end

            natural_focus_pos_m = trans_pos_m + [0, 0, matrix_tp.curved.curv_radius_mm / 1000]';

            % Apply Clover setup if requested
            if transducer_pars.is_clover_setup
                elem_pos_m = create_clover_array(elem_pos_m, transducer_pars, trans_pos_m);
            end

            transducer_pars.n_elements = size(elem_pos_m, 2);

            % Convert to grid units
            elem_pos_grid = round(elem_pos_m * 1e3 / grid_step_mm);

            natural_focus_pos_grid = round(natural_focus_pos_m * 1e3 / grid_step_mm);

            % Loop through each transducer element to create its geometry
            for el_i = 1:transducer_pars.n_elements
                el_pos_grid_i = elem_pos_grid(:, el_i);

                r_c = matrix_tp.curved.curv_radius_mm / 1e3;
                A_target = matrix_tp.elem_height_mm * matrix_tp.elem_width_mm;
                a_min = 0.001; % Lower bound of aperture radius (in m)
                a_max = r_c * 0.999; % Slightly less than full bowl radius

                % Define the function to solve: A_cap(a) - A_target = 0
                area_diff = @(a) 2*pi*r_c*(r_c - sqrt(r_c^2 - a^2)) - A_target;

                % Check if the function changes sign in the interval
                if sign(area_diff(a_min)) == sign(area_diff(a_max))
                    error('No sign change in interval: cannot find bowl diameter. Possibly A_target is out of bounds.');
                end

                a_solution = fzero(area_diff, [a_min, a_max]);
                diameter = 2 * a_solution;

                % Convert to grid dimensions
                r_c_grid = r_c * 1e3 / grid_step_mm;
                diameter_grid = diameter * 1e3 / grid_step_mm;

                bowl = makeBowl(grid_dims, el_pos_grid_i, r_c_grid, diameter_grid, natural_focus_pos_grid);

                % Add the current element's bowl geometry to the binary mask
                transducer_mask = transducer_mask + bowl;

                % Assign a unique label to this element in the source label matrix
                source_label = source_label + el_i * bowl;
            end
        otherwise
            error('Array type %s is unknown or not implemented.', transducer_pars.array_shape.type)
    end

end