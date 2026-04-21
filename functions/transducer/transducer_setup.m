function [transducer_mask, source_label, tr] = ...
    transducer_setup(tr, trans_pos, focus_pos, grid_dims, grid_res_mm)

% TRANSDUCER_SETUP  Create a transducer mask and label matrix for a k-Wave grid
%
% Generates a binary transducer mask and element-labelled source matrix
% for a multi-element ultrasound transducer, and updates the transducer
% struct with grid-unit equivalents of the geometric parameters.
%
% Use as:
%   [transducer_mask, source_label, tr] = ...
%       transducer_setup(tr, trans_pos, focus_pos, grid_dims, grid_res_mm)
%
% Input:
%   tr          - (1,1) transducer parameter struct
%   trans_pos   - [1x2/3] transducer position in grid coordinates
%   focus_pos   - [1x2/3] focus position in grid coordinates
%   grid_dims   - [1x2/3] computational grid dimensions [voxels]
%   grid_res_mm - grid step size [mm]
%
% Output:
%   transducer_mask - grid array with non-zero values at active element voxels
%   source_label    - grid array with unique integer label per element
%   tr              - updated transducer struct with grid-unit geometry fields
%
% See also: FOCAL_DISTANCE_CALCULATION, CONVERT_TO_ELEMENT_POS

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
    switch tr.type
        case 'annular'
            % Convert element diameters from millimeters to grid points and ensure they are odd integers
			tr.annular.Elements_OD = ...
                2*floor(tr.annular.elem_od_mm / grid_res_mm / 2) + 1; % Outer diameter in grid points
			tr.annular.Elements_ID = ...
                2*floor(tr.annular.elem_id_mm / grid_res_mm / 2) + 1; % Inner diameter in grid points

			% Handle cases where inner diameter is zero (e.g., for flat elements)
			tr.annular.Elements_ID(tr.annular.elem_id_mm == 0) = 0;

			% Convert the curvature radius from millimeters to grid points
			tr.annular.radius_grid = round(tr.annular.curv_radius_mm / grid_res_mm); % Radius in grid points

            % Loop through each transducer element to create its geometry
            for el_i = 1:tr.annular.elem_n
                % Create the outer bowl geometry for the current element
                bowl = makeBowl(...
                    grid_dims, ...
                    trans_pos, ...
                    tr.annular.radius_grid,...
                    tr.annular.Elements_OD(el_i), ...
                    focus_pos);

                % If the inner diameter is greater than zero, subtract the inner bowl geometry
                if tr.annular.Elements_ID(el_i) > 0
                    bowl = bowl - makeBowl(...
                        grid_dims, ...
                        trans_pos, ...
                        tr.annular.radius_grid, ...
                        tr.annular.Elements_ID(el_i), ...
                        focus_pos);
                end

                % Add the current element's bowl geometry to the binary mask
                transducer_mask = transducer_mask + bowl;

                % Assign a unique label to this element in the source label matrix
                source_label = source_label + el_i * bowl;
            end
        case 'matrix'
            trans_pos_m = trans_pos * grid_step_mm / 1e3;
            focus_pos_m = focus_pos * grid_step_mm / 1e3;

            switch tr.matrix.matrix_shape.type
                case 'define_here'
                    [elem_pos_m, tr] = convert_to_element_pos(parameters, tr,  trans_pos_m, focus_pos_m);
                case 'extract_from_file'
                    elem_pos_m = extract_element_pos(parameters, tr,  trans_pos_m);
                otherwise
                    error('Matrix shape %s is unknown or not implemented.', tr.matrix.matrix_shape.type)
            end

            natural_focus_pos_m = trans_pos_m + [0, 0, tr.matrix.curv_radius_mm / 1000]';

            % Apply Clover setup if requested
            if tr.is_clover_setup
                elem_pos_m = create_clover_array(elem_pos_m, tr, trans_pos_m);
            end

            tr.elem_n = size(elem_pos_m, 2);

            % Convert to grid units
            elem_pos_grid = round(elem_pos_m * 1e3 / grid_step_mm);

            natural_focus_pos_grid = round(natural_focus_pos_m * 1e3 / grid_step_mm);

            % Loop through each transducer element to create its geometry
            for el_i = 1:tr.elem_n
                el_pos_grid_i = elem_pos_grid(:, el_i);

                r_c = tr.matrix.curv_radius_mm / 1e3;
                A_target = tr.matrix.elem_height_mm * tr.matrix.elem_width_mm;
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
            error('Array type %s is unknown or not implemented.', tr.type)
    end

end