function [transducer_mask, source_label, transducer_pars] = ...
    transducer_setup(transducer_pars, trans_pos, focus_pos, grid_dims, grid_step_mm, perform_focus_rotation)

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
%   transducer_pars         - Struct containing geometric parameters of the transducer:
%                            * n_elements: Number of elements in the transducer.
%                            * Elements_OD_mm: Outer diameter of each element (in mm).
%                            * Elements_ID_mm: Inner diameter of each element (in mm).
%                            * curv_radius_mm: Curvature radius of the elements (in mm).
%
%   trans_pos                - [1x2/3] or [2/3x1] array specifying the Cartesian coordinates 
%                              of the transducer position.
%
%   focus_pos                - [1x2/3] or [2/3x1] array specifying the Cartesian coordinates 
%                              of the geometric focus position.
%
%   grid_dims                - [Nx, Ny, Nz] array defining the dimensions of the computational grid.
%
%   grid_step_mm             - Scalar specifying the step size of the computational grid (in mm).
%
%   perform_focus_rotation  - Boolean flag:
%                             true  -> rotate transducer position toward the focus (2D steering)
%                             false -> keep transducer position fixed and adjust phases to reach the focus (3D steering)
% 
% Output:
%   transducer_mask          - Binary matrix of size `grid_dims`. Non-zero values represent 
%                              active regions of specific elements in the transducer.
%
%   source_label             - Labeled matrix of size `grid_dims`. Each non-zero value uniquely 
%                              identifies an individual element in the transducer.
%
%   transducer_pars          - Updated struct with additional fields:
%                              * Elements_OD: Outer diameter in grid points.
%                              * Elements_ID: Inner diameter in grid points.
%                              * radius_grid: Curvature radius in grid points.

    % Check if perform_focus_rotation was given, otherwise assume default
    % behavior (2D steering)
    if nargin < 6 || isempty(perform_focus_rotation)
        perform_focus_rotation = true;
    end

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

    if strcmp(transducer_pars.element_shape, 'annular')
        % Convert element diameters from millimeters to grid points and ensure they are odd integers
        transducer_pars.Elements_OD = 2*floor(transducer_pars.Elements_OD_mm / grid_step_mm / 2) + 1; % Outer diameter in grid points
        transducer_pars.Elements_ID = 2*floor(transducer_pars.Elements_ID_mm / grid_step_mm / 2) + 1; % Inner diameter in grid points
    
        % Handle cases where inner diameter is zero (e.g., for flat elements)
        transducer_pars.Elements_ID(transducer_pars.Elements_ID_mm == 0) = 0;
    
        % Convert the curvature radius from millimeters to grid points
        transducer_pars.radius_grid = round(transducer_pars.curv_radius_mm / grid_step_mm); % Radius in grid points
        
        % Visualiation for flat transducer is not implemented. When
        % transducer is flat handle it as if the radius is very small.
        % MakeBowl does not accept a radius of zero.
        if transducer_pars.radius_grid ==  0
            transducer_pars.radius_grid = 1;
        end

        % Initialize computational grids for the transducer mask and source label matrix
        transducer_mask = zeros(grid_dims); % Binary mask representing active transducer regions
        source_label = zeros(grid_dims);   % Label matrix for identifying individual elements
    
        % Loop through each transducer element to create its geometry
        for el_i = 1:transducer_pars.n_elements
            if ~perform_focus_rotation
                focus_pos = [trans_pos(1), trans_pos(2), trans_pos(3) + transducer_pars.curv_radius_mm / grid_step_mm];
            end
            % Create the outer bowl geometry for the current element
            bowl = makeBowl(grid_dims, trans_pos, transducer_pars.radius_grid, transducer_pars.Elements_OD(el_i), focus_pos);
            
            % If the inner diameter is greater than zero, subtract the inner bowl geometry
            if transducer_pars.Elements_ID(el_i) > 0
                bowl = bowl - makeBowl(grid_dims, trans_pos, transducer_pars.radius_grid, transducer_pars.Elements_ID(el_i), focus_pos);
            end
    
            % Add the current element's bowl geometry to the binary mask
            transducer_mask = transducer_mask + bowl;
    
            % Assign a unique label to this element in the source label matrix
            source_label = source_label + el_i * bowl;
        end
    else
        [totalMask] = setup_dolphin_mask(transducer_pars);
        [transducer_mask, source_label, transducer_pars] = setup_dolphin_transducer_mask(grid_dims, transducer_pars, trans_pos, grid_step_mm, totalMask);
    end

end