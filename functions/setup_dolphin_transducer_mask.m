function [transducer_mask, source_label, transducer_pars] = setup_dolphin_transducer_mask(grid_dims, transducer_pars, trans_pos, grid_step_mm, totalMask)
% This function generates a binary mask (`transducer_mask`) and a labeled matrix 
% (`source_label`) for the Dolphin transducer. It computes these 
% based on the transducer's geometric parameters and its position relative to 
% a computational grid. Additionally, it updates the `transducer_pars` structure 
% with converted values (e.g., diameters in grid points).
% Input:
%   grid_dims         - [Nx, Ny, Nz] array defining the dimensions of the computational grid.
%   transducer_pars   - Struct containing geometric parameters of the transducer:
%                       * n_elements: Number of elements in the transducer.
%                       * Elements_OD_mm: Outer diameter of each element (in mm).
%                       * Elements_ID_mm: Inner diameter of each element (in mm).
%                       * curv_radius_mm: Curvature radius of the elements (in mm).
%   trans_pos         - [1x3] or [3x1] array specifying the Cartesian coordinates 
%                       of the transducer position.
%   grid_step_mm      - Scalar specifying the step size of the computational grid (in mm).
%   totalMask         - binary mask representation of Dolphin transducer.
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

    % Convert physical dimensions to grid points
    grid_step_m = grid_step_mm / 1000;
    elem_width_grid = round(transducer_pars.elem_width / grid_step_m);  
    elem_height_grid = round(transducer_pars.elem_height / grid_step_m);
    elem_spacing_width_grid = round(transducer_pars.elem_spacing_width / grid_step_m);
    elem_spacing_height_grid = round(transducer_pars.elem_spacing_height / grid_step_m);  

    totalMaskCols = size(totalMask, 2);  % Get the number of columns in the new mask 
    totalMaskRows = size(totalMask, 1);

    % Calculate the updated center X position for the totalMask
    totalMaskCenterX = (totalMaskCols + 1) / 2;
    totalMaskCenterY = (totalMaskRows + 1) / 2;

    % Element counter for labeling
    element_counter = 0;
    
    % Initialize computational grids for the transducer mask and source label matrix
    transducer_mask = zeros(grid_dims); % Binary mask representing active transducer regions
    source_label = zeros(grid_dims);   % Label matrix for identifying individual elements
    
    % Process each element in the grid using the total mask dimensions
    for row = 1:transducer_pars.n_elem_row
        for col = 1:totalMaskCols
            % Check if this element is active according to the mask
            if totalMask(row, col)
                element_counter = element_counter + 1;
                
                % Calculate grid position for this element
                grid_x = trans_pos(1) + (col - totalMaskCenterX) * (elem_width_grid + elem_spacing_width_grid);
                grid_y = trans_pos(2) + (row - totalMaskCenterY) * (elem_height_grid + elem_spacing_height_grid);
                grid_z = trans_pos(3);
                
                % Convert to integer grid positions
                x_start = round(grid_x - elem_width_grid/2);
                y_start = round(grid_y - elem_height_grid/2);
                z_start = round(grid_z);
                
                % Create element using rectangular shape
                for x = x_start:x_start+elem_width_grid
                    for y = y_start:y_start+elem_height_grid
                        % Check if position is within grid bounds
                        if x >= 1 && x <= grid_dims(1) && y >= 1 && y <= grid_dims(2) && z_start >= 1 && z_start <= grid_dims(3)
                            transducer_mask(x, y, z_start) = 1;
                            source_label(x, y, z_start) = element_counter;
                        end
                    end
                end
            end
        end
    end
    
    % Update transducer_pars with number of active elements
    transducer_pars.n_elements = element_counter;
end