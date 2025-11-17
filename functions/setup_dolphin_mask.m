function [totalMask] = setup_dolphin_mask(transducer_pars)
% This function creates the binary mask representing the Dolphin transducer
% design.
% Input:
%   transducer_pars   - Struct containing geometric parameters of the transducer:
%                       * n_elements: Number of elements in the transducer.
%                       * Elements_OD_mm: Outer diameter of each element (in mm).
%                       * Elements_ID_mm: Inner diameter of each element (in mm).
% Output:
%   totalMask         - binary mask representation of Dolphin transducer.

    % Define geometry
    % Create coordinate grids for the transducer elements (not the simulation grid)
    [X, Y] = meshgrid(1:transducer_pars.n_elem_col, 1:transducer_pars.n_elem_row);
    
    % Define circle parameters for the mask
    centerX = (transducer_pars.n_elem_col + 1) / 2;  % X center of the circle
    centerY = (transducer_pars.n_elem_row + 1) / 2;  % Y center of the circle

    % Create circular mask using the equation of a circle
    mask = (X - centerX).^2 + (Y - centerY).^2 <= (transducer_pars.n_elem_col/2+0.1)^2;
    mask(:,1) = [];

    % Create a symmetric mask by flipping and concatenating
    totalMask = [fliplr(mask), mask];
        
end