function [id, od] = calc_elements_id_od_mm(end_elem, n_elem)
%CALC_ELEMENTS_ID_OD_MM Calculates inner and outer diameters for transducer elements
%   Inputs:
%       END_ELEM - The outer diameter of the entire transducer in mm
%       N_ELEM - Number of elements in the transducer
%
%   Outputs:
%       ID - 1×N_ELEM array of inner diameters for each element (in mm)
%       OD - 1×N_ELEM array of outer diameters for each element (in mm)
%
    % Calculate total distance between the inner and outer diameters
    dist = end_elem - 0;
    
    % Calculate the homogeneous distance per virtual element, considering the element's percentage of the gap
    elem_dist = (dist / n_elem);
    
    % Calculate the homogeneous space between the virtual elements
    space_dist = (dist - elem_dist * n_elem) / (n_elem - 1);
    
    % Pre-allocate arrays for storing inner and outer diameters of virtual elements
    id = zeros(1, n_elem); % Inner diameters
    od = zeros(1, n_elem); % Outer diameters
    
    % Set the inner and outer diameters for the first virtual element
    id(1) = 0;       % Inner diameter of the first virtual element
    od(1) = 0 + elem_dist; % Outer diameter of the first virtual element
    
    % Loop to calculate inner and outer diameters for the remaining virtual elements
    for i = 1:n_elem - 1
        id(i + 1) = od(i) + space_dist;        % Inner diameter of the next element is after the space
        od(i + 1) = id(i + 1) + elem_dist;    % Outer diameter is calculated by adding the element distance
    end
end