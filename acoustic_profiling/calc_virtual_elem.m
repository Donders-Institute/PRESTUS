clear;

start_elem = 0; % Inner diameter of first original element
end_elem = 45.4914; % Outer diameter of last original element
n_v_elem = 10; % Number of desired virtual elements
perc_elem = 90; % Percentage of gap used for element itself, rest is space between elements


dist = end_elem - start_elem; % total distance between original inner and outer element
elem_dist = (dist/n_v_elem)*(perc_elem/100); % homogeneous distance per virtual element
space_dist = (dist - elem_dist*n_v_elem)/(n_v_elem-1); % homogenous distance inbetween elements

id = zeros(1, n_v_elem);
od = zeros(1, n_v_elem);

id(1) = start_elem; % Inner diameter of each virtual element
od(1) = start_elem+elem_dist; % Outer diameter of each virtual element 
for i = 1:n_v_elem-1
    id(i+1) = od(i)+space_dist;
    od(i+1) = id(i+1)+elem_dist;
end

fprintf('Inner diameter of each virtual element: %s \n', mat2str(id))
fprintf('Outer diameter of each virtual element: %s \n', mat2str(od))