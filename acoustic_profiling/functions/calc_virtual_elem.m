clear;

% Define input parameters
start_elem = 0;            % Inner diameter of the first original element (initial value)
end_elem = 45.4914;        % Outer diameter of the last original element (final value)
n_v_elem = 10;             % Number of desired virtual elements
perc_elem = 90;            % Percentage of the gap allocated for the element itself, the rest is space between elements

% Calculate total distance between the inner and outer diameters
dist = end_elem - start_elem;

% Calculate the homogeneous distance per virtual element, considering the element's percentage of the gap
elem_dist = (dist / n_v_elem) * (perc_elem / 100);

% Calculate the homogeneous space between the virtual elements
space_dist = (dist - elem_dist * n_v_elem) / (n_v_elem - 1);

% Pre-allocate arrays for storing inner and outer diameters of virtual elements
id = zeros(1, n_v_elem); % Inner diameters
od = zeros(1, n_v_elem); % Outer diameters

% Set the inner and outer diameters for the first virtual element
id(1) = start_elem;       % Inner diameter of the first virtual element
od(1) = start_elem + elem_dist; % Outer diameter of the first virtual element

% Loop to calculate inner and outer diameters for the remaining virtual elements
for i = 1:n_v_elem - 1
    id(i + 1) = od(i) + space_dist;        % Inner diameter of the next element is after the space
    od(i + 1) = id(i + 1) + elem_dist;    % Outer diameter is calculated by adding the element distance
end

% Print the results
fprintf('Inner diameter of each virtual element: [%s%s]\n', sprintf('%.3f, ', id(1:end-1)), sprintf('%.3f', id(end))); % Display the array of inner diameters
fprintf('Outer diameter of each virtual element: [%s%s]\n', sprintf('%.3f, ', od(1:end-1)), sprintf('%.3f', od(end))); % Display the array of outer diameters