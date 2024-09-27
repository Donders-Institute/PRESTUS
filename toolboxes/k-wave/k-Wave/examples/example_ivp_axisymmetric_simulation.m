% Simulations In An Axisymmetric Coordinate System Example
%
% This example provides a simple demonstration of using k-Wave for the
% simulation and detection of the pressure field generated by an initial
% pressure distribution within an axisymmetric coordinate system. It builds
% on the Homogeneous Propagation Medium and Heterogeneous Propagation
% Medium examples.
%
% author: Bradley Treeby
% date: 24th April 2018
% last update: 25th April 2018
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

clearvars;

% =========================================================================
% SIMULATION
% =========================================================================

% create the computational grid
Nx = 128;           % number of grid points in the axial (x) direction
Ny = 64;            % number of grid points in the radial (y) direction
dx = 0.1e-3;        % grid point spacing in the axial (x) direction [m]
dy = 0.1e-3;        % grid point spacing in the radial (y) direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1500 * ones(Nx, Ny);       % [m/s]
medium.sound_speed(Nx/2:end, :) = 1800;         % [m/s]
medium.density = 1000 * ones(Nx, Ny);           % [kg/m^3]
medium.density(Nx/2:end, :) = 1200;             % [kg/m^3]

% create initial pressure distribution in the shape of a disc - this is
% generated on a 2D grid that is doubled in size in the radial (y)
% direction, and then trimmed so that only half the disc is retained
source.p0 = 10 * makeDisc(Nx, 2 * Ny, Nx/4 + 8, Ny + 1, 5);
source.p0 = source.p0(:, Ny + 1:end);

% define a Cartesian sensor mask with points in the shape of a circle
sensor.mask = makeCartCircle(40 * dx, 50);

% remove points from sensor mask where y < 0
sensor.mask(:, sensor.mask(2, :) < 0) = [];

% run the simulation using the axisymmetric code
sensor_data = kspaceFirstOrderAS(kgrid, medium, source, sensor);

% =========================================================================
% VISUALISATION
% =========================================================================

% create plot axis
x_vec = 1e3 * kgrid.x_vec;
y_vec = 1e3 * (kgrid.y_vec - kgrid.y_vec(1));

% create the simulation layout, removing the PML
pml_size = 20;
sim_layout = double(source.p0 | cart2grid(kgrid, sensor.mask, true));
sim_layout = sim_layout(1 + pml_size:end - pml_size, 1:end - pml_size); 

% plot the simulation layout
figure;
imagesc(y_vec, x_vec, sim_layout, [0, 1]);
axis image;
colormap(flipud(gray));
xlabel('y (radial) position [mm]');
ylabel('x (axial) position [mm]');

% plot the simulated sensor data
figure;
imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;

% plot a trace from the simulated sensor data
figure;
plot(sensor_data(5, :));
xlabel('Time Step');
ylabel('Pressure [Pa]');