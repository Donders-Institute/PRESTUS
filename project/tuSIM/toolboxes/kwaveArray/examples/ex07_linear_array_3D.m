% Script to show how to setup a linear array transducer using the
% kWaveArray class.
%
% NOTE: setting the position does not correctly rotate the angle of the
% elements if an element angle is set. As a work around, set the element
% angle in addRectElement to the same angle as the array rotation.
%
% author: Bradley Treeby
% date: 1st November 2020
% last update: 1st November 2020
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2020 Bradley Treeby

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
% DEFINE LITERALS
% =========================================================================
    
% select which k-Wave code to run
%   1: MATLAB CPU code
%   2: MATLAB GPU code
%   3: C++ code
%   4: CUDA code
model           = 1;

% medium parameters
c0              = 1500;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]

% source parameters
source_f0       = 1e6;      % source frequency [Hz]
source_amp      = 1e6;      % source pressure [Pa]
source_cycles   = 5;        % number of toneburst cycles
source_focus    = 20e-3;    % focal length [m]
element_num     = 15;       % number of elements
element_width   = 1e-3;     % width [m]
element_length  = 10e-3;    % elevation height [m]
element_pitch   = 2e-3;     % pitch [m]

% source movement
translation     = [5e-3, 0, 8e-3];
rotation        = [0, 20, 0];

% grid parameters
grid_axial      = 40e-3;    % total grid size in the z dimension [m]
grid_width      = 40e-3;    % total grid size in the x dimension [m]
grid_height     = 20e-3;    % total grid size in the y dimension [m]

% computational parameters
ppw             = 3;        % number of points per wavelength
t_end           = 30e-6;    % total compute time [s]
cfl             = 0.5;      % CFL number

% =========================================================================
% RUN SIMULATION
% =========================================================================

% --------------------
% GRID
% --------------------

% calculate the grid spacing based on the PPW and F0
dx = c0 / (ppw * source_f0);   % [m]

% compute the size of the grid
Nx = roundEven(grid_width / dx);
Ny = roundEven(grid_height / dx);
Nz = roundEven(grid_axial / dx);

% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dx, Nz, dx);

% compute points per temporal period
PPP = round(ppw / cfl);

% compute corresponding time spacing
dt = 1 / (PPP * source_f0);

% create the time array using an integer number of points per period
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
disp(['PPW = ' num2str(c0 / (dx * source_f0))]);
disp(['CFL = ' num2str(c0 * dt / dx)]);

% --------------------
% SOURCE
% --------------------

% set indices for each element
if rem(element_num, 2)
    ids = (1:element_num) - ceil(element_num/2);
else
    ids = (1:element_num) - (element_num + 1)/2;
end

% set time delays for each element to focus at source_focus
time_delays = -(sqrt((ids .* element_pitch).^2 + source_focus.^2) - source_focus) ./ c0;
time_delays = time_delays - min(time_delays);

% create time varying source signals (one for each physical element)
source_sig = source_amp .* toneBurst(1/kgrid.dt, source_f0, source_cycles, 'SignalOffset', round(time_delays / kgrid.dt));

% create empty kWaveArray
karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

% add rectangular elements facing in the z-direction
for ind = 1:element_num
    
    % set element y position
    x_pos = 0 - (element_num * element_pitch / 2 - element_pitch / 2) + (ind - 1) * element_pitch;
    
    % add element (see note in header)
    karray.addRectElement([x_pos, 0, kgrid.z_vec(1)], element_width, element_length, rotation);
    
end

% move the array (see note in header)
karray.setArrayPosition(translation, rotation)

% assign binary mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% plot the off-grid source mask
voxelPlot(single(source.p_mask));
title('Off-grid source mask');

% assign source signals
source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
    
% --------------------
% MEDIUM
% --------------------

% assign medium properties
medium.sound_speed = c0;
medium.density = rho0;

% --------------------
% SENSOR
% --------------------

% set sensor mask to record central plane, not including the source point
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(:, Ny/2, 2:end) = 1;

% record the pressure
sensor.record = {'p_max'};

% --------------------
% SIMULATION
% --------------------

% set input options
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};

% run code
switch model
    case 1
        
        % MATLAB CPU
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
            input_args{:}, ...
            'DataCast', 'single', ...
            'PlotScale', [-1, 1] * source_amp);
        
    case 2
        
        % MATLAB GPU
        sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
            input_args{:}, ...
            'DataCast', 'gpuArray-single', ...
            'PlotScale', [-1, 1] * source_amp);
        
    case 3
        
        % C++
        sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, sensor, input_args{:});
        
    case 4
        
        % C++/CUDA GPU
        sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});
        
end

% reshape data
p_max = reshape(sensor_data.p_max, Nx, Nz - 1);

% =========================================================================
% VISUALISATION
% =========================================================================

% define axis vectors for plotting
z_vec = kgrid.z_vec(2:end, :) - kgrid.z_vec(1);
x_vec = kgrid.x_vec;

% plot the pressure field 
figure;
imagesc(1e3 * z_vec, 1e3 * x_vec, p_max);
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
axis image;
title('Pressure Field');