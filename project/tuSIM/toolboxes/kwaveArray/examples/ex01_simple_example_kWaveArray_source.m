% Script to demonstrate using kWaveArray to create an array with three
% arc-shaped transducer elements. All of the elements point towards the
% same focal position. Each element is driven with a different time-varying
% signal.
%
% author: Bradley Treeby
% date: 4th September 2018
% last update: 6th July 2019
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2018-2019 Bradley Treeby

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
% DEFINE KWAVEARRAY
% =========================================================================

% create empty array
karray = kWaveArray;

% add arc shaped element
elem_pos  = [10e-3, -40e-3];
radius    = 50e-3;
diameter  = 30e-3;
focus_pos = [-20e-3, 0];
karray.addArcElement(elem_pos, radius, diameter, focus_pos);

% add arc shaped element
elem_pos  = [20e-3, 0];
radius    = 50e-3;
diameter  = 30e-3;
focus_pos = [-20e-3, 0];
karray.addArcElement(elem_pos, radius, diameter, focus_pos);

% add arc shaped element
elem_pos  = [10e-3, 40e-3];
radius    = 50e-3;
diameter  = 30e-3;
focus_pos = [-20e-3, 0];
karray.addArcElement(elem_pos, radius, diameter, focus_pos);

% move the array down, and rotate (this moves all the elements together)
karray.setArrayPosition([10e-3, 0], 10);

% =========================================================================
% DEFINE GRID PROPERTIES
% =========================================================================

% grid properties
Nx = 256;
dx = 0.5e-3;
Ny = 256;
dy = 0.5e-3;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% medium properties
medium.sound_speed = 1500;

% time array
kgrid.makeTime(medium.sound_speed);

% =========================================================================
% SIMULATION
% =========================================================================

% assign binary mask from karray to the source mask
source.p_mask = karray.getArrayBinaryMask(kgrid);

% set source signals, one for each physical array element
f1 = 100e3;
f2 = 200e3;
f3 = 500e3;
sig1 = toneBurst(1/kgrid.dt, f1, 3);
sig2 = toneBurst(1/kgrid.dt, f2, 5);
sig3 = toneBurst(1/kgrid.dt, f3, 5);

% combine source signals into one array
source_signal = zeros(2, max(length(sig1), length(sig2)));
source_signal(1, 1:length(sig1)) = sig1;
source_signal(2, 1:length(sig2)) = sig2;
source_signal(3, 1:length(sig3)) = sig3;

% get distributed source signals (this automatically returns a weighted
% source signal for each grid point that forms part of the source)
source.p = karray.getDistributedSourceSignal(kgrid, source_signal);

% run k-Wave simulation (no sensor is used for this example)
kspaceFirstOrder2D(kgrid, medium, source, []);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot binary source mask (note it is non-local)
figure;
imagesc(kgrid.y_vec, kgrid.x_vec, source.p_mask);
axis image;
colormap(flipud(gray));

% overlay the physical source positions
hold on;
karray.plotArray(false);