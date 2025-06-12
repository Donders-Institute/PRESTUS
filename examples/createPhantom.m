% Create a 2D nifti image according to the desired tissue specification

%% path management

% the following detects the path of this script (when run), which we can
% use to set up relative paths. in this example, PRESTUS is located in a
% directory called 'project/tools', and we would like to depoit the benchmark
% images into parallel folders project/data/simnibs/m2m-sub<xxx>

currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
cd(fullfile(pathstr,'..', '..'))
rootpath = pwd;

pn.benchmark1 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-001'); % specify the path to a SimNIBS m2m folder (optional)
pn.benchmark2 = fullfile(rootpath, '..', 'data', 'simnibs', 'm2m_sub-002'); % specify second benchmark as a different subject (optional)

% if no path is provided, the folder containing PRESTUS is specified by default
% if paths are specified and do not yet exist, create the folders
if isempty(pn.benchmark1)
    pn.benchmark1 = rootpath;
else
    if ~exist(pn.benchmark1)
        mkdir(pn.benchmark1);
    end
end
if isempty(pn.benchmark2)
    pn.benchmark2 = rootpath;
else
    if ~exist(pn.benchmark2)
        mkdir(pn.benchmark2);
    end
end

%% Benchmark 1

% This benchmark models a water / soft tissue interface.

pml = 20; % include a pml layer (on the axial dimension only)

% Initialize the matrix with global value 0
data = zeros(120+pml, 70);
% Generate coordinate grids for each pixel
[X, Y] = meshgrid(1:70, 1:120+pml);
% Fill region with white matter (soft tissue) labels
data([30:120]+pml,:) = 1;
% Reshape data to 3D for NIfTI compatibility
data_3d = reshape(data, [120+pml, 70, 1]);
% Invert image to match what we would like
% data_3d = data_3d(end:-1:1,:,:);
% Write the data to a NIfTI file with 1x1x1 mm voxel size
niftiwrite(data_3d, fullfile(pn.benchmark1, 'benchmark1.nii'));
gzip(fullfile(pn.benchmark1, 'benchmark1.nii')) % Note: PRESTUS expects compressed images
delete(fullfile(pn.benchmark1, 'benchmark1.nii')) % remove uncompressed image
movefile(fullfile(pn.benchmark1, 'benchmark1.nii.gz'), ...
    fullfile(pn.benchmark1, 'final_tissues.nii.gz')) % rename to SimNIBS convention

% Note: the output file should be called final_tissues.nii.gz and be 
% located in a "subject-"/benchmark specific simnibs m2m_sub-xxx folder

%% Benchmark 2

% This benchmark models a water / skull / soft tissue interface.

% Step 1: Initialize the matrix with global value 0 (dimensions: 120 rows (y), 70 columns (x))
data = zeros(120+pml, 70);

% Step 2: Generate coordinate grids for each pixel (assuming 1 mm per pixel)
[X, Y] = meshgrid(1:70, 1:120+pml);

% Step 3: Compute the center of curvature
center_x = 35;                  % x-coordinate (column)
center_y = 90 - 75;      % y-coordinate (row) => pml + 15

% Step 4: Calculate the Euclidean distance from the center
distance = sqrt((X - center_x).^2 + (Y - center_y).^2);

% Step 5: Assign values based on distance from center and y position
% Only fill below the center (Y >= center_y)
annulus_mask = (distance > 68.5) & (distance <= 75) & (Y >= center_y);
inner_mask   = (distance <= 68.5);

data(annulus_mask) = 4;    % Annular layer
data(inner_mask)   = 1;    % Inner region

% Step 6: Reshape data to 3D for NIfTI compatibility (singleton z-dimension)
data_3d = reshape(data, [pml + 120, 70, 1]);

% Invert image to match what we would like
data_3d = data_3d(end:-1:1,:,:);

% Step 7: Write the data to a NIfTI file with 1x1x1 mm voxel size
niftiwrite(data_3d, fullfile(pn.benchmark2, 'benchmark2.nii'));
gzip(fullfile(pn.benchmark2, 'benchmark2.nii'))
delete(fullfile(pn.benchmark2, 'benchmark2.nii'))
movefile(fullfile(pn.benchmark2, 'benchmark2.nii.gz'), ...
    fullfile(pn.benchmark2, 'final_tissues.nii.gz')) % rename to SimNIBS convention

%% Create a 0.1 mm version

% % Step 1: Define new image size (1200 rows (y), 700 columns (x)) for 0.1 mm resolution
% data = zeros(1200, 700);
% % Step 2: Generate coordinate grids (now in 0.1 mm steps)
% [X, Y] = meshgrid(0.1:0.1:70, 0.1:0.1:120); % X: 0.1 to 70 mm, Y: 0.1 to 120 mm
% ...
