% Create a 2D nifti image according to the desired tissue specification

%% Benchmark 1: Create a 1 mm version

% Initialize the matrix with global value 0
data = zeros(120, 70);
% Generate coordinate grids for each pixel
[X, Y] = meshgrid(1:70, 1:120);
% Fill region with white matter (soft tissue) labels
data(30:120,:) = 1;
% Reshape data to 3D for NIfTI compatibility
data_3d = reshape(data, [120, 70, 1]);
% Invert image to match what we would like
% data_3d = data_3d(end:-1:1,:,:);
% Write the data to a NIfTI file with 1x1x1 mm voxel size
niftiwrite(data_3d, 'benchmark1.nii');

% Note: copy this to a file called final_tissues.nii.gz in m2m_sub-xxx folder in the simnibs folder.

%% Benchmark 2: Create a 1 mm version

% Step 1: Initialize the matrix with global value 0 (dimensions: 120 rows (y), 70 columns (x))
data = zeros(120, 70);

% Step 2: Generate coordinate grids for each pixel (assuming 1 mm per pixel)
[X, Y] = meshgrid(1:70, 1:120);

% Step 3: Compute the center of curvature
center_x = 35;           % x-coordinate (column)
center_y = 90 - 75;      % y-coordinate (row) => 15

% Step 4: Calculate the Euclidean distance from the center
distance = sqrt((X - center_x).^2 + (Y - center_y).^2);

% Step 5: Assign values based on distance from center and y position
% Only fill below the center (Y >= center_y)
annulus_mask = (distance > 68.5) & (distance <= 75) & (Y >= center_y);
inner_mask   = (distance <= 68.5);

data(annulus_mask) = 4;    % Annular layer
data(inner_mask)   = 1;    % Inner region

% Step 6: Reshape data to 3D for NIfTI compatibility (singleton z-dimension)
data_3d = reshape(data, [120, 70, 1]);

% Invert image to match what we would like
data_3d = data_3d(end:-1:1,:,:);

% Step 7: Write the data to a NIfTI file with 1x1x1 mm voxel size
niftiwrite(data_3d, 'benchmark2.nii');

%% Create a 0.1 mm version

% % Step 1: Define new image size (1200 rows (y), 700 columns (x)) for 0.1 mm resolution
% data = zeros(1200, 700);
% % Step 2: Generate coordinate grids (now in 0.1 mm steps)
% [X, Y] = meshgrid(0.1:0.1:70, 0.1:0.1:120); % X: 0.1 to 70 mm, Y: 0.1 to 120 mm
% ...
