function plot_coronal_slices(img, coord, options)

% PLOT_CORONAL_SLICES Visualizes coronal slices of a 3D image.
%
% This function displays coronal, sagittal, and axial slices of a 3D image (`img`) 
% at specified coordinates (`coord`). It supports integer-labeled images (e.g., 
% segmented masks) and applies a colormap with labels if provided. The function 
% also rescales intensity values if required.
%
% Input:
%   img    - [Nx x Ny x Nz] 3D matrix representing the image to visualize.
%   coord  - [3x1] array specifying the slice coordinates for coronal, sagittal, 
%            and axial views (default: center of the image).
%   options:
%     * cmap    - Custom colormap for integer-labeled images (default: SimNIBS-based colors).
%     * rescale - Boolean flag to rescale intensity values to the range [1, 255] (default: 0).
%     * labels  - Cell array of strings specifying labels for integer-labeled images 
%                 (default: SimNIBS tissue labels).
%
% Output:
%   None. The function displays the slices in a montage with optional legends.

    arguments
        img (:,:,:)
        coord (3,1) = round(size(img)/2) % Default: center of the image
        options.cmap = [] % Custom colormap
        options.rescale = 0 % Rescale intensity values
        options.labels = ["Background", "White Matter", "Gray Matter", "CSF", "Bone", ...
                          "Scalp", "Eye balls", "Compact bone", "Spongy bone", ...
                          "Blood", "Muscle", "Electrode", "Saline or gel"] % Default SimNIBS labels
    end

    %% Define default colors for SimNIBS tissue labels
    default_colors = [230,230,230; % Background
                      129,129,129; % White Matter
                      104,163,255; % Gray Matter
                      255,239,179; % CSF
                      255,166,133; % Bone
                      255,240,0;   % Scalp
                      255,239,179; % Eye balls
                      255,138,57;  % Compact bone
                      0,65,142;    % Spongy bone
                      0,118,14;    % Blood
                      37,79,255;   % Muscle
                      103,255,226];% Electrode

    %% Set colormap for integer-labeled images if not provided
    if isinteger(img) && isempty(options.cmap)
        N = length(unique(img(:))); % Number of unique labels in the image
        options.cmap = [[0,0,0]; default_colors(1:N,:) / 255]; % Add black for background and normalize colors
    end

    %% Rescale intensity values if requested
    if options.rescale
        img = round(rescale(img, 1, 255)); % Rescale to range [1, 255]
    end

    %% Visualize slices in a montage
    figure;
    montage({imrotate(squeeze(img(coord(1),:,:)),90), ... % Coronal slice
             imrotate(squeeze(img(:,coord(2),:)),90), ... % Sagittal slice
             squeeze(img(:,:,coord(3)))'}, ...           % Axial slice
             options.cmap, 'Size', [1 3]);               % Apply colormap and arrange slices in a row

    %% Add legend for integer-labeled images (optional)
    if isinteger(img)
        hold on;
        N = length(unique(img(:))); % Number of unique labels in the image
        
        % Generate lines with colors matching the colormap for legend entries
        L = line(ones(N), ones(N), 'LineWidth', 2); 
        set(L, {'color'}, mat2cell(options.cmap(1:N,:), ones(1,N), 3)); 
        
        % Add legend entries based on unique labels in the image
        legend(options.labels(unique(img(:)) + 1)); 
    end

end