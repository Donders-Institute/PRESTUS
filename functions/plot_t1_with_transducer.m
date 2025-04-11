function [res_image, transducer_pars] = plot_t1_with_transducer(t1_image, voxel_size_mm, trans_pos_grid, focus_pos_grid, parameters, plot_options)

% PLOT_T1_WITH_TRANSDUCER Creates a plot of a T1 slice oriented along the transducer's axis.
%
% This function generates a visualization of a T1-weighted image slice overlaid with 
% the transducer mask and focus position. The slice is oriented along the transducer's axis.
%
% Input:
%   t1_image        - [Nx x Ny x Nz] matrix representing the T1-weighted image.
%   voxel_size_mm   - Scalar specifying the voxel size in mm.
%   trans_pos_grid  - [1x3] array specifying the transducer position in grid coordinates.
%   focus_pos_grid  - [1x3] array specifying the focus position in grid coordinates.
%   parameters      - Struct containing transducer configuration parameters.
%   plot_options    - Struct containing optional plotting settings:
%                     * slice_dim: Dimension for slicing (default: 2).
%                     * slice_ind: Index of the slice along the specified dimension (default: transducer position).
%
% Output:
%   res_image       - [Mx x My x 3] RGB image matrix representing the T1 slice with overlays.
%   transducer_pars - Struct containing parameters of the transducer setup.

    % Checks if all data is in the right format
    arguments
        t1_image (:,:,:) double
        voxel_size_mm (1,1) double
        trans_pos_grid (1,3) double
        focus_pos_grid (1,3) double
        parameters struct
        plot_options.slice_dim (1,1) double = 2 % Default slicing dimension (y-axis)
        plot_options.slice_ind (1,1) double = 0 % Default slice index (transducer position)
    end

    %% Determine slice index if not provided
    if isempty(plot_options.slice_ind) || plot_options.slice_ind == 0
        plot_options.slice_ind = trans_pos_grid(plot_options.slice_dim); % Use transducer position as default slice index
    end

    %% Pad T1 image to ensure it includes transducer position
    im_size = max(size(t1_image), trans_pos_grid); % Determine required image size based on transducer position
    if ~isequal(im_size, size(t1_image))
        t1_image = padarray(t1_image, im_size - size(t1_image), 'post'); % Pad image to include transducer position
    end

    %% Create transducer mask and setup parameters
    [transducer_bowl, ~, transducer_pars] = transducer_setup(parameters.transducer, ...
                                                             trans_pos_grid, focus_pos_grid, ...
                                                             im_size, voxel_size_mm);

    %% Define slicing indices based on specified dimension and index
    slice_pointer = repmat({':'}, 1, 3); % Initialize slicing indices for all dimensions
    slice_pointer{plot_options.slice_dim} = plot_options.slice_ind; % Set slicing index for specified dimension

    %% Extract T1 slice along specified axis
    t1_slice = double(t1_image(slice_pointer{:})); % Extract slice using slicing indices
    res_image = repmat(mat2gray(squeeze(t1_slice)), [1, 1, 3]); % Normalize and convert to RGB format

    %% Overlay transducer mask on T1 slice
    transducer_bowl = squeeze(transducer_bowl(slice_pointer{:})); % Extract corresponding transducer mask slice
    res_image(:, :, 1) = max(res_image(:, :, 1), transducer_bowl); % Overlay mask on red channel

    %% Highlight focus position on T1 slice
    focus_pos_ind = find(~cellfun(@(x) ~strcmp(x, ':'), slice_pointer)); % Identify dimensions used for slicing
    res_image(focus_pos_grid(focus_pos_ind(1)) + (-2:2), ... % Draw square around focus position in green channel
              focus_pos_grid(focus_pos_ind(2)) + (-2:2), ...
              2) = 1;

end