function [res_image, tr] = plot_t1_with_transducer(t1_image, voxel_size_mm, trans_pos_grid, focus_pos_grid, parameters, plot_options)

% PLOT_T1_WITH_TRANSDUCER  Plot a T1 slice with transducer mask overlaid
%
% Generates a visualisation of a T1 image slice oriented along the
% transducer axis, with the transducer bowl mask and focus position
% overlaid as an RGB image.
%
% Use as:
%   [res_image, tr] = plot_t1_with_transducer(t1_image, voxel_size_mm, ...
%       trans_pos_grid, focus_pos_grid, parameters)
%   [res_image, tr] = plot_t1_with_transducer(..., plot_options)
%
% Input:
%   t1_image       - [Nx x Ny x Nz] T1-weighted image volume
%   voxel_size_mm  - voxel size [mm]
%   trans_pos_grid - [1x3] transducer position in grid coordinates
%   focus_pos_grid - [1x3] focus position in grid coordinates
%   parameters     - (1,1) simulation parameters struct
%   plot_options   - name-value plotting options (slice_dim, slice_ind)
%
% Output:
%   res_image - [Mx x My x 3] RGB image of the T1 slice with overlays
%   tr        - transducer parameter struct from transducer_setup
%
% See also: PLOT_OVERLAY, TRANSDUCER_SETUP

    % Checks if all data is in the right format
    arguments
        t1_image (:,:,:) double
        voxel_size_mm (1,1) double
        trans_pos_grid (1,:) double
        focus_pos_grid (1,:) double
        parameters (1,1) struct
        plot_options.slice_dim (1,1) double = 2 % Default slicing dimension (y-axis)
        plot_options.slice_ind (1,1) double = 0 % Default slice index (transducer position)
    end

    %% Pad T1 image to ensure it includes transducer/focus position (Negative/Pre-padding)
    min_coords = min([1 1 1; trans_pos_grid; focus_pos_grid]);
    if any(min_coords < 1)
        pad_amount = ceil(1 - min_coords); % Amount to shift to make min 1
        % Pad at the beginning (pre)
        t1_image = cat(1, zeros([pad_amount(1), size(t1_image,2), size(t1_image,3)], class(t1_image)), t1_image);
        t1_image = cat(2, zeros([size(t1_image,1), pad_amount(2), size(t1_image,3)], class(t1_image)), t1_image);
        t1_image = cat(3, zeros([size(t1_image,1), size(t1_image,2), pad_amount(3)], class(t1_image)), t1_image);
        
        % Shift coordinates
        trans_pos_grid = trans_pos_grid + pad_amount;
        focus_pos_grid = focus_pos_grid + pad_amount;
        
        % Shift slice index if provided
         if ~isempty(plot_options.slice_ind) && plot_options.slice_ind ~= 0
             plot_options.slice_ind = plot_options.slice_ind + pad_amount(plot_options.slice_dim);
         end
    end

    %% Determine slice index if not provided
    if isempty(plot_options.slice_ind) || plot_options.slice_ind == 0
        plot_options.slice_ind = trans_pos_grid(plot_options.slice_dim); % Use transducer position as default slice index
    end

    %% Pad T1 image to ensure it includes transducer position (Positive/Post-padding)
    im_size = max([size(t1_image); trans_pos_grid; focus_pos_grid]);
    if ~isequal(im_size, size(t1_image))
        pad_post = im_size - size(t1_image);
        t1_image = cat(1, t1_image, zeros([pad_post(1), size(t1_image,2), size(t1_image,3)], class(t1_image)));
        t1_image = cat(2, t1_image, zeros([size(t1_image,1), pad_post(2), size(t1_image,3)], class(t1_image)));
        t1_image = cat(3, t1_image, zeros([size(t1_image,1), size(t1_image,2), pad_post(3)], class(t1_image)));
    end

    %% Create transducer mask and setup parameters
    [transducer_bowl, ~, tr] = transducer_setup(parameters.transducer(1), ...
                                                             trans_pos_grid, focus_pos_grid, ...
                                                             im_size, voxel_size_mm, parameters);

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
    try
    res_image(focus_pos_grid(focus_pos_ind(1)) + (-2:2), ... % Draw square around focus position in green channel
              focus_pos_grid(focus_pos_ind(2)) + (-2:2), ...
              2) = 1;
    catch warn('Boxes may be out of bounds...');
    end

end