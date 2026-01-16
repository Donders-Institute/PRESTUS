function out_img = get_slice_by_label(img, slice_label, slice_n)

% GET_SLICE_BY_LABEL Extracts a specific slice from a 3D image based on axis label and slice number.
%
% This function extracts a 2D slice from a 3D image (`img`) along the specified 
% axis (`slice_label`) at the given slice number (`slice_n`). The axis can be 
% 'x', 'y', or 'z'.
%
% Input:
%   img         - [Nx x Ny x Nz] matrix representing the 3D image.
%   slice_label - String specifying the axis along which to extract the slice ('x', 'y', or 'z').
%   slice_n     - Integer specifying the slice number along the specified axis.
%
% Output:
%   out_img     - [Mx x My] matrix representing the extracted 2D slice.
%
% Notes:
%   - If an invalid `slice_label` is provided, an error is raised.
%   - The function uses MATLAB's `squeeze` to remove singleton dimensions.

    % Initialize ranges for slicing along each axis
    slice_x = 1:size(img,1); % Range for x-axis
    slice_y = 1:size(img,2); % Range for y-axis
    slice_z = 1:size(img,3); % Range for z-axis

    % Update slicing range based on the specified axis
    if strcmp(slice_label, 'x')
        slice_x = slice_n; % Extract a single x-slice
    elseif strcmp(slice_label, 'y')
        slice_y = slice_n; % Extract a single y-slice
    elseif strcmp(slice_label, 'z')
        slice_z = slice_n; % Extract a single z-slice
    else
        % Raise an error for invalid axis labels
        error("slice_label must be one of 'x', 'y', or 'z', and slice_n must be a valid slice number.");
    end

    % Extract the specified 2D slice and remove singleton dimensions
    out_img = squeeze(img(slice_x, slice_y, slice_z));
end
