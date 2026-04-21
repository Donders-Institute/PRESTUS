function out_img = get_slice_by_label(img, slice_label, slice_n)
% GET_SLICE_BY_LABEL  Extract a 2D slice from a 3D volume along a named axis
%
%   Returns the slice at position slice_n along the axis specified by
%   slice_label ('x', 'y', or 'z'). Singleton dimensions are removed with squeeze.
%
% Use as:
%   out_img = get_slice_by_label(img, slice_label, slice_n)
%
% Input:
%   img         - [Nx x Ny x Nz] numeric, 3D image volume
%   slice_label - [1x1] char, axis to slice along: 'x', 'y', or 'z'
%   slice_n     - [1x1] positive integer, slice index along slice_label axis
%
% Output:
%   out_img     - [Mx x My] numeric, extracted 2D slice
%
% See also: SQUEEZE, NDGRID

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
