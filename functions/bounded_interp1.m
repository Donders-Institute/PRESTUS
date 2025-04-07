function new_y = bounded_interp1(x, y, new_x)

% BOUNDED_INTERP1 Performs linear interpolation with bounds handling.
%
% This function interpolates values using `interp1` but ensures that values 
% outside the bounds of the input `x` are assigned the corresponding boundary 
% values of `y`. It prevents extrapolation by bounding the results to the 
% minimum and maximum values of `x`.
%
% Input:
%   x      - [1xN] array of input x-coordinates.
%   y      - [1xN] array of input y-values corresponding to `x`.
%   new_x  - [1xM] array of new x-coordinates for interpolation.
%
% Output:
%   new_y  - [1xM] array of interpolated y-values, bounded by the minimum and 
%            maximum values of `x`.

    %% Step 1: Perform linear interpolation
    % Use MATLAB's `interp1` to interpolate values at `new_x` based on `x` and `y`.
    % Points outside the range of `x` are assigned NaN by default.
    new_y = interp1(x, y, new_x, 'linear', NaN);

    %% Step 2: Handle values below the minimum bound
    new_y(new_x < min(x)) = y(x == min(x));

    %% Step 3: Handle values above the maximum bound
    new_y(new_x > max(x)) = y(x == max(x));
end
