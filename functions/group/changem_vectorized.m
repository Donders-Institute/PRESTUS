function B = changem_vectorized(A, newval, oldval)

% CHANGEM_VECTORIZED Replace values in an array using vectorized operations.
%
% This function replaces elements in array `A` that match values in `oldval` 
% with corresponding values from `newval`. It uses vectorized operations for 
% efficient computation.
%
% Input:
%   A       - [MxN] array of input values.
%   newval  - [1xP] array of new values to replace matching elements in `A`.
%   oldval  - [1xP] array of old values to be replaced in `A`.
%
% Output:
%   B       - [MxN] array with replaced values.

    %% Step 1: Initialize output array
    % Start by copying the input array `A` to the output array `B`.
    B = A;

    %% Step 2: Find matching elements
    % Use `bsxfun` to compare each element of `A` with all elements of `oldval`.
    % The result is a logical matrix where each row corresponds to an element in `A`
    % and each column corresponds to a value in `oldval`. The maximum value along 
    % each row indicates whether a match was found.
    [valid, id] = max(bsxfun(@eq, A(:), oldval(:).'), [], 2); %//'

    %% Step 3: Replace matched elements
    % Replace elements in `B` where matches were found (`valid == true`) with 
    % the corresponding values from `newval`. The index `id(valid)` maps the 
    % matched positions in `oldval` to their replacements in `newval`.
    B(valid) = newval(id(valid));

end