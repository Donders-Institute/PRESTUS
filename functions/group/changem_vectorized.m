function B = changem_vectorized(A, newval, oldval)

% CHANGEM_VECTORIZED  Replace specified values in an array using vectorised operations
%
% Replaces elements in A that match oldval with the corresponding elements
% of newval. Equivalent to the Mapping Toolbox changem but vectorised.
%
% Use as:
%   B = changem_vectorized(A, newval, oldval)
%
% Input:
%   A      - input array
%   newval - [1xP] replacement values
%   oldval - [1xP] values to replace
%
% Output:
%   B - array with matching elements replaced
%
% See also: CREATE_GROUP_MNI_PLOTS

arguments
    A      {mustBeNumeric}
    newval (1,:) {mustBeNumeric}
    oldval (1,:) {mustBeNumeric}
end

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