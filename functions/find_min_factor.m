function min_factor = find_min_factor(min_number, max_number)

% FIND_MIN_FACTOR Finds the number with the smallest maximum factor in a range.
%
% This function computes the number within a specified range `[min_number, max_number]` 
% that has the smallest maximum factor. It uses MATLAB's `factor` function to 
% determine the factors of each number in the range.
%
% Input:
%   min_number - Integer specifying the start of the range (inclusive).
%   max_number - Integer specifying the end of the range (inclusive).
%
% Output:
%   min_factor - The number within the range `[min_number, max_number]` that 
%                has the smallest maximum factor.

    % Initialize arrays to store the number of factors and maximum factor for each number
    facs = zeros(1, max_number - min_number); % Number of factors for each number
    fac_max = facs;                           % Maximum factor for each number

    % Loop through each number in the specified range
    for index = min_number:max_number
        % Compute all factors of the current number
        current_factors = factor(index);

        % Store the number of factors and the maximum factor
        facs(index - min_number + 1) = length(current_factors);
        fac_max(index - min_number + 1) = max(current_factors);
    end

    % Find the number(s) with the smallest maximum factor
    % Add `min_number` to adjust for zero-based indexing
    min_factor = min(min_number + find(fac_max == min(fac_max)) - 1);

end
