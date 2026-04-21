function min_factor = find_min_factor(min_number, max_number)
% FIND_MIN_FACTOR  Find the integer in a range whose largest prime factor is smallest
%
%   Searches [min_number, max_number] and returns the value whose maximum
%   prime factor is minimised. Used to choose FFT-friendly grid dimensions for
%   k-Wave simulations; k-Wave FFTs are most efficient when the grid dimension
%   is a product of small primes.
%
% Use as:
%   min_factor = find_min_factor(min_number, max_number)
%
% Input:
%   min_number - [1x1] lower bound of the search range (inclusive)
%   max_number - [1x1] upper bound of the search range (inclusive)
%
% Output:
%   min_factor - [1x1] integer in [min_number, max_number] with the smallest maximum prime factor
%
% See also: FACTOR, PREPROC_CROP_GRID

arguments
    min_number (1,1) {mustBeInteger, mustBePositive}
    max_number (1,1) {mustBeInteger, mustBePositive}
end

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
