function x = round_if_integer(x, err_msg)

% ROUND_IF_INTEGER Rounds a value if it is sufficiently close to an integer.
%
% This function checks whether the input value `x` is sufficiently close to an integer 
% (within a tolerance of 1e-6). If the condition is met, it rounds `x` to the nearest integer. 
% Otherwise, it raises an error with the provided error message.
%
% Input:
%   x        - Numeric value or array to be checked and rounded.
%   err_msg  - String specifying the error message to display if `x` is not close to an integer.
%
% Output:
%   x        - Rounded numeric value or array.

    % Assert that the input value is sufficiently close to an integer
    assert(abs(round(x) - x) < 1e-6, err_msg);

    % Round the input value to the nearest integer
    x = round(x);
end