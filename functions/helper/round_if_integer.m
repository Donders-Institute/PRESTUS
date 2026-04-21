function x = round_if_integer(x, err_msg)
% ROUND_IF_INTEGER  Round a value to the nearest integer or error if not close enough
%
% Asserts that x is within 1e-6 of an integer, then rounds. Used to
% validate time-step count calculations where non-integer counts indicate
% an invalid protocol configuration.
%
% Use as:
%   x = round_if_integer(x, err_msg)
%
% Input:
%   x       - value or array to check and round
%   err_msg - error message to display if x is not close to an integer
%
% Output:
%   x - rounded value
%
% See also: THERMAL_PARAMETERS

arguments
    x       {mustBeNumeric}
    err_msg (1,:) char
end

    % Assert that the input value is sufficiently close to an integer
    assert(abs(round(x) - x) < 1e-6, err_msg);

    % Round the input value to the nearest integer
    x = round(x);
end