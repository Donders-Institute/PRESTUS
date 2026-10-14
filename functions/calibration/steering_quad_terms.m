function A = steering_quad_terms(Z)
% STEERING_QUAD_TERMS  Design matrix of a full quadratic surface in three coordinates
%
% Use as:
%   A = steering_quad_terms(Z)
%
% Input:
%   Z - [N x 3] (standardised) coordinates
%
% Output:
%   A - [N x 10] columns 1, x, y, z, x^2, y^2, z^2, xy, xz, yz
%
% See also: FIT_STEERING_CORRECTION, STEERING_CORRECTION_FACTOR

    A = [ones(size(Z, 1), 1), Z, Z.^2, Z(:, 1) .* Z(:, 2), Z(:, 1) .* Z(:, 3), Z(:, 2) .* Z(:, 3)];
end
