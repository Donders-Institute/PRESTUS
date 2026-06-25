function [flhm_center, z_first, z_last] = compute_flhm_center(I, z)
% COMPUTE_FLHM_CENTER  Full-length-half-maximum centre of an axial profile.
%
%   [flhm_center, z_first, z_last] = compute_flhm_center(I, z)
%
%   Finds the contiguous region around the global maximum where I >= max(I)/2,
%   and returns the midpoint of that region's axial extent.
%
%   Falls back to the peak position when no half-maximum crossing is found
%   (e.g. profile is entirely above half-max, or has only one sample).
%
% Input:
%   I  — intensity vector (any non-negative units)
%   z  — corresponding axial position vector [same length as I, any unit]
%
% Output:
%   flhm_center — midpoint of the FLHM region [same unit as z]
%   z_first     — proximal edge of the FLHM region
%   z_last      — distal edge of the FLHM region

[pk, pk_idx] = max(I);
half_max = pk / 2;

% Walk proximal from peak to find first crossing
i_prox = pk_idx;
while i_prox > 1 && I(i_prox - 1) >= half_max
    i_prox = i_prox - 1;
end

% Walk distal from peak to find last crossing
i_dist = pk_idx;
while i_dist < numel(I) && I(i_dist + 1) >= half_max
    i_dist = i_dist + 1;
end

z_first     = z(i_prox);
z_last      = z(i_dist);
flhm_center = (z_first + z_last) / 2;
