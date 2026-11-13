function [flhm_center, z_first, z_last] = compute_flhm_center(I, z, focal_reference)
% COMPUTE_FLHM_CENTER  Focal length at half maximum (FLHM) centre of an axial profile.
%
%   [flhm_center, z_first, z_last] = compute_flhm_center(I, z)
%   [flhm_center, z_first, z_last] = compute_flhm_center(I, z, focal_reference)
%
%   Finds the contiguous region around the global maximum where
%   I >= max(I)*frac, and returns the midpoint of that region's axial
%   extent, where frac is the threshold fraction for the requested
%   focal_reference convention (see FOCAL_REFERENCE_TO_FRAC).
%
%   Falls back to the peak position when no threshold crossing is found
%   (e.g. profile is entirely above threshold, or has only one sample).
%
%   The threshold convention varies by manufacturer/equipment: some define
%   "focal depth" as the -3 dB intensity point ('flhm_3db', I >= peak/2,
%   ITRUSST/Martin et al. 2025 default), others as the -6 dB intensity
%   point ('flhm_6db', I >= peak/4, i.e. half of peak PRESSURE amplitude,
%   as in Imasonic/IGT steering tables).
%
% Input:
%   I               — intensity vector (any non-negative units)
%   z               — corresponding axial position vector [same length as I, any unit]
%   focal_reference — convention for the FLHM threshold: 'flhm_3db'
%                     (default, I >= peak/2) | 'flhm_6db' (I >= peak/4) |
%                     'peak' (centre = peak position; edges at -3 dB).
%                     Should match calibration.focal_reference for the
%                     equipment being analysed.
%
% Output:
%   flhm_center — midpoint of the FLHM region [same unit as z]
%   z_first     — proximal edge of the FLHM region
%   z_last      — distal edge of the FLHM region
%
% See also: FOCAL_REFERENCE_TO_FRAC, COMPUTE_FOCAL_METRICS

arguments
    I               {mustBeNumeric}
    z               {mustBeNumeric}
    focal_reference (1,:) char = 'flhm_3db'
end

% 'peak' has no FLHM threshold of its own: the centre is the peak itself, and
% the returned edges fall back to the ITRUSST -3 dB threshold.
use_peak = strcmp(focal_reference, 'peak');
if use_peak
    threshold_frac = focal_reference_to_frac('flhm_3db');
else
    threshold_frac = focal_reference_to_frac(focal_reference);
end

[pk, pk_idx] = max(I);
half_max = pk * threshold_frac;

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
if use_peak
    flhm_center = z(pk_idx);
else
    flhm_center = (z_first + z_last) / 2;
end
