function phase_deg = annular_geometric_phases(geom, apex_mm, aim_mm, freq_hz, c_water)
% ANNULAR_GEOMETRIC_PHASES  Ring phases that steer a curved annular array to a point, in water
%
% phi_e = k * (area-weighted mean distance of ring e to the aim point), so
% that all rings arrive at the aim point in phase (PRESTUS sources are
% A*sin(wt + phi): arrival phase phi - k*r). This is the drive that
% RECIPROCITY_SCORE_POSES assumes when it scores a pose.
%
% Use as:
%   phase_deg = annular_geometric_phases(tr.annular, apex_mm, aim_mm, freq_hz, c_water)
%
% Input:
%   geom     - annular geometry (curv_radius_mm, elem_id_mm, elem_od_mm)
%   apex_mm  - [1x3] bowl apex [mm] (any Cartesian frame, e.g. RAS)
%   aim_mm   - [1x3] aim point [mm] (same frame)
%   freq_hz  - frequency [Hz]
%   c_water  - steering sound speed [m/s]
%
% Output:
%   phase_deg - [1 x n_rings] phases [deg], 0-360, ring 1 at 0
%
% See also: ANNULAR_APERTURE_POINTS, RECIPROCITY_SCORE_POSES,
%           TIME_REVERSAL_POSE_PARAMETERS

arguments
    geom    (1,1) struct
    apex_mm (1,3) {mustBeNumeric}
    aim_mm  (1,3) {mustBeNumeric}
    freq_hz (1,1) {mustBeNumeric, mustBePositive}
    c_water (1,1) {mustBeNumeric, mustBePositive}
end

    apex_m = apex_mm * 1e-3;
    aim_m  = aim_mm * 1e-3;
    [pts, el, area] = annular_aperture_points(geom, apex_m, aim_m - apex_m, 0.25 * c_water / freq_hz);
    r   = sqrt(sum((pts - aim_m).^2, 2));
    n   = numel(geom.elem_od_mm);
    r_e = accumarray(el, area .* r, [n 1]) ./ accumarray(el, area, [n 1]);
    phi = 2 * pi * freq_hz / c_water * r_e;
    phase_deg = mod(rad2deg(phi - phi(1)), 360)';
end
