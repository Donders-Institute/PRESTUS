function [pts_m, elem, area_m2] = annular_aperture_points(geom, apex_m, axis_u, spacing_m)
% ANNULAR_APERTURE_POINTS  Sample points on the radiating surface of a curved annular array
%
% The bowl is a spherical cap of radius curv_radius_mm with its apex (the
% centre of the rear surface, PRESTUS trans_pos) at apex_m and its axis along
% axis_u; the geometric focus is apex_m + R*axis_u. Ring e covers the polar
% angles asin(elem_id_mm/2/R) to asin(elem_od_mm/2/R), seen from the
% geometric focus. Points are spaced ~spacing_m along both surface
% directions, each with its surface area, so sum(area_m2(elem == e)) is the
% ring area.
%
% Used to score transducer poses on a reverse field (RECIPROCITY_SCORE_POSES),
% where building kWaveArray weights per pose would be too slow. Flat elements
% are not supported.
%
% Use as:
%   [pts_m, elem, area_m2] = annular_aperture_points(tr.annular, apex_m, axis_u, spacing_m)
%
% Input:
%   geom      - annular geometry: curv_radius_mm, elem_id_mm, elem_od_mm [mm]
%   apex_m    - [1 x 3] bowl apex [m]
%   axis_u    - [1 x 3] beam direction (apex -> focus)
%   spacing_m - point spacing on the surface [m] (e.g. lambda/4)
%
% Output:
%   pts_m   - [N x 3] point positions [m]
%   elem    - [N x 1] ring index of each point
%   area_m2 - [N x 1] surface area per point [m^2]
%
% See also: RECIPROCITY_SCORE_POSES, SOURCE_CREATE

arguments
    geom      (1,1) struct
    apex_m    (1,3) {mustBeNumeric}
    axis_u    (1,3) {mustBeNumeric}
    spacing_m (1,1) {mustBeNumeric, mustBePositive}
end

    R = geom.curv_radius_mm * 1e-3;
    if ~isfinite(R)
        error('PRESTUS:annularAperturePoints:flat', 'Flat annular elements are not supported.');
    end
    u = axis_u / norm(axis_u);
    % Orthonormal pair perpendicular to u
    [~, imin] = min(abs(u));
    e1 = zeros(1, 3); e1(imin) = 1;
    e1 = e1 - dot(e1, u) * u; e1 = e1 / norm(e1);
    e2 = cross(u, e1);
    centre = apex_m + R * u;

    id = geom.elem_id_mm(:)' * 1e-3;
    od = geom.elem_od_mm(:)' * 1e-3;
    pts_c = {}; elem_c = {}; area_c = {};
    for e = 1:numel(od)
        th0 = asin(min(id(e) / 2 / R, 1));
        th1 = asin(min(od(e) / 2 / R, 1));
        n_th = max(1, ceil((th1 - th0) * R / spacing_m));
        dth  = (th1 - th0) / n_th;
        for it = 1:n_th
            th   = th0 + (it - 0.5) * dth;
            circ = 2 * pi * R * sin(th);
            n_al = max(1, ceil(circ / spacing_m));
            al   = (0:n_al-1)' * 2 * pi / n_al;
            dirs = cos(th) * u + sin(th) * (cos(al) * e1 + sin(al) * e2);
            pts_c{end+1}  = centre - R * dirs;                          %#ok<AGROW>
            elem_c{end+1} = repmat(e, n_al, 1);                         %#ok<AGROW>
            area_c{end+1} = repmat(R^2 * sin(th) * dth * 2 * pi / n_al, n_al, 1); %#ok<AGROW>
        end
    end
    pts_m   = vertcat(pts_c{:});
    elem    = vertcat(elem_c{:});
    area_m2 = vertcat(area_c{:});
end
