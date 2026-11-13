function depth_mm = exit_plane_depth_mm(geom)
% EXIT_PLANE_DEPTH_MM  Distance from the bowl apex to the exit plane used for placement
%
% The exit plane is the library one, curv_radius_mm - dist_geom_ep_mm
% (ANNULAR_BOWL_DEPTH_MM; BabelBrain's FocalLength - NaturalOutPlaneDistance).
% A transducer without dist_geom_ep_mm falls back to the bowl rim, the
% sagitta of the largest element, with a warning.
%
% Use as:
%   depth_mm = exit_plane_depth_mm(tr.(tr.type))
%
% Input:
%   geom     - transducer geometry struct with curv_radius_mm, elem_od_mm and
%              (preferably) dist_geom_ep_mm [mm]
%
% Output:
%   depth_mm - apex-to-exit-plane distance along the beam axis [mm]
%
% See also: ANNULAR_BOWL_DEPTH_MM, TRANSDUCER_SCALP_GEOMETRY

    if isfield(geom, 'dist_geom_ep_mm') && ~isempty(geom.dist_geom_ep_mm)
        depth_mm = annular_bowl_depth_mm(geom);
        return;
    end
    depth_mm = 0;
    if isfield(geom, 'curv_radius_mm') && ~isempty(geom.curv_radius_mm) && isfinite(geom.curv_radius_mm)
        R = geom.curv_radius_mm;
        depth_mm = R - 0.5 * sqrt(4*R^2 - max(geom.elem_od_mm)^2);
        warn('prestus:placement:rimExitPlane', ...
            'transducer has no dist_geom_ep_mm; placement uses the bowl rim (%.2f mm from the apex) as exit plane.', depth_mm);
    end
end
