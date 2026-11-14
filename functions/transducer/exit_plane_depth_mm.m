function depth_mm = exit_plane_depth_mm(geom)
% EXIT_PLANE_DEPTH_MM  Distance from the bowl apex to the exit plane used for placement
%
% The exit plane is the library one, curv_radius_mm - dist_geom_ep_mm
% (ANNULAR_BOWL_DEPTH_MM; BabelBrain's FocalLength - NaturalOutPlaneDistance).
% LOAD_TRANSDUCER_PARAMETERS sets a missing dist_geom_ep_mm from the bowl rim
% (VALIDATE_ANNULAR_TRANSDUCER, VALIDATE_MATRIX_TRANSDUCER). A struct that did
% not pass through it falls back to the same rim, with a warning: the sagitta
% of outer_diameter_mm (matrix) or of the largest element (annular).
%
% Use as:
%   depth_mm = exit_plane_depth_mm(tr.(tr.type))
%
% Input:
%   geom     - transducer geometry struct with curv_radius_mm, dist_geom_ep_mm
%              and, for the fallback, outer_diameter_mm or elem_od_mm [mm]
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
        if isfield(geom, 'outer_diameter_mm') && ~isempty(geom.outer_diameter_mm)
            D = geom.outer_diameter_mm;            % matrix aperture
        else
            D = max(geom.elem_od_mm);              % annular: largest element
        end
        depth_mm = R - sqrt(R^2 - (D/2)^2);
        warn('prestus:placement:rimExitPlane', ...
            'transducer has no dist_geom_ep_mm; placement uses the bowl rim (%.2f mm from the apex) as exit plane.', depth_mm);
    end
end
