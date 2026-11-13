function depth_mm = annular_bowl_depth_mm(geom)
% ANNULAR_BOWL_DEPTH_MM  Distance from the transducer position (bowl apex) to the exit plane
%
% The depth is curv_radius_mm - dist_geom_ep_mm (focal_distance_offset). It
% equals the rim sagitta only when dist_geom_ep_mm was derived from the
% aperture (VALIDATE_ANNULAR_TRANSDUCER); library entries can place the exit
% plane in front of the rim (CTX_250_001: 10.56 vs 7.06 mm). A flat element (curv_radius_mm = Inf) has no
% depth: its face is the exit plane.
%
% Use as:
%   depth_mm = annular_bowl_depth_mm(tr.annular)
%
% Input:
%   geom     - transducer geometry struct with curv_radius_mm and dist_geom_ep_mm [mm]
%
% Output:
%   depth_mm - bowl depth [mm]; 0 for flat elements or when either field is missing
%
% See also: VALIDATE_ANNULAR_TRANSDUCER, FOCAL_DISTANCE_CALCULATION

    depth_mm = 0;
    if isfield(geom, 'curv_radius_mm') && isfield(geom, 'dist_geom_ep_mm') && ...
            ~isempty(geom.curv_radius_mm) && ~isempty(geom.dist_geom_ep_mm) && ...
            isfinite(geom.curv_radius_mm)
        depth_mm = geom.curv_radius_mm - geom.dist_geom_ep_mm;
    end
end
