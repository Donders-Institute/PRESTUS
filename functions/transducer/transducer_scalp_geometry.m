function [trans_pos, geom_focus, ex_plane, norm_v] = ...
    transducer_scalp_geometry(scalp_pts, target_vox, tr, pixel_size, skin_gap_mm)
% TRANSDUCER_SCALP_GEOMETRY  Place a transducer on the scalp with a standoff gap
%
% Given one or more scalp/outer-boundary points and a target, returns the
% transducer position offset outward along the focal axis by a skin gap plus
% the apex-to-exit-plane distance (EXIT_PLANE_DEPTH_MM), together with the
% geometric focus and exit plane. This is the exact geometry the heuristic
% placement applies to every candidate position (TP_CANDIDATE_MESH): the
% transducer axis points from the target to the scalp point, and the exit
% plane sits skin_gap_mm off the scalp (placement.skin_gap_mm).
%
% Use as:
%   trans_pos = transducer_scalp_geometry(scalp_pts, target_vox, tr, pixel_size)
%   [trans_pos, geom_focus, ex_plane, norm_v] = ...
%       transducer_scalp_geometry(scalp_pts, target_vox, tr, pixel_size, skin_gap_mm)
%
% Input:
%   scalp_pts   - [N x 3] scalp / outer-boundary voxel coordinates
%   target_vox  - [1 x 3] focus / target voxel coordinates
%   tr          - one transducer struct (parameters.transducer(1)); uses
%                 tr.(tr.type).curv_radius_mm, dist_geom_ep_mm and elem_od_mm
%   pixel_size  - voxel size [mm]
%   skin_gap_mm - gap from the scalp to the exit plane along the axis [mm]
%                 (optional, default 0)
%
% Output:
%   trans_pos  - [N x 3] transducer position(s) in voxel space
%   geom_focus - [N x 3] geometric focus position(s) in voxel space
%   ex_plane   - [N x 3] exit-plane position(s) in voxel space
%   norm_v     - [N x 3] unit transducer axis (target -> scalp), per point
%
% See also: TP_CANDIDATE_MESH, TP_SCALP_BOUNDARY, EXIT_PLANE_DEPTH_MM

arguments
    scalp_pts               % [N x 3] (double or gpuArray; untyped to preserve gpuArray)
    target_vox              % [1 x 3] (double or gpuArray)
    tr          struct
    pixel_size  (1,1) double
    skin_gap_mm (1,1) double = 5
end

    % Unit transducer axis from target to each scalp point
    d      = scalp_pts - target_vox;
    norm_v = d ./ sqrt(sum(d.^2, 2));

    % Transducer geometry (annular-style fields, as in tp_candidate_mesh)
    g                = tr.(tr.type);
    dist_tp_to_ep_mm = exit_plane_depth_mm(g);

    % Offset off the scalp along the axis: skin gap + apex-to-exit-plane distance
    pos_shift_mm = skin_gap_mm + dist_tp_to_ep_mm;
    trans_pos    = scalp_pts + norm_v * (pos_shift_mm / pixel_size);

    % Exit plane from the apex (finite also for flat elements, whose
    % geometric focus lies at infinity)
    geom_focus = trans_pos - norm_v * (g.curv_radius_mm / pixel_size);
    ex_plane   = trans_pos - norm_v * (dist_tp_to_ep_mm / pixel_size);
end
