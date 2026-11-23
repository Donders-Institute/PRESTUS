function d = steering_path_lengths(tran, depths_ep_mm)
% STEERING_PATH_LENGTHS  Ring-centre to target distances used for geometric steering
%
% Same geometry as COMPUTE_PHASES, but unwrapped: the geometric steering
% phase of ring i for depth d at steering sound speed c is
%   phi_i(d) = 2*pi*f * d(d,i) / c      (COMPUTE_PHASES returns it mod 360°)
%
% Use as:
%   d = steering_path_lengths(tran, depths_ep_mm)
%
% Input:
%   tran         - transducer struct (uses .transducer.annular)
%   depths_ep_mm - focal depths from the exit plane [mm]
%
% Output:
%   d - [N_depths x N_elem] distances [m]
%
% See also: COMPUTE_PHASES, FIT_TRANSDUCER_MODEL

arguments
    tran (1,1) struct
    depths_ep_mm {mustBeNumeric}
end

    ann   = tran.transducer.annular;
    R_c   = ann.curv_radius_mm;
    r_mid = (ann.elem_id_mm(:)' / 2 + ann.elem_od_mm(:)' / 2) / 2;
    if isfield(ann, 'steer_r_mid_mm') && ~isempty(ann.steer_r_mid_mm)
        r_mid = ann.steer_r_mid_mm(:)';   % sub-rings steer with their ring (RING_PROFILE_TRANSDUCER)
    end
    z_t   = (depths_ep_mm(:) - ann.dist_geom_ep_mm) / 1000;   % aim w.r.t. natural focus [m]
    ex    = r_mid / 1000;
    ez    = -sqrt(R_c^2 - r_mid.^2) / 1000;
    d     = sqrt(ex.^2 + (ez - z_t).^2);
end
