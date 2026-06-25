function phases = compute_phases(SOUND_SPEED_WATER, tran, desired_focal_distance_ep)
% COMPUTE_PHASES  Compute element steering phases for a target focal depth
%
% Calculates the per-element phase delays required to focus an annular
% transducer at a specified depth by computing geometric path lengths from
% element positions to the target point.
%
% Element positions are derived from ring geometry (elem_id_mm, elem_od_mm,
% curv_radius_mm) stored in the transducer struct — no external file required.
%
% Use as:
%   phases = compute_phases(SOUND_SPEED_WATER, tran, desired_focal_distance_ep)
%
% Input:
%   SOUND_SPEED_WATER         - sound speed in water [m/s]
%   tran                      - struct with transducer parameters; requires
%                               tran.transducer.freq_hz,
%                               tran.transducer.annular.elem_id_mm,
%                               tran.transducer.annular.elem_od_mm,
%                               tran.transducer.annular.curv_radius_mm,
%                               tran.transducer.annular.dist_geom_ep_mm
%   desired_focal_distance_ep - focal depth relative to the exit plane [mm]
%
% Output:
%   phases - [1xN] computed phase values per element [°], one per annular ring
%
% See also: SET_REAL_PHASES

arguments
    SOUND_SPEED_WATER     (1,1) {mustBeNumeric}
    tran                  (1,1) struct
    desired_focal_distance_ep (1,1) {mustBeNumeric}
end

    wavelen = SOUND_SPEED_WATER / tran.transducer.freq_hz; % [m]

    R_c        = tran.transducer.annular.curv_radius_mm;
    elem_id_mm = tran.transducer.annular.elem_id_mm(:)';
    elem_od_mm = tran.transducer.annular.elem_od_mm(:)';
    n_pos      = numel(elem_id_mm);

    % Mid-radius of each ring [mm]
    r_mid = (elem_id_mm / 2 + elem_od_mm / 2) / 2;

    % Convert focal distance relative to exit plane to relative to natural focus.
    % Bowl apex is at z = -R_c; target at z = aim_wrt_natural_focus.
    dist_to_mid_bowl      = R_c - tran.transducer.annular.dist_geom_ep_mm;
    focus_wrt_mid_bowl    = desired_focal_distance_ep + dist_to_mid_bowl;
    aim_wrt_natural_focus = focus_wrt_mid_bowl - R_c;

    % Target point [m]
    x_t = 0;
    y_t = 0;
    z_t = aim_wrt_natural_focus / 1000;

    phases = zeros(1, n_pos);
    for i = 1:n_pos
        r = r_mid(i);
        % Element sits on the bowl surface behind the natural focus (z < 0)
        elem_x =  r / 1000;
        elem_y =  0;
        elem_z = -sqrt(R_c^2 - r^2) / 1000;

        dist = sqrt((elem_x - x_t)^2 + (elem_y - y_t)^2 + (elem_z - z_t)^2);
        phases(i) = mod(dist / wavelen, 1) * 360;
    end

end
