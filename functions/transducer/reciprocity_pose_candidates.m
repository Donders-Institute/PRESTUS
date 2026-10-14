function cand = reciprocity_pose_candidates(parameters, medium_masks, tr, trans_vox, target_vox, options)
% RECIPROCITY_POSE_CANDIDATES  Transducer poses around the planned one, for scoring on one reverse field
%
% Generates poses of a curved annular transducer that all aim at the target:
% the beam direction is tilted away from the planned one (target -> transducer)
% in steps of options.AngleStep up to options.MaxAngle, i.e. the transducer
% moves over the scalp on an orbit around the target. For each direction the
% bowl is placed against the scalp with the planned stand-off: the apex sits
% at (skin distance along that direction) + stand-off + bowl depth from the
% target, so the required focal depth changes with the direction. A single
% element cannot steer: its apex is placed at the curvature radius from the
% target instead, so the geometric focus is on the target. Poses whose
% radiating surface would enter the head or leave the grid are dropped.
%
% The returned sensor mask covers the radiating surfaces of all poses,
% widened by options.MarginMm so that the aim correction in
% RECIPROCITY_SCORE_POSES (small tilts about the apex) stays inside it. One
% REVERSE_POINT_SOURCE_FIELD run on this mask then scores every pose.
%
% Use as:
%   cand = reciprocity_pose_candidates(parameters, medium_masks, tr, trans_vox, target_vox)
%   cand = reciprocity_pose_candidates(..., 'MaxAngle', 20, 'AngleStep', 2.5, 'StandoffMm', [])
%
% Input:
%   parameters   - PRESTUS config (grid.dims, grid.resolution_mm, medium_properties)
%   medium_masks - grid label volume; 0 = water (coupling), > 0 = head
%   tr           - transducer struct (type 'annular', curved)
%   trans_vox    - [1 x 3] planned bowl apex [grid voxels]
%   target_vox   - [1 x 3] target [grid voxels]
%   options      - 'MaxAngle' [deg] (20), 'AngleStep' [deg] (2.5),
%                  'StandoffMm' scalp-to-exit-plane distance [mm]
%                  ([] = that of the planned pose), 'MarginMm' (4),
%                  'SpacingFraction' aperture point spacing in water
%                  wavelengths (0.25), 'JitterMm' (e.g. 0.5): also the
%                  planned pose with its apex moved by +-JitterMm along each
%                  grid axis, still aimed at the target (kind 'jitter')
%
% Output:
%   cand - struct: poses (struct array: kind ('orbit' | 'jitter'), apex_m,
%          axis_u, polar_deg, azimuth_deg, skin_mm, exit_to_target_mm), sensor_mask, dx_m,
%          target_m, spacing_m, standoff_mm, bowl_depth_mm, n_rejected,
%          fixed_focus (single element)
%
% See also: RECIPROCITY_SCORE_POSES, ANNULAR_APERTURE_POINTS, REVERSE_POINT_SOURCE_FIELD

arguments
    parameters   (1,1) struct
    medium_masks {mustBeNumeric}
    tr           (1,1) struct
    trans_vox    (1,3) {mustBeNumeric}
    target_vox   (1,3) {mustBeNumeric}
    options.MaxAngle        (1,1) double = 20
    options.AngleStep       (1,1) double = 2.5
    options.StandoffMm      double = []
    options.MarginMm        (1,1) double = 4
    options.SpacingFraction (1,1) double = 0.25
    options.JitterMm        double = []
end

    if ~strcmp(tr.type, 'annular') || ~isfinite(tr.annular.curv_radius_mm)
        error('PRESTUS:reciprocityPlacement:type', ...
            'Pose search supports curved annular transducers only (got %s).', tr.type);
    end
    dims = size(medium_masks);
    dx   = parameters.grid.resolution_mm * 1e-3;
    to_m = @(v) (v - 1) * dx;
    target_m = to_m(target_vox);
    apex0_m  = to_m(trans_vox);
    d0 = (apex0_m - target_m) / norm(apex0_m - target_m);

    depth_mm = annular_bowl_depth_mm(tr.annular);
    lambda_w = parameters.medium_properties.water.sound_speed / tr.freq_hz;
    spacing  = options.SpacingFraction * lambda_w;

    skin0 = skin_distance(medium_masks, target_m, d0, dx);
    fixed_focus = numel(tr.annular.elem_od_mm) == 1;     % single element: no steering
    if isempty(options.StandoffMm)
        standoff_mm = norm(apex0_m - target_m) * 1e3 - skin0 * 1e3 - depth_mm;
    else
        standoff_mm = options.StandoffMm;
    end

    % Directions: rings of polar angle around d0
    [~, imin] = min(abs(d0));
    e1 = zeros(1, 3); e1(imin) = 1;
    e1 = e1 - dot(e1, d0) * d0; e1 = e1 / norm(e1);
    e2 = cross(d0, e1);
    dirs = d0; pol = 0; azi = 0;
    for b = options.AngleStep:options.AngleStep:options.MaxAngle
        n_az = max(1, ceil(360 * sind(b) / options.AngleStep));
        a = (0:n_az-1)' * 360 / n_az;
        dirs = [dirs; cosd(b) * d0 + sind(b) * (cosd(a) * e1 + sind(a) * e2)]; %#ok<AGROW>
        pol  = [pol; repmat(b, n_az, 1)];                                      %#ok<AGROW>
        azi  = [azi; a];                                                       %#ok<AGROW>
    end

    poses = struct('kind', {}, 'apex_m', {}, 'axis_u', {}, 'polar_deg', {}, 'azimuth_deg', {}, ...
                   'skin_mm', {}, 'exit_to_target_mm', {});
    sensor_mask = false(dims);
    n_rej = 0;
    % Jitter poses: the planned apex moved along each grid axis, aimed at the target
    jit = zeros(0, 3);
    for j = options.JitterMm(:)'
        jit = [jit; j * [eye(3); -eye(3)]]; %#ok<AGROW>
    end
    n_orbit = size(dirs, 1);
    for i = 1:n_orbit + size(jit, 1)
        if i <= n_orbit
            kind = 'orbit';
            d = dirs(i, :) / norm(dirs(i, :));
            skin = skin_distance(medium_masks, target_m, d, dx);
            if ~isfinite(skin), n_rej = n_rej + 1; continue; end
        else
            kind = 'jitter';
            apex_j = apex0_m + jit(i - n_orbit, :) * 1e-3;
            d = (apex_j - target_m) / norm(apex_j - target_m);
            skin = skin_distance(medium_masks, target_m, d, dx);
        end
        if strcmp(kind, 'jitter')
            apex_m = apex_j;
        elseif fixed_focus
            apex_m = target_m + tr.annular.curv_radius_mm * 1e-3 * d;
        else
            apex_m = target_m + (skin + (standoff_mm + depth_mm) * 1e-3) * d;
        end
        pts = annular_aperture_points(tr.annular, apex_m, -d, spacing);
        vox = round(pts / dx) + 1;
        if any(vox(:) < 2) || any(vox(:, 1) > dims(1) - 1) || any(vox(:, 2) > dims(2) - 1) || any(vox(:, 3) > dims(3) - 1)
            n_rej = n_rej + 1; continue;
        end
        lin = sub2ind(dims, vox(:, 1), vox(:, 2), vox(:, 3));
        if any(medium_masks(lin) > 0)
            n_rej = n_rej + 1; continue;            % surface would cut into the head
        end
        sensor_mask(lin) = true;
        if i <= n_orbit, pa = [pol(i), azi(i)]; else, pa = [acosd(min(1, dot(d, d0))), NaN]; end
        poses(end+1) = struct('kind', kind, 'apex_m', apex_m, 'axis_u', -d, 'polar_deg', pa(1), ...
            'azimuth_deg', pa(2), 'skin_mm', skin * 1e3, ...
            'exit_to_target_mm', norm(apex_m - target_m) * 1e3 - depth_mm); %#ok<AGROW>
    end
    if isempty(poses)
        error('PRESTUS:reciprocityPlacement:noPoses', 'No candidate pose fits in the grid outside the head.');
    end

    % Widen to cover small re-aiming tilts and the interpolation stencil
    sensor_mask = bwdist(sensor_mask) <= options.MarginMm / parameters.grid.resolution_mm + 1;

    cand = struct('poses', poses, 'sensor_mask', sensor_mask, 'dx_m', dx, 'fixed_focus', fixed_focus, ...
        'target_m', target_m, 'spacing_m', spacing, 'standoff_mm', standoff_mm, ...
        'bowl_depth_mm', depth_mm, 'n_rejected', n_rej);
    fprintf('[time_reversal] pose search: %d poses (%d rejected), stand-off %.1f mm, %d sensor points\n', ...
        numel(poses), n_rej, standoff_mm, nnz(sensor_mask));
end

function s = skin_distance(medium_masks, target_m, d, dx)
% Distance from the target to the last head voxel along d [m]; Inf if the
% ray leaves the grid inside the head
    dims = size(medium_masks);
    step = dx / 2;
    s_last = 0; s = 0;
    while true
        v = round((target_m + s * d) / dx) + 1;
        if any(v < 1) || any(v > dims)
            if medium_masks(clamp_idx(target_m + (s - step) * d, dx, dims)) > 0
                s = Inf; return
            end
            break
        end
        if medium_masks(v(1), v(2), v(3)) > 0
            s_last = s;
        end
        s = s + step;
    end
    s = s_last + dx / 2;
end

function lin = clamp_idx(p, dx, dims)
    v = min(max(round(p / dx) + 1, 1), dims);
    lin = sub2ind(dims, v(1), v(2), v(3));
end
