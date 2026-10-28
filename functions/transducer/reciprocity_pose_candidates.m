function cand = reciprocity_pose_candidates(parameters, medium_masks, tr, trans_vox, target_vox, options)
% RECIPROCITY_POSE_CANDIDATES  Transducer poses around the planned one, for scoring on one reverse field
%
% Generates poses of a curved annular transducer that all aim at the target:
% the beam direction is tilted away from the planned one (target -> transducer)
% in steps of options.AngleStep up to options.MaxAngle, i.e. the transducer
% moves over the scalp on an orbit around the target. For each direction the
% bowl is placed against the scalp with the planned stand-off: the apex sits
% at (skin distance along that direction) + stand-off + bowl depth from the
% target, so the required focal depth changes with the direction; where the
% scalp curves, the bowl is then moved outward until its radiating surface
% overlaps the head no more than the planned pose's surface does, plus
% options.MaxInsideFraction (at most options.MaxPushMm), so that all poses
% are seated alike. A single element cannot steer:
% its apex is placed at the curvature radius from the target instead, so the
% geometric focus is on the target. Poses whose radiating surface still
% overlaps the head by more than that, or leaves the grid, are dropped; the
% planned pose is always kept.
%
% With options.AllowOutside, poses that leave the grid are kept and
% cand.pad_lo / cand.pad_hi give the voxels the grid would need on each side
% (TIME_REVERSAL_STAGE pads the reverse-run grid with water by that much and
% calls this function again on the padded grid).
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
%   medium_masks - grid label volume; 0 and the water label (MEDIUM_SETUP
%                  relabels 0 to the index of 'water' in medium_properties)
%                  are coupling medium, every other label is head
%   tr           - transducer struct (type 'annular', curved)
%   trans_vox    - [1 x 3] planned bowl apex [grid voxels]
%   target_vox   - [1 x 3] target [grid voxels]
%   options      - 'MaxAngle' [deg] (20), 'AngleStep' [deg] (2.5),
%                  'StandoffMm' scalp-to-exit-plane distance [mm]
%                  ([] = that of the planned pose), 'MarginMm' (4),
%                  'SpacingFraction' aperture point spacing in water
%                  wavelengths (0.25), 'JitterMm' (one or more distances,
%                  e.g. [0.5 1 2]): also the planned pose with its apex
%                  moved by each distance, still aimed at the target (kind
%                  'jitter'), 'JitterMode' 'grid' (default: +- along each
%                  grid axis, as when editing voxel coordinates; 6 poses per
%                  distance) or 'scalp' (slid over the scalp perpendicular to
%                  the planned beam in 8 directions and seated against the
%                  scalp like an orbit pose; 8 poses per distance),
%                  'MaxInsideFraction' (0.02), 'MaxPushMm' (30),
%                  'AllowOutside' (false)
%
% Output:
%   cand - struct: poses (struct array: kind ('orbit' | 'jitter'), apex_m,
%          axis_u, polar_deg, azimuth_deg, skin_mm, exit_to_target_mm), sensor_mask, dx_m,
%          target_m, spacing_m, standoff_mm, bowl_depth_mm, n_rejected,
%          fixed_focus (single element), pad_lo, pad_hi ([1 x 3] voxels the
%          grid lacks below 1 / above dims; zeros unless AllowOutside).
%          Each pose also carries push_mm (outward move to clear the head)
%          and inside_fraction (share of its surface points in the head).
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
    options.JitterMode      (1,:) char {mustBeMember(options.JitterMode, {'grid', 'scalp'})} = 'grid'
    options.MaxInsideFraction (1,1) double = 0.02
    options.MaxPushMm       (1,1) double = 30
    options.AllowOutside    (1,1) logical = false
end

    if ~strcmp(tr.type, 'annular') || ~isfinite(tr.annular.curv_radius_mm)
        error('PRESTUS:reciprocityPlacement:type', ...
            'Pose search supports curved annular transducers only (got %s).', tr.type);
    end
    dims = size(medium_masks);
    dx   = parameters.grid.resolution_mm * 1e-3;
    head = head_mask(parameters, medium_masks);
    to_m = @(v) (v - 1) * dx;
    target_m = to_m(target_vox);
    apex0_m  = to_m(trans_vox);
    d0 = (apex0_m - target_m) / norm(apex0_m - target_m);

    depth_mm = annular_bowl_depth_mm(tr.annular);
    lambda_w = parameters.medium_properties.water.sound_speed / tr.freq_hz;
    spacing  = options.SpacingFraction * lambda_w;

    skin0 = skin_distance(head, target_m, d0, dx);
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
                   'skin_mm', {}, 'exit_to_target_mm', {}, 'push_mm', {}, 'inside_fraction', {}, ...
                   'jitter_mm', {});
    all_vox = cell(0, 1);
    n_rej = 0;
    margin_vox = ceil(options.MarginMm / parameters.grid.resolution_mm) + 2;
    planned_frac = 0;
    % Jitter poses: the planned apex moved by each distance, aimed at the target
    scalp_jitter = strcmp(options.JitterMode, 'scalp');
    jit = zeros(0, 3); jit_d = zeros(0, 1); jit_az = zeros(0, 1);
    for j = options.JitterMm(:)'
        if scalp_jitter
            a = (0:45:315)';
            jit = [jit; j * (cosd(a) * e1 + sind(a) * e2)]; %#ok<AGROW>
            jit_az = [jit_az; a];                                %#ok<AGROW>
        else
            jit = [jit; j * [eye(3); -eye(3)]];               %#ok<AGROW>
            jit_az = [jit_az; nan(6, 1)];                       %#ok<AGROW>
        end
        jit_d = [jit_d; repmat(j, size(jit, 1) - numel(jit_d), 1)]; %#ok<AGROW>
    end
    n_orbit = size(dirs, 1);
    for i = 1:n_orbit + size(jit, 1)
        if i <= n_orbit
            kind = 'orbit';
            d = dirs(i, :) / norm(dirs(i, :));
            skin = skin_distance(head, target_m, d, dx);
            if ~isfinite(skin), n_rej = n_rej + 1; continue; end
        else
            kind = 'jitter';
            apex_j = apex0_m + jit(i - n_orbit, :) * 1e-3;
            d = (apex_j - target_m) / norm(apex_j - target_m);
            skin = skin_distance(head, target_m, d, dx);
            if scalp_jitter && ~isfinite(skin), n_rej = n_rej + 1; continue; end
        end
        if strcmp(kind, 'jitter') && ~scalp_jitter
            apex_m = apex_j;                        % moved as given, contact not re-seated
        elseif fixed_focus
            apex_m = target_m + tr.annular.curv_radius_mm * 1e-3 * d;
        else
            apex_m = target_m + (skin + (standoff_mm + depth_mm) * 1e-3) * d;
        end
        planned = i == 1;
        % Every pose is seated as the planned one: its surface may overlap the
        % head as much as the planned surface does (+ MaxInsideFraction). A
        % steerable bowl is moved outward until it does no more than that.
        limit = options.MaxInsideFraction + planned_frac;
        push = 0;
        movable = (strcmp(kind, 'orbit') || scalp_jitter) && ~fixed_focus && ~planned;
        while true
            pts = annular_aperture_points(tr.annular, apex_m + push * d, -d, spacing);
            vox = round(pts / dx) + 1;
            frac = mean(in_head(head, vox));
            if frac <= limit || ~movable || push * 1e3 >= options.MaxPushMm, break; end
            push = push + dx;
        end
        apex_m = apex_m + push * d;
        outside = any(vox < 2, 2) | vox(:, 1) > dims(1) - 1 | vox(:, 2) > dims(2) - 1 | vox(:, 3) > dims(3) - 1;
        if planned, planned_frac = frac; end
        if ~planned && frac > limit
            n_rej = n_rej + 1; continue;            % surface would cut into the head
        end
        if any(outside) && ~options.AllowOutside
            n_rej = n_rej + 1; continue;
        end
        if planned && frac > 0
            fprintf('[time_reversal] pose search: %.1f%% of the planned surface lies in the head\n', 100 * frac);
        end
        all_vox{end+1, 1} = vox; %#ok<AGROW>
        if i <= n_orbit
            pa = [pol(i), azi(i)]; jd = 0;
        else
            pa = [acosd(min(1, dot(d, d0))), jit_az(i - n_orbit)]; jd = jit_d(i - n_orbit);
        end
        poses(end+1) = struct('kind', kind, 'apex_m', apex_m, 'axis_u', -d, 'polar_deg', pa(1), ...
            'azimuth_deg', pa(2), 'skin_mm', skin * 1e3, ...
            'exit_to_target_mm', norm(apex_m - target_m) * 1e3 - depth_mm, ...
            'push_mm', push * 1e3, 'inside_fraction', frac, 'jitter_mm', jd); %#ok<AGROW>
    end
    if isempty(poses)
        error('PRESTUS:reciprocityPlacement:noPoses', 'No candidate pose fits in the grid outside the head.');
    end

    % Voxels the grid lacks (with the margin), and the sensor mask inside it
    V = vertcat(all_vox{:});
    pad_lo = max(0, 1 + margin_vox - min(V, [], 1));
    pad_hi = max(0, max(V, [], 1) + margin_vox - dims);
    if ~options.AllowOutside, pad_lo = zeros(1, 3); pad_hi = zeros(1, 3); end
    sensor_mask = false(dims);
    keep = all(V >= 1 & V <= dims, 2);
    sensor_mask(sub2ind(dims, V(keep, 1), V(keep, 2), V(keep, 3))) = true;
    % Widen to cover small re-aiming tilts and the interpolation stencil
    sensor_mask = bwdist(sensor_mask) <= options.MarginMm / parameters.grid.resolution_mm + 1;

    cand = struct('poses', poses, 'sensor_mask', sensor_mask, 'dx_m', dx, 'fixed_focus', fixed_focus, ...
        'target_m', target_m, 'spacing_m', spacing, 'standoff_mm', standoff_mm, ...
        'bowl_depth_mm', depth_mm, 'n_rejected', n_rej, 'pad_lo', pad_lo, 'pad_hi', pad_hi);
    fprintf(['[time_reversal] pose search: %d poses (%d rejected), stand-off %.1f mm, %d sensor points, ' ...
        'median push %.1f mm%s\n'], numel(poses), n_rej, standoff_mm, nnz(sensor_mask), median([poses.push_mm]), ...
        ternary(any(pad_lo + pad_hi > 0), sprintf(', grid lacks [%s] / [%s] voxels', num2str(pad_lo), num2str(pad_hi)), ''));
end

function h = in_head(head, vox)
% Head membership of voxel subscripts; outside the grid counts as coupling medium
    dims = size(head);
    h = false(size(vox, 1), 1);
    ok = all(vox >= 1 & vox <= dims, 2);
    h(ok) = head(sub2ind(dims, vox(ok, 1), vox(ok, 2), vox(ok, 3)));
end

function out = ternary(c, a, b)
    if c, out = a; else, out = b; end
end

function head = head_mask(parameters, medium_masks)
% Head voxels: labels other than 0 and the water label of MEDIUM_SETUP
    head = medium_masks > 0;
    if isfield(parameters, 'medium_properties') && isstruct(parameters.medium_properties)
        w = find(strcmp(fieldnames(parameters.medium_properties), 'water'), 1);
        if ~isempty(w), head = head & medium_masks ~= w; end
    end
end

function s = skin_distance(head, target_m, d, dx)
% Distance from the target to the last head voxel along d [m]; Inf if the
% ray leaves the grid inside the head
    dims = size(head);
    step = dx / 2;
    s_last = 0; s = 0;
    while true
        v = round((target_m + s * d) / dx) + 1;
        if any(v < 1) || any(v > dims)
            if head(clamp_idx(target_m + (s - step) * d, dx, dims))
                s = Inf; return
            end
            break
        end
        if head(v(1), v(2), v(3))
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
