function T = reciprocity_score_poses(cand, P_vol, tr, k_water, k_target, options)
% RECIPROCITY_SCORE_POSES  Score transducer poses on one reverse field (acoustic reciprocity)
%
% P_vol holds the reverse field of a point source at the target
% (REVERSE_POINT_SOURCE_FIELD) on the sensor mask of RECIPROCITY_POSE_CANDIDATES
% (NaN elsewhere). For a pose with radiating-surface points x_j (area a_j,
% ring e(j)) driven with A_e*exp(1i*phi_e), the pressure at the target is,
% up to one constant common to all poses,
%
%   f = sum_j A_e(j) * exp(1i*phi_e(j)) * a_j * P(x_j)
%
% so every pose costs one interpolation of P instead of one forward run.
% Per pose:
%   1. geometric aim: beam axis through the target, rings steered to the
%      target in water (as PRESTUS steers without calibration offsets)
%   2. predicted focal shift of that drive (PREDICT_FOCAL_SHIFT)
%   3. aim correction (options.AimIterations): aim and steer at
%      target - predicted shift (lateral part only for a single element,
%      which cannot steer), tilting the bowl about its apex, and predict
%      again. This is the aberration-corrected aim: where to point and steer
%      so that the focus lands on the target. Poses whose corrected surface
%      leaves the recorded region are flagged (complete = false).
%   4. ring refocusing: target pressure with time-reversal ring phases
%      (what a ring array could add on top; 1 for a single element)
% Pressures are relative to the planned pose (polar angle 0) with geometric
% aim, so p_corrected_rel = 1.3 means 30 % more pressure at the target.
%
% Use as:
%   T = reciprocity_score_poses(cand, P_vol, tr, k_water, k_target)
%   T = reciprocity_score_poses(..., 'AimIterations', 2)
%
% Input:
%   cand     - from RECIPROCITY_POSE_CANDIDATES
%   P_vol    - complex grid volume of the reverse field (NaN outside the sensor mask)
%   tr       - transducer struct (annular); elem_p_amp sets the ring amplitudes
%   k_water  - wavenumber in water [rad/m] (geometric steering)
%   k_target - wavenumber in the tissue at the target [rad/m]
%   options  - 'AimIterations' (default 2), 'PredictShift' (default true)
%
% Output:
%   T - table, one row per pose, sorted by p_corrected_rel (descending)
%
% See also: RECIPROCITY_POSE_CANDIDATES, PREDICT_FOCAL_SHIFT, TIME_REVERSAL_STAGE

arguments
    cand     (1,1) struct
    P_vol
    tr       (1,1) struct
    k_water  (1,1) double
    k_target (1,1) double
    options.AimIterations (1,1) double = 2
    options.PredictShift  (1,1) logical = true
end

    dims = size(P_vol);
    [g1, g2, g3] = deal(1:dims(1), 1:dims(2), 1:dims(3));
    Fr = griddedInterpolant({g1, g2, g3}, real(P_vol), 'cubic', 'none');
    Fi = griddedInterpolant({g1, g2, g3}, imag(P_vol), 'cubic', 'none');
    dx = cand.dx_m;
    sample = @(pts) Fr(pts / dx + 1) + 1i * Fi(pts / dx + 1);

    amp = tr.annular.elem_p_amp(:);
    n_el = numel(tr.annular.elem_od_mm);
    if isscalar(amp), amp = repmat(amp, n_el, 1); end
    t_m = cand.target_m;

    n = numel(cand.poses);
    rows = cell(n, 1);
    for i = 1:n
        ps = cand.poses(i);
        r = struct('pose', i, 'kind', ps.kind, 'polar_deg', ps.polar_deg, 'azimuth_deg', ps.azimuth_deg, ...
            'apex_vox_x', ps.apex_m(1) / dx + 1, 'apex_vox_y', ps.apex_m(2) / dx + 1, ...
            'apex_vox_z', ps.apex_m(3) / dx + 1, 'skin_mm', ps.skin_mm, ...
            'exit_to_target_mm', ps.exit_to_target_mm, 'push_mm', field_or_nan(ps, 'push_mm'), ...
            'inside_fraction', field_or_nan(ps, 'inside_fraction'), 'jitter_mm', field_or_nan(ps, 'jitter_mm'));

        % 1-2. geometric aim at the target
        aim = t_m;
        [f, coh, q, pts, ring_gain] = pose_field(ps.apex_m, aim, tr, amp, k_water, cand.spacing_m, sample);
        r.p_geometric = abs(f);
        r.coherence   = coh;
        r.ring_refocus_gain = ring_gain;
        complete = isfinite(f);
        sh = nan_shift();
        if complete && options.PredictShift
            sh = predict_focal_shift(pts, q, t_m, aim - ps.apex_m, k_target);
        end
        r.shift_axial_mm   = sh.shift_axial_mm;
        r.shift_lateral_mm = sh.shift_lateral_mm;
        r.shift_at_edge    = sh.at_edge;

        % 3. aberration-corrected aim
        f_c = f; sh_c = sh;
        for it = 1:options.AimIterations
            if ~complete || ~options.PredictShift || sh_c.at_edge || ~isfinite(sh_c.shift_mm), break; end
            corr = sh_c.shift_m;
            if isfield(cand, 'fixed_focus') && cand.fixed_focus
                % single element: the focal depth is fixed, only tilt (lateral) helps
                u_ax = (aim - ps.apex_m) / norm(aim - ps.apex_m);
                corr = corr - dot(corr, u_ax) * u_ax;
            end
            aim = aim - corr;
            [f_c, ~, q, pts] = pose_field(ps.apex_m, aim, tr, amp, k_water, cand.spacing_m, sample);
            if ~isfinite(f_c), complete = false; break; end
            sh_c = predict_focal_shift(pts, q, t_m, aim - ps.apex_m, k_target);
        end
        r.p_corrected       = abs(f_c);
        r.residual_shift_mm = sh_c.shift_mm;
        r.aim_vox_x = aim(1) / dx + 1;
        r.aim_vox_y = aim(2) / dx + 1;
        r.aim_vox_z = aim(3) / dx + 1;
        r.aim_offset_mm = norm(aim - t_m) * 1e3;
        r.tilt_deg = acosd(min(1, dot((aim - ps.apex_m) / norm(aim - ps.apex_m), ...
                                      (t_m - ps.apex_m) / norm(t_m - ps.apex_m))));
        r.complete = complete;
        rows{i} = r;
    end
    T = struct2table([rows{:}]);

    ref = T.p_geometric(strcmp(T.kind, 'orbit') & T.polar_deg == 0);
    if isempty(ref) || ~isfinite(ref(1)), ref = max(T.p_geometric); else, ref = ref(1); end
    T.p_geometric_rel = T.p_geometric / ref;
    T.p_corrected_rel = T.p_corrected / ref;
    T = sortrows(T, 'p_corrected_rel', 'descend', 'MissingPlacement', 'last');
end

function [f, coh, q, pts, ring_gain] = pose_field(apex_m, aim_m, tr, amp, k_water, spacing, sample)
% Target pressure of one pose (axis apex -> aim, rings steered to aim in water)
    u = (aim_m - apex_m) / norm(aim_m - apex_m);
    [pts, el, area] = annular_aperture_points(tr.annular, apex_m, u, spacing);
    P = sample(pts);
    n_el = numel(amp);
    r_aim = sqrt(sum((pts - aim_m).^2, 2));
    r_e   = accumarray(el, area .* r_aim, [n_el 1]) ./ accumarray(el, area, [n_el 1]);
    phi   = k_water * r_e;                          % arrive in phase at the aim point (water)
    q     = amp(el) .* exp(1i * phi(el)) .* area .* P;
    f     = sum(q);
    coh   = abs(f) / sum(abs(q));
    G     = accumarray(el, area .* P, [n_el 1]);
    ring_gain = sum(amp .* abs(G)) / abs(f);
end

function sh = nan_shift()
    sh = struct('shift_m', [NaN NaN NaN], 'shift_axial_mm', NaN, 'shift_lateral_mm', NaN, ...
        'shift_mm', NaN, 'peak_rel', NaN, 'at_edge', true, 'f_target', NaN);
end

function v = field_or_nan(s, f)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = NaN; end
end
