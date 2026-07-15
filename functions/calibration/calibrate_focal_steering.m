function [D_star_ep_mm, diag] = calibrate_focal_steering(tr, target_ep_mm, parameters, opts)
% CALIBRATE_FOCAL_STEERING  Analytical steer that lands the FWHM/peak at a target
%
% PRESTUS specifies a focal depth (focal_distance_ep) as the FWHM/FLHM-centre of
% the achieved field.  Steering, however, is geometric: focusing geometrically
% at depth D produces an intensity peak proximal to D (diffraction focal shift)
% and an FLHM-centre near D.  To place the achieved FLHM-centre (or peak) exactly
% at a requested depth this routine root-finds the geometric steer D* whose
% analytically predicted FLHM-centre (or peak) equals the target.
%
% The forward model is the in-repo Rayleigh–Sommerfeld solver
% (RAYLEIGH_AXIAL_INTENSITY); geometric phases for a candidate steer come from
% SET_REAL_PHASES / COMPUTE_PHASES.  focusedAnnulusONeil is intentionally NOT
% used (it is resolved from the MATLAB path at runtime and is not guaranteed to
% be present).
%
% Use as:
%   [D_star_ep_mm, diag] = calibrate_focal_steering(tr, target_ep_mm, parameters)
%   [D_star_ep_mm, diag] = calibrate_focal_steering(tr, target_ep_mm, parameters, opts)
%
% Input:
%   tr            - single transducer struct (parameters.transducer(1) style):
%                   uses tr.freq_hz and tr.annular.{elem_n, elem_id_mm,
%                   elem_od_mm, curv_radius_mm, dist_geom_ep_mm, elem_p_amp}
%   target_ep_mm  - desired FWHM-centre (or peak) depth from the exit plane [mm]
%   parameters    - PRESTUS config; uses medium_properties.water.{sound_speed,density}
%   opts.target   - 'flhm_center' (default) | 'peak'
%   opts.phases_rad - optional [1xN] fixed drive phases.  When supplied, the
%                   routine does NOT re-steer: it evaluates the achieved focus of
%                   these phases once and returns it (locate-only mode, for
%                   measured/library phases that must not be re-derived).
%   opts.eval_max_mm - far end of the axial evaluation window [mm from bowl]
%                   (default: 1.6 * (target + bowl depth), min 130).
%
% Output:
%   D_star_ep_mm  - geometric steer depth [mm EP] that places the achieved
%                   FLHM-centre (or peak) at target_ep_mm.  In locate-only mode
%                   this is the achieved focus of the supplied phases.
%   diag          - struct: achieved_flhm_ep_mm, achieved_peak_ep_mm,
%                   flhm_span_mm, flhm_truncated (logical), iters, mode
%
% See also: SET_REAL_PHASES, COMPUTE_PHASES, RAYLEIGH_AXIAL_INTENSITY,
%           COMPUTE_FOCAL_METRICS

arguments
    tr            (1,1) struct
    target_ep_mm  (1,1) {mustBeNumeric}
    parameters    (1,1) struct
    opts          (1,1) struct = struct()
end

    if ~isfield(opts, 'target'); opts.target = 'flhm_center'; end
    ann = tr.annular;

    % Bowl <-> exit-plane conversion (bowl surface / apex = 0):
    %   ep = bowl + (dist_geom_ep_mm - curv_radius_mm)   [ep_offset < 0]
    ep_offset_mm = ann.dist_geom_ep_mm - ann.curv_radius_mm;

    % Axial evaluation grid (from bowl apex), fine enough to resolve the −3 dB span.
    if isfield(opts, 'eval_max_mm') && ~isempty(opts.eval_max_mm)
        eval_max_mm = opts.eval_max_mm;
    else
        eval_max_mm = max(130, 1.6 * (target_ep_mm - ep_offset_mm));
    end

    % Context bundle shared by the forward-model evaluations (local functions).
    ctx = struct();
    ctx.dist_bowl_mm  = (0.5:0.25:eval_max_mm)';
    ctx.ep_offset_mm  = ep_offset_mm;
    ctx.target        = opts.target;
    rho = parameters.medium_properties.water.density;
    c   = parameters.medium_properties.water.sound_speed;
    ctx.velocity = ann.elem_p_amp(1) / (rho * c);      % particle velocity [m/s]
    ctx.p_ray = struct();
    ctx.p_ray.transducer.annular  = ann;
    ctx.p_ray.transducer.freq_hz  = tr.freq_hz;
    ctx.p_ray.medium_properties.water = parameters.medium_properties.water;
    % Adapter for set_real_phases / compute_phases.
    ctx.tran_ph = struct('n_elem', ann.elem_n);
    ctx.tran_ph.transducer.freq_hz = tr.freq_hz;
    ctx.tran_ph.transducer.annular = ann;
    ctx.parameters = parameters;

    % ---- Locate-only mode: fixed phases supplied, do not re-steer ----------
    if isfield(opts, 'phases_rad') && ~isempty(opts.phases_rad)
        [flhm_ep, peak_ep, span_mm, trunc] = achieved_focus(opts.phases_rad(:)', ctx);
        D_star_ep_mm = pick_metric(flhm_ep, peak_ep, opts.target);
        diag = pack_diag(flhm_ep, peak_ep, span_mm, trunc, 0, 'locate');
        return
    end

    % ---- Re-steer mode: root-find D* so achieved focus == target -----------
    g = @(D) focus_metric_for_steer(D, ctx) - target_ep_mm;   % monotonic increasing

    lo = max(1, target_ep_mm - 20);
    hi = target_ep_mm + 30;
    glo = g(lo); ghi = g(hi);
    iters = 2;

    if isnan(glo) || isnan(ghi) || glo * ghi > 0
        warning('calibrate_focal_steering:noBracket', ...
            ['Could not bracket the steer for target %.1f mm EP ' ...
             '(g(%.1f)=%.2f, g(%.1f)=%.2f). Falling back to D*=target.'], ...
            target_ep_mm, lo, glo, hi, ghi);
        D_star_ep_mm = target_ep_mm;
        ph = deg2rad(set_real_phases(ctx.tran_ph, D_star_ep_mm, parameters));
        [flhm_ep, peak_ep, span_mm, trunc] = achieved_focus(ph, ctx);
        diag = pack_diag(flhm_ep, peak_ep, span_mm, trunc, iters, 'naive_fallback');
        return
    end

    optim = optimset('TolX', 0.05);   % 0.05 mm steer tolerance
    [D_star_ep_mm, ~, ~, out] = fzero(g, [lo hi], optim);
    iters = iters + out.iterations;

    ph = deg2rad(set_real_phases(ctx.tran_ph, D_star_ep_mm, parameters));
    [flhm_ep, peak_ep, span_mm, trunc] = achieved_focus(ph, ctx);
    diag = pack_diag(flhm_ep, peak_ep, span_mm, trunc, iters, 're-steer');
    if trunc
        warning('calibrate_focal_steering:flhmTruncated', ...
            ['−3 dB focal lobe reached the evaluation window edge; ' ...
             'D* may be unreliable (increase opts.eval_max_mm).']);
    end
end

% -------------------------------------------------------------------------
function m = focus_metric_for_steer(D, ctx)
% Achieved focus metric (EP) for a geometric steer to depth D (EP).
    ph = deg2rad(set_real_phases(ctx.tran_ph, D, ctx.parameters));
    [flhm_ep, peak_ep] = achieved_focus(ph, ctx);
    m = pick_metric(flhm_ep, peak_ep, ctx.target);
end

% -------------------------------------------------------------------------
function [flhm_ep, peak_ep, span_mm, trunc] = achieved_focus(phases_rad, ctx)
    I = rayleigh_axial_intensity(phases_rad, ctx.velocity, ctx.p_ray, ctx.dist_bowl_mm);
    [flhm_bowl, peak_bowl, span_mm, trunc] = flhm_from_profile(I, ctx.dist_bowl_mm);
    flhm_ep = flhm_bowl + ctx.ep_offset_mm;
    peak_ep = peak_bowl + ctx.ep_offset_mm;
end

% -------------------------------------------------------------------------
function m = pick_metric(flhm_ep, peak_ep, target)
    if strcmpi(target, 'peak'); m = peak_ep; else; m = flhm_ep; end
end

% -------------------------------------------------------------------------
function d = pack_diag(flhm_ep, peak_ep, span_mm, trunc, iters, mode)
    d = struct('achieved_flhm_ep_mm', flhm_ep, 'achieved_peak_ep_mm', peak_ep, ...
               'flhm_span_mm', span_mm, 'flhm_truncated', trunc, ...
               'iters', iters, 'mode', mode);
end

% -------------------------------------------------------------------------
function [flhm_center, peak_pos, span, trunc] = flhm_from_profile(I, axis_pos)
% Contiguous −3 dB (half-max) span bracketing the peak; midpoint = FLHM centre.
% Matches the logic in compute_focal_metrics.m.
    I = I(:); axis_pos = axis_pos(:);
    [pk, k] = max(I);
    thr = pk / 2;
    lo = k; while lo > 1        && I(lo-1) >= thr; lo = lo - 1; end
    hi = k; while hi < numel(I) && I(hi+1) >= thr; hi = hi + 1; end
    flhm_center = 0.5 * (axis_pos(lo) + axis_pos(hi));
    peak_pos    = axis_pos(k);
    span        = axis_pos(hi) - axis_pos(lo);
    trunc       = (lo == 1) || (hi == numel(I));
end
