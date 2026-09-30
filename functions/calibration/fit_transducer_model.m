function model = fit_transducer_model(tran, prof, options)
% FIT_TRANSDUCER_MODEL  Fit a static annular-array model to hydrophone axial profiles
%
% Fits ONE depth-independent source model to measured on-axis profiles at
% one or more focal depths:
%
%   w_i(d) = g_i * exp(1i * (phi_geo,i(d; c_steer) + phi_i))
%   p_d(z) = H(z + dz, :) * w(d)
%
%   phi_i    per-element phase offset           ('phase')
%   g_i      per-element amplitude weight        ('amp')      hardware gain / sensitivity
%   c_steer  sound speed of the steering law     ('c_steer')  driver's delay calculation
%   dz       axial registration of the profiles  ('z_shift')  hydrophone / exit-plane offset
%
% By default each measured depth gets its own closed-form scale, so only
% profile SHAPE is fitted (hydrophone profiles are recorded at a matched
% focal ISPPA per depth and carry no cross-depth amplitude information).
% The residual is a pseudo-Huber loss on the peak-normalised profile
% (pressure amplitude or intensity, see .Domain) beyond the near field,
% weighted so every depth counts equally — secondary lobes and split foci
% therefore carry weight, while isolated artefacts do not dominate. The
% optimum is taken over several starts (lsqnonlin).
%
% With a single depth and a Prior this is the depth-specific refinement:
% fit phases (and gains), regularised towards the static model.
%
% Advanced options reproduce PERFORM_GLOBAL_SEARCH (the pre-2026-09
% per-depth calibration): phase precession constraints (.PhaseConstraint),
% an absolute-intensity fit with a fitted drive velocity (.AmplitudeFit),
% global optimisers (.Optimizer) and focal weighting / a fit window
% (.FocalWeight, .FitWindowMm). All results are still returned as offsets
% relative to geometric steering, so the library and loader are unchanged.
%
% Use as:
%   model = fit_transducer_model(tran, prof)
%   model = fit_transducer_model(tran, prof, options)
%
% Input:
%   tran    - transducer struct from load_equipment_config
%   prof    - struct: .depths_ep_mm [1 x N_d], .z_ep_mm [N_z x 1],
%             .I [N_z x N_d] measured intensity (any scale; NaN = missing;
%             W/cm² for .AmplitudeFit 'absolute')
%   options - struct, all optional:
%     .FreeParams      - cellstr subset of {'phase','amp','c_steer','z_shift'}
%                        (default {'phase','amp','z_shift'})
%     .ForwardModel    - 'rayleigh' (default) | 'oneil'
%     .NearFieldSkipMm - ignore z < this [mm from EP] (default 10)
%     .HuberDelta      - pseudo-Huber scale on the normalised profile (default 0.05)
%     .Domain          - 'amplitude' (default) | 'intensity': quantity the loss
%                        compares. 'amplitude' weights low-level lobes and tails
%                        more; 'intensity' matches the reported R² and keeps the
%                        tail from outweighing a dominant (e.g. near-field) peak
%     .Lambda          - weight pulling phase offsets to the prior (default 1e-3)
%     .LambdaAmp       - weight pulling log amplitude weights to the prior (default 1e-2)
%     .AmpBounds       - bounds on g_i relative to the prior (default [0.1 10])
%     .PlausibleGain   - gains outside this range (relative to the mean) are
%                        flagged in model.flags as a possible hardware fault
%                        (default [0.5 2])
%     .NearFieldPenalty - weight of a soft penalty on model near-field
%                        intensity exceeding the measured near field
%                        (NearFieldMinMm..NearFieldSkipMm) by more than
%                        NearFieldMargin of the far-field peak; 0 disables
%                        (default 1, margin 0.1, min 2 mm)
%     .ZShiftBounds    - bounds on dz [mm] (default [-3 3])
%     .CSteerBounds    - bounds on c_steer [m/s] (default [1440 1560])
%     .CWater, .Rho    - propagation medium (default 1500 m/s, 994 kg/m³)
%     .NStarts         - number of lsqnonlin starts incl. the prior (default 16)
%     .ZShiftStarts    - 'full' (default): extra starts draw dz stratified over
%                        the whole ZShiftBounds range | 'central': uniform over
%                        the middle half (previous behaviour)
%     .Seed            - RNG seed for the extra starts (default 1)
%     .Prior           - model struct to start from / regularise towards
%                        (default: geometric, unit gains, c_steer = CWater, dz = 0)
%     .ChannelGroups   - [1 x N_elem] channel index per ring; rings in one
%                        group share phase and gain (default: contiguous
%                        groups when tran.n_elem < elem_n, else one per ring)
%     .Display         - lsqnonlin display (default 'off')
%     .EvaluateOnly    - true: no fitting, score the Prior (default false)
%   Advanced (PERFORM_GLOBAL_SEARCH behaviour):
%     .PhaseConstraint - 'none' (default): free per-ring offsets |
%                        'linear': absolute phases psi_g = a + (g-1)*b |
%                        'monotonic': psi_g = a + cumsum(delta), delta within
%                        PhaseStepBounds. Constraints act on the absolute drive
%                        phases at the fitted depth, so they need one depth.
%     .PhaseStepBounds - bounds on b ('linear', default [-2*pi 2*pi]) or on each
%                        delta ('monotonic', default [0 2*pi], as
%                        PERFORM_GLOBAL_SEARCH; with mod 2*pi this does not
%                        restrict the phases — use e.g. [0 pi] to enforce order)
%     .AmplitudeFit    - 'shape' (default): closed-form scale per depth |
%                        'absolute': drive velocity fitted to the measured
%                        W/cm² (loss relative to the largest measured peak)
%     .SharedVelocity  - with 'absolute': one velocity for all depths, i.e. the
%                        profiles were recorded at the same drive (default false)
%     .Optimizer       - 'lsqnonlin' (default, NStarts starts) |
%                        'fexminimize' (external/FEX-minimize multi-start search,
%                        as PERFORM_GLOBAL_SEARCH's default) |
%                        'globalsearch' (MATLAB GlobalSearch + fmincon; needs
%                        the Global Optimization Toolbox)
%     .PopSize, .MaxFunEvals - FEX-minimize settings (default 5000, 2e5)
%     .Polish          - refine a global optimiser's result with lsqnonlin
%                        (default true)
%     .FocalWeight     - 0 (default): uniform weights | w >= 1: Gaussian weights
%                        around the measured FLHM centre, sigma = centre / w
%                        (distance from the bowl apex), as PERFORM_GLOBAL_SEARCH
%     .FitWindowMm     - [z_min z_max] mm from EP compared in the loss and the
%                        metrics, in addition to NearFieldSkipMm (default [])
%
% Output (model):
%   .elem_phase_correction_deg - [1 x N_elem] phase offsets [deg, 0..360)
%   .elem_amp_weight           - [1 x N_elem] gains, mean 1
%   .steering_sound_speed      - c_steer [m/s]
%   .z_shift_mm                - fitted registration offset [mm] (diagnostic;
%                                 not applied at run time)
%   .forward_model, .free_params, .fit_depths_ep_mm
%   .phase_constraint, .phase_params - constraint and its raw parameters [rad]
%   .elem_phase_deg            - [N_d x N_elem] absolute drive phases [deg, 0..360)
%   .amplitude_fit, .optimizer
%   .drive_velocity            - [1 x N_d] mean element velocity [m/s] that gives
%                                 the measured peak intensity at each depth
%                                 ('absolute': the fitted velocity)
%                                 (meaningful when prof.I is in W/cm²)
%   .measured_peak_wcm2        - [1 x N_d] measured peak beyond the near field
%   .metrics                   - table per depth: r2, r2_pressure, focal_err_mm,
%                                 peak_err_mm (near field excluded),
%                                 nearfield_excess = max model − max measured
%                                 near-field intensity, as a fraction of the
%                                 far-field peak (> 0: model hotter than measured),
%                                 and peak_ratio = model / measured peak intensity
%                                 ('absolute' only, else NaN)
%   .flags                     - cellstr: 'gain_implausible', 'gain_at_bound',
%                                 'nearfield_excess' (empty = none)
%   .z_ep_mm, .I_meas, .I_model - normalised profiles for plotting
%   .cost, .start_costs, .exitflag
%
% See also: AXIAL_TRANSFER_MATRIX, REFINE_TRANSDUCER_LIBRARY,
%           LOAD_TRANSDUCER_FROM_LIBRARY, CALIBRATION_GUI, PERFORM_GLOBAL_SEARCH

    if nargin < 3; options = struct(); end
    o.FreeParams      = {'phase', 'amp', 'z_shift'};
    o.ForwardModel    = 'rayleigh';
    o.NearFieldSkipMm = 10;
    o.HuberDelta      = 0.05;
    o.Lambda          = 1e-3;
    o.LambdaAmp       = 1e-2;
    o.AmpBounds       = [0.1 10];
    o.PlausibleGain   = [0.5 2];
    o.NearFieldPenalty = 1;
    o.NearFieldMinMm  = 2;
    o.NearFieldMargin = 0.1;
    o.ZShiftBounds    = [-3 3];
    o.CSteerBounds    = [1440 1560];
    o.CWater          = 1500;
    o.Rho             = 994;
    o.NStarts         = 16;
    o.ZShiftStarts    = 'full';
    o.Seed            = 1;
    o.Prior           = [];
    o.ChannelGroups   = [];
    o.Display         = 'off';
    o.EvaluateOnly    = false;
    o.Domain          = 'amplitude';
    o.PhaseConstraint = 'none';
    o.PhaseStepBounds = [];
    o.AmplitudeFit    = 'shape';
    o.SharedVelocity  = false;
    o.Optimizer       = 'lsqnonlin';
    o.PopSize         = 5000;
    o.MaxFunEvals     = 2e5;
    o.Polish          = true;
    o.FocalWeight     = 0;
    o.FitWindowMm     = [];
    fn = fieldnames(options);
    for k = 1:numel(fn); o.(fn{k}) = options.(fn{k}); end
    free = @(name) any(strcmp(o.FreeParams, name));
    switch o.Domain
        case 'amplitude'; o.pw = 1;
        case 'intensity'; o.pw = 2;
        otherwise; error('fit_transducer_model: unknown Domain ''%s''.', o.Domain);
    end
    o.PhaseConstraint = char(o.PhaseConstraint);
    if ~any(strcmp(o.PhaseConstraint, {'none', 'linear', 'monotonic'}))
        error('fit_transducer_model: unknown PhaseConstraint ''%s''.', o.PhaseConstraint);
    end
    if ~any(strcmp(o.AmplitudeFit, {'shape', 'absolute'}))
        error('fit_transducer_model: unknown AmplitudeFit ''%s''.', o.AmplitudeFit);
    end
    if ~any(strcmp(o.Optimizer, {'lsqnonlin', 'fexminimize', 'globalsearch'}))
        error('fit_transducer_model: unknown Optimizer ''%s''.', o.Optimizer);
    end
    constrained = ~strcmp(o.PhaseConstraint, 'none') && free('phase') && ~o.EvaluateOnly;
    absolute    = strcmp(o.AmplitudeFit, 'absolute');
    if constrained && numel(prof.depths_ep_mm) ~= 1
        error(['fit_transducer_model: PhaseConstraint ''%s'' constrains the absolute phases ' ...
               'at one depth; fit one depth at a time.'], o.PhaseConstraint);
    end
    if isempty(o.PhaseStepBounds)
        o.PhaseStepBounds = ternary(strcmp(o.PhaseConstraint, 'linear'), [-2*pi 2*pi], [0 2*pi]);
    end

    ann    = tran.transducer.annular;
    n_elem = ann.elem_n;
    f_hz   = tran.transducer.freq_hz;
    bowl_offset = ann.curv_radius_mm - ann.dist_geom_ep_mm;

    %% Channel groups -------------------------------------------------------
    grp = o.ChannelGroups;
    if isempty(grp)
        n_ch = n_elem;
        if isfield(tran, 'n_elem') && tran.n_elem < n_elem && mod(n_elem, tran.n_elem) == 0
            n_ch = tran.n_elem;
        end
        grp = kron(1:n_ch, ones(1, n_elem / n_ch));
    end
    grp  = grp(:)';
    n_g  = max(grp);

    %% Prior ----------------------------------------------------------------
    pr = o.Prior;
    if isempty(pr); pr = struct(); end
    pr_phase = get_or(pr, 'elem_phase_correction_deg', zeros(1, n_elem)) * pi/180;
    pr_gain  = get_or(pr, 'elem_amp_weight', ones(1, n_elem));
    pr_c     = get_or(pr, 'steering_sound_speed', o.CWater);
    pr_dz    = get_or(pr, 'z_shift_mm', 0);
    pr_phase_g = group_mean_angle(pr_phase, grp, n_g);
    pr_lgain_g = accumarray(grp', log(pr_gain(:)), [n_g 1], @mean)';

    %% Data -----------------------------------------------------------------
    depths = prof.depths_ep_mm(:)';
    n_d    = numel(depths);
    z_ep   = prof.z_ep_mm(:);
    I      = prof.I;
    I(~isfinite(I) | I < 0) = 0;
    mask   = z_ep >= o.NearFieldSkipMm;
    if ~isempty(o.FitWindowMm)
        mask = mask & z_ep >= o.FitWindowMm(1) & z_ep <= o.FitWindowMm(2);
    end
    if nnz(mask) < 3; error('fit_transducer_model: fewer than 3 profile samples in the fit window.'); end
    far_peak = max(sqrt(I(mask, :)), [], 1);              % per-depth peak amplitude
    % Shape fit: each depth normalised to its own peak. Absolute fit: all depths
    % relative to the largest measured peak, so the fitted drive carries the level.
    nrm = far_peak;
    if absolute; nrm(:) = max(far_peak); end
    M      = (sqrt(I(mask, :)) ./ max(nrm, eps)).^o.pw;    % measured profile, normalised
    zm_bowl = z_ep(mask) + bowl_offset;
    zm_ep   = z_ep(mask);
    nf_mask  = z_ep >= o.NearFieldMinMm & z_ep < o.NearFieldSkipMm;
    NF.zb     = z_ep(nf_mask) + bowl_offset;
    NF.M      = (sqrt(I(nf_mask, :)) ./ max(nrm, eps)).^o.pw;   % same units as M
    NF.lambda = o.NearFieldPenalty * ~isempty(NF.zb);
    NF.margin = o.NearFieldMargin;

    % Residual weights: uniform, or Gaussian around the measured FLHM centre
    % (sigma = centre / FocalWeight from the bowl apex, as PERFORM_GLOBAL_SEARCH),
    % normalised to mean 1 so the loss keeps its scale
    W = ones(size(M));
    if o.FocalWeight >= 1
        for j = 1:n_d
            Ij  = M(:, j).^(2/o.pw);                        % FLHM of the measured intensity
            c_b = compute_flhm_center(Ij / max(Ij), zm_ep) + bowl_offset;
            sg  = c_b / o.FocalWeight;
            W(:, j) = exp(-(zm_bowl - c_b).^2 / (2 * sg^2));
            W(:, j) = W(:, j) / max(mean(W(:, j)), eps);
        end
    end
    o.W = W;

    %% Transfer matrix on a fine grid covering all shifts ---------------------
    dzb = max(abs(o.ZShiftBounds)) * free('z_shift') + abs(pr_dz) + 1;
    zf  = (min([zm_bowl; NF.zb]) - dzb : 0.1 : max(zm_bowl) + dzb)';
    Hf  = axial_transfer_matrix(tran, zf, o.CWater, o.Rho, o.ForwardModel);

    % Unwrapped geometric steering path lengths (same geometry as COMPUTE_PHASES)
    dist_m = steering_path_lengths(tran, depths);      % [n_d x n_elem] m
    geo_g  = @(c) group_mean_angle(2*pi*f_hz/c * dist_m(1, :), grp, n_g);   % absolute geometric phase per group

    % Absolute-fit units: intensity [W/cm²] = |v p|^2 / (2 rho c) * 1e-4
    o.abs_k = 1 / sqrt(2 * o.Rho * o.CWater * 1e4) / max(nrm(1), eps);

    %% Parameter layout -----------------------------------------------------
    idx = struct('phase', [], 'amp', [], 'c', [], 'dz', [], 'vel', []);
    x0 = []; lb = []; ub = [];
    if free('phase')
        idx.phase = numel(x0) + (1:n_g);
        if constrained
            psi0 = unwrap(geo_g(pr_c) + pr_phase_g);            % prior absolute phases
            switch o.PhaseConstraint
                case 'linear'
                    cf  = [ones(n_g, 1), (0:n_g-1)'] \ psi0(:);
                    xp  = [angle(exp(1i*cf(1))), min(max(cf(2), o.PhaseStepBounds(1)), o.PhaseStepBounds(2))];
                    idx.phase = numel(x0) + (1:2);
                    lbp = [-2*pi, o.PhaseStepBounds(1)]; ubp = [2*pi, o.PhaseStepBounds(2)];
                case 'monotonic'
                    dl  = mod(diff(psi0), 2*pi);
                    xp  = [angle(exp(1i*psi0(1))), min(max(dl, o.PhaseStepBounds(1)), o.PhaseStepBounds(2))];
                    lbp = [-2*pi, repmat(o.PhaseStepBounds(1), 1, n_g-1)];
                    ubp = [ 2*pi, repmat(o.PhaseStepBounds(2), 1, n_g-1)];
            end
            x0 = [x0, xp]; lb = [lb, lbp]; ub = [ub, ubp];
        else
            x0 = [x0, pr_phase_g];            lb = [lb, -inf(1, n_g)]; ub = [ub, inf(1, n_g)];
        end
    end
    if free('amp') && n_g > 1
        idx.amp = numel(x0) + (1:n_g-1);  % group 1 is the reference gain
        rel = pr_lgain_g(2:end) - pr_lgain_g(1);
        x0 = [x0, rel];
        lb = [lb, rel + log(o.AmpBounds(1))]; ub = [ub, rel + log(o.AmpBounds(2))];
    end
    if free('c_steer')
        idx.c = numel(x0) + 1;
        x0 = [x0, pr_c]; lb = [lb, o.CSteerBounds(1)]; ub = [ub, o.CSteerBounds(2)];
    end
    if free('z_shift')
        idx.dz = numel(x0) + 1;
        x0 = [x0, pr_dz]; lb = [lb, o.ZShiftBounds(1)]; ub = [ub, o.ZShiftBounds(2)];
    end
    if o.EvaluateOnly
        idx = struct('phase', [], 'amp', [], 'c', [], 'dz', [], 'vel', []);
        x0 = []; lb = []; ub = [];
    elseif isempty(x0)
        error('fit_transducer_model: no free parameters.');
    end
    P = struct('idx', idx, 'grp', grp, 'n_g', n_g, 'pr_phase_g', pr_phase_g, ...
               'pr_lgain_g', pr_lgain_g, 'pr_c', pr_c, 'pr_dz', pr_dz, ...
               'constraint', ternary(constrained, o.PhaseConstraint, 'none'), 'geo_g', geo_g);
    if absolute
        % log drive velocity per depth (or shared), started at the closed-form scale
        P0 = P; P0.idx.vel = [];
        [~, ~, ~, amax0] = residuals(x0, P0, Hf, zf, zm_bowl, M, dist_m, f_hz, o, NF);
        v0 = log(max(far_peak, eps) * sqrt(2 * o.Rho * o.CWater * 1e4) ./ max(amax0, eps));
        if o.SharedVelocity; v0 = log(mean(exp(v0))); end
        P.idx.vel = numel(x0) + (1:numel(v0));
        x0 = [x0, v0]; lb = [lb, v0 - log(100)]; ub = [ub, v0 + log(100)];
    end
    idx = P.idx;

    resid = @(x) residuals(x, P, Hf, zf, zm_bowl, M, dist_m, f_hz, o, NF);

    %% Optimise -------------------------------------------------------------
    if o.EvaluateOnly
        r0 = resid(x0);
        best = struct('x', x0, 'cost', sum(r0.^2), 'flag', NaN);
        start_costs = best.cost;
    else
        lso = optimoptions('lsqnonlin', 'Display', o.Display, 'MaxIterations', 400, ...
            'FunctionTolerance', 1e-10, 'StepTolerance', 1e-10);
        switch o.Optimizer
            case 'lsqnonlin'
                rs  = RandStream('mt19937ar', 'Seed', o.Seed);
                best = struct('cost', inf);
                start_costs = nan(1, o.NStarts);
                for s = 1:o.NStarts
                    xs = x0;
                    if s > 1; xs = random_start(x0, lb, ub, idx, P, o, rs, s); end
                    [x, cost, ~, flag] = lsqnonlin(resid, xs, lb, ub, lso);
                    start_costs(s) = cost;
                    if cost < best.cost
                        best = struct('x', x, 'cost', cost, 'flag', flag);
                    end
                end
            otherwise
                % Global optimisers need finite bounds: free offsets are periodic
                glb = lb; gub = ub;
                if ~constrained && ~isempty(idx.phase)
                    glb(idx.phase) = x0(idx.phase) - 2*pi; gub(idx.phase) = x0(idx.phase) + 2*pi;
                end
                fobj = @(x) sum(resid(x(:)').^2);
                rng(o.Seed, 'twister');
                if strcmp(o.Optimizer, 'fexminimize')
                    if isempty(which('minimize')) || isempty(which('setoptimoptions'))
                        error('fit_transducer_model: FEX-minimize (external/FEX-minimize) is not on the path.');
                    end
                    opts = setoptimoptions('popsize', o.PopSize, 'FinDiffType', 'central', ...
                        'TolCon', 1e-8, 'MaxFunEvals', o.MaxFunEvals);
                    [xg, fg] = minimize(fobj, x0, [], [], [], [], glb, gub, [], opts);
                else
                    if ~license('test', 'GADS_Toolbox') || isempty(which('GlobalSearch'))
                        error('fit_transducer_model: Optimizer ''globalsearch'' needs the Global Optimization Toolbox.');
                    end
                    problem = createOptimProblem('fmincon', 'x0', x0, 'objective', fobj, ...
                        'lb', glb, 'ub', gub, 'options', optimoptions('fmincon', ...
                        'OptimalityTolerance', 1e-8, 'Display', 'off'));
                    [xg, fg] = run(GlobalSearch('Display', 'off'), problem);
                end
                xg = xg(:)'; flag = NaN;
                if o.Polish
                    [xp, fp, ~, flag] = lsqnonlin(resid, min(max(xg, lb), ub), lb, ub, lso);
                    if fp < fg; xg = xp; fg = fp; end
                end
                best = struct('x', xg, 'cost', fg, 'flag', flag);
                start_costs = fg;
        end
    end

    %% Package ----------------------------------------------------------------
    [phase, gain, c_st, dz, psi_g] = unpack_params(best.x, P);
    [~, A, ~, amax, Anf, vel] = residuals(best.x, P, Hf, zf, zm_bowl, M, dist_m, f_hz, o, NF);
    gm = mean(gain);
    model.elem_phase_correction_deg = mod(phase * 180/pi, 360);
    model.elem_amp_weight           = gain / gm;
    model.steering_sound_speed      = c_st;
    model.z_shift_mm                = dz;
    model.forward_model             = o.ForwardModel;
    model.free_params               = o.FreeParams;
    model.fit_depths_ep_mm          = depths;
    model.phase_constraint          = P.constraint;
    model.phase_params              = [];
    if constrained; model.phase_params = best.x(idx.phase); end
    k_st = 2 * pi * f_hz / c_st;
    model.elem_phase_deg            = mod((k_st * dist_m + phase) * 180/pi, 360);
    if constrained; model.elem_phase_deg = mod(psi_g(grp) * 180/pi, 360); end
    model.amplitude_fit             = o.AmplitudeFit;
    model.optimizer                 = o.Optimizer;
    model.cost                      = best.cost;
    model.start_costs               = start_costs;
    model.exitflag                  = best.flag;

    % Mean element velocity [m/s] that reproduces the measured peak intensity
    % at each depth (gains normalised to mean 1): I = |p|^2/(2 rho c) * 1e-4
    I_peak = max(I(mask, :), [], 1);                              % [W/cm²]
    model.drive_velocity = sqrt(I_peak * 1e4 * 2 * o.Rho * o.CWater) .* gm ./ amax;
    peak_ratio = nan(n_d, 1);
    if absolute
        model.drive_velocity = vel .* gm;
        peak_ratio = ((vel .* amax).^2 / (2 * o.Rho * o.CWater) * 1e-4 ./ max(I_peak, eps))';
    end
    model.measured_peak_wcm2 = I_peak;

    % Metrics (intensity domain, peak-normalised, fit window)
    A = A.^(1/o.pw); M = M.^(1/o.pw); Anf = Anf.^(1/o.pw); NF.M = NF.M.^(1/o.pw);   % back to amplitude
    Imod = A.^2; Imeas = M.^2;
    r2 = nan(n_d, 1); r2p = r2; fe = r2; pe = r2;
    for j = 1:n_d
        ym = Imod(:, j) / max(max(Imod(:, j)), eps);
        yt = Imeas(:, j) / max(max(Imeas(:, j)), eps);
        mt = M(:, j) / max(max(M(:, j)), eps);
        r2(j)  = 1 - sum((yt - ym).^2) / max(sum((yt - mean(yt)).^2), eps);
        r2p(j) = 1 - sum((mt - A(:,j)/max(A(:,j))).^2) / max(sum((mt - mean(mt)).^2), eps);
        fe(j)  = compute_flhm_center(ym, zm_ep) - compute_flhm_center(yt, zm_ep);
        [~, a1] = max(ym); [~, a2] = max(yt);
        pe(j)  = zm_ep(a1) - zm_ep(a2);
    end
    nfx = nan(n_d, 1);
    if ~isempty(NF.zb)
        nfx = (max(Anf.^2, [], 1) - max(NF.M.^2, [], 1))' ./ max(max(A.^2, [], 1)', eps);
    end
    model.metrics = table(depths(:), r2, r2p, fe, pe, nfx, peak_ratio, 'VariableNames', ...
        {'depth_ep_mm', 'r2', 'r2_pressure', 'focal_err_mm', 'peak_err_mm', 'nearfield_excess', 'peak_ratio'});

    % Plausibility flags (reported, not enforced)
    g_rel = gain / gm;
    flags = {};
    if any(g_rel < o.PlausibleGain(1) | g_rel > o.PlausibleGain(2)); flags{end+1} = 'gain_implausible'; end
    if ~isempty(idx.amp)
        xa = best.x(idx.amp);
        if any(abs(xa - lb(idx.amp)) < 1e-3 | abs(xa - ub(idx.amp)) < 1e-3); flags{end+1} = 'gain_at_bound'; end
    end
    if any(nfx > NF.margin); flags{end+1} = 'nearfield_excess'; end
    model.flags = flags;
    model.z_ep_mm = zm_ep;
    model.I_meas  = Imeas ./ max(Imeas, [], 1);
    model.I_model = Imod ./ max(Imod, [], 1);
end

%% ========================================================================
function [phase, gain, c_st, dz, psi_g] = unpack_params(x, P)
% phase: per-ring offset relative to geometric steering [rad]
    idx = P.idx;
    lg_g = P.pr_lgain_g - P.pr_lgain_g(1);
    if ~isempty(idx.amp); lg_g = [0, x(idx.amp)]; end
    c_st = P.pr_c;  if ~isempty(idx.c);  c_st = x(idx.c);  end
    dz   = P.pr_dz; if ~isempty(idx.dz); dz   = x(idx.dz); end
    psi_g = [];
    ph_g = P.pr_phase_g;
    if ~isempty(idx.phase)
        xp = x(idx.phase);
        switch P.constraint
            case 'linear';    psi_g = xp(1) + (0:P.n_g-1) * xp(2);
            case 'monotonic'; psi_g = xp(1) + [0, cumsum(xp(2:end))];
            otherwise;        ph_g  = xp;
        end
        if ~isempty(psi_g); ph_g = psi_g - P.geo_g(c_st); end   % absolute → offset
    end
    phase = ph_g(P.grp);
    gain  = exp(lg_g(P.grp));
    if P.n_g < 1; phase = []; end
end

function [r, A, scale, amax, Anf, vel] = residuals(x, P, Hf, zf, zm_bowl, M, dist_m, f_hz, o, NF)
    [phase, gain, c_st, dz] = unpack_params(x, P);
    idx = P.idx;
    n_d = size(M, 2);
    k_st = 2 * pi * f_hz / c_st;
    A = zeros(size(M)); scale = zeros(1, n_d); amax = zeros(1, n_d); vel = nan(1, n_d);
    Anf = zeros(size(NF.M)); r_nf = zeros(numel(NF.M), 1); n_nf = size(NF.M, 1);
    r = zeros(numel(M), 1);
    n_z = size(M, 1);
    for j = 1:n_d
        w  = gain(:) .* exp(1i * (k_st * dist_m(j, :)' + phase(:)));
        pf = Hf * w;
        p  = interp1(zf, pf, zm_bowl + dz, 'linear', 0);
        amax(j) = max(abs(p));
        a  = abs(p).^o.pw;
        if ~isempty(idx.vel)
            vel(j) = exp(x(idx.vel(min(j, numel(idx.vel)))));
            s = (vel(j) * o.abs_k)^o.pw;                  % fitted drive → measured units
        else
            s = (a' * M(:, j)) / max(a' * a, eps);        % closed-form shape scale
        end
        A(:, j) = s * a;
        scale(j) = s;
        e  = (A(:, j) - M(:, j)) / o.HuberDelta;
        rr = o.HuberDelta * sign(e) .* sqrt(2 * (sqrt(1 + e.^2) - 1));   % pseudo-Huber
        r((j-1)*n_z + (1:n_z)) = rr .* sqrt(o.W(:, j)) / sqrt(n_z);
        if n_nf > 0
            Anf(:, j) = s * abs(interp1(zf, pf, NF.zb + dz, 'linear', 0)).^o.pw;
            if NF.lambda > 0
                excess = max(0, Anf(:, j) - (NF.M(:, j) + NF.margin));
                r_nf((j-1)*n_nf + (1:n_nf)) = sqrt(NF.lambda) * excess / sqrt(n_nf);
            end
        end
    end
    if NF.lambda > 0; r = [r; r_nf]; end
    if ~isempty(idx.phase)
        if strcmp(P.constraint, 'none')
            dphi = angle(exp(1i * (x(idx.phase) - P.pr_phase_g)));
        else
            ph_g = accumarray(P.grp', exp(1i * phase(:)), [P.n_g 1])';   % implied offsets per group
            dphi = angle(ph_g .* exp(-1i * P.pr_phase_g));
        end
        r = [r; sqrt(o.Lambda) * 2 * sin(dphi(:) / 2)];
    end
    if ~isempty(idx.amp)
        r = [r; sqrt(o.LambdaAmp) * (x(idx.amp)' - (P.pr_lgain_g(2:end)' - P.pr_lgain_g(1)))];
    end
end

function xs = random_start(x0, lb, ub, idx, P, o, rs, s)
    xs = x0;
    if ~isempty(idx.phase)
        switch P.constraint
            case 'none'; xs(idx.phase) = 2*pi*rand(rs, 1, numel(idx.phase)) - pi;
            otherwise
                ip = idx.phase;
                xs(ip(1)) = 2*pi*rand(rs) - pi;
                xs(ip(2:end)) = lb(ip(2:end)) + (ub(ip(2:end)) - lb(ip(2:end))) .* rand(rs, 1, numel(ip) - 1);
        end
    end
    if ~isempty(idx.amp); xs(idx.amp) = x0(idx.amp) + 0.4*rand(rs, 1, numel(idx.amp)) - 0.2; end
    if ~isempty(idx.c);   xs(idx.c)   = lb(idx.c) + (ub(idx.c) - lb(idx.c)) * rand(rs); end
    % 'full': stratified over the whole range, since split foci /
    % near-field peaks give separate basins in dz, often near the bounds
    if ~isempty(idx.dz)
        if strcmp(o.ZShiftStarts, 'central')
            xs(idx.dz) = lb(idx.dz)/2 + (ub(idx.dz) - lb(idx.dz))/2 * rand(rs);
        else
            xs(idx.dz) = lb(idx.dz) + (ub(idx.dz) - lb(idx.dz)) * (mod(s - 2, o.NStarts - 1) + rand(rs)) / (o.NStarts - 1);
        end
    end
    if ~isempty(idx.vel); xs(idx.vel) = x0(idx.vel) + 0.4*rand(rs, 1, numel(idx.vel)) - 0.2; end
    xs = min(max(xs, lb), ub);
end

function m = group_mean_angle(phi, grp, n_g)
    m = zeros(1, n_g);
    for g = 1:n_g; m(g) = angle(mean(exp(1i * phi(grp == g)))); end
end

function v = get_or(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = double(s.(f)(:)'); else; v = d; end
end

function out = ternary(c, a, b)
    if c; out = a; else; out = b; end
end
