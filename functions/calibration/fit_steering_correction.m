function sc = fit_steering_correction(targets_mm, p_model, p_measured, opts)
% FIT_STEERING_CORRECTION  Fit measured / modelled focal peak over a matrix array's steering range
%
% A steered array delivers a focal peak that differs from the model by a
% factor that varies with the steering target (element directivity,
% crosstalk, manufacturing spread). This fits that factor,
%
%   k(x) = p_measured(x) / p_model(x),
%
% from hydrophone measurements at a set of steering targets x, as either
%   'scale'     - one constant (the mean ratio), or
%   'quadratic' - a full quadratic surface in x, y, z (10 coefficients:
%                 1, x, y, z, x^2, y^2, z^2, xy, xz, yz) on coordinates
%                 standardised over the measured targets, by linear least
%                 squares.
% Both are scored by leave-one-out cross-validation: each target is
% predicted from all the others, error = |k_pred p_model - p_measured| / p_measured.
% With model 'auto' (default) the quadratic surface is kept only when there
% are at least min_targets targets and its median cross-validated error is
% at most min_gain times that of the constant; otherwise the constant is
% used. The fitted map is applied by STEERING_CORRECTION_FACTOR.
%
% Use as:
%   sc = fit_steering_correction(targets_mm, p_model, p_measured)
%   sc = fit_steering_correction(targets_mm, p_model, p_measured, opts)
%
% Input:
%   targets_mm - [N x 3] steering targets [mm], transducer frame
%   p_model    - [N x 1] modelled focal peak (same drive as the measurement)
%   p_measured - [N x 1] measured focal peak (same unit as p_model)
%   opts       - optional: .model ('auto' | 'scale' | 'quadratic', default
%                'auto'), .min_targets (default 15), .min_gain (default 0.8)
%
% Output:
%   sc - struct: model, n_targets, scale, coef [10x1] (quadratic), mean_mm,
%        sd_mm [1x3], bounds_mm [2x3] (min; max of the targets),
%        cv (struct: scale_med_pct, scale_p95_pct, quadratic_med_pct,
%        quadratic_p95_pct; NaN when not computed), ratio_range [min max]
%
% See also: STEERING_CORRECTION_FACTOR, BUILD_STEERING_CORRECTION, MATRIX_FOCAL_PEAK

    if nargin < 4; opts = struct(); end
    model = char(get_opt(opts, 'model', 'auto'));
    min_n = get_opt(opts, 'min_targets', 15);
    min_gain = get_opt(opts, 'min_gain', 0.8);

    X = targets_mm; r = p_measured(:) ./ p_model(:);
    assert(size(X, 2) == 3 && size(X, 1) == numel(r), ...
        'fit_steering_correction:size', 'targets_mm must be N x 3 with one row per measured peak.');
    assert(all(isfinite(r)) && all(r > 0), 'fit_steering_correction:ratio', ...
        'Measured and modelled peaks must be positive and finite.');
    N = numel(r);

    sc = struct('model', '', 'n_targets', N, 'scale', mean(r), 'coef', [], ...
        'mean_mm', mean(X, 1), 'sd_mm', std(X, 0, 1), 'bounds_mm', [min(X, [], 1); max(X, [], 1)], ...
        'cv', struct('scale_med_pct', NaN, 'scale_p95_pct', NaN, 'quadratic_med_pct', NaN, 'quadratic_p95_pct', NaN), ...
        'ratio_range', [min(r) max(r)]);
    sc.sd_mm(sc.sd_mm == 0) = 1;   % a coordinate not varied: no standardisation

    % Leave-one-out errors of both models [%]
    e_scale = NaN(N, 1); e_quad = NaN(N, 1);
    A = steering_quad_terms((X - sc.mean_mm) ./ sc.sd_mm);
    do_quad = N > size(A, 2) && rank(A) == size(A, 2);
    for i = 1:N
        o = [1:i-1, i+1:N];
        e_scale(i) = mean(r(o)) / r(i) - 1;
        if do_quad
            e_quad(i) = (A(i, :) * (A(o, :) \ r(o))) / r(i) - 1;
        end
    end
    sc.cv.scale_med_pct = 100 * median(abs(e_scale));
    sc.cv.scale_p95_pct = 100 * prctile(abs(e_scale), 95);
    if do_quad
        sc.cv.quadratic_med_pct = 100 * median(abs(e_quad));
        sc.cv.quadratic_p95_pct = 100 * prctile(abs(e_quad), 95);
    end

    switch model
        case 'scale'
            use_quad = false;
        case 'quadratic'
            if ~do_quad
                error('fit_steering_correction:tooFew', ...
                    'A quadratic surface needs more than 10 targets spread in all three directions (got %d).', N);
            end
            use_quad = true;
        case 'auto'
            use_quad = do_quad && N >= min_n && sc.cv.quadratic_med_pct <= min_gain * sc.cv.scale_med_pct;
        otherwise
            error('fit_steering_correction:model', 'Unknown model ''%s''.', model);
    end
    if use_quad
        sc.model = 'quadratic';
        sc.coef = A \ r;
    else
        sc.model = 'scale';
    end
end

function v = get_opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
