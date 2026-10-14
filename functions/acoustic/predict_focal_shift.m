function pred = predict_focal_shift(pts_m, q, target_m, beam_u, k, options)
% PREDICT_FOCAL_SHIFT  Focus near the target predicted from the reverse field over an aperture
%
% Given, for every aperture point j, the drive c_j and the reverse field P_j
% (point source at the target, REVERSE_POINT_SOURCE_FIELD), the forward field
% at the target is sum_j c_j P_j. Near the target it is approximated by
%
%   F(y) = sum_j q_j * exp(-1i*k*(|x_j - y| - |x_j - t|)),   q_j = c_j * P_j
%
% i.e. the skull aberration seen from the target is taken to be the same for
% points y within a few millimetres of it (isoplanatic patch); only the
% geometric path difference in the tissue around the target (wavenumber k)
% changes. The peak of |F| then gives the focal shift for that drive. Beyond
% the isoplanatic patch (a few mm; it shrinks with stronger aberration) the
% prediction degrades, so large predicted shifts are a ranking signal, not
% a position: check them with a forward simulation.
%
% |F| is evaluated on a coarse cube around the target, then on a fine cube
% around the coarse peak. Aperture points are first merged into cells of
% options.CellSize (default lambda/6) to bound the cost.
%
% Use as:
%   pred = predict_focal_shift(pts_m, q, target_m, beam_u, k)
%   pred = predict_focal_shift(..., 'HalfWidth', 0.01, 'Step', 1e-3, 'FineStep', 2.5e-4)
%
% Input:
%   pts_m    - [N x 3] (or N x 2) aperture point positions [m]
%   q        - [N x 1] complex contributions c_j * P_j
%   target_m - [1 x 3] target position [m] (same frame as pts_m)
%   beam_u   - [1 x 3] unit beam direction (transducer -> target)
%   k        - wavenumber in the tissue around the target [rad/m]
%   options  - 'HalfWidth' [m] (default 10 mm), 'Step' [m] (1 mm),
%              'FineStep' [m] (0.25 mm), 'CellSize' [m] (default 2*pi/k/6)
%
% Output:
%   pred - struct:
%     .shift_m        [1 x 3] predicted peak - target [m]
%     .shift_axial_mm component along beam_u (> 0: deeper than the target)
%     .shift_lateral_mm component perpendicular to beam_u
%     .shift_mm       |shift|
%     .peak_rel       |F(peak)| / |F(target)|
%     .at_edge        true when the peak lies on the coarse cube boundary
%                     (shift beyond HalfWidth: prediction unreliable)
%     .f_target       F(target) (= sum(q))
%
% See also: TIME_REVERSAL_REFOCUS, RECIPROCITY_SCORE_POSES

arguments
    pts_m    (:,:) {mustBeNumeric}
    q        (:,1)
    target_m (1,:) {mustBeNumeric}
    beam_u   (1,:) {mustBeNumeric}
    k        (1,1) {mustBeNumeric, mustBePositive}
    options.HalfWidth (1,1) double = 10e-3
    options.Step      (1,1) double = 1e-3
    options.FineStep  (1,1) double = 0.25e-3
    options.CellSize  (1,1) double = NaN
end

    nd = size(pts_m, 2);
    target_m = target_m(1:nd);
    beam_u   = beam_u(1:nd) / norm(beam_u(1:nd));
    keep  = isfinite(q) & all(isfinite(pts_m), 2) & q ~= 0;
    pts_m = pts_m(keep, :);
    q     = q(keep);

    % Merge points into cells (sum of contributions, |q|-weighted centroid)
    cell_size = options.CellSize;
    if ~isfinite(cell_size), cell_size = 2*pi/k/6; end
    [~, ~, cid] = unique(round(pts_m / cell_size), 'rows');
    w     = abs(q);
    q_c   = accumarray(cid, q);
    w_c   = accumarray(cid, w);
    pts_c = zeros(numel(q_c), nd);
    for d = 1:nd
        pts_c(:, d) = accumarray(cid, pts_m(:, d) .* w) ./ max(w_c, eps);
    end
    r_t = sqrt(sum((pts_c - target_m).^2, 2));

    % Coarse cube
    hw = options.HalfWidth;
    ax = -hw:options.Step:hw;
    [Fc, Yc] = eval_cube(pts_c, q_c, r_t, target_m, ax, nd, k);
    [~, ic] = max(abs(Fc));
    peak_c  = Yc(ic, :);
    at_edge = any(abs(peak_c - target_m) >= hw - options.Step/2);

    % Fine cube around the coarse peak
    ax_f = -options.Step:options.FineStep:options.Step;
    [Ff, Yf] = eval_cube(pts_c, q_c, r_t, peak_c, ax_f, nd, k);
    [fmax, jf] = max(abs(Ff));
    peak = Yf(jf, :);

    f_target = sum(q_c);
    shift = peak - target_m;
    ax_comp = dot(shift, beam_u);
    pred = struct( ...
        'shift_m',          shift, ...
        'shift_axial_mm',   ax_comp * 1e3, ...
        'shift_lateral_mm', norm(shift - ax_comp * beam_u) * 1e3, ...
        'shift_mm',         norm(shift) * 1e3, ...
        'peak_rel',         fmax / max(abs(f_target), eps), ...
        'at_edge',          at_edge, ...
        'f_target',         f_target);
end

function [F, Y] = eval_cube(pts, q, r_t, centre, ax, nd, k)
% |F| on a cube (square in 2D) of offsets ax around centre, in chunks
    if nd == 3
        [a, b, c] = ndgrid(ax, ax, ax);
        Y = centre + [a(:), b(:), c(:)];
    else
        [a, b] = ndgrid(ax, ax);
        Y = centre + [a(:), b(:)];
    end
    F = complex(zeros(size(Y, 1), 1));
    chunk = max(1, floor(2e7 / numel(q)));
    for i0 = 1:chunk:size(Y, 1)
        ii = i0:min(i0 + chunk - 1, size(Y, 1));
        d2 = zeros(numel(q), numel(ii));
        for d = 1:nd
            d2 = d2 + (pts(:, d) - Y(ii, d)').^2;
        end
        F(ii) = (q.' * exp(-1i * k * (sqrt(d2) - r_t))).';
    end
end
