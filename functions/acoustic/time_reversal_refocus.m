function res = time_reversal_refocus(P, W, C_elem, pts_m, target_m, beam_u, k_target, tr_type, options)
% TIME_REVERSAL_REFOCUS  Element transfers, refocusing phases and aberration metrics of one transducer
%
% From the reverse field P at the source points (REVERSE_POINT_SOURCE_FIELD)
% and the grid weight of every element at those points (SOURCE_CREATE), the
% transfer of element e to the target is, by reciprocity,
%
%   G_e = sum_i W(i,e) * P_i
%
% and the field at the target for element drives C_e = A_e*exp(1i*phi_e) is
% sum_e C_e*G_e (up to one constant common to all elements and drives).
%
% Time reversal (phase conjugation) sets phi_e = -angle(G_e) (plus one
% common offset, chosen here to keep the phases close to the current ones),
% which makes every element arrive at the target in phase: the largest
% target pressure reachable with the given element amplitudes. What a
% transducer can use depends on how its phases are driven:
%   matrix                    per-element phases: full refocusing
%   annular, elem_n > 1       per-ring phases: the refocusing of a ring
%                             array (axially symmetric part of the aberration)
%   annular, elem_n == 1      one phase: no refocusing; the metrics still
%                             say how much of the aperture adds up in phase
%                             and where the focus is predicted to land
%
% Metrics (all from the same reverse run; no forward simulation):
%   coherence_elements  |sum_e C_e G_e| / sum_e |C_e G_e|; 1 = elements in phase
%   coherence_aperture  the same over all source points (cells of ~lambda/6),
%                       i.e. including phase spread within elements, which
%                       no element phase can correct
%   refocus_gain        target pressure with time-reversal phases / with the
%                       current phases (= 1/coherence_elements); intensity
%                       gain is its square
%   predicted focal shift and relative peak for the current and the
%                       time-reversal drive (PREDICT_FOCAL_SHIFT)
%
% Use as:
%   res = time_reversal_refocus(P, W, C_elem, pts_m, target_m, beam_u, k_target, tr_type)
%   res = time_reversal_refocus(..., 'PredictShift', false)
%
% Input:
%   P        - [num_src x 1] complex reverse field at the source points
%   W        - [num_src x n_elem] element weights (columns of this transducer)
%   C_elem   - [n_elem x 1] current complex drive A_e*exp(1i*phi_e)
%   pts_m    - [num_src x n_dims] source point positions [m]
%   target_m - [1 x n_dims] target position [m]
%   beam_u   - [1 x n_dims] unit beam direction (transducer -> target)
%   k_target - wavenumber in the tissue at the target [rad/m]
%   tr_type  - 'annular' | 'matrix'
%   options  - 'PredictShift' (default true): run PREDICT_FOCAL_SHIFT
%
% Output:
%   res - struct: capability ('full' | 'rings' | 'none'), G, phase_tr_rad,
%         phase_error_deg (per element, current drive, relative to the
%         amplitude-weighted mean arrival phase), coherence_elements,
%         coherence_aperture, refocus_gain, pred_current, pred_refocused
%
% See also: REVERSE_POINT_SOURCE_FIELD, PREDICT_FOCAL_SHIFT, TIME_REVERSAL_STAGE

arguments
    P        (:,1)
    W        (:,:)
    C_elem   (:,1)
    pts_m    (:,:) {mustBeNumeric}
    target_m (1,:) {mustBeNumeric}
    beam_u   (1,:) {mustBeNumeric}
    k_target (1,1) {mustBeNumeric, mustBePositive}
    tr_type  {mustBeTextScalar}
    options.PredictShift (1,1) logical = true
end

    n_elem = numel(C_elem);
    if size(W, 2) ~= n_elem
        error('PRESTUS:timeReversal:weights', 'W has %d columns for %d elements.', size(W, 2), n_elem);
    end
    W = double(W);
    P = double(P);

    % Only the points of this transducer
    rows  = find(any(W ~= 0, 2));
    W     = W(rows, :);
    P     = P(rows);
    pts_m = pts_m(rows, :);

    G = full(W.' * P);                                   % element -> target transfer

    if strcmp(tr_type, 'matrix')
        capability = 'full';
    elseif n_elem > 1
        capability = 'rings';
    else
        capability = 'none';
    end

    arrive   = C_elem .* G;                              % per-element contribution at the target
    f_cur    = sum(arrive);
    coh_elem = abs(f_cur) / max(sum(abs(arrive)), eps);

    % Time-reversal phases, offset to the amplitude-weighted mean of the current ones
    phase_tr = -angle(G);
    offset   = angle(sum(abs(arrive) .* exp(1i * (angle(C_elem) - phase_tr))));
    phase_tr = mod(phase_tr + offset, 2*pi);
    C_tr     = abs(C_elem) .* exp(1i * phase_tr);
    f_tr     = sum(C_tr .* G);

    mean_arrival    = angle(f_cur);
    phase_error_deg = rad2deg(angle(arrive .* exp(-1i * mean_arrival)));

    % Per-point contributions (current and time-reversal drive)
    q_cur = full(W * C_elem) .* P;
    q_tr  = full(W * C_tr)   .* P;
    coh_ap = point_coherence(pts_m, q_cur, k_target);

    res = struct( ...
        'capability',         capability, ...
        'G',                  G, ...
        'phase_current_rad',  angle(C_elem), ...
        'phase_tr_rad',       phase_tr, ...
        'phase_error_deg',    phase_error_deg, ...
        'coherence_elements', coh_elem, ...
        'coherence_aperture', coh_ap, ...
        'refocus_gain',       abs(f_tr) / max(abs(f_cur), eps), ...
        'f_current',          f_cur, ...
        'f_refocused',        f_tr, ...
        'pred_current',       [], ...
        'pred_refocused',     []);

    if options.PredictShift
        res.pred_current = predict_shift_widening(pts_m, q_cur, target_m, beam_u, k_target);
        if ~strcmp(capability, 'none')
            res.pred_refocused = predict_shift_widening(pts_m, q_tr, target_m, beam_u, k_target);
        end
    end
end

function pred = predict_shift_widening(pts_m, q, target_m, beam_u, k)
% Focal-shift prediction in a +-10 mm cube; when the peak lies on its boundary
% (shift beyond the cube), search once more in a +-60 mm cube at 2 mm steps, so
% large shifts (e.g. a geometric drive that cannot reach a shallow target) are
% predicted rather than clipped at the cube edge. at_edge stays true if the peak
% is still on the wider boundary.
pred = predict_focal_shift(pts_m, q, target_m, beam_u, k);
if pred.at_edge
    pred = predict_focal_shift(pts_m, q, target_m, beam_u, k, 'HalfWidth', 60e-3, 'Step', 2e-3, 'FineStep', 0.5e-3);
end
end

function c = point_coherence(pts_m, q, k)
% |sum q| / sum |q| over cells of lambda/6 (cells merge the BLI side lobes of one surface point)
    cell_size = 2*pi/k/6;
    [~, ~, cid] = unique(round(pts_m / cell_size), 'rows');
    q_c = accumarray(cid, q);
    c = abs(sum(q_c)) / max(sum(abs(q_c)), eps);
end
