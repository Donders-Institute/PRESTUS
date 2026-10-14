function p = matrix_rayleigh_pressure(L, freq_hz, c, amp_pa, phase_rad, pts_m, opts)
% MATRIX_RAYLEIGH_PRESSURE  Free-field complex pressure of a matrix array (Rayleigh integral)
%
% Sums the Rayleigh integral over every element of a matrix array, each
% element split into small flat patches:
%
%   p(x) = sum_elements a_e exp(i phi_e) (i k / 2 pi) sum_patches dA exp(-i k r) / r
%
% with a_e the element surface pressure [Pa] (velocity a_e / (rho c)), phi_e
% its drive phase, k = 2 pi f / c and r the patch-to-point distance; time
% convention exp(+i w t), so the steering phase phi_e = k |element - target|
% (CREATE_MATRIX_KARRAY) focuses at the target. Each element radiates as a
% baffled piston: its directivity comes from the patch sum. Bowl elements
% are integrated as flat discs (their sag over a few mm is far below a
% wavelength).
%
% Use as:
%   p = matrix_rayleigh_pressure(L, freq_hz, c, amp_pa, phase_rad, pts_m)
%   p = matrix_rayleigh_pressure(L, freq_hz, c, amp_pa, phase_rad, pts_m, opts)
%
% Input:
%   L         - element layout (MATRIX_ELEMENT_LAYOUT)
%   freq_hz   - drive frequency [Hz]
%   c         - sound speed [m/s]
%   amp_pa    - [1xN] (or scalar) element surface pressure [Pa]
%   phase_rad - [1xN] element phases [rad]
%   pts_m     - [3xM] evaluation points [m], transducer frame
%   opts      - optional: .patch_per_lambda (default 6), .chunk (points per
%               block, default 1000)
%
% Output:
%   p - [Mx1] complex pressure [Pa]
%
% See also: MATRIX_ELEMENT_LAYOUT, MATRIX_FOCAL_PEAK, RAYLEIGH_AXIAL_INTENSITY

    if nargin < 7; opts = struct(); end
    ppl = get_opt(opts, 'patch_per_lambda', 6);
    chunk = get_opt(opts, 'chunk', 1000);
    k = 2 * pi * freq_hz / c;
    dmax = (c / freq_hz) / ppl;
    amp_pa = amp_pa(:)' .* ones(1, L.n);
    phase_rad = phase_rad(:)';

    % Patches of one element in its local frame: positions [2 x P] and areas [1 x P]
    switch L.shape
        case 'rect'
            nx = max(1, ceil(L.width_m / dmax)); ny = max(1, ceil(L.height_m / dmax));
            u = ((1:nx) - 0.5) / nx * L.width_m - L.width_m / 2;
            v = ((1:ny) - 0.5) / ny * L.height_m - L.height_m / 2;
            [U, V] = ndgrid(u, v);
            loc = [U(:)'; V(:)']; dA = repmat(L.width_m * L.height_m / numel(U), 1, numel(U));
        case {'disc', 'bowl'}
            a = L.diameter_m / 2; nr = max(1, ceil(a / dmax));
            loc = zeros(2, 0); dA = zeros(1, 0);
            for ir = 1:nr
                r0 = (ir - 1) / nr * a; r1 = ir / nr * a; rm = (r0 + r1) / 2;
                nt = max(1, ceil(2 * pi * rm / dmax));
                th = ((1:nt) - 0.5) / nt * 2 * pi;
                loc = [loc, [rm * cos(th); rm * sin(th)]]; %#ok<AGROW>
                dA = [dA, repmat(pi * (r1^2 - r0^2) / nt, 1, nt)]; %#ok<AGROW>
            end
        otherwise
            error('matrix_rayleigh_pressure:shape', 'Unknown element shape %s.', L.shape);
    end
    P = size(loc, 2);

    % All patches of all elements [3 x N*P] and their complex weights
    src = zeros(3, L.n * P); w = zeros(1, L.n * P);
    for e = 1:L.n
        idx = (e - 1) * P + (1:P);
        src(:, idx) = L.pos_m(:, e) + L.ex(:, e) * loc(1, :) + L.ey(:, e) * loc(2, :);
        w(idx) = amp_pa(e) * exp(1i * phase_rad(e)) * (1i * k / (2 * pi)) * dA;
    end

    M = size(pts_m, 2);
    p = complex(zeros(M, 1));
    for s0 = 1:chunk:M
        id = s0:min(M, s0 + chunk - 1);
        r = sqrt((pts_m(1, id)' - src(1, :)).^2 + (pts_m(2, id)' - src(2, :)).^2 + (pts_m(3, id)' - src(3, :)).^2);
        p(id) = (exp(-1i * k * r) ./ r) * w.';
    end
end

function v = get_opt(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
