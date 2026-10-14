function H = field_transfer_matrix(tran, r_mm, z_bowl_mm, c, rho, n_phi)
% FIELD_TRANSFER_MATRIX  Complex pressure per unit element drive at points off the axis
%
% Generalises AXIAL_TRANSFER_MATRIX ('rayleigh') to points at a radial
% distance r from the beam axis. The field of an annular array is
% axisymmetric, so a point is given by (r, z): the Rayleigh integral over
% each ring runs over the bowl radius ρ (Gauss–Legendre) and the azimuth φ
% (trapezoidal rule, exact for the periodic integrand once n_phi resolves it):
%
%   p(r, z) = H(r, z, :) * w,   d = sqrt(r² + ρ² - 2 r ρ cos φ + (z - z_s(ρ))²)
%
% At r = 0 it equals AXIAL_TRANSFER_MATRIX. A within-ring velocity profile
% (annular.elem_weight, RING_PROFILE_TRANSDUCER) multiplies each column.
%
% Use as:
%   H = field_transfer_matrix(tran, r_mm, z_bowl_mm, c, rho)
%   H = field_transfer_matrix(tran, r_mm, z_bowl_mm, c, rho, n_phi)
%
% Input:
%   tran      - annular transducer struct (load_equipment_config / RING_PROFILE_TRANSDUCER)
%   r_mm      - [N x 1] radial distance from the axis [mm]
%   z_bowl_mm - [N x 1] axial position from the bowl apex [mm] (same size as r_mm)
%   c, rho    - sound speed [m/s] and density [kg/m³]
%   n_phi     - azimuth samples (default: per block of points, enough for the
%               largest phase variation k·2rρ/d around the ring there, at least 64)
%
% Output:
%   H - [N x N_elem] complex transfer matrix (pressure [Pa] per element velocity [m/s])
%
% See also: AXIAL_TRANSFER_MATRIX, FIT_RING_PROFILE, RING_PROFILE_TRANSDUCER

    ann = tran.transducer.annular;
    n   = ann.elem_n;
    f   = tran.transducer.freq_hz;
    k   = 2 * pi * f / c;
    R   = ann.curv_radius_mm * 1e-3;
    r   = r_mm(:) * 1e-3;
    z   = z_bowl_mm(:) * 1e-3;
    assert(numel(r) == numel(z), 'field_transfer_matrix:size', 'r_mm and z_bowl_mm must have the same size.');
    fixed_phi = nargin >= 6 && ~isempty(n_phi);
    rho_max = max(ann.elem_od_mm) / 2 * 1e-3;
    z_rim   = R - sqrt(max(R^2 - rho_max^2, 0));
    [xi, wi] = gauss_legendre(32);
    H = complex(zeros(numel(r), n));
    % Blocks of points ordered by distance to the bowl, so each block gets the
    % azimuth sampling it needs (points near the bowl need the most)
    dpt = max(sqrt(max(r - rho_max, 0).^2 + max(z - z_rim, 0).^2), 0.5e-3);
    [~, order] = sort(dpt);
    blk = 512;
    for p0 = 1:blk:numel(r)
        idx = order(p0:min(numel(r), p0 + blk - 1));
        if fixed_phi
            nph = n_phi;
        else
            nph = max(64, 2 * ceil(2 * k * max(r(idx)) * rho_max / min(dpt(idx))) + 16);
        end
        phi = (0:nph - 1) * 2 * pi / nph;  dphi = 2 * pi / nph;  cphi = cos(phi);
        sub = max(1, floor(2e6 / (32 * nph)));          % points per pass (memory)
        for q0 = 1:sub:numel(idx)
          jdx = idx(q0:min(numel(idx), q0 + sub - 1));
          rp = reshape(r(jdx), 1, 1, []); zp = reshape(z(jdx), 1, 1, []);
          for i = 1:n
            a = ann.elem_id_mm(i) * 1e-3 / 2; b = ann.elem_od_mm(i) * 1e-3 / 2;
            s  = 0.5 * (b - a) * xi + 0.5 * (b + a);   % [nq x 1] ring radius samples
            ws = 0.5 * (b - a) * wi;
            zs = R - sqrt(max(R^2 - s.^2, 0));
            d  = sqrt(rp.^2 + s.^2 - 2 * rp .* s .* cphi + (zp - zs).^2);   % [nq x n_phi x np]
            kern = exp(-1i * k * d) ./ d;
            I = squeeze(sum(sum(kern .* (ws .* s), 1) * dphi, 2));         % ∫∫ e^{-ikd}/d ρ dρ dφ
            H(jdx, i) = rho * c * (1i * k / (2 * pi)) * I(:);
          end
        end
    end
    if isfield(ann, 'elem_weight') && ~isempty(ann.elem_weight)
        H = H .* reshape(ann.elem_weight, 1, []);
    end
end

function [x, w] = gauss_legendre(n)
    x = cos(pi * ((1:n)' - 0.25) / (n + 0.5));
    for iter = 1:100
        P0 = ones(n, 1); P1 = x;
        for kk = 2:n
            P2 = ((2*kk - 1) * x .* P1 - (kk - 1) * P0) / kk;
            P0 = P1; P1 = P2;
        end
        dP = n * (P0 - x .* P1) ./ (1 - x.^2);
        dx = P1 ./ dP; x = x - dx;
        if max(abs(dx)) < 1e-14; break; end
    end
    w = 2 ./ ((1 - x.^2) .* dP.^2);
end
