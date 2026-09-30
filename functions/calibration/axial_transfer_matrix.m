function H = axial_transfer_matrix(tran, z_bowl_mm, c, rho, forward_model)
% AXIAL_TRANSFER_MATRIX  On-axis complex pressure per unit element drive
%
% The on-axis field of an annular array is linear in the complex element
% drives w_i = v_i * exp(1i*phi_i):
%
%   p(z) = H(z,:) * w        [Pa],  w in [m/s]
%
% H depends only on geometry, frequency and medium, so it is computed once
% and every forward evaluation during fitting is a matrix product.
%
% Use as:
%   H = axial_transfer_matrix(tran, z_bowl_mm, c, rho)
%   H = axial_transfer_matrix(tran, z_bowl_mm, c, rho, forward_model)
%
% Input:
%   tran          - transducer struct from load_equipment_config (uses .transducer)
%   z_bowl_mm     - [N_z x 1] axial positions from the bowl apex [mm]
%   c, rho        - sound speed [m/s] and density [kg/m³] of the medium
%   forward_model - 'rayleigh' (default) | 'oneil'
%
% Output:
%   H - [N_z x N_elem] complex transfer matrix
%
% Notes:
%   'rayleigh' reproduces RAYLEIGH_AXIAL_INTENSITY (same quadrature and bowl
%   sagitta) but vectorised: |H*w|.^2/(2*rho*c)*1e-4 equals
%   rayleigh_axial_intensity(angle(w), v, ...) for uniform |w| = v.
%   'oneil' builds each column from the complex per-ring O'Neil field
%   (FOCUSEDBOWLONEIL outer minus inner bowl), i.e. the terms that
%   FOCUSEDANNULUSONEIL sums before taking abs(); same (z - 0.5 mm)
%   evaluation offset as elsewhere in PRESTUS calibration.
%
% See also: RAYLEIGH_AXIAL_INTENSITY, FIT_TRANSDUCER_MODEL

    if nargin < 5 || isempty(forward_model); forward_model = 'rayleigh'; end
    ann  = tran.transducer.annular;
    n    = ann.elem_n;
    f    = tran.transducer.freq_hz;
    z    = z_bowl_mm(:) * 1e-3;
    H    = complex(zeros(numel(z), n));

    switch forward_model
        case 'rayleigh'
            R_c = ann.curv_radius_mm * 1e-3;
            a   = ann.elem_id_mm(:)' * 1e-3 / 2;
            b   = ann.elem_od_mm(:)' * 1e-3 / 2;
            k   = 2 * pi * f / c;
            [xi, wi] = gauss_legendre(32);
            for i = 1:n
                rho_s = (0.5 * (b(i) - a(i)) * xi + 0.5 * (b(i) + a(i)))';   % [1 x nq]
                dw    = (0.5 * (b(i) - a(i)) * wi)';
                z_s   = R_c - sqrt(max(R_c^2 - rho_s.^2, 0));
                r     = sqrt(rho_s.^2 + (z - z_s).^2);                     % [N_z x nq]
                kern  = (1i * k / (2*pi)) .* exp(-1i * k .* r) ./ r .* (2*pi * rho_s);
                H(:, i) = rho * c * (kern * dw');
            end
        case 'oneil'
            R_c = ann.curv_radius_mm * 1e-3;
            for i = 1:n
                [~, ~, p_out] = focusedBowlONeil(R_c, ann.elem_od_mm(i) * 1e-3, 1, f, c, rho, z - 0.5e-3, []);
                p_in = 0;
                if ann.elem_id_mm(i) > 0
                    [~, ~, p_in] = focusedBowlONeil(R_c, ann.elem_id_mm(i) * 1e-3, 1, f, c, rho, z - 0.5e-3, []);
                end
                H(:, i) = p_out(:) - p_in(:);
            end
        otherwise
            error('axial_transfer_matrix: unknown forward model ''%s''.', forward_model);
    end
end

function [x, w] = gauss_legendre(n)
    x = cos(pi * ((1:n)' - 0.25) / (n + 0.5));
    for iter = 1:100
        P0 = ones(n, 1); P1 = x;
        for k = 2:n
            P2 = ((2*k - 1) * x .* P1 - (k - 1) * P0) / k;
            P0 = P1; P1 = P2;
        end
        dP = n * (P0 - x .* P1) ./ (1 - x.^2);
        dx = P1 ./ dP;
        x  = x - dx;
        if max(abs(dx)) < 1e-14; break; end
    end
    w = 2 ./ ((1 - x.^2) .* dP.^2);
end
