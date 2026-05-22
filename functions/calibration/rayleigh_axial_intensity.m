function I_axial = rayleigh_axial_intensity(phases_rad, velocity, parameters, dist_bowl_mm)
% RAYLEIGH_AXIAL_INTENSITY  On-axis intensity via the Rayleigh–Sommerfeld integral
%
% Computes the on-axis complex pressure by numerically integrating the
% Rayleigh–Sommerfeld kernel over the surface of each annular ring element,
% then converts to intensity. This is the forward model used by BabelBrain's
% calibration pipeline (TxCalibration.py, RayleighCoeff).
%
% For element i with inner radius a_i, outer radius b_i, amplitude u0,
% phase phi_i:
%
%   p(z) = (i k rho c) / (2 pi) * integral_over_element
%              exp(-i k r) / r  dS
%
% where r = sqrt(rho_s^2 + z^2) for a source point at lateral position rho_s.
% For on-axis evaluation this reduces to a 1-D integral over the radial
% coordinate of each ring (azimuthal integral is 2pi * rho_s * d_rho_s).
%
% This model is more accurate than O'Neil near the natural focus and for
% configurations where phase steering moves the focus well away from the
% geometric focal point, because O'Neil assumes uniform piston motion on a
% continuous spherical bowl, whereas the Rayleigh integral respects ring gaps.
%
% Use as:
%   I_axial = rayleigh_axial_intensity(phases_rad, velocity, parameters, dist_bowl_mm)
%
% Input:
%   phases_rad    - [1 x N_elem] element phases [rad]
%   velocity      - scalar particle velocity [m/s]
%   parameters    - PRESTUS config with transducer.annular and
%                   medium_properties.water
%   dist_bowl_mm  - [N_pts x 1] axial evaluation positions measured from
%                   bowl face (same convention as O'Neil: bowl surface = 0) [mm]
%
% Output:
%   I_axial  - [N_pts x 1] on-axis intensity [W/cm²]
%
% Notes:
%   - Integration uses Gauss–Legendre quadrature over each ring.
%   - n_quad points per ring; increase for very narrow rings or high freq.
%   - For moderate frequencies and ring widths, 32 points is sufficient.
%
% See also: PERFORM_JOINT_DEPTH_FIT, PHASE_OPTIMIZATION_ANNULUS_FULL_CURVE

arguments
    phases_rad   (1,:) {mustBeNumeric}
    velocity     (1,1) {mustBeNumeric}
    parameters   (1,1) struct
    dist_bowl_mm (:,1) {mustBeNumeric}
end

% Transducer parameters
n_elem   = parameters.transducer.annular.elem_n;
R_c      = parameters.transducer.annular.curv_radius_mm * 1e-3;   % [m]
id_mm    = parameters.transducer.annular.elem_id_mm(:)' * 1e-3;   % inner radii [m]
od_mm    = parameters.transducer.annular.elem_od_mm(:)' * 1e-3;   % outer radii [m]
a_vec    = id_mm / 2;   % inner radius per element [m]
b_vec    = od_mm / 2;   % outer radius per element [m]

freq     = parameters.transducer.freq_hz;
c        = parameters.medium_properties.water.sound_speed;
rho      = parameters.medium_properties.water.density;
k        = 2 * pi * freq / c;   % wavenumber [rad/m]
omega    = 2 * pi * freq;

% Evaluation z-positions (measured from bowl surface, positive toward target)
z_eval   = dist_bowl_mm(:) * 1e-3;   % [m]
N_pts    = numel(z_eval);

% Gauss-Legendre quadrature nodes and weights on [-1, 1]
n_quad   = 32;
[xi, wi] = gauss_legendre_nodes(n_quad);

% Accumulate pressure contribution from each element
p_total  = zeros(N_pts, 1);

for i = 1:n_elem
    a  = a_vec(i);
    b  = b_vec(i);
    ph = phases_rad(i);

    % Map quadrature nodes from [-1,1] to [a, b]
    rho_s = 0.5 * (b - a) * xi + 0.5 * (b + a);   % [n_quad x 1]
    dw    = 0.5 * (b - a) * wi;                     % scaled weights

    % Element source phasor: velocity * exp(1i * phase)
    u0 = velocity * exp(1i * ph);

    % For each axial point, integrate kernel over ring
    for iz = 1:N_pts
        z = z_eval(iz);

        % Distance from each ring point to axial evaluation point
        % Bowl surface offset: source points on curved bowl at z_s = -(R_c - sqrt(R_c^2 - rho_s^2))
        % For simplicity (flat-piston approximation per ring, valid when ring width << R_c):
        %   r = sqrt(rho_s^2 + z^2)
        % For full spherical bowl: add z_s offset per source point
        z_s = R_c - sqrt(max(R_c^2 - rho_s.^2, 0));   % axial position of source on bowl [m]
        r   = sqrt(rho_s.^2 + (z - z_s).^2);          % source-to-field distance [m]

        % Rayleigh–Sommerfeld kernel (first Rayleigh integral):
        %   (i k / 2 pi) * exp(-i k r) / r * (2 pi * rho_s)  per unit rho_s
        kernel = (1i * k / (2*pi)) .* exp(-1i * k .* r) ./ r .* (2*pi * rho_s);

        p_total(iz) = p_total(iz) + u0 * rho * c * sum(dw .* kernel);
    end
end

% Pressure to intensity [W/cm²]
I_axial = abs(p_total).^2 / (2 * rho * c) * 1e-4;

end

% --- Gauss-Legendre nodes and weights on [-1, 1] -----------------------
function [x, w] = gauss_legendre_nodes(n)
% Simple iterative computation of Gauss-Legendre nodes/weights.
    x = cos(pi * ((1:n)' - 0.25) / (n + 0.5));
    P0 = ones(n, 1);  P1 = x;
    for iter = 1:100
        P0 = ones(n, 1);  P1 = x;
        for k = 2:n
            P2 = ((2*k-1) * x .* P1 - (k-1) * P0) / k;
            P0 = P1;  P1 = P2;
        end
        dP = n * (P0 - x .* P1) ./ (1 - x.^2);
        dx = P1 ./ dP;
        x  = x - dx;
        if max(abs(dx)) < 1e-14; break; end
    end
    w = 2 ./ ((1 - x.^2) .* dP.^2);
end
