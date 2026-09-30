function [elem_pos_m, tr] = convert_to_element_pos(parameters, tr, trans_pos_m, focus_pos_m) %#ok<INUSD>
% CONVERT_TO_ELEMENT_POS  Generate element positions for matrix transducers
%
% Generates 3-D element coordinates of a matrix array transducer from the
% configuration in the parameter structure. Supported grid layouts are
% rectangular (columns along x at elem_width_mm + elem_spacing_width_mm,
% rows along y at elem_height_mm + elem_spacing_height_mm, clipped to
% outer_diameter_mm), Fibonacci spiral, and Fermat spiral (makeCartBowl /
% makeCartDisc). Elements lie on the spherical cap when is_curved, with the
% apex at trans_pos_m and the beam along +z (as CREATE_MATRIX_KARRAY
% assumes). rect.sparsity_factor is not applied.
%
% Use as:
%   [elem_pos_m, tr] = convert_to_element_pos(parameters, tr, trans_pos_m, focus_pos_m)
%
% Input:
%   parameters  - (1,1) simulation parameters struct (unused; kept for the call signature)
%   tr          - (1,1) transducer parameter struct
%   trans_pos_m - [3x1] transducer position [m]
%   focus_pos_m - [3x1] focus position [m] (unused; steering is applied via phases)
%
% Output:
%   elem_pos_m - [3xN] element center coordinates [m]
%   tr         - transducer struct with tr.matrix.elem_n set
%
% See also: EXTRACT_ELEMENT_POS, TRANSDUCER_SETUP

    % Layout in mm in the transducer frame: apex at the origin, beam along +z
    % (as CREATE_MATRIX_KARRAY assumes). focus_pos_m is not used: elements sit
    % on the natural bowl and steering is applied through element phases.
    m  = tr.matrix;
    gs = m.matrix_shape.define_here.grid_shape;
    a  = m.outer_diameter_mm / 2;        % aperture radius [mm]
    R  = m.curv_radius_mm;               % radius of curvature [mm]
    trans_pos_m = trans_pos_m(:);

    switch gs.type
        case 'fibonacci'
            % Sunflower (golden-angle) spiral, uniform density over the aperture
            N  = double(gs.fibonacci.elem_n);
            i  = (1:N)';
            r  = a * sqrt((i - 0.5) / N);
            th = pi * (3 - sqrt(5)) * (i - 1);
            elem_pos_mm = curve([r .* cos(th), r .* sin(th)], m);

            % Report pitch and the element size the kerf leaves room for
            nn = inf(N, 1);
            for k = 1:N
                d = hypot(elem_pos_mm(:, 1) - elem_pos_mm(k, 1), elem_pos_mm(:, 2) - elem_pos_mm(k, 2));
                d(k) = inf;
                nn(k) = min(d);
            end
            fprintf('Estimated avg pitch: %.3f mm\n', mean(nn));
            fprintf('Recommended element size (kerf %.3f mm): %.3f mm\n', ...
                gs.fibonacci.kerf_mm, mean(nn) - gs.fibonacci.kerf_mm);

        case 'rect'
            % Columns along x (element width), rows along y (element height),
            % centred on the axis; centres outside the aperture are dropped
            g = gs.rect;
            n_col = double(g.elem_n_col); n_row = double(g.elem_n_row);
            x = ((1:n_col) - (n_col + 1) / 2) * (m.elem_width_mm  + g.elem_spacing_width_mm);
            y = ((1:n_row) - (n_row + 1) / 2) * (m.elem_height_mm + g.elem_spacing_height_mm);
            [X, Y] = meshgrid(x, y);
            keep = hypot(X(:), Y(:)) <= a;
            elem_pos_mm = curve([X(keep), Y(keep)], m);

        case 'fermat'
            % k-Wave Fermat-spiral sampling of the bowl (or the flat disc)
            N = double(gs.fermat.elem_n);
            if m.is_curved
                p = makeCartBowl([0 0 0], R * 1e-3, 2 * a * 1e-3, [0 0 R * 1e-3], N, false);
            else
                p = makeCartDisc([0 0 0], a * 1e-3, [0 0 1], N, false, true);
            end
            elem_pos_mm = p' * 1e3;

        otherwise
            error('Grid shape %s is unknown or not implemented.', gs.type)
    end

    if isempty(elem_pos_mm)
        error('convert_to_element_pos: the %s grid places no elements inside outer_diameter_mm = %g.', ...
              gs.type, m.outer_diameter_mm);
    end
    tr.matrix.elem_n = size(elem_pos_mm, 1);
    elem_pos_m = elem_pos_mm' * 1e-3 + trans_pos_m;
end

function pos = curve(xy, m)
% Place 2-D element centres [N x 2, mm] on the spherical cap (or the plane)
    if m.is_curved
        R = m.curv_radius_mm;
        s = R^2 - sum(xy.^2, 2);
        if any(s < 0)
            error('Some elements fall outside the spherical cap (curv_radius_mm %g).', R);
        end
        z = R - sqrt(s);
    else
        z = zeros(size(xy, 1), 1);
    end
    pos = [xy, z];
end
