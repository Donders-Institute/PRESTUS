function diameter_m = matrix_elem_diameter(tr_matrix)
% MATRIX_ELEM_DIAMETER  Aperture diameter of a circular (disc/bowl) matrix element
%
% Disc and bowl elements are circular, while the element size is configured
% as width x height. By default the diameter is chosen so that the element
% has the same radiating area as the width x height rectangle: a flat disc of
% area w*h for 'disc', a spherical cap of area w*h on a sphere of radius
% curv_radius_mm for 'bowl'. Set tr_matrix.elem_diameter_mm to use an explicit
% aperture diameter instead (e.g. the manufacturer's element diameter).
%
% Use as:
%   diameter_m = matrix_elem_diameter(tr_matrix)
%
% Input:
%   tr_matrix - matrix transducer struct: elem_shape ('disc' | 'bowl'),
%               elem_width_mm, elem_height_mm [mm], curv_radius_mm [mm] (bowl),
%               optional elem_diameter_mm [mm]
%
% Output:
%   diameter_m - element aperture diameter [m]
%
% See also: CREATE_MATRIX_KARRAY, TRANSDUCER_SETUP

    if isfield(tr_matrix, 'elem_diameter_mm') && ~isempty(tr_matrix.elem_diameter_mm)
        diameter_m = tr_matrix.elem_diameter_mm * 1e-3;
        return
    end

    A_target = tr_matrix.elem_width_mm * tr_matrix.elem_height_mm * 1e-6;   % [m²]

    switch lower(tr_matrix.elem_shape)
        case 'disc'
            diameter_m = sqrt(A_target * 4 / pi);

        case 'bowl'
            r_c = tr_matrix.curv_radius_mm / 1e3;
            a_min = 0.001;        % Lower bound of aperture radius [m]
            a_max = r_c * 0.999;  % Slightly less than full bowl radius

            % Cap area A(a) = 2*pi*r_c*(r_c - sqrt(r_c^2 - a^2)) for aperture radius a
            area_diff = @(a) 2*pi*r_c*(r_c - sqrt(r_c^2 - a^2)) - A_target;
            if sign(area_diff(a_min)) == sign(area_diff(a_max))
                error('No sign change in interval: cannot find bowl diameter. Possibly A_target is out of bounds.');
            end
            diameter_m = 2 * fzero(area_diff, [a_min, a_max]);

        otherwise
            error('matrix_elem_diameter: elem_shape "%s" has no circular aperture.', tr_matrix.elem_shape);
    end
end
