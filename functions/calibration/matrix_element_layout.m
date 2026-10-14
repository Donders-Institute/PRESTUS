function L = matrix_element_layout(parameters, tr)
% MATRIX_ELEMENT_LAYOUT  Matrix-array elements in the transducer frame (no grid)
%
% Element centres, normals and in-plane axes of a matrix array in the
% transducer frame: origin at the transducer reference point (apex of a
% curved array, centre of a flat one), beam along +z. Uses the same layout
% code as the simulation (CONVERT_TO_ELEMENT_POS / EXTRACT_ELEMENT_POS) and
% the same element orientation as CREATE_MATRIX_KARRAY, so an analytical
% model built on it (MATRIX_RAYLEIGH_PRESSURE) sees the simulated array.
%
% Use as:
%   L = matrix_element_layout(parameters, tr)
%
% Input:
%   parameters - PRESTUS config (passed to the layout functions)
%   tr         - matrix transducer struct (tr.matrix with matrix_shape,
%                elem_shape, elem_width_mm, elem_height_mm, is_curved,
%                curv_radius_mm; optional elem_diameter_mm)
%
% Output:
%   L - struct: pos_m [3xN], ex, ey, ez [3xN] element axes (ez = normal),
%       shape ('rect' | 'disc' | 'bowl'), width_m, height_m, diameter_m, n
%
% Clover setups are not supported.
%
% See also: MATRIX_RAYLEIGH_PRESSURE, MATRIX_FOCAL_PEAK, CREATE_MATRIX_KARRAY

    m = tr.matrix;
    if isfield(m, 'is_clover_setup') && isequal(m.is_clover_setup, true)
        error('matrix_element_layout:clover', 'Clover setups are not supported.');
    end
    if is_flat_matrix(m); m.curv_radius_mm = inf; tr.matrix = m; end
    origin = [0; 0; 0];
    switch m.matrix_shape.type
        case 'define_here'
            [pos_m, tr] = convert_to_element_pos(parameters, tr, origin, origin);
        case 'extract_from_file'
            [pos_m, tr] = extract_element_pos(parameters, tr, origin);
        otherwise
            error('matrix_element_layout:shape', 'Matrix shape %s is unknown.', m.matrix_shape.type);
    end
    m = tr.matrix;
    n = size(pos_m, 2);
    aim = matrix_element_aim(m, pos_m, origin);
    rot = zeros(1, n);
    if isfield(m, 'elem_rot_deg') && ~isempty(m.elem_rot_deg); rot = m.elem_rot_deg(:)'; end

    ex = zeros(3, n); ey = ex; ez = ex;
    for i = 1:n
        nv = (aim(:, i) - pos_m(:, i)) / norm(aim(:, i) - pos_m(:, i));
        % Same rotation as CREATE_MATRIX_KARRAY: [0;0;1] onto the normal, then
        % the element's own rotation about it
        v = cross([0; 0; 1], nv); s = norm(v); c = dot([0; 0; 1], nv);
        vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
        R = eye(3) + vx + vx^2 * ((1 - c) / (s^2 + eps));
        R = R * [cosd(rot(i)) -sind(rot(i)) 0; sind(rot(i)) cosd(rot(i)) 0; 0 0 1];
        ex(:, i) = R(:, 1); ey(:, i) = R(:, 2); ez(:, i) = R(:, 3);
    end

    L = struct('pos_m', pos_m, 'ex', ex, 'ey', ey, 'ez', ez, 'shape', lower(char(m.elem_shape)), ...
        'width_m', m.elem_width_mm * 1e-3, 'height_m', m.elem_height_mm * 1e-3, 'diameter_m', [], 'n', n);
    if ismember(L.shape, {'disc', 'bowl'})
        if strcmp(L.shape, 'bowl') && is_flat_matrix(m); m.elem_shape = 'disc'; end
        L.diameter_m = matrix_elem_diameter(m);
    end
end
