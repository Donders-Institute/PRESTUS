function [karray, tr] = create_matrix_karray(kgrid, karray, parameters, tr, elem_pos_m, trans_pos, focus_pos)
%CREATE_MATRIX_KARRAY Adds elements of a matrix transducer to a kWave array.
%
% Inputs:
%   kgrid            - kWave grid object
%   karray           - kWaveArray object to which elements are added
%   parameters       - Simulation parameters struct
%   tr  - Struct containing transducer geometry, type, curvature, and element properties
%   elem_pos_m          - Nx3 matrix of element positions in meters
%   trans_pos        - 1x3 transducer reference position (indices in kgrid)
%   focus_pos        - 1x3 focus position (indices in kgrid)
%
% Output:
%   karray           - Updated kWaveArray object with elements added

    % Convert positions from kgrid indices to meters
    trans_pos_m = [kgrid.x_vec(trans_pos(1)), ...
                   kgrid.y_vec(trans_pos(2)), ...
                   kgrid.z_vec(trans_pos(3))]';
               
    focus_pos_m = [kgrid.x_vec(focus_pos(1)), ...
                   kgrid.y_vec(focus_pos(2)), ...
                   kgrid.z_vec(focus_pos(3))]';

    natural_focus_pos_m = trans_pos_m + [0, 0, tr.matrix.curv_radius_mm / 1000]';

    % Apply Clover setup if requested
    if tr.matrix.is_clover_setup    
        elem_pos_m = create_clover_array(parameters, tr.matrix, elem_pos_m, trans_pos_m, focus_pos_m);
    end
    
    tr.matrix.elem_n = size(elem_pos_m, 2);

    % Initialize source amplitudes (uniform)
    tr.matrix.elem_amp = tr.matrix.elem_amp(1) * ones(1, tr.matrix.elem_n);

    % Wavelength and wavenumber for phase calculation
    lambda = parameters.medium_properties.water.sound_speed / tr.freq_hz;
    k = 2 * pi / lambda;

    % Initialize source phases, scaled vectors, tx, ty, tz
    elem_phase_rad = zeros(1, tr.matrix.elem_n);
    scaled_vectors = zeros(3, tr.matrix.elem_n);
    tx = zeros(1, tr.matrix.elem_n);
    ty = zeros(1, tr.matrix.elem_n);
    tz = zeros(1, tr.matrix.elem_n);

    % Loop over each element and add it to karray
    for ind = 1:tr.matrix.elem_n
        el_pos_m_i = elem_pos_m(:, ind);

        % Vector from element to natural focus to position elements to
        % natural focus
        vec_to_nat_focus = natural_focus_pos_m - el_pos_m_i;
        norm_vec = vec_to_nat_focus / norm(vec_to_nat_focus);
        scaled_vectors(:, ind) = vec_to_nat_focus;

        % Compute rotation matrix to align default normal [0;0;1] to vec_to_focus
        z0 = [0;0;1];
        v = cross(z0, norm_vec);
        s = norm(v);
        c = dot(z0, norm_vec);
        vx = [ 0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0 ];
        R = eye(3) + vx + vx^2*((1-c)/(s^2+eps));


        % Convert rotation matrix to ZYX Euler angles (extrinsic)
        yaw   = atan2d(R(2,1), R(1,1));  % around z
        pitch = atan2d(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));  % around y
        roll  = atan2d(R(3,2), R(3,3));  % around x

        tx(ind) = roll;
        ty(ind) = pitch;
        tz(ind) = yaw;

        % Determine element type
        switch lower(tr.matrix.elem_shape)
            case 'rect'
                karray.addRectElement(el_pos_m_i, tr.matrix.elem_height_mm, tr.matrix.elem_width_mm, [roll, pitch, yaw]);

            case 'disc'
                % Disc with same area as rectangular element
                diameter = sqrt(tr.matrix.elem_height_mm * tr.matrix.elem_width_mm * 4 / pi);
                karray.addDiscElement(el_pos_m_i, diameter, natural_focus_pos_m);

            case 'bowl'
                r_c = tr.matrix.curv_radius_mm / 1e3;
                A_target = tr.matrix.elem_height_mm * tr.matrix.elem_width_mm;
                a_min = 0.001; % Lower bound of aperture radius (in m)
                a_max = r_c * 0.999; % Slightly less than full bowl radius

                
                % Define the function to solve: A_cap(a) - A_target = 0
                area_diff = @(a) 2*pi*r_c*(r_c - sqrt(r_c^2 - a^2)) - A_target;

                % Check if the function changes sign in the interval
                if sign(area_diff(a_min)) == sign(area_diff(a_max))
                    error('No sign change in interval: cannot find bowl diameter. Possibly A_target is out of bounds.');
                end

                a_solution = fzero(area_diff, [a_min, a_max]);
                diameter = 2 * a_solution;
                karray.addBowlElement(el_pos_m_i, r_c, diameter, natural_focus_pos_m);
        end

        % Calculate phase delay based on set focus
        distance = sqrt((el_pos_m_i(1) - focus_pos_m(1))^2 + (el_pos_m_i(2) - focus_pos_m(2))^2 + (el_pos_m_i(3) - focus_pos_m(3))^2);

        elem_phase_rad(ind) = mod(k * distance, 2 * pi);
    end

    tr.elem_phase_rad = elem_phase_rad;
    tr.elem_phase_deg = rad2deg(elem_phase_rad);

    % [DEBUG] visualize matrix element orientation
    if parameters.simulation.debug == 1
        h = figure;
        hold on;
        axis equal;
    
        scatter3(elem_pos_m(1,:), elem_pos_m(2,:), elem_pos_m(3,:), 'b.');
        scatter3(natural_focus_pos_m(1,:), natural_focus_pos_m(2,:), natural_focus_pos_m(3,:), 'r.');
        plot3(focus_pos_m(1), focus_pos_m(2), focus_pos_m(3), 'go', 'MarkerFaceColor', 'g');
    
        for ind = 1:tr.elem_n
            el_pos = elem_pos_m(:, ind);
    
            % --- Ground truth direction (element -> focus)
            vec_truth = focus_pos_m - el_pos;
            vec_truth = vec_truth / norm(vec_truth);
    
            % --- Euler angles (degrees -> radians)
            roll  = deg2rad(tx(ind));
            pitch = deg2rad(ty(ind));
            yaw   = deg2rad(tz(ind));
    
            % --- Rotation matrix (extrinsic ZYX)
            Rz = [cos(yaw) -sin(yaw) 0;
                sin(yaw)  cos(yaw) 0;
                0         0        1];
    
            Ry = [cos(pitch) 0 sin(pitch);
                0          1 0;
                -sin(pitch) 0 cos(pitch)];
    
            Rx = [1 0 0;
                0 cos(roll) -sin(roll);
                0 sin(roll)  cos(roll)];
    
            R = Rz * Ry * Rx;
    
            % --- Predicted element normal
            vec_pred = R * [0; 0; 1];
    
            % --- Scale vectors for visualization
            scale = norm(focus_pos_m - el_pos);
            vec_pred  = vec_pred  * scale;
            vec_truth = vec_truth * scale;
    
            % --- Plot vectors
            quiver3(el_pos(1), el_pos(2), el_pos(3), ...
                vec_pred(1), vec_pred(2), vec_pred(3), ...
                0, 'r', 'LineWidth', 1.5);
    
            quiver3(el_pos(1), el_pos(2), el_pos(3), ...
                vec_truth(1), vec_truth(2), vec_truth(3), ...
                0, 'g--', 'LineWidth', 1.2);
    
        end
    
        view(0, 0);
        xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
        title('Element Orientation Validation (Red = Euler, Green = Ground Truth)')
        grid on;
    
        legend({'Element Position','Natural Focus','Focus','From Euler Angles','Ground Truth'})
    
        base_name = sprintf('sub-%03d_%s_transducer_element_orientation%s', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix);
    
        saveas(h, fullfile(parameters.io.debug_dir, [base_name '.fig']))
        saveas(h, fullfile(parameters.io.debug_dir, [base_name '.png']))
    
        close(h)
    end
end