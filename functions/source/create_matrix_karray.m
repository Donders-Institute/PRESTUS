function [karray, tr] = create_matrix_karray(kgrid, karray, parameters, tr, elem_pos_m, trans_pos, focus_pos)
% CREATE_MATRIX_KARRAY  Add matrix transducer elements to a kWaveArray object
%
% For each element in elem_pos_m, computes the rotation matrix that aligns
% the default element normal [0;0;1] with the vector from the element to
% the natural focus (defined by the radius of curvature along Z; for a
% Clover setup, the natural focus of the element's leaf). The per-element
% Euler angles (ZYX) are passed to kWaveArray.addRectElement,
% addDiscElement, or addBowlElement depending on tr.matrix.elem_shape.
% Phase delays for steering to the requested focus_pos are also computed
% and stored in tr.matrix.elem_phase_rad, unless calibrated phases are
% given in tr.matrix.elem_phase_override_deg (one value per element), which
% then replace the geometric steering phases.
%
% Flat arrays (is_curved false) have every element facing +z. Rectangular
% elements are turned about their own normal by tr.matrix.elem_rot_deg
% (one value per element, from EXTRACT_ELEMENT_POS) when it is set.
%
% Steering correction: when tr.matrix.steering_correction holds a fitted
% map (FIT_STEERING_CORRECTION, usually from the transducer library) and
% tr.matrix.use_steering_correction is not false, every element amplitude
% is multiplied by k(target) = measured / modelled focal peak at the
% steering target (STEERING_CORRECTION_FACTOR). elem_p_amp then stays the
% nominal drive setting, and the simulated field is what the device
% delivers at that setting. k is stored in tr.matrix.steering_k.
%
% Use as:
%   [karray, tr] = create_matrix_karray(kgrid, karray, parameters, tr, elem_pos_m, trans_pos, focus_pos)
%
% Input:
%   kgrid      - kWaveGrid object
%   karray     - kWaveArray object to which elements are added
%   parameters - PRESTUS config; must contain medium_properties.water.sound_speed [m/s]
%                and the debug settings (DEBUG_ON)
%   tr         - transducer geometry: tr.matrix.elem_shape, tr.matrix.curv_radius_mm [mm],
%                tr.matrix.elem_height_mm [mm], tr.matrix.elem_width_mm [mm],
%                tr.matrix.elem_p_amp, tr.freq_hz [Hz]; optional
%                tr.matrix.elem_diameter_mm [mm] (disc/bowl aperture, see
%                MATRIX_ELEM_DIAMETER) and tr.matrix.elem_phase_override_deg [deg]
%   elem_pos_m - [3xN] element positions [m]
%   trans_pos  - [1x3] transducer reference position (grid indices)
%   focus_pos  - [1x3] desired focus position (grid indices)
%
% Output:
%   karray - kWaveArray with all elements added
%   tr     - updated with tr.matrix.elem_phase_rad, tr.matrix.elem_phase_deg,
%            tr.matrix.elem_n, tr.matrix.elem_p_amp, tr.matrix.steering_k
%
% See also: CREATE_CLOVER_ARRAY, SOURCE_CREATE, MATRIX_ELEMENT_AIM,
%           STEERING_CORRECTION_FACTOR

arguments
    kgrid      (1,1)
    karray     (1,1)
    parameters (1,1) struct
    tr         (1,1) struct
    elem_pos_m (3,:) {mustBeNumeric}
    trans_pos  (1,3) {mustBeNumeric}
    focus_pos  (1,3) {mustBeNumeric}
end

    % Convert positions from kgrid indices to meters
    trans_pos_m = [kgrid.x_vec(trans_pos(1)), ...
                   kgrid.y_vec(trans_pos(2)), ...
                   kgrid.z_vec(trans_pos(3))]';
               
    focus_pos_m = [kgrid.x_vec(focus_pos(1)), ...
                   kgrid.y_vec(focus_pos(2)), ...
                   kgrid.z_vec(focus_pos(3))]';

    % Assumes transducer beam axis is +Z. For transducers oriented along other
    % axes (e.g. a second transducer firing in -Z in a bilateral setup) element
    % orientations and phase delays will be wrong. Annular transducers are
    % unaffected (this function is only called for matrix type).
    % Elements face the natural focus (curved) or along +z (flat)
    aim_pos_m = matrix_element_aim(tr.matrix, elem_pos_m, trans_pos_m);
    natural_focus_pos_m = mean(aim_pos_m, 2);   % for the debug figure only

    % Apply Clover setup if requested
    if isequal(tr.matrix.is_clover_setup, true)    
        [elem_pos_m, aim_pos_m] = create_clover_array(parameters, tr.matrix, elem_pos_m, trans_pos_m, focus_pos_m);
    end
    
    tr.matrix.elem_n = size(elem_pos_m, 2);

    % Initialize source amplitudes (uniform)
    tr.matrix.elem_p_amp = tr.matrix.elem_p_amp(1) * ones(1, tr.matrix.elem_n);

    % Steering correction: measured / modelled focal peak at this target
    tr.matrix.steering_k = 1;
    use_sc = ~isfield(tr.matrix, 'use_steering_correction') || isempty(tr.matrix.use_steering_correction) ...
        || isequal(tr.matrix.use_steering_correction, true);
    if use_sc && isfield(tr.matrix, 'steering_correction') && isstruct(tr.matrix.steering_correction)
        target_mm = (focus_pos_m - trans_pos_m)' * 1e3;   % transducer frame, beam along +z
        tr.matrix.steering_k = steering_correction_factor(tr.matrix.steering_correction, target_mm);
        tr.matrix.elem_p_amp = tr.matrix.elem_p_amp * tr.matrix.steering_k;
        fprintf('Matrix transducer: steering correction k = %.3f at target [%.1f %.1f %.1f] mm.\n', ...
            tr.matrix.steering_k, target_mm);
    end

    % Per-element rotation about the element normal (rectangular elements)
    elem_rot_deg = zeros(1, tr.matrix.elem_n);
    if isfield(tr.matrix, 'elem_rot_deg') && ~isempty(tr.matrix.elem_rot_deg)
        if numel(tr.matrix.elem_rot_deg) ~= tr.matrix.elem_n
            error('create_matrix_karray:rotSize', ...
                'elem_rot_deg has %d values but the matrix array has %d elements.', ...
                numel(tr.matrix.elem_rot_deg), tr.matrix.elem_n);
        end
        elem_rot_deg = tr.matrix.elem_rot_deg(:)';
    end
    if is_flat_matrix(tr.matrix) && strcmpi(tr.matrix.elem_shape, 'bowl')
        error('create_matrix_karray:flatBowl', ...
            'A flat matrix array (is_curved false) needs elem_shape ''rect'' or ''disc''.');
    end

    % Wavelength and wavenumber for phase calculation
    lambda = parameters.medium_properties.water.sound_speed / tr.freq_hz;
    k = 2 * pi / lambda;

    % Initialize source phases, scaled vectors, tx, ty, tz
    elem_phase_rad = zeros(1, tr.matrix.elem_n);
    scaled_vectors = zeros(3, tr.matrix.elem_n);
    tx = zeros(1, tr.matrix.elem_n);
    ty = zeros(1, tr.matrix.elem_n);
    tz = zeros(1, tr.matrix.elem_n);

    % Aperture diameter of circular elements (same for every element)
    if ismember(lower(tr.matrix.elem_shape), {'disc', 'bowl'})
        elem_diameter_m = matrix_elem_diameter(tr.matrix);
    end

    % Loop over each element and add it to karray
    for ind = 1:tr.matrix.elem_n
        el_pos_m_i = elem_pos_m(:, ind);

        % Vector from element to natural focus to position elements to
        % natural focus
        aim_i = aim_pos_m(:, ind);
        vec_to_nat_focus = aim_i - el_pos_m_i;
        norm_vec = vec_to_nat_focus / norm(vec_to_nat_focus);
        scaled_vectors(:, ind) = vec_to_nat_focus;

        % Compute rotation matrix to align default normal [0;0;1] to vec_to_focus
        z0 = [0;0;1];
        v = cross(z0, norm_vec);
        s = norm(v);
        c = dot(z0, norm_vec);
        vx = [ 0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0 ];
        R = eye(3) + vx + vx^2*((1-c)/(s^2+eps));

        % Turn the element about its own normal (applied before aiming)
        if elem_rot_deg(ind) ~= 0
            R = R * [cosd(elem_rot_deg(ind)) -sind(elem_rot_deg(ind)) 0;
                     sind(elem_rot_deg(ind))  cosd(elem_rot_deg(ind)) 0;
                     0 0 1];
        end


        % Convert rotation matrix to ZYX Euler angles (extrinsic)
        yaw   = atan2d(R(2,1), R(1,1));  % around z
        pitch = atan2d(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));  % around y
        roll  = atan2d(R(3,2), R(3,3));  % around x

        tx(ind) = roll;
        ty(ind) = pitch;
        tz(ind) = yaw;

        % Determine element type (kWaveArray sizes are in m; width along x,
        % height along y, matching CONVERT_TO_ELEMENT_POS)
        elem_w_m = tr.matrix.elem_width_mm  * 1e-3;
        elem_h_m = tr.matrix.elem_height_mm * 1e-3;
        switch lower(tr.matrix.elem_shape)
            case 'rect'
                karray.addRectElement(el_pos_m_i, elem_w_m, elem_h_m, [roll, pitch, yaw]);

            case 'disc'
                % Disc with same area as rectangular element (or elem_diameter_mm)
                karray.addDiscElement(el_pos_m_i, elem_diameter_m, aim_i);

            case 'bowl'
                % Spherical cap with same area as rectangular element (or elem_diameter_mm)
                r_c = tr.matrix.curv_radius_mm / 1e3;
                karray.addBowlElement(el_pos_m_i, r_c, elem_diameter_m, aim_i);
        end

        % Calculate phase delay based on set focus
        distance = sqrt((el_pos_m_i(1) - focus_pos_m(1))^2 + (el_pos_m_i(2) - focus_pos_m(2))^2 + (el_pos_m_i(3) - focus_pos_m(3))^2);

        elem_phase_rad(ind) = mod(k * distance, 2 * pi);
    end

    % Calibrated phases (e.g. from hydrophone characterisation) replace the
    % geometric steering phases when provided
    if isfield(tr.matrix, 'elem_phase_override_deg') && ~isempty(tr.matrix.elem_phase_override_deg)
        override_deg = tr.matrix.elem_phase_override_deg(:)';
        if numel(override_deg) ~= tr.matrix.elem_n
            error('create_matrix_karray:phaseOverrideSize', ...
                'elem_phase_override_deg has %d values but the matrix array has %d elements.', ...
                numel(override_deg), tr.matrix.elem_n);
        end
        fprintf('Matrix transducer: using elem_phase_override_deg instead of geometric steering phases.\n');
        elem_phase_rad = mod(deg2rad(override_deg), 2 * pi);
    end

    tr.matrix.elem_phase_rad = elem_phase_rad;
    tr.matrix.elem_phase_deg = rad2deg(elem_phase_rad);

    % [DEBUG] visualize matrix element orientation
    if debug_on(parameters, 'source', 'figures')
        h = figure;
        hold on;
        axis equal;
    
        scatter3(elem_pos_m(1,:), elem_pos_m(2,:), elem_pos_m(3,:), 'b.');
        scatter3(natural_focus_pos_m(1,:), natural_focus_pos_m(2,:), natural_focus_pos_m(3,:), 'r.');
        plot3(focus_pos_m(1), focus_pos_m(2), focus_pos_m(3), 'go', 'MarkerFaceColor', 'g');
    
        for ind = 1:tr.matrix.elem_n
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
    
        xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
        title('Element Orientation Validation (Red = Euler, Green = Ground Truth)')
        grid on;
    
        legend({'Element Position','Natural Focus','Focus','From Euler Angles','Ground Truth'})
    
        base_name = sprintf('sub-%03d_%s_transducer_element_orientation%s', ...
            parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix);
    
        debug_save_figure(parameters, 'source', h, base_name, ...
            'Matrix elements: orientation from Euler angles (red) against ground truth (green)', 'Views', {'beam', 'side', 'oblique'});
    end
end