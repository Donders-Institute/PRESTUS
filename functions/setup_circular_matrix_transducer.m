function [transducer_pars, karray] = setup_circular_matrix_transducer(kgrid, trans_pos, transducer_pars, focus_pos, parameters, karray, is_rect_grid, is_curved, is_cut)
    if is_rect_grid
        % Step 1: Compute total width and height of the array
        tran_width = (transducer_pars.n_elem_row * transducer_pars.elem_width) + ...
                     ((transducer_pars.n_elem_row - 1) * transducer_pars.elem_spacing_width);
        
        tran_height = (transducer_pars.n_elem_col * transducer_pars.elem_height) + ...
                      ((transducer_pars.n_elem_col - 1) * transducer_pars.elem_spacing_height);
        
        % Step 2: Define x and y coordinates of element centers
        x_vec = linspace(-tran_width/2 + transducer_pars.elem_width/2, ...
                          tran_width/2 - transducer_pars.elem_width/2, transducer_pars.n_elem_row);
        
        y_vec = linspace(-tran_height/2 + transducer_pars.elem_height/2, ...
                          tran_height/2 - transducer_pars.elem_height/2, transducer_pars.n_elem_col);
        
        [X, Y] = meshgrid(x_vec, y_vec);
        
        % Step 4: Apply translation to position the transducer in the simulation grid
        X = X + kgrid.x_vec(trans_pos(1));
        Y = Y + kgrid.y_vec(trans_pos(2));
        Z = kgrid.z_vec(trans_pos(3)) * ones(size(X)); % For 2D array in XY plane 
    
        h = figure;
        scatter(X(:), Y(:), 60, 'filled'); % 60 is marker size
        axis equal;
        xlabel('X [m]');
        ylabel('Y [m]');
        title('Grid - Element Center Positions');
        grid on;
    
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_transducer_grid%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);
    
        dx = X(:) - kgrid.x_vec(trans_pos(1));
        dy = Y(:) - kgrid.y_vec(trans_pos(2));
        distances = sqrt(dx.^2 + dy.^2);
        
        % Logical mask for elements within aperture
        radius = transducer_pars.Elements_OD_mm(end) * 1e-3 / 2 ;
        mask = distances <= radius;
        
        % Select only elements within the circular aperture
        elem_pos = [X(:), Y(:), Z(:)];
        elem_pos = elem_pos(mask, :);  % Efficient filtering
    
        transducer_pars.n_elements = size(elem_pos, 1);
        transducer_pars.source_amp = transducer_pars.source_amp(1) * ones(1, transducer_pars.n_elements); 
        
        % Calculate inner and outer diameters for elements based on total diameter and count
        [id, od] = calc_elements_id_od_mm(transducer_pars.Elements_OD_mm(end), transducer_pars.n_elements);
        
        % Store element dimensions in standard location for visualization compatibility
        transducer_pars.Elements_ID_mm = id;
        transducer_pars.Elements_OD_mm = od;
    
        h = figure;
        scatter(elem_pos(:, 1)*1e3, elem_pos(:, 2)*1e3, 60, 'filled'); % 60 is marker size
        axis equal;
        xlabel('X [mm]');
        ylabel('Y [mm]');
        title('Circular Aperture - Element Center Positions');
        grid on;
    
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_transducer_grid%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);
        
        if is_curved
            % Compute Z coordinates based on spherical surface (assuming a spherical cap)
            ROC = transducer_pars.curv_radius_mm / 1000;
            R2 = ROC^2;
            
            % Compute the sag (z-offset from flat plane)
            sagitta_term = R2 - elem_pos(:, 1).^2 - elem_pos(:, 2).^2;
            
            % Safety check to prevent imaginary numbers
            if any(sagitta_term < 0)
                error('Some elements fall outside the spherical cap surface.');
            end
            
            % Apply curvature (spherical cap)
            elem_pos(:, 3) = elem_pos(:, 3) + ROC - sqrt(sagitta_term);
        end
    else
        elem_pos = fibonacci_spherical_cap(transducer_pars.n_elements, transducer_pars.curv_radius_mm, transducer_pars.Elements_OD_mm(end));
        elem_pos = elem_pos / 1e3; % convert from mm to m

        % Apply translation to position the transducer in the simulation grid
        elem_pos(: ,1) = elem_pos(: ,1) + kgrid.x_vec(trans_pos(1));
        elem_pos(: ,2) = elem_pos(: ,2) + kgrid.y_vec(trans_pos(2));
        elem_pos(: ,3) = elem_pos(: ,3) + kgrid.z_vec(trans_pos(3)); % For 2D array in XY plane
    end

    elem_pos = elem_pos';

    h = figure;
    scatter3(elem_pos(1, :)*1e3, elem_pos(2, :)*1e3, elem_pos(3, :)*1e3, 60, 'filled');
    xlabel('X [mm]');
    ylabel('Y [mm]');
    zlabel('Z [mm]');
    title('3D Circular Aperture');
    axis equal;
    grid on;

    output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_circular_transducer_grid%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    saveas(h, output_plot_filename, 'png')
    close(h);
    
    if is_cut
    % Backup before cutting
        elem_pos_full = elem_pos;

        % Remove two parts of module to be able to combine modules
        x = elem_pos(1, :);
        y = elem_pos(2, :);
        
        % Remove left lobe
        cut1 = x < kgrid.x_vec(trans_pos(1)) - 23.695/1000;
        
        % Remove lower wedge-shaped region
        cx = kgrid.x_vec(trans_pos(1));
        cy = kgrid.y_vec(trans_pos(2));

        chord_len = 27.915e-3;
        center_dist = 23.695e-3;
        start_dist = 27.5e-3;
        chord_start = [cx, cy + start_dist];
        
        angle = asin(center_dist/start_dist);
        
        x_dist = sin(angle) * chord_len;
        y_dist = cos(angle) * chord_len;

        chord_end = chord_start + [-x_dist, y_dist];
        
        % Define the line vector
        chord_vec = chord_end - chord_start;
        
        % For each point, determine if it's on the "keep" side of the line
        % The cross product of (point - chord_start) and chord_vec will be positive
        % if the point is on one side and negative if on the other side
        cut2 = zeros(size(x));
        for i = 1:length(x)
            point = [x(i), y(i)];
            point_vec = point + chord_start;
            % In 2D, the cross product is essentially: 
            cross_prod = point_vec(1)*chord_vec(2) - point_vec(2)*chord_vec(1);
            
            % Negative cross product means the point is on the "remove" side (bottom left)
            cut2(i) = (cross_prod < 0);
        end
        
        % Keep only the elements outside the clover cuts
        cut = cut1 | cut2;
        elem_pos_cut = elem_pos(:, ~cut);

        % === Plot Before and After ===
        h = figure;
        hold on;
        scatter3(elem_pos_full(1, :)*1e3, elem_pos_full(2, :)*1e3, elem_pos_full(3, :)*1e3, ...
                 60, [0.7 0.7 0.7], 'filled'); % full set in gray
        scatter3(elem_pos_cut(1, :)*1e3, elem_pos_cut(2, :)*1e3, elem_pos_cut(3, :)*1e3, ...
                 60, 'r', 'filled');           % cut set in red
        xlabel('X [mm]');
        ylabel('Y [mm]');
        zlabel('Z [mm]');
        legend({'Before Cut', 'After Cut'});
        title('Comparison: Full vs Cut 3D Curved Aperture');
        axis equal;
        grid on;
        hold off;
        
        output_plot_filename = fullfile(parameters.debug_dir, ...
            sprintf('sub-%03d_%s_compare_transducer_grid%s.png', ...
            parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png');
        close(h);

        elem_pos = elem_pos_cut;
        
    end

    if is_rect_grid
        circular_elements = elem_pos;
        n_circular = size(circular_elements, 2);
        n_selected = round(n_circular * transducer_pars.sparsity_factor);
    
        % Ensure we don't try to select more elements than available
        n_selected = min(n_selected, n_circular);
        
        % Random selection from circular elements
        if n_circular > 0
            selected_indices = randperm(n_circular, n_selected);
            elem_pos = circular_elements(:, selected_indices);
        else
            elem_pos = [];
            warning('No elements found within circular aperture!');
        end

            h = figure;
            scatter3(elem_pos(1, :)*1e3, elem_pos(2, :)*1e3, elem_pos(3, :)*1e3, 60, 'filled');
            xlabel('X [mm]');
            ylabel('Y [mm]');
            zlabel('Z [mm]');
            title('3D Circular Aperture');
            axis equal;
            grid on;
        
            output_plot_filename = fullfile(parameters.debug_dir, ...
                    sprintf('sub-%03d_%s_circular_transducer_sparse_grid%s.png', ...
                    parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
    end

    % Define focus in meters
    focus_pos_m = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))]';
    trans_pos_m = [kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))]';

    natural_focus_pos_m = trans_pos_m + [0, 0, transducer_pars.curv_radius_mm / 1000]';

    if transducer_pars.is_clover
        parent_mid_bowl_m = trans_pos_m;

        x0 = elem_pos(1, :)'*1e3;  % element x positions in mm
        y0 = elem_pos(2, :)'*1e3;
        z0 = elem_pos(3, :)'*1e3;
        
        % --- Parameters ---
        parent_mid_bowl_m(3) = trans_pos_m(3) + transducer_pars.ROC_bowl*1e-3; % middle of parent bowl
        theta_az = 2 * pi / transducer_pars.n_leaves;  % 120° between leaves

        % Compute radius of circle where centers lie to avoid overlap
        radius_circle_of_centers = transducer_pars.Elements_OD_mm(end) / (2 * transducer_pars.ROC_bowl * sin(theta_az / 2));
        
        % Compute elevation angle based on ROC
        elevation_angle = asin(radius_circle_of_centers);

        % --- Focus position ---
        parent_mid_bowl_mm = parent_mid_bowl_m * 1e3;

        % --- First transducer center ---
        % Place first transducer at desired elevation on ROC sphere (along X-axis azimuthally)
        x_c = transducer_pars.ROC_bowl * cos(elevation_angle);
        y_c = 0;
        z_c = transducer_pars.ROC_bowl * sin(elevation_angle);
        center0 = parent_mid_bowl_mm + [x_c; y_c; z_c];
        
        % --- Positions relative to center ---
        positions_local = [x0, y0, z0] - center0';
        
        % --- Storage ---
        x_all = [];
        y_all = [];
        z_all = [];
        
        % --- Plot ---
        h = figure; hold on; axis equal;
        colors = lines(transducer_pars.n_leaves);
        
        legend_labels = cell(1, length(transducer_pars.n_leaves)); % Pre-allocate legend labels
        for i = 0:transducer_pars.n_leaves-1
            % --- Step 1: rotate around Z to spread evenly
            angle_z = i * theta_az;
            Rz = [cos(angle_z), -sin(angle_z), 0;
                  sin(angle_z),  cos(angle_z), 0;
                  0,             0,            1];
        
            % --- Step 2: rotate around Y to tilt downward (toward focus)
            Ry = [cos(elevation_angle), 0, sin(elevation_angle);
                  0,                    1, 0;
                 -sin(elevation_angle), 0, cos(elevation_angle)];
        
            % --- Combined rotation
            R_total = Rz * Ry;
        
            % --- Rotate center and elements
            % --- Compute rotated center
            center_rot = (R_total * (center0 - parent_mid_bowl_mm)) + parent_mid_bowl_mm;
        
            % --- Rotate elements
            elems_rot = (R_total * positions_local')' + center_rot';

            % Residual function for sphere fitting
            residuals = @(c) sqrt(sum((elems_rot - c).^2, 2)) - transducer_pars.curv_radius_mm;
            
            % Initial guess: mean of points
            initial_guess = mean(elems_rot);
            
            % Least squares optimization
            center = lsqnonlin(residuals, initial_guess);
            
            % Compute apex (middle point on the bowl)
            mean_point = mean(elems_rot);
            direction = (mean_point - center) / norm(mean_point - center);
            apex = center + transducer_pars.curv_radius_mm * direction;
        
            % Vector from true center to focus
            vector_to_focus = parent_mid_bowl_mm' - apex;
            
            % Euclidean distance
            distance_to_focus = norm(vector_to_focus);
            
            % Create legend label for this transducer
            legend_labels{i+1} = sprintf('Transducer %d: Dist. to focus = %.2f mm', ...
                i+1, distance_to_focus);
        
            % --- Store
            x_all = [x_all; elems_rot(:,1)];
            y_all = [y_all; elems_rot(:,2)];
            z_all = [z_all; elems_rot(:,3)];
                
            % --- Plot
            scatter3(elems_rot(:,1), elems_rot(:,2), elems_rot(:,3), 15, colors(i+1,:), 'filled');
            h0(i+1) = plot3([apex(1), parent_mid_bowl_mm(1)], [apex(2), parent_mid_bowl_mm(2)], [apex(3), parent_mid_bowl_mm(3)], 'k--');               
        end
        
        % Plot focus and positions
        h1 = scatter3(parent_mid_bowl_m(1)*1e3, parent_mid_bowl_m(2)*1e3, parent_mid_bowl_m(3)*1e3, 100, 'r', 'filled');
        h2 = scatter3(focus_pos_m(1)*1e3, focus_pos_m(2)*1e3, focus_pos_m(3)*1e3, 100, 'b', 'filled');
        h3 = scatter3(trans_pos_m(1)*1e3, trans_pos_m(2)*1e3, trans_pos_m(3)*1e3, 100, 'g', 'filled');
        
        % Combine legend handles: all transducer bowls + fixed points
        legend_handles = [h0, h1, h2, h3];
        
        % Create combined legend labels: all transducers + fixed points
        fixed_labels = { "Middle of parent bowl, ROC " + num2str(transducer_pars.ROC_bowl, '%.2f') + " mm", "Focus", "Transducer middle" };
        legend(legend_handles, [legend_labels, fixed_labels]);
        
        xlabel('X [mm]');
        ylabel('Y [mm]');
        zlabel('Z [mm]');
        title(strcat('3D Clover Transducer Layout, ROC sub-array ', {' '}, num2str(transducer_pars.curv_radius_mm), {' '}, 'mm'));
        grid on;
        view([20 25 30]);

        output_plot_filename = fullfile(parameters.debug_dir, ...
        sprintf('sub-%03d_%s_clover_formation%s.png', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
        saveas(h, output_plot_filename, 'png')
        close(h);

        elem_pos = [x_all, y_all, z_all]'/1e3;
    end

    % Update the main variable to the potential cut and/or sparse version
    transducer_pars.n_elements = size(elem_pos, 2);
    transducer_pars.source_amp = transducer_pars.source_amp(1) * ones(1, transducer_pars.n_elements); 

    % Calculate inner and outer diameters for elements based on total diameter and count
    [id, od] = calc_elements_id_od_mm(transducer_pars.Elements_OD_mm(end), transducer_pars.n_elements);
    
    % Store element dimensions in standard location for visualization compatibility
    transducer_pars.Elements_ID_mm = id;
    transducer_pars.Elements_OD_mm = od;
   
    % Calculate wavelength and wavenumber for phase calculation in
    % water medium
    lambda = parameters.medium.water.sound_speed / transducer_pars.source_freq_hz; % [m]
    k = 2 * pi / lambda; % [rad/m]

    source_phase_rad = zeros(1, transducer_pars.n_elements);

    % Initialize tx, ty, tz
    scaled_vectors = zeros(3, size(elem_pos, 2));
    tx = zeros(1, size(elem_pos, 2));
    ty = zeros(1, size(elem_pos, 2));
    tz = zeros(1, size(elem_pos, 2));

    % Add elements to the array
    for ind = 1:transducer_pars.n_elements
        el_pos = elem_pos(:, ind);
        scaled_vec = natural_focus_pos_m - el_pos;  % vector from element to focus
        norm_vec = scaled_vec / norm(scaled_vec);            % normalize to unit vector

        % Convert direction vector to Euler angles
        % Assumes initial normal of element is along +z
        % Compute rotation that aligns [0; 0; 1] to vec
        z0 = [0; 0; 1];
        v = cross(z0, norm_vec);
        s = norm(v);
        c = dot(z0, norm_vec);
        vx = [  0    -v(3)  v(2);
               v(3)   0    -v(1);
              -v(2)  v(1)   0  ];
        R = eye(3) + vx + vx^2 * ((1 - c) / (s^2 + eps));  % rotation matrix
        
        % Convert rotation matrix to ZYX Euler angles (extrinsic)
        yaw   = atan2d(R(2,1), R(1,1));  % around z
        pitch = atan2d(-R(3,1), sqrt(R(3,2)^2 + R(3,3)^2));  % around y
        roll  = atan2d(R(3,2), R(3,3));  % around x
      
        tx(ind) = roll;
        ty(ind) = pitch;
        tz(ind) = yaw;

        if strcmp(transducer_pars.element, 'rect')
            karray.addRectElement(el_pos, transducer_pars.elem_height, transducer_pars.elem_width, [roll, pitch, yaw])
        elseif strcmp(transducer_pars.element, 'disc')
            % Diameter of disc with equal area as rectangular element
            diameter = (2 * transducer_pars.elem_height) / sqrt(pi); 
            karray.addDiscElement(el_pos, diameter, natural_focus_pos_m);
        elseif strcmp(transducer_pars.element, 'bowl')
            % Diameter of bowl with equal area as rectangular element
            r_c = transducer_pars.curv_radius_mm / 1e3;
            % Define the area to match
            A_target = transducer_pars.elem_height^2;
            
            % Define the function to solve: A_cap(a) - A_target = 0
            a_min = 0.001;           % Lower bound of aperture radius (in m)
            a_max = r_c * 0.999;     % Slightly less than full bowl radius
            area_diff = @(a) 2*pi*r_c*(r_c - sqrt(r_c^2 - a^2)) - A_target;

            % Check if the function changes sign in the interval
            if sign(area_diff(a_min)) == sign(area_diff(a_max))
                error('No sign change in interval: cannot find bowl diameter. Possibly A_target is out of bounds.');
            end
            
            % Solve for 'a' (base radius of spherical cap)
            a_solution = fzero(area_diff, [a_min, a_max]);  % base radius in m
            diameter = 2 * a_solution;  % diameter of bowl element
            karray.addBowlElement(el_pos, r_c, diameter, natural_focus_pos_m);
        end

        distance = sqrt((elem_pos(1, ind) - focus_pos_m(1))^2 + (elem_pos(2, ind) - focus_pos_m(2))^2 + (elem_pos(3, ind) - focus_pos_m(3))^2);
                    
        source_phase_rad(ind) = mod(k * distance, 2 * pi);

    end

    transducer_pars.source_phase_rad = source_phase_rad;
    transducer_pars.source_phase_deg = rad2deg(source_phase_rad);

    h = figure;
    hold on;
    axis equal;
    
    scatter3(elem_pos(1,:), elem_pos(2,:), elem_pos(3,:), 'b.');
    scatter3(natural_focus_pos_m(1,:), natural_focus_pos_m(2,:), natural_focus_pos_m(3,:), 'r.');
    plot3(focus_pos_m(1), focus_pos_m(2), focus_pos_m(3), 'go', 'MarkerFaceColor', 'g');
    
    for ind = 1:transducer_pars.n_elements
        el_pos = elem_pos(:, ind);
    
        % Ground truth direction (focus - element)
        vec_truth = focus_pos_m - el_pos;
        vec_truth = vec_truth / norm(vec_truth);
    
        % Get Euler angles (in degrees)
        roll  = tx(ind);
        pitch = ty(ind);
        yaw   = tz(ind);
    
        % Convert to radians
        r = deg2rad(roll);
        p = deg2rad(pitch);
        y = deg2rad(yaw);
    
        % Compute rotation matrix (extrinsic ZYX)
        Rz = [cos(y) -sin(y) 0; sin(y) cos(y) 0; 0 0 1];
        Ry = [cos(p) 0 sin(p); 0 1 0; -sin(p) 0 cos(p)];
        Rx = [1 0 0; 0 cos(r) -sin(r); 0 sin(r) cos(r)];
    
        R = Rz * Ry * Rx;
    
        % Rotate default normal to get predicted direction
        vec_rot = R * [0; 0; 1];

        scale_factor = norm(focus_pos_m - el_pos);

        % Scale the vectors for visualization
        vec_rot_scaled = vec_rot * scale_factor;
        vec_truth_scaled = vec_truth * scale_factor;
    
        % Plot predicted orientation (from Euler angles)
        quiver3(el_pos(1), el_pos(2), el_pos(3), ...
                vec_rot_scaled(1), vec_rot_scaled(2), vec_rot_scaled(3), 0, 'r', 'LineWidth', 1.5);
    
        % Plot ground truth orientation (focus - element)
        quiver3(el_pos(1), el_pos(2), el_pos(3), ...
                vec_truth_scaled(1), vec_truth_scaled(2), vec_truth_scaled(3), 0, 'g--', 'LineWidth', 1.2);
    end
    
    view(90, 0);
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title('Validation of Element Orientations: Red = Computed, Green = Ground Truth');
    legend('Element Position', 'Focus', 'From Euler Angles', 'Ground Truth');

    output_plot_filename = fullfile(parameters.debug_dir, ...
    sprintf('sub-%03d_%s_transducer_element_orientation%s.fig', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    saveas(h, output_plot_filename, 'fig')
    
    output_plot_filename = fullfile(parameters.debug_dir, ...
    sprintf('sub-%03d_%s_transducer_element_orientation%s.png', ...
        parameters.subject_id, parameters.simulation_medium, parameters.results_filename_affix));
    saveas(h, output_plot_filename, 'png')

    close(h);
end