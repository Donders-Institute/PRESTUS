function [transducer_pars, karray] = setup_stanford_transducer(kgrid, trans_pos, transducer_pars, focus_pos, parameters, karray)
            tran_info = readtable('TR64pcd1d65r75_Stanford_fromED.xlsx');
    
            phys_positions = table2array(tran_info(1:64, 2:4));

            % Convert from mm to m, and recenter to [0,0,0] by subtracting the actual center
            phys_positions_m = (phys_positions - [0, 0, 75]) / 1000;
            phys_positions_m(:,3) = -phys_positions_m(:,3);  % flip Z curvature
           
            elem_pos = phys_positions_m + [kgrid.x_vec(trans_pos(1)), kgrid.y_vec(trans_pos(2)), kgrid.z_vec(trans_pos(3))];
            elem_pos = elem_pos';

            % Define focus in meters
            focus_pos_m = [kgrid.x_vec(focus_pos(1)), kgrid.y_vec(focus_pos(2)), kgrid.z_vec(focus_pos(3))]';

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
                scaled_vec = focus_pos_m - el_pos;  % vector from element to focus
                norm_vec = scaled_vec / norm(scaled_vec);            % normalize to unit vector

                scaled_vectors(:, ind) = scaled_vec;

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

                karray.addRectElement(el_pos, transducer_pars.elem_height, transducer_pars.elem_width, [roll, pitch, yaw])

                distance = sqrt((el_pos(1) - focus_pos_m(1))^2 + (el_pos(2) - focus_pos_m(2))^2 + (el_pos(3) - focus_pos_m(3))^2);
                            
                source_phase_rad(ind) = mod(k * distance, 2 * pi);

            end

            transducer_pars.source_phase_rad = source_phase_rad;
            transducer_pars.source_phase_deg = rad2deg(source_phase_rad);

            h = figure;
            hold on;
            axis equal;
            
            scatter3(elem_pos(1,:), elem_pos(2,:), elem_pos(3,:), 'b.');
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
            
            view(0, 0);
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