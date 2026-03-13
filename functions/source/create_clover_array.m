          parent_mid_bowl_m = trans_pos_m;
        
                x0 = elem_pos_m(1, :)'*1e3;  % element x positions in mm
                y0 = elem_pos_m(2, :)'*1e3;
                z0 = elem_pos_m(3, :)'*1e3;
                
                % --- Parameters ---
                parent_mid_bowl_m(3) = trans_pos_m(3) + transducer_pars.ROC_bowl*1e-3; % middle of parent bowl
                theta_az = 2 * pi / 3;  % 120° between leaves
        
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
        
                elem_pos_m = [x_all, y_all, z_all]'/1e3;