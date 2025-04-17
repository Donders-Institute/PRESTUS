function show_3d_head(segmented_img, target_xyz, trans_xyz, parameters, pixel_size, coord_mesh_xyz, crop_at_target, view_angle, open_figure)

% SHOW_3D_HEAD Visualizes a 3D segmented brain image with transducer and target positions.
%
% This function generates a 3D visualization of a segmented brain image. It includes:
%   - Skin surface rendering.
%   - Visualization of GM (gray matter), WM (white matter), CSF (cerebrospinal fluid), and bone.
%   - Transducer placement and exit plane visualization.
%   - Target location as a sphere.
%
% Input:
%   segmented_img    - [Nx x Ny x Nz] matrix representing the segmented brain image.
%   target_xyz       - [1x3] array specifying the target position in voxel coordinates.
%   trans_xyz        - [1x3] array specifying the transducer position in voxel coordinates.
%   parameters       - Struct containing transducer parameters (e.g., curvature radius, element diameters).
%   pixel_size       - Scalar specifying the voxel size in mm.
%   coord_mesh_xyz   - [Mx x 3] array representing the voxel coordinates of the computational grid.
%   crop_at_target   - [1x3] array specifying cropping directions relative to the target (default: [0,0,0]).
%                      * 1: Crop along positive direction of axis.
%                      * -1: Crop along negative direction of axis.
%                      * 0: No cropping along axis.
%   view_angle       - [1x2] array specifying azimuth and elevation angles for the 3D view (default: [0,0]).
%   open_figure      - Boolean flag to open a new figure for visualization (default: true).
%
% Output:
%   None. The function generates a 3D visualization.

    arguments
        segmented_img (:,:,:) 
        target_xyz (1,3)
        trans_xyz (1,3)
        parameters struct
        pixel_size (1,1)
        coord_mesh_xyz (:,3) 
        crop_at_target (1,3) = [0,0,0]
        view_angle (1,2) = [0,0]
        open_figure (1,1) = 1
    end

    %% Open figure if required
    if open_figure
        figure;
    end
    colormap(gray(80));

    %% Plot skin surface
    segmented_img_to_plot = segmented_img(1:2:end, 1:2:end, 1:2:end); % Downsample for visualization

    % Apply cropping if specified
    origin_shift = [0 0 0];
    if any(crop_at_target)
        original_size = size(segmented_img_to_plot);
        for i = 1:3
            cropping_indices = {':', ':', ':'};
            if crop_at_target(i) == 1
                cropping_indices{i} = floor(target_xyz(i)/2):size(segmented_img_to_plot,i);
            elseif crop_at_target(i) == -1
                cropping_indices{i} = 1:(floor(target_xyz(i)/2)+1);
                origin_shift(i) = floor(target_xyz(i)/2);
            else
                continue;
            end
            S = struct;
            S.type = '()';
            S.subs = cropping_indices;
            segmented_img_to_plot = subsasgn(segmented_img_to_plot, S, []);
        end
        origin_shift = origin_shift([2, 1, 3]);
    end

    % Smooth and render skin surface
    if gpuDeviceCount > 0 
        Ds = smooth3(gpuArray(double(segmented_img_to_plot > 0)));
    else 
        Ds = smooth3(double(segmented_img_to_plot > 0));
    end

    skin_isosurface = isosurface(Ds, 0.5);
    skin_isosurface.vertices = skin_isosurface.vertices - origin_shift;
    
    hiso = patch(skin_isosurface,...
       'FaceColor', [1,.75,.65],...
       'EdgeColor', 'none',...
       'facealpha', 0.8);
    
    isonormals(Ds,hiso);
    hiso.SpecularColorReflectance = 0;
    hiso.SpecularExponent = 50;

    %% Plot GM, WM, CSF, and bone surfaces
    for_caps = segmented_img_to_plot;
    for_caps(for_caps > 4) = 0; % Only GM, WM, CSF, bone
    
    end_slices = isocaps(for_caps, 4);
    %end_slices.vertices = end_slices.vertices - origin_shift;
    
    hcap = patch(end_slices,...
       'FaceColor', 'interp',...
       'EdgeColor', 'none',...
       'facealpha', 0.8);
    
    lighting gouraud;
    hcap.AmbientStrength = 0.6;

    %% Plot transducer slice at exit plane (if not cropped)
    if ~any(crop_at_target)

        if parameters.unique_tran_design == 1
            max_od_mm = max(parameters.transducer.Elements_OD_mm);
            max_od_grid = max_od_mm / pixel_size;
            
            norm_vector = (trans_xyz - target_xyz) / norm(trans_xyz - target_xyz);
    
            % Location of geometric focus
            geom_focus_xyz = trans_xyz - norm_vector * (parameters.transducer.curv_radius_mm) / pixel_size;
    
            % Distance from geometric focus to exit plane
            dist_gf_to_ep_mm = 0.5 * sqrt(4 * parameters.transducer.curv_radius_mm^2 - max_od_mm^2);
    
            % Coordinates of exit plane
            ex_plane_xyz = geom_focus_xyz + norm_vector * dist_gf_to_ep_mm / pixel_size;
    
            % Calculate distance metric and identify exit plane disk
            d = sum(norm_vector .* ex_plane_xyz);
            orth_plane_disk = find(abs(sum(coord_mesh_xyz .* norm_vector, 2) - d) < 0.5);
            orth_plane_disk = orth_plane_disk(sqrt(sum((coord_mesh_xyz(orth_plane_disk,:) - ex_plane_xyz).^2,2)) < max_od_grid / 2);
            orth_plane_3d = zeros(size(segmented_img));
            orth_plane_3d(orth_plane_disk) = 1;        
            orth_plane_3d_surf = isosurface(smooth3(orth_plane_3d(1:2:end,1:2:end,1:2:end)));
            orth_plane_3d_surf.vertices = orth_plane_3d_surf.vertices - origin_shift;

            patch(orth_plane_3d_surf,...
               'FaceColor', 'blue',...
               'EdgeColor', 'none');
  
        elseif parameters.unique_tran_design == 2
                % Initialize computational grids for the transducer mask and source label matrix
                transducer_mask = zeros(size(segmented_img)); % Binary mask representing active transducer regions

                % Define geometry
                % Create coordinate grids for the transducer elements (not the simulation grid)
                [X, Y] = meshgrid(1:parameters.transducer.n_elem_col, 1:parameters.transducer.n_elem_row);
                
                % Define circle parameters for the mask
                centerX = (parameters.transducer.n_elem_col + 1) / 2;  % X center of the circle
                centerY = (parameters.transducer.n_elem_row + 1) / 2;  % Y center of the circle

                % Create circular mask using the equation of a circle
                mask = (X - centerX).^2 + (Y - centerY).^2 <= parameters.transducer.modular_radius^2;
                mask(:,1) = [];
        
                % Create a symmetric mask by flipping and concatenating
                totalMask = [fliplr(mask), mask];
                totalMaskCols = size(totalMask, 2);  % Get the number of columns in the new mask     
            
                % Calculate the updated center X position for the totalMask
                totalMaskCenterX = (totalMaskCols + 1) / 2;
                
                % Convert physical dimensions to grid points
                elem_width_grid = round(parameters.transducer.elem_width / parameters.grid_step_m);  
                elem_height_grid = round(parameters.transducer.elem_height / parameters.grid_step_m);
                elem_spacing_width_grid = round(parameters.transducer.elem_spacing_width / parameters.grid_step_m);
                elem_spacing_height_grid = round(parameters.transducer.elem_spacing_height / parameters.grid_step_m);
                
                
                % Process each element in the grid using the total mask dimensions
                for row = 1:parameters.transducer.n_elem_row
                    for col = 1:totalMaskCols
                        % Check if this element is active according to the mask
                        if totalMask(row, col)
                            
                            % Calculate grid position for this element
                            grid_x = trans_xyz(1) + (col - totalMaskCenterX) * (elem_width_grid + elem_spacing_width_grid);
                            grid_y = trans_xyz(2) + (row - centerY) * (elem_height_grid + elem_spacing_height_grid);
                            grid_z = trans_xyz(3);
                            
                            % Convert to integer grid positions
                            x_start = round(grid_x - elem_width_grid/2);
                            y_start = round(grid_y - elem_height_grid/2);
                            z_start = round(grid_z); % Assuming 1-point thickness in z direction
                            
                            % Create element using rectangular shape
                            for x = x_start:x_start+elem_width_grid
                                for y = y_start:y_start+elem_height_grid
                                    % Check if position is within grid bounds
                                    if x >= 1 && x <= size(transducer_mask, 1) && y >= 1 && y <= size(transducer_mask, 2) && z_start >= 1 && z_start <= size(transducer_mask, 3)
                                        transducer_mask(x, y, z_start) = 1;
                                    end
                                end
                            end
                        end
                    end
                end
                
                % Add visualization code similar to design 1
                % Create an isosurface from the transducer mask
                transducer_surf = isosurface(smooth3(transducer_mask(1:2:end,1:2:end,1:2:end)));
                
                % Apply the same origin shift as in design 1
                transducer_surf.vertices = transducer_surf.vertices - origin_shift;
                
                % Display the isosurface with patch
                patch(transducer_surf,...
                   'FaceColor', 'blue',...
                   'EdgeColor', 'none');
        else
            fprintf("Transducer design number %i hasn't been implemented. \n", parameters.unique_tran_design)
        end
    end

    %% Plot target sphere
    sphere_3d = zeros(size(segmented_img));
    sphere_3d(pdist2(coord_mesh_xyz,target_xyz)<3) = 1;
    sphere_surf = isosurface(smooth3(sphere_3d(1:2:end,1:2:end,1:2:end)));
    sphere_surf.vertices = sphere_surf.vertices - origin_shift;
    patch(sphere_surf,'FaceColor','red','EdgeColor','red','facealpha',0.8);

    %% Adjust viewing angles and labels
    if ~isempty(view_angle)
        view(view_angle);
    end
    xlabel('X');
    ylabel('Y');
    set(gca,'YDir','reverse');
    zlabel('Z');
    if any(crop_at_target)
        original_size = original_size([2, 1, 3]);
        
        if any(origin_shift)
            axis(reshape([-origin_shift ; -origin_shift + original_size], [1 6]));
        end
    end
 
    if view_angle(1) < 0
        camlight('right');
    else
        camlight('left');
    end
    view(view_angle);

    axis equal tight;

end
