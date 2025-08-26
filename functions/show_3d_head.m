function show_3d_head(segmented_img, target_xyz, trans_xyz, parameters, pixel_size, coord_mesh_xyz, crop_at_target, view_angle, open_figure)
% SHOW_3D_HEAD Visualizes a 3D segmented brain image with multiple transducer and target positions.
%
% Visualizes a 3D segmented brain, with:
%   - Skin, tissue surfaces
%   - One or more transducer positions and target spheres (in unique colors)
%   - Transducer exit-plane disks, using an anti-striped direct-downsampling approach
%
% INPUTS:
%   segmented_img    - [Nx x Ny x Nz] volume, segmented tissue labels
%   target_xyz       - [K x 3] array, each row is a target position (voxel indices)
%   trans_xyz        - [K x 3] array, each row is a transducer position (voxel indices)
%   parameters       - struct, must include transducer geometry fields
%   pixel_size       - scalar, voxel size in mm
%   coord_mesh_xyz   - [M x 3] grid of all voxels (voxel indices)
%   crop_at_target   - [1x3], crop along axes (default: [0,0,0])
%   view_angle       - [1x2], [az, el] for MATLAB view (default: [0,0])
%   open_figure      - bool, whether to create new figure (default: true)
%
% OUTPUT: None—shows a 3D plot with all transducers/targets overlaid.
%
% NOTE:
% - All positions must be in voxel (image index) space. Convert before calling!
% - No stripes: exit-plane masks are generated in the downsampled space.
% - Up to K pairs of transducers/targets supported, each with its unique color.

    arguments
        segmented_img (:,:,:) 
        target_xyz (:,3)
        trans_xyz (:,3)
        parameters struct
        pixel_size (1,1)
        coord_mesh_xyz (:,3) 
        crop_at_target (1,3) = [0,0,0]
        view_angle (1,2) = [0,0]
        open_figure (1,1) = 1
    end

    % Ensure both target_xyz and trans_xyz are of same number of points
    nPairs = size(target_xyz,1);
    if size(trans_xyz,1) ~= nPairs
        error('target_xyz and trans_xyz must have the same number of rows.');
    end

    %% Color palette for pairs
    color_list = lines(max(7, nPairs)); % MATLAB's distinct color set

    %% Open figure if requested
    if open_figure
        figure;
    end
    colormap(gray(80));

    %% Downsample anatomy for visualization
    ds_factor = 2;
    segmented_img_ds = segmented_img(1:ds_factor:end, 1:ds_factor:end, 1:ds_factor:end);

    % Crop if requested
    origin_shift = [0 0 0];
    if any(crop_at_target)
        original_size = size(segmented_img_ds);
        for i = 1:3
            cropping_indices = {':', ':', ':'};
            if crop_at_target(i) == 1
                cropping_indices{i} = floor(target_xyz(1,i)/ds_factor):size(segmented_img_ds,i);
            elseif crop_at_target(i) == -1
                cropping_indices{i} = 1:(floor(target_xyz(1,i)/ds_factor) + 1);
                origin_shift(i) = floor(target_xyz(1,i)/ds_factor);
            else
                continue;
            end
            % Crop using direct indexing, not subsasgn
            segmented_img_ds = segmented_img_ds(cropping_indices{:});
        end
        origin_shift = origin_shift([2, 1, 3]); % (xyz → yxz legacy)
    end

    % Smooth and render skin
    if gpuDeviceCount > 0 
        Ds = smooth3(gpuArray(double(segmented_img_ds > 0)));
    else 
        Ds = smooth3(double(segmented_img_ds > 0));
    end
    skin_iso = isosurface(Ds, 0.5);
    skin_iso.vertices = skin_iso.vertices - origin_shift;
    hiso = patch(skin_iso, 'FaceColor', [1,.75,.65], 'EdgeColor', 'none', 'facealpha', 0.8);

    isonormals(Ds, hiso);
    hiso.SpecularColorReflectance = 0;
    hiso.SpecularExponent = 50;

    % Plot caps for GM/WM/CSF/bone
    for_caps = segmented_img_ds;
    for_caps(for_caps > 4) = 0;
    cap_surf = isocaps(for_caps*10,4);
    cap_surf.vertices = cap_surf.vertices - origin_shift;
    hcap = patch(cap_surf, 'FaceColor', 'interp', 'EdgeColor', 'none', 'facealpha', 0.8);
    lighting gouraud;
    hcap.AmbientStrength = 0.6;

    %% Plot for each transducer/target pair (K pairs)
    for k = 1:nPairs
        thisTrans = trans_xyz(k,:);
        thisTarg = target_xyz(k,:);
        c = color_list(k,:);
        % Plot transducer exit plane if not cropped
        if ~any(crop_at_target)
            % All shapes in downsampled space
            max_od_mm = max(parameters.transducer.Elements_OD_mm);
            max_od_grid = max_od_mm / pixel_size;
            norm_vec = (thisTrans - thisTarg) / norm(thisTrans - thisTarg);

            % Geometric focus in grid space
            geom_focus = thisTrans - norm_vec * (parameters.transducer.curv_radius_mm) / pixel_size;
            dist_gf_to_ep_mm = 0.5 * sqrt(4 * parameters.transducer.curv_radius_mm^2 - max_od_mm^2);
            ex_plane = geom_focus + norm_vec * dist_gf_to_ep_mm / pixel_size;

            % Now get full-res (voxel) coordinates, then downsample
            d = sum(norm_vec .* ex_plane);
            % This finds full-res indices close to the disk
            dists = abs(sum(coord_mesh_xyz .* norm_vec, 2) - d);
            radius2 = sum((coord_mesh_xyz - ex_plane).^2,2);
            orth_plane_idx = find(dists < 0.5 & sqrt(radius2) < max_od_grid / 2);

            % Downsample coordinates to match image
            [x_full, y_full, z_full] = ind2sub(size(segmented_img), orth_plane_idx);
            x_ds = round((x_full + (ds_factor-1))/ds_factor);
            y_ds = round((y_full + (ds_factor-1))/ds_factor);
            z_ds = round((z_full + (ds_factor-1))/ds_factor);

            % Remove out-of-bounds indices
            inside_vol = x_ds >= 1 & y_ds >= 1 & z_ds >= 1 & ...
                x_ds <= size(segmented_img_ds,1) & ...
                y_ds <= size(segmented_img_ds,2) & ...
                z_ds <= size(segmented_img_ds,3);

            x_ds = x_ds(inside_vol); y_ds = y_ds(inside_vol); z_ds = z_ds(inside_vol);

            orth_plane_mask = zeros(size(segmented_img_ds));
            orth_idx_ds = sub2ind(size(orth_plane_mask), x_ds, y_ds, z_ds);
            orth_plane_mask(orth_idx_ds) = 1;

            ep_surf = isosurface(smooth3(orth_plane_mask));
            ep_surf.vertices = ep_surf.vertices - origin_shift;
            patch(ep_surf, 'FaceColor', c, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
        end

        % Plot target sphere
        sphere_mask = zeros(size(segmented_img_ds));
        % Downsampled target for rounded indexing
        for k = 1:nPairs
            targ_ds = round((target_xyz(k,:) + (ds_factor-1))/ds_factor);
            % Get valid points within mask
            valid_mask = ...
               coord_mesh_xyz(:,1) >= 1 & coord_mesh_xyz(:,1) <= size(sphere_mask,1) & ...
               coord_mesh_xyz(:,2) >= 1 & coord_mesh_xyz(:,2) <= size(sphere_mask,2) & ...
               coord_mesh_xyz(:,3) >= 1 & coord_mesh_xyz(:,3) <= size(sphere_mask,3);
            mesh_valid = coord_mesh_xyz(valid_mask, :);
        
            dists = sqrt(sum((mesh_valid - targ_ds).^2, 2));
            flag = dists < 3;
        
            idx_x = mesh_valid(flag,1);
            idx_y = mesh_valid(flag,2);
            idx_z = mesh_valid(flag,3);
        
            lin_idx = sub2ind(size(sphere_mask), idx_x, idx_y, idx_z);
            sphere_mask(lin_idx) = 1;
        
            sphere_surf = isosurface(smooth3(sphere_mask));
            sphere_surf.vertices = sphere_surf.vertices - origin_shift;
            patch(sphere_surf, 'FaceColor', c, 'EdgeColor', c, 'facealpha', 0.9);
    
            % Optionally, plot the actual transducer voxel as a small sphere
            hold on
            [X,Y,Z] = sphere(8);
            surf(targ_ds(2) + 1.2*X, targ_ds(1) + 1.2*Y, targ_ds(3) + 1.2*Z, ...
                'FaceColor', c, 'EdgeColor', 'none', 'FaceAlpha', 1);
    
            trans_ds = round((thisTrans +(ds_factor-1)) / ds_factor);
            surf(trans_ds(2) + 1.2*X, trans_ds(1) + 1.2*Y, trans_ds(3) + 1.2*Z, ...
                'FaceColor', c*.7, 'EdgeColor', 'none', 'FaceAlpha', 1);
        end
        
    end

    %% Adjust viewing angles and axis labels
    if ~isempty(view_angle)
        view(view_angle);
    end
    xlabel('X');
    ylabel('Y');
    set(gca,'YDir','reverse'); % For radiological convention
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

end
