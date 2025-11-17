close all; clear;

coord_dirs = ["X", "Y", "Z"];
steps = ["-20", "-10", "ref", "+10", "+20"]; % ["-40", "-30", "-20", "-10", "ref", "+10", "+20", "+30", "+40"];
%max_isppa_for_scaling = 90;
medium = 'layered';

path = 'C:\Users\marge\OneDrive - Radboud Universiteit\Documenten\Software\repositories\prestus\Research\New transducer development\CLOVER\simulations_after_toolkit\p75_s75_f75\steering\Ernie_noRot\';
fig_title = 'Clover 2 leaves p75 s75 64e bowl fibo ernie';
pre_name = 'sub-001_layered_results_p75_s75_f75_2_leaves_clover_';
post_name = '_noRot_64e_bowl_fibo_ernie.mat';

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 

for d = 1:length(coord_dirs)
    coord_dir = coord_dirs(d);

    for s = 1:length(steps)
        step = steps(s);

        if strcmp(step, "ref")
            filename = strcat(path, pre_name, step, post_name);
            dir_step = step;
        else
            filename = strcat(path, pre_name, coord_dir, step, post_name);
            dir_step = strcat(coord_dir, step);
        end
        
        load(filename, 'parameters', 'sensor_data', 'kwave_medium');
        
        acoustic_pressure = gather(sensor_data.p_max_all);
        acoustic_isppa = acoustic_pressure.^2./...
            (2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-6;
        max_Isppa = max(acoustic_isppa(:));
        
        slice_i_x = squeeze(acoustic_isppa(parameters.focus_pos_grid(1), :,:));
        slice_i_y = squeeze(acoustic_isppa(:, parameters.focus_pos_grid(2),:));
        slice_i_z = squeeze(acoustic_isppa(:, :, parameters.focus_pos_grid(3)));
        
        % Creates the foundation for a mask before the exit plane to calculate max values outside of it
        comp_grid_size = size(sensor_data.p_max_all);
        after_exit_plane_mask = ones(comp_grid_size);
        bowl_depth_grid = round((parameters.transducer.curv_radius_mm-...
            parameters.transducer.dist_to_plane_mm)/parameters.grid_step_mm);
        % Places the exit plane mask in the grid, adjusted to the amount of dimensions
        if parameters.n_sim_dims == 3
            if parameters.transducer.pos_grid(3) > comp_grid_size(3)/2
                after_exit_plane_mask(:,:,(parameters.transducer.pos_grid(parameters.n_sim_dims)-...
                    bowl_depth_grid):end) = 0;
            else
                after_exit_plane_mask(:,:,1:(parameters.transducer.pos_grid(parameters.n_sim_dims)+...
                    bowl_depth_grid)) = 0;
            end
        else
            if parameters.transducer.pos_grid(2) > comp_grid_size(2)/2
                after_exit_plane_mask(:,(parameters.transducer.pos_grid(parameters.n_sim_dims)-...
                    bowl_depth_grid):end) = 0;
            else
                after_exit_plane_mask(:,1:(parameters.transducer.pos_grid(parameters.n_sim_dims)+...
                    bowl_depth_grid)) = 0;
            end
        end
        
        % Calculates the X, Y and Z coordinates of the max. intensity
        [max_Isppa_after_exit_plane, Ix_eplane, Iy_eplane, Iz_eplane] = ...
            masked_max_3d(acoustic_isppa, after_exit_plane_mask);
        
        % Combines these coordinates into a point of max. intensity in the grid
        if parameters.n_sim_dims==3
            max_isppa_eplane_pos = [Ix_eplane, Iy_eplane, Iz_eplane];
        else 
            max_isppa_eplane_pos = [Ix_eplane, Iy_eplane];
        end
        
        % Calculates the average Isppa within a circle around the target
        real_focal_distance = norm(max_isppa_eplane_pos-parameters.transducer.pos_grid)*...
            parameters.grid_step_mm;
        distance_target_real_maximum = norm(max_isppa_eplane_pos-parameters.focus_pos_grid)*...
            parameters.grid_step_mm;
        avg_radius = round(parameters.focus_area_radius/parameters.grid_step_mm); %grid
        avg_isppa_around_target = acoustic_isppa(...
            (parameters.focus_pos_grid(1)-avg_radius):(parameters.focus_pos_grid(1)+avg_radius),...
            (parameters.focus_pos_grid(2)-avg_radius):(parameters.focus_pos_grid(2)+avg_radius),...
            (parameters.focus_pos_grid(3)-avg_radius):(parameters.focus_pos_grid(3)+avg_radius));
        avg_isppa_around_target = mean(avg_isppa_around_target(:));
        
        % Reports the Isppa within the original stimulation target
        isppa_at_target = acoustic_isppa(parameters.focus_pos_grid(1),parameters.focus_pos_grid(2),parameters.focus_pos_grid(3));
        
        if strcmp(medium, 'water')
            % Add annotation textbox outside the plots
            info_str = sprintf(['Max ISPPA: %.2f W/cm^2\n', ...
                                'ISPPA at target: %.2f W/cm^2\n', ...
                                'Avg ISPPA around target: %.2f W/cm^2\n', ...
                                'Amplitude per element: %.2e Pa\n', ...
                                'Focal gain max.: %.2f\n', ...
                                'Focal gain target.: %.2f'], ...
                                max_Isppa, isppa_at_target, avg_isppa_around_target, ...
                                amp_per_element, focal_gain_max, focal_gain_target);
        else
            % Creates a logical skull mask and register skull_ids
            labels = fieldnames(parameters.layer_labels);
            skull_i = find(strcmp(labels, 'skull_cortical'));
            trabecular_i = find(strcmp(labels, 'skull_trabecular'));
            all_skull_ids = [skull_i, trabecular_i];
            mask_skull = ismember(medium_masks,all_skull_ids);
            brain_i = find(strcmp(labels, 'brain'));
            mask_brain = ismember(medium_masks,brain_i);
            skin_i = find(strcmp(labels, 'skin'));
            mask_skin = ismember(medium_masks,skin_i);
            
            [max_Isppa_brain, Ix_brain, Iy_brain, Iz_brain] = masked_max_3d(acoustic_isppa, mask_brain);
            [min_Isppa_brain] = min(acoustic_isppa(mask_brain));
            half_max = acoustic_isppa >= max_Isppa_brain/2 & mask_brain;
            half_max_ISPPA_volume_brain = sum(half_max(:))*(parameters.grid_step_mm^3);
            
            [max_Isppa_skull, Ix_skull, Iy_skull, Iz_skull] = masked_max_3d(acoustic_isppa, mask_skull);
            
            [max_Isppa_skin, Ix_skin, Iy_skin, Iz_skin] = masked_max_3d(acoustic_isppa, mask_skin);

            max_pressure_brain = sqrt(max_Isppa_brain * 2 * kwave_medium.sound_speed(1, 1, 1) * kwave_medium.density(1, 1, 1) * 1e6);
            focal_gain_brain =  max_pressure_brain / amp_per_element;

            info_str = sprintf(['Max ISPPA: %.2f W/cm^2\n', ...
                                'ISPPA at target: %.2f W/cm^2\n', ...
                                'Avg ISPPA around target: %.2f W/cm^2\n', ...
                                'Max ISPPA skin: %.2f W/cm^2\n', ...
                                'Max ISPPA skull: %.2f W/cm^2\n', ...
                                'Max ISPPA brain: %.2f W/cm^2\n', ...
                                'Amplitude per element: %.2e Pa\n', ...
                                'Focal gain max.: %.2f\n', ...
                                'Focal gain target.: %.2f', ...
                                'Focal gain brain.: %.2f'], ...
                                max_Isppa, isppa_at_target, avg_isppa_around_target, ...
                                max_Isppa_skin, max_Isppa_skull, max_Isppa_brain,...
                                amp_per_element, focal_gain_max, focal_gain_target, ...
                                focal_gain_brain);
        
        end

        amp_per_element = parameters.transducer.source_amp(1);
        max_pressure = sqrt(max_Isppa * 2 * kwave_medium.sound_speed(1, 1, 1) * kwave_medium.density(1, 1, 1) * 1e6);
        focal_gain_max = max_pressure / amp_per_element;
        pressure_at_target = sqrt(isppa_at_target * 2 * kwave_medium.sound_speed(1, 1, 1) * kwave_medium.density(1, 1, 1) * 1e6);
        focal_gain_target =  pressure_at_target / amp_per_element;
        
        figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]); 
        
        tiledlayout(4,1, 'Padding', 'compact', 'TileSpacing', 'compact');
        
        axis image
        axis off;
        
        % First subplot
        ax1 = nexttile;
        imagesc(slice_i_x);
        colormap('viridis');
        title('Slice in X direction');
        xlabel('Y'); ylabel('Z');
        
        % Second subplot
        ax2 = nexttile;
        imagesc(slice_i_y);
        title('Slice in Y direction');
        xlabel('X'); ylabel('Z');
        
        % Third subplot
        ax3 = nexttile;
        imagesc(slice_i_z);
        title('Slice in Z direction');
        xlabel('X'); ylabel('Y');
        
        % Set shared color limits across all subplots (optional but good for comparison)
        clim = [0, max_Isppa];  % or manually set a fixed range
        set([ax1, ax2, ax3], 'CLim', clim);
        
        % Shared colorbar
        cb = colorbar(ax3, 'Location', 'eastoutside');
        ylabel(cb, 'ISPPA (W/cm^2)');
        
        ax3 = nexttile;
        
        axis off; % no axes visible
        text(0, 1, info_str, ...
            'Units', 'normalized', ...
            'FontSize', 10, ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', 'white', ...
            'EdgeColor', 'black', ...
            'Margin', 5);
        
        updated_title = strcat(fig_title, {' - '}, dir_step);
        sgtitle(updated_title);
    
        % Convert figure to image
        frame = getframe(gcf);            % capture current figure
        img = frame2im(frame);            % convert to image (RGB)

        % Convert to indexed image for TIFF saving
        [imind, cm] = rgb2ind(img, 256);
        
        output_dir = fullfile(path, 'summaries');
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end

        % Define output file
        out_file = fullfile(output_dir, strcat('summary_', coord_dir,'_', fig_title, '.tiff'));
        
        % If first image, create the file. Otherwise, append.
        if s == 1
            imwrite(imind, cm, out_file, 'tiff', 'WriteMode', 'overwrite', 'Compression', 'none');
        else
            imwrite(imind, cm, out_file, 'tiff', 'WriteMode', 'append', 'Compression', 'none');
        end
        
        close(gcf); % close the figure to save memory
    end
end

