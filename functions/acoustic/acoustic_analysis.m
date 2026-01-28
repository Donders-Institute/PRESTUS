function [results, acoustic_isppa, acoustic_MI, acoustic_pressure, highlighted_pos] = ...
            acoustic_analysis(parameters, kwave_medium, medium_masks, sensor_data, trans_pos, focus_pos, ...
                                segmentation, source_labels)
    
    disp('Processing the results of acoustic simulations...')
    
    % intialize output structure
    results = struct();

    % What is the highest pressure level for every gridpoint
    acoustic_pressure = gather(sensor_data.p_max_all); % gather is used since it could be a GPU array
    results.max_pressure = max(acoustic_pressure(:));

    % Calculates the Isppa for every gridpoint
    acoustic_isppa = acoustic_pressure.^2./...
        (2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4;
    % Calculates the max Isppa
    results.max_Isppa = max(acoustic_isppa(:));

    % Calculates the Mechanical Index for every gridpoint
    acoustic_MI = (acoustic_pressure/10^6)/sqrt((parameters.transducer(1).source_freq_hz/10^6));

    % Creates the foundation for a mask before the exit plane to calculate max values outside of it
    if numel(parameters.transducer) > 1
        warning(['Multi-transducer: exit-plane related metrics (after_exit_plane_mask, ' ...
                    'max_Isppa_after_exit_plane, real_focal_distance, etc.) are computed ' ...
                    'only w.r.t. the canonical first transducer.']);
    end
    comp_grid_size = size(sensor_data.p_max_all);
    after_exit_plane_mask = ones(comp_grid_size);
    bowl_depth_grid = round((parameters.transducer(1).curv_radius_mm-...
        parameters.transducer(1).dist_to_plane_mm)/parameters.grid_step_mm);
    % Places the exit plane mask in the grid, adjusted to the amount of dimensions
    if parameters.n_sim_dims == 3
        if trans_pos(3) > comp_grid_size(3)/2
            after_exit_plane_mask(:,:,(trans_pos(parameters.n_sim_dims)-...
                bowl_depth_grid):end) = 0;
        else
            after_exit_plane_mask(:,:,1:(trans_pos(parameters.n_sim_dims)+...
                bowl_depth_grid)) = 0;
        end
    elseif parameters.n_sim_dims == 2
        if trans_pos(2) > comp_grid_size(2)/2
            after_exit_plane_mask(:,(trans_pos(parameters.n_sim_dims)-...
                bowl_depth_grid):end) = 0;
        else
            after_exit_plane_mask(:,1:(trans_pos(parameters.n_sim_dims)+...
                bowl_depth_grid)) = 0;
        end
    end

    % [1] average Isppa within a circle around the target
    % convert the radius from mm to voxels
    avg_radius = round(parameters.focus_area_radius/parameters.grid_step_mm); % [voxel]
    idx = arrayfun(@(d) max(1,focus_pos(d)-avg_radius):...
        min(size(acoustic_isppa,d),focus_pos(d)+avg_radius), ...
            1:ndims(acoustic_isppa), 'UniformOutput', false);
    results.avg_isppa_around_target = acoustic_isppa(idx{:});
    results.avg_isppa_around_target = mean(results.avg_isppa_around_target(:));
    
    % [2] Isppa within the original stimulation target
    idx = num2cell(focus_pos);
    results.isppa_at_target = acoustic_isppa(idx{:});

    % Get tissue-specific masks
    mask = tissuemask_binary(parameters, medium_masks);

    % calculate X, Y and Z coordinates of max. intensity 
    % if medium includes brain, consider only intensities within the brain
    % otherwise: consider all media beyond the exit plane
    if ~isempty(mask.brain)
        [results.max_Isppa_after_exit_plane, Ix_eplane, Iy_eplane,Iz_eplane] = ...
            masked_max_3d(acoustic_isppa, mask.brain);
    else
        [results.max_Isppa_after_exit_plane, Ix_eplane, Iy_eplane, Iz_eplane] = ...
            masked_max_3d(acoustic_isppa, after_exit_plane_mask);
    end
    % combine coordinates into one point of max. intensity in the grid
    if parameters.n_sim_dims==3
        results.max_isppa_eplane_pos = [Ix_eplane, Iy_eplane, Iz_eplane];
    elseif parameters.n_sim_dims==2
        results.max_isppa_eplane_pos = [Ix_eplane, Iy_eplane];
    end
    disp('Final transducer, expected focus, and max ISPPA positions')
    
    % [3] calculate the realized focal distance
    real_focal_distance = norm(results.max_isppa_eplane_pos-trans_pos)*parameters.grid_step_mm; % [mm]

    % Layer-specific outcomes (in case a layered simulation)
    if contains(parameters.simulation_medium, {'layered'; 'phantom'})

        % calculate max. isppa and location across full space
        [~, Ix, Iy, Iz] = masked_max_3d(acoustic_isppa, ones(size(medium_masks)));

        % highlight this position in the future
        highlighted_pos = [Ix, Iy, Iz];
        highlighted_pos = highlighted_pos(1:numel(trans_pos));
        
        % extract indices in brain medium
        [results.max_Isppa_brain, Ix_brain, Iy_brain, Iz_brain] = ...
            masked_max_3d(acoustic_isppa, mask.brain);
        [results.min_Isppa_brain] = min(acoustic_isppa(mask.brain));
        half_max = acoustic_isppa >= results.max_Isppa_brain/2 & mask.brain;
        half_max_ISPPA_volume_brain = sum(half_max(:))*(parameters.grid_step_mm^3);
        [results.max_pressure_brain] = masked_max_3d(acoustic_pressure, mask.brain);
        [results.max_MI_brain] = masked_max_3d(acoustic_MI, mask.brain);
        
        % extract indices in skull medium
        [results.max_Isppa_skull] = masked_max_3d(acoustic_isppa, mask.skull);
        [results.max_pressure_skull] = masked_max_3d(acoustic_pressure, mask.skull);
        [results.max_MI_skull] = masked_max_3d(acoustic_MI, mask.skull);
        
        % extract indices in skin medium
        [results.max_Isppa_skin] = masked_max_3d(acoustic_isppa, mask.skin);
        [results.max_pressure_skin] = masked_max_3d(acoustic_pressure, mask.skin);
        [results.max_MI_skin] = masked_max_3d(acoustic_MI, mask.skin);

        writetable(table(parameters.subject_id, results.max_Isppa, results.max_Isppa_after_exit_plane, ...
            real_focal_distance, results.max_Isppa_skin, results.max_Isppa_skull, ...
            results.max_Isppa_brain, results.max_pressure_skin, results.max_pressure_skull, ...
            results.max_pressure_brain, results.max_MI_skin, results.max_MI_skull, ...
            results.max_MI_brain, Ix_brain, Iy_brain, Iz_brain, trans_pos, focus_pos, ...
            results.isppa_at_target, results.avg_isppa_around_target, half_max_ISPPA_volume_brain, ...
            'VariableNames', { 'subject_id', 'max_Isppa', 'max_Isppa_after_exitplane', 'real_focal_distance_mm', ...
            'max_Isppa_skin', 'max_Isppa_skull', 'max_Isppa_brain', ...
            'max_pressure_skin_Pa', 'max_pressure_skull_Pa', 'max_pressure_brain_Pa', ...
            'max_MI_skin', 'max_MI_skull', 'max_MI_brain', ...
            'Ix_brain_vox', 'Iy_brain_vox', 'Iz_brain_vox', ...
            'trans_pos_vox', 'focus_pos_vox', ...
            'isppa_at_target', 'avg_isppa_around_target', ...
            'halfmax_ISPPA_volume_brain_mm3'}), ...
            parameters.filename_output_table);
    else 
        % If no layered tissue was selected, the max Isppa is highlighted on the plane and written in a table.
        highlighted_pos = results.max_isppa_eplane_pos;
        writetable(table(parameters.subject_id, results.max_Isppa, results.max_Isppa_after_exit_plane, results.max_pressure, real_focal_distance, trans_pos, focus_pos, results.isppa_at_target, results.avg_isppa_around_target, ...
        'VariableNames', {'subject_id','max_Isppa', 'max_Isppa_after_exitplane', 'max_pressure_Pa', 'real_focal_distance_mm', 'trans_pos_vox', 'focus_pos_vox', 'isppa_at_target', 'avg_isppa_around_target'}), ...
        parameters.filename_output_table);
    end

    % Plot intensity on the segmented image (up to 2 transducers)
    n_plots = min(2, numel(parameters.transducer));
    if numel(parameters.transducer) > n_plots
        warning('More than two transducers: ISPPA plots on segmentation will be created only for the first 2 transducers');
    end
    for ti = 1:n_plots
        if ti == 1
            tpos_sim = trans_pos;
            fpos_sim = focus_pos;
        else
            tpos_sim = parameters.transducer(ti).trans_pos;
            fpos_sim = parameters.transducer(ti).focus_pos;
        end
    
        if parameters.n_sim_dims==3
            [~,~,~,~,~,~,~,h]=plot_overlay(...
                acoustic_isppa, ...
                segmentation, ...
                source_labels, ...
                parameters, ...
                {'y', fpos_sim(2)}, ...
                tpos_sim, ...
                fpos_sim, ...
                highlighted_pos);
        else
            h = plot_overlay_2d(...
                acoustic_isppa, ...
                segmentation, ...
                source_labels, ...
                after_exit_plane_mask, ...
                tpos_sim, ...
                fpos_sim, ...
                highlighted_pos);
        end
    
        output_plot = fullfile(parameters.output_dir, ...
            sprintf('sub-%03d_%s_isppa_T%02d%s.png', ...
            parameters.subject_id, parameters.simulation_medium, ti, parameters.results_filename_affix));
        set(h, 'InvertHardcopy', 'off'); % keep original colours
        saveas(h, output_plot, 'png')
        close(h);
    end