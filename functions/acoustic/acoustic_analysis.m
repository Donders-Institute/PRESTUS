function [results, acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos] = ...
            acoustic_analysis(parameters, kwave_medium, medium_masks, sensor_data, ...
                                segmentation, source_labels)
% ACOUSTIC_ANALYSIS  Compute acoustic intensity metrics and generate output plots and tables
%
% Derives spatial pressure, intensity (Isppa), and mechanical index (MI) maps
% from k-Wave sensor output. For layered and phantom media, computes per-tissue
% peak values (skin, skull, brain) and transcranial MI. Writes a summary CSV and
% saves overlay plots of the intensity on the segmentation.
%
% Use as:
%   [results, acoustic_Ipa, acoustic_MI, acoustic_pressure, highlighted_pos] = ...
%       acoustic_analysis(parameters, kwave_medium, medium_masks, sensor_data, ...
%                         segmentation, source_labels)
%
% Input:
%   parameters   - PRESTUS config with transducer, grid, analysis, io fields
%   kwave_medium - spatial medium property maps (sound_speed, density, etc.)
%   medium_masks - voxel-wise tissue layer label map
%   sensor_data  - k-Wave output; requires p_max_all [Pa]
%   segmentation - CHARM tissue label map
%   source_labels- transducer element label map (for overlay plots)
%
% Output:
%   results          - struct with fields: Isppa, Isppa_after_exit_plane, Psptp,
%                      Ipa_target, Ipa_target_radius, real_focal_distance, and per-tissue
%                      Isppa/Psptp/MI fields for layered/phantom media [W/cm², Pa, -]
%   acoustic_Ipa     - spatial Isppa map [W/cm²]
%   acoustic_MI      - spatial mechanical index map [-]
%   acoustic_pressure- temporal peak pressure map [Pa]
%   highlighted_pos  - [1x3] grid coordinates of peak intensity for plots [voxels]
%
% See also: ACOUSTIC_WRAPPER, ACOUSTIC_SIMULATION, SIMULATION_NIFTI, MASKED_MAX_3D

arguments
    parameters   (1,1) struct
    kwave_medium (1,1) struct
    medium_masks {mustBeNumericOrLogical}
    sensor_data  (1,1) struct
    segmentation {mustBeNumericOrLogical}
    source_labels{mustBeNumericOrLogical}
end

    disp('Processing the results of acoustic simulations...')

    % select transducer info
    tr = parameters.transducer(1);
    trans_pos = tr.trans_pos;
    focus_pos = tr.focus_pos;
    
    % intialize output structure
    results = struct();

    % Temporal peak pressure at every gridpoint (p_max_all = peak over last steady-state cycles)
    acoustic_pressure = gather(sensor_data.p_max_all); % gather is used since it could be a GPU array
    results.Psptp = max(acoustic_pressure(:)); % spatial peak temporal peak pressure

    % Calculates the Isppa for every gridpoint
    acoustic_Ipa = acoustic_pressure.^2./...
        (2*(kwave_medium.sound_speed.*kwave_medium.density)).*1e-4;
    % Calculates the ISPPA (spatial peak pulse-average intensity)
    results.Isppa = max(acoustic_Ipa(:));

    % Calculates the Mechanical Index for every gridpoint

    freq_Hz = tr.freq_hz;
    freq_MHz = freq_Hz/10^6;
    acoustic_MI = (acoustic_pressure/10^6)/sqrt(freq_MHz);

    % Creates the foundation for a mask before the exit plane to calculate max values outside of it
    if numel(parameters.transducer) > 1
        warn(['Multi-transducer: exit-plane related metrics (after_exit_plane_mask, ' ...
                    'max_Isppa_after_exit_plane, real_focal_distance, etc.) are computed ' ...
                    'only w.r.t. the canonical first transducer.']);
    end
    comp_grid_size = size(sensor_data.p_max_all);
    after_exit_plane_mask = ones(comp_grid_size);
    bowl_depth_grid = round((tr.(tr.type).curv_radius_mm-...
        tr.(tr.type).dist_geom_ep_mm)/parameters.grid.resolution_mm);
    % Places the exit plane mask in the grid, adjusted to the amount of dimensions
    if numel(parameters.grid.dims) == 3
        if trans_pos(3) > comp_grid_size(3)/2
            after_exit_plane_mask(:,:,(trans_pos(numel(parameters.grid.dims))-...
                bowl_depth_grid):end) = 0;
        else
            after_exit_plane_mask(:,:,1:(trans_pos(numel(parameters.grid.dims))+...
                bowl_depth_grid)) = 0;
        end
    elseif numel(parameters.grid.dims) == 2
        if trans_pos(2) > comp_grid_size(2)/2
            after_exit_plane_mask(:,(trans_pos(numel(parameters.grid.dims))-...
                bowl_depth_grid):end) = 0;
        else
            after_exit_plane_mask(:,1:(trans_pos(numel(parameters.grid.dims))+...
                bowl_depth_grid)) = 0;
        end
    end

    % [1] average Isppa within a circle around the target
    % convert the radius from mm to voxels
    avg_radius = round(parameters.analysis.focus_area_radius/parameters.grid.resolution_mm); % [voxel]
    idx = arrayfun(@(d) max(1,focus_pos(d)-avg_radius):...
        min(size(acoustic_Ipa,d),focus_pos(d)+avg_radius), ...
            1:ndims(acoustic_Ipa), 'UniformOutput', false);
    results.Ipa_target_radius = acoustic_Ipa(idx{:});
    results.Ipa_target_radius = mean(results.Ipa_target_radius(:));
    
    % [2] IPA and temporal peak pressure at the target
    idx = num2cell(focus_pos);
    results.Ipa_target  = acoustic_Ipa(idx{:});
    results.Ptp_target  = acoustic_pressure(idx{:}); % temporal peak pressure at target

    % Get tissue-specific masks
    mask = tissuemask_binary(parameters, medium_masks);

    % calculate X, Y and Z coordinates of max. intensity 
    % if medium includes brain, consider only intensities within the brain
    % otherwise: consider all media beyond the exit plane
    if any(mask.brain(:))
        [results.Isppa_after_exit_plane, Ix_eplane, Iy_eplane,Iz_eplane] = ...
            masked_max_3d(acoustic_Ipa, mask.brain);
    else
        [results.Isppa_after_exit_plane, Ix_eplane, Iy_eplane, Iz_eplane] = ...
            masked_max_3d(acoustic_Ipa, after_exit_plane_mask);
    end
    % combine coordinates into one point of max. intensity in the grid
    if numel(parameters.grid.dims)==3
        results.max_isppa_eplane_pos = [Ix_eplane, Iy_eplane, Iz_eplane];
    elseif numel(parameters.grid.dims)==2
        results.max_isppa_eplane_pos = [Ix_eplane, Iy_eplane];
    end
    disp('Final transducer, expected focus, and max ISPPA positions')
    
    % [3] calculate the realized focal distance
    real_focal_distance = norm(results.max_isppa_eplane_pos-trans_pos)*parameters.grid.resolution_mm; % [mm]

    % Layer-specific outcomes (in case a layered simulation)
    if contains(parameters.simulation.medium, {'layered'; 'phantom'})

        % calculate max. isppa and location across full space
        [~, Ix, Iy, Iz] = masked_max_3d(acoustic_Ipa, ones(size(medium_masks)));

        % highlight this position in the future
        highlighted_pos = [Ix, Iy, Iz];
        highlighted_pos = highlighted_pos(1:numel(trans_pos));
        
        % extract indices in brain medium
        [results.Isppa_brain, Ix_brain, Iy_brain, Iz_brain] = ...
            masked_max_3d(acoustic_Ipa, mask.brain);
        [results.min_Isppa_brain] = min(acoustic_Ipa(mask.brain));
        half_max = acoustic_Ipa >= results.Isppa_brain/2 & mask.brain;
        half_max_ISPPA_volume_brain = sum(half_max(:))*(parameters.grid.resolution_mm^3);
        [results.Psptp_brain] = masked_max_3d(acoustic_pressure, mask.brain);
        [results.MI_brain] = masked_max_3d(acoustic_MI, mask.brain);

        % extract indices in skull medium
        [results.Isppa_skull] = masked_max_3d(acoustic_Ipa, mask.skull);
        [results.Psptp_skull] = masked_max_3d(acoustic_pressure, mask.skull);
        [results.MI_skull] = masked_max_3d(acoustic_MI, mask.skull);

        % extract indices in skin medium
        [results.Isppa_skin] = masked_max_3d(acoustic_Ipa, mask.skin);
        [results.Psptp_skin] = masked_max_3d(acoustic_pressure, mask.skin);
        [results.MI_skin] = masked_max_3d(acoustic_MI, mask.skin);

        % MItc: transcranial MI — max MI across intracranial voxels (WM, GM, CSF, blood)
        intracranial_mask = ismember(segmentation, charm_seg_labels().intracranial);
        if any(intracranial_mask(:))
            [results.MI_tc] = masked_max_3d(acoustic_MI, intracranial_mask);
        else
            results.MI_tc = NaN; % intracranial tissues not modelled
        end

        writetable(table(parameters.subject_id, freq_Hz, ...
            results.Isppa, results.Isppa_after_exit_plane, ...
            real_focal_distance, results.Isppa_skin, results.Isppa_skull, ...
            results.Isppa_brain, results.Psptp_skin, results.Psptp_skull, ...
            results.Psptp_brain, results.Ptp_target, results.MI_skin, results.MI_skull, ...
            results.MI_brain, results.MI_tc, Ix_brain, Iy_brain, Iz_brain, trans_pos, focus_pos, ...
            results.Ipa_target, results.Ipa_target_radius, half_max_ISPPA_volume_brain, ...
            'VariableNames', { 'subject_id', 'freq_Hz', ...
            'Isppa', 'Isppa_after_exitplane', 'real_focal_distance_mm', ...
            'Isppa_skin', 'Isppa_skull', 'Isppa_brain', ...
            'Psptp_skin', 'Psptp_skull', 'Psptp_brain', 'Ptp_target', ...
            'MI_skin', 'MI_skull', 'MI_brain', 'MI_tc', ...
            'Ix_brain_vox', 'Iy_brain_vox', 'Iz_brain_vox', ...
            'trans_pos_vox', 'focus_pos_vox', ...
            'Ipa_target', 'Ipa_target_radius', ...
            'halfmax_ISPPA_volume_brain_mm3'}), ...
            parameters.io.filename_table);
    else
        % If no layered tissue was selected, the max Isppa is highlighted on the plane and written in a table.
        highlighted_pos = results.max_isppa_eplane_pos;
        writetable(table(parameters.subject_id, freq_Hz, results.Isppa, results.Isppa_after_exit_plane, ...
            results.Psptp, results.Ptp_target, real_focal_distance, trans_pos, focus_pos, ...
            results.Ipa_target, results.Ipa_target_radius, ...
            'VariableNames', {'subject_id', 'freq_Hz', 'Isppa', 'Isppa_after_exitplane', ...
            'Psptp', 'Ptp_target', 'real_focal_distance_mm', 'trans_pos_vox', 'focus_pos_vox', ...
            'Ipa_target', 'Ipa_target_radius'}), ...
        parameters.io.filename_table);
    end

    % Plot intensity on the segmented image (up to 2 transducers)
    n_plots = min(2, numel(parameters.transducer));
    if numel(parameters.transducer) > n_plots
        warn('More than two transducers: intensity plots on segmentation will be created only for the first 2 transducers');
    end
    for ti = 1:n_plots

        tpos_sim = parameters.transducer(ti).trans_pos;
        fpos_sim = parameters.transducer(ti).focus_pos;
		
		trans_suffix = '';
		if n_plots > 1; trans_suffix = sprintf('_T%02d', ti); end
        
        if numel(parameters.grid.dims)==3
            slices = struct( ...
                'dim', {'x', 'y', 'z'}, ...
                'pos', {fpos_sim(1), fpos_sim(2), fpos_sim(3)} ...
                );

            for i_slice = 1:numel(slices)
                slice = slices(i_slice);

                [~,~,~,~,~,~,~,h]=plot_overlay(...
                    acoustic_Ipa, ...
                    segmentation, ...
                    source_labels, ...
                    parameters, ...
                    {slice.dim, slice.pos}, ...
                    tpos_sim, ...
                    fpos_sim, ...
                    highlighted_pos);

                % Construct output filename
                output_plot = fullfile( ...
                    parameters.io.dir_img, ...
                    sprintf('sub-%03d_%s_intensity_%s%s%s.png', ...
                    parameters.subject_id, ...
                    parameters.simulation.medium, ...
                    slice.dim, ...
                    trans_suffix, ...
                    parameters.io.output_affix) ...
                    );

                % Keep original colors and save
                set(h, 'InvertHardcopy', 'off');
                saveas(h, output_plot, 'png');
                close(h);
            end
        else
            h = plot_overlay_2d(...
                acoustic_Ipa, ...
                segmentation, ...
                source_labels, ...
                after_exit_plane_mask, ...
                tpos_sim, ...
                fpos_sim, ...
                highlighted_pos);

            output_plot = fullfile(parameters.io.dir_img, ...
                sprintf('sub-%03d_%s_intensity%s%s.png', ...
                parameters.subject_id, ...
                parameters.simulation.medium, ...
                trans_suffix, ...
                parameters.io.output_affix));
            set(h, 'InvertHardcopy', 'off'); % keep original colours
            saveas(h, output_plot, 'png')
            close(h);
        end
	end
end
