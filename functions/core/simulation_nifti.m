function simulation_nifti(parameters, planimg, results_acoustic, acoustic_isppa, acoustic_MI, acoustic_pressure, ...
                         medium_masks, results_heating, kwave_medium, trans_pos, ...
                         focus_pos, highlighted_pos)

% SIMULATION_NIFTI - Export k-Wave simulation results to NIfTI (orig + MNI space).
%
%   simulation_nifti(PARAMETERS, PLANIMG, ACOUSTIC_ISPPA, ACOUSTIC_MI, ACOUSTIC_PRESSURE, ...
%                   MEDIUM_MASKS, HEATING_MAXT, HEATING_CEM43, KWAVE_MEDIUM, ...
%                   TRANS_POS, FOCUS_POS, HIGHLIGHTED_POS)
%
% Inputs:
%   - parameters (struct) - Simulation config: output_dir, simulation_medium, 
%                           results_filename_affix, run_heating_sims, etc.
%   - planimg (struct) - SimNibs planning: t1_image_orig, inv_transf, t1_header
%   - acoustic_isppa (array) - Peak spatial-average intensity [W/cm²], from acoustic_analysis
%   - acoustic_MI (array) - Mechanical Index grid
%   - acoustic_pressure (array) - Peak pressure [Pa]
%   - medium_masks (array) - Layer label mask (uint8)
%   - results_heating.maxT (array) - Max temperature from thermal sim [°C] (if run_heating_sims)
%   - results_heating.CEM43 (array) - Cumulative Equivalent Minutes at 43°C
%   - kwave_medium (struct) - k-Wave medium (for temp_0)
%   - trans_pos, focus_pos, highlighted_pos (arrays) - Transducer/focus/max-Isppa voxel coords
%
% Outputs:
%   - None (side effects): Saves NIfTIs (*.nii.gz) to parameters.output_dir:
%       * orig coord: sub-XXX_final_{isppa|MI|...}_orig_coord.nii.gz
%       * MNI space:  sub-XXX_final_{isppa|MI|...}_MNI.nii.gz
%     Plots: Isppa overlays on T1 for first 2 transducers.

    if contains(parameters.simulation_medium, {'layered'; 'phantom'})

        data_types = "medium_masks";
        if parameters.acoustics_available == 1 
            data_types  = [data_types, "isppa","MI","pressure"];
        end
        if isfield(parameters, 'run_heating_sims') && parameters.run_heating_sims 
            data_types  = [data_types, "heating", "heatrise", "CEM43"];
        end
        for data_type = data_types
            orig_file = fullfile(parameters.output_dir, sprintf('sub-%03d_final_%s_orig_coord%s',...
                parameters.subject_id, data_type, parameters.results_filename_affix));
            mni_file  = fullfile(parameters.output_dir, sprintf('sub-%03d_final_%s_MNI%s.nii.gz',...
                parameters.subject_id, data_type, parameters.results_filename_affix));

            if strcmp(data_type, "isppa")
                data = single(acoustic_isppa);
            elseif strcmp(data_type, "MI")
                data = single(acoustic_MI);
            elseif strcmp(data_type, "pressure")
                data = single(acoustic_pressure);
            elseif strcmp(data_type, "medium_masks")
                data = medium_masks;
            elseif strcmp(data_type, "heating")
                data = single(results_heating.maxT);
            elseif strcmp(data_type, "heatrise")
                data = single(results_heating.maxT-kwave_medium.temp_0);
            elseif strcmp(data_type, "CEM43")
                data = single(results_heating.CEM43);
            end
            orig_file_with_ext = strcat(orig_file, '.nii.gz');
    
            if confirm_overwriting(orig_file_with_ext, parameters)
                if ~strcmp(parameters.simulation_medium, 'phantom')
                    % Transforms the data to original T1 image dimensions and orientation
                    orig_hdr = planimg.t1_header;
    
                    if strcmp(data_type, "medium_masks")
                        data_backtransf = tformarray(uint8(data), planimg.inv_transf, ...
                                                makeresampler('nearest', 'fill'), [1 2 3], [1 2 3], size(planimg.t1_image_orig), [], 0) ;
    
                        orig_hdr.Datatype = 'uint8';
                        orig_hdr.BitsPerPixel = 8;
                    else
                        data_backtransf = tformarray(data, planimg.inv_transf, ...
                                                makeresampler('cubic', 'fill'), [1 2 3], [1 2 3], size(planimg.t1_image_orig), [], 0) ;
                        
                        orig_hdr.Datatype = 'single';

                        if strcmp(data_type, "heating")
                            % Removes edge artifacts
                            heatmap_mask = data_backtransf > 0;
                            % Removes the edge that often contains a mix of values
                            heatmap_mask_eroded = imerode(heatmap_mask, strel('sphere', 3));
                            % Takes a second edge within the new mask to determine what
                            % values to replace the 0's with
                            heatmap_mask_eroded_edge = imdilate(heatmap_mask_eroded, strel('sphere', 1)) & heatmap_mask & ~heatmap_mask_eroded;
                            heatmap_band_values = data_backtransf(heatmap_mask_eroded_edge);
                            if ~isempty(heatmap_band_values)
                                heatmap_band_value = median(heatmap_band_values);
                            else
                                heatmap_band_value = parameters.thermal.temp_0.water;
                            end
                            data_backtransf(~heatmap_mask_eroded) = heatmap_band_value;
                            data_backtransf(1:2,:,:)     = parameters.thermal.temp_0.water;  % First two planes in the 1st dimension
                            data_backtransf(end-1:end,:,:) = parameters.thermal.temp_0.water;  % Last two planes in the 1st dimension
                            
                            data_backtransf(:,1:2,:)     = parameters.thermal.temp_0.water;  % First two planes in the 2nd dimension
                            data_backtransf(:,end-1:end,:) = parameters.thermal.temp_0.water;  % Last two planes in the 2nd dimension
                            
                            data_backtransf(:,:,1:2)     = parameters.thermal.temp_0.water;  % First two planes in the 3rd dimension
                            data_backtransf(:,:,end-1:end) = parameters.thermal.temp_0.water;  % Last two planes in the 3rd dimension
                        elseif strcmp(data_type, "CEM43")
                            % Removes edge artifacts
                            data_backtransf(data_backtransf <= 0) = 0.0000001;
                        end
                    end
                    niftiwrite(data_backtransf, orig_file, orig_hdr, 'Compressed', true);
                else
                    niftiwrite(data, orig_file, 'Compressed', true);
                end
            else 
                data_backtransf = niftiread(orig_file_with_ext);
            end

            if strcmp(data_type, "isppa") && ~strcmp(parameters.simulation_medium, 'phantom')
            
                max_plots = min(2, numel(parameters.transducer));
                if numel(parameters.transducer) > max_plots
                    warning('More than two transducers: ISPPA-over-T1 plots will be created only for the first 2 transducers');
                end

                % define the maximum value to plot
                if isfield(results_acoustic, 'max_Isppa_brain') & ~isempty(results_acoustic.max_Isppa_brain)
                    max_val = results_acoustic.max_Isppa_brain;
                else
                    max_val = results_acoustic.max_Isppa;
                end
            
                for ti = 1:max_plots
                    if ti == 1
                        tpos_sim = trans_pos;
                        fpos_sim = focus_pos;
                    else
                        tpos_sim = parameters.transducer(ti).trans_pos;
                        fpos_sim = parameters.transducer(ti).focus_pos;
                    end
            
                    % map canonical / per-T positions + highlighted_pos back to T1 space
                    backtransf_coordinates = round(tformfwd(...
                        [tpos_sim; fpos_sim; highlighted_pos], ...
                        planimg.inv_transf));
            
                    % Creates a visual overlay of this transducer
                    [~, source_labels] = transducer_setup(...
                        parameters.transducer(1), ...
                        backtransf_coordinates(1,:), ...
                        backtransf_coordinates(2,:), ...
                        size(planimg.t1_image_orig), ...
                        planimg.t1_header.PixelDimensions(1));
            
                    % Plots the Isppa over the untransformed image
                    [~,~,~,~,~,~,~,h]=plot_overlay(...
                        data_backtransf, ...
                        planimg.t1_image_orig, ...
                        source_labels, ...
                        parameters, ...
                        {'y', backtransf_coordinates(2,2)}, ...
                        backtransf_coordinates(1,:), ...
                        backtransf_coordinates(2,:), ...
                        backtransf_coordinates(3,:), ...
                        'show_rectangles', 0, ...
                        'grid_step', planimg.t1_header.PixelDimensions(1), ...
                        'overlay_color_range', [0 max_val], ...
                        'overlay_threshold_low', 0, ...
                        'overlay_threshold_high', max_val, ...
                        'rotation', 0); % rotation = 90 not implemented for transducer overlay
            
                    output_plot_filename = fullfile(parameters.output_dir, ...
                        sprintf('sub-%03d_%s_isppa_t1_T%02d%s.png', ...
                        parameters.subject_id, parameters.simulation_medium, ti, ...
                        parameters.results_filename_affix));
                    saveas(h, output_plot_filename, 'png')
                    close(h);
                end
            end
            
            m2m_folder= fullfile(parameters.seg_path, sprintf('m2m_sub-%03d', parameters.subject_id));
            
            % transform outputs to MNI space (using SimNibs or applying transformation matrix)
            if ~confirm_overwriting(mni_file, parameters) || strcmp(parameters.simulation_medium, 'phantom')
                continue
            end
            if strcmp(parameters.segmentation_software, 'headreco')
                if strcmp(data_type, "medium_masks")
                   convert_final_to_MNI_matlab(data, m2m_folder, planimg.inv_transf, parameters, 'nifti_filename', mni_file,  'nifti_data_type', 'uint8', 'BitsPerPixel', 8);
                else
                   convert_final_to_MNI_matlab(data, m2m_folder, planimg.inv_transf, parameters, 'nifti_filename', mni_file);
                end
            elseif strcmp(parameters.segmentation_software, 'charm')
                convert_final_to_MNI_simnibs(orig_file_with_ext , m2m_folder, mni_file, parameters, 'interpolation_order', 0);
            end

            clear data;
        end
        
        % Since charm does not transform the T1 into MNI space, one is manually created here
        if strcmp(parameters.segmentation_software, 'charm') && ~strcmp(parameters.simulation_medium, 'phantom')
            path_to_input_img = fullfile(m2m_folder,'T1.nii.gz');
            path_to_output_img = fullfile(m2m_folder,'toMNI','T1_to_MNI_post-hoc.nii.gz');

            if ~exist(path_to_output_img,'file')
                convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters)
            end
        end
    end