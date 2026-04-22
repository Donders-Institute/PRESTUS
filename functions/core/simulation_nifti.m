function simulation_nifti(parameters, planimg, results_acoustic, ...
    acoustic_Ipa, acoustic_MI, acoustic_pressure, medium_masks, ...
    results_heating, kwave_medium, highlighted_pos)
% SIMULATION_NIFTI  Export k-Wave simulation results to NIfTI in native and MNI space
%
% For layered and phantom media, writes compressed NIfTI files for each
% requested data type (medium_masks, intensity, MI, pressure, and — if heating
% was run — heating, heatrise, CEM43, and their end-of-train variants). For
% layered (non-phantom) media each volume is back-transformed to the original
% T1 image space via tformarray before being converted to MNI space with
% SimNIBS subject2mni. Intensity overlays on the T1 are saved as PNG. Thermal
% back-transforms include edge-artifact removal using morphological erosion.
%
% Use as:
%   simulation_nifti(parameters, planimg, results_acoustic, acoustic_Ipa, ...
%                    acoustic_MI, acoustic_pressure, medium_masks, ...
%                    results_heating, kwave_medium, highlighted_pos)
%
% Input:
%   parameters      - (1,1) struct, PRESTUS config with io.output_dir,
%                     simulation.medium, io.output_affix, modules flags, etc.
%   planimg         - (1,1) struct, SimNIBS planning data: t1_image_orig,
%                     inv_transf (affine tform), t1_header
%   results_acoustic- (1,1) struct, acoustic analysis results (Isppa_brain, Isppa, etc.)
%   acoustic_Ipa    - (:,:,:) single, Isppa map [W/cm²]
%   acoustic_MI     - (:,:,:) single, mechanical index map [-]
%   acoustic_pressure- (:,:,:) single, temporal peak pressure map [Pa]
%   medium_masks    - (:,:,:) uint8, voxel-wise tissue label map
%   results_heating - (1,1) struct, heating results: maxT, heating_endT, CEM43,
%                     CEM43_end, CEM43_iso, CEM43_iso_end [°C / CEM43 units]
%   kwave_medium    - (1,1) struct, k-Wave medium properties; used for temp_0 [°C]
%   highlighted_pos - [1x3] numeric, grid-space position of peak intensity [voxels]
%
% See also: ACOUSTIC_ANALYSIS, CONVERT_FINAL_TO_MNI_SIMNIBS, PLOT_OVERLAY

arguments
    parameters       (1,1) struct
    planimg          (1,1) struct
    results_acoustic (1,1) struct
    acoustic_Ipa     {mustBeNumeric}
    acoustic_MI      {mustBeNumeric}
    acoustic_pressure{mustBeNumeric}
    medium_masks     {mustBeNumericOrLogical}
    results_heating  struct
    kwave_medium     (1,1) struct
    highlighted_pos  (1,:) {mustBeNumeric}
end

    if contains(parameters.simulation.medium, {'layered'; 'phantom'})

        data_types = "medium_masks";
        if parameters.state.acoustics_available == 1 
            data_types  = [data_types, "intensity","MI","pressure"];
        end
        if parameters.state.heating_available == 1
            data_types  = [data_types, "heating", "heating_end", "heatrise", "heatrise_end", "CEM43", "CEM43_end", "CEM43_iso", "CEM43_iso_end"];
        end
        for data_type = data_types
            orig_file = fullfile(parameters.io.output_dir, sprintf('sub-%03d_final_%s_orig_coord%s',...
                parameters.subject_id, data_type, parameters.io.output_affix));
            mni_file  = fullfile(parameters.io.output_dir, sprintf('sub-%03d_final_%s_MNI%s.nii.gz',...
                parameters.subject_id, data_type, parameters.io.output_affix));

            if strcmp(data_type, "medium_masks")
                data = medium_masks;
            elseif strcmp(data_type, "intensity")
                data = single(acoustic_Ipa);
            elseif strcmp(data_type, "MI")
                data = single(acoustic_MI);
            elseif strcmp(data_type, "pressure")
                data = single(acoustic_pressure);
            elseif strcmp(data_type, "heating")
                data = single(results_heating.maxT);
            elseif strcmp(data_type, "heating_end")
                data = single(results_heating.heating_endT);
            elseif strcmp(data_type, "heatrise")
                data = single(results_heating.maxT-kwave_medium.temp_0);
            elseif strcmp(data_type, "heatrise_end")
                data = single(results_heating.heating_endT-kwave_medium.temp_0);
            elseif strcmp(data_type, "CEM43")
                data = single(results_heating.CEM43);
            elseif strcmp(data_type, "CEM43_end")
                data = single(results_heating.CEM43_end);
            elseif strcmp(data_type, "CEM43_iso")
                data = single(results_heating.CEM43_iso);
            elseif strcmp(data_type, "CEM43_iso_end")
                data = single(results_heating.CEM43_iso_end);
            end
            orig_file_with_ext = strcat(orig_file, '.nii.gz');
    
            if confirm_overwriting(orig_file_with_ext, parameters)
                if ~strcmp(parameters.simulation.medium, 'phantom')
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
                        elseif strcmp(data_type, "CEM43") || strcmp(data_type, "CEM43_end") || ...
                               strcmp(data_type, "CEM43_iso") || strcmp(data_type, "CEM43_iso_end")
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

            if strcmp(data_type, "intensity") && ~strcmp(parameters.simulation.medium, 'phantom')
            
                max_plots = min(2, numel(parameters.transducer));
                if numel(parameters.transducer) > max_plots
                    warning('More than two transducers: intensity-over-T1 plots will be created only for the first 2 transducers');
                end

                % define the maximum value to plot
                if isfield(results_acoustic, 'Isppa_brain') && ~isempty(results_acoustic.Isppa_brain) && ~isnan(results_acoustic.Isppa_brain)
                    max_val = results_acoustic.Isppa_brain;
                else
                    max_val = results_acoustic.Isppa;
                end
            
                for ti = 1:max_plots
                    tpos_sim = parameters.transducer(ti).trans_pos;
                    fpos_sim = parameters.transducer(ti).focus_pos;
            
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
            
                    % Plots the intensity over the untransformed image
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
            
                    trans_suffix = '';
                    if max_plots > 1; trans_suffix = sprintf('_T%02d', ti); end
                    output_plot_filename = fullfile(parameters.io.output_dir, ...
                        sprintf('sub-%03d_%s_intensity_t1%s%s.png', ...
                        parameters.subject_id, parameters.simulation.medium, trans_suffix, ...
                        parameters.io.output_affix));
                    saveas(h, output_plot_filename, 'png')
                    close(h);
                end
            end
            
            m2m_folder= fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id));
            
            % transform outputs to MNI space (using SimNibs or applying transformation matrix)
            if ~confirm_overwriting(mni_file, parameters) || strcmp(parameters.simulation.medium, 'phantom')
                continue
            end
            convert_final_to_MNI_simnibs(orig_file_with_ext , m2m_folder, mni_file, parameters, 'interpolation_order', 0);

            clear data;
        end
        
        % Since charm does not transform the T1 into MNI space, one is manually created here
        if ~strcmp(parameters.simulation.medium, 'phantom')
            path_to_input_img = fullfile(m2m_folder,'T1.nii.gz');
            path_to_output_img = fullfile(m2m_folder,'toMNI','T1_to_MNI_post-hoc.nii.gz');

            if ~exist(path_to_output_img,'file')
                convert_final_to_MNI_simnibs(path_to_input_img, m2m_folder, path_to_output_img, parameters)
            end
        end
    end