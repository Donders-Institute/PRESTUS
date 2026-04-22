function kwave_medium = medium_setup(parameters, medium_masks, planimg, pseudoCT)
% MEDIUM_SETUP  Build the k-Wave medium property struct from tissue masks and config
%
% Starting from uniform NaN grids, assigns acoustic and thermal properties
% layer-by-layer based on medium_masks. Skull bone properties are
% optionally overridden with pseudo-CT-derived values using the algorithms
% specified in parameters.pct (density, sound speed, attenuation). Alpha
% coefficients are remapped via fitPowerLawParamsMulti to correct for the
% fractional Laplacian second-order dispersion at the source frequency. For
% layered simulations, property maps are also smoothed if
% headmodel.smooth_properties is set. Medium property NIfTIs (in T1 space)
% are saved unconditionally for layered simulations; raw grid matrices are
% saved only in debug mode.
%
% Use as:
%   kwave_medium = medium_setup(parameters, medium_masks, planimg)
%   kwave_medium = medium_setup(parameters, medium_masks, planimg, pseudoCT)
%
% Input:
%   parameters   - PRESTUS config; must contain medium_properties, layers,
%                  grid.dims, grid.resolution_mm [mm], pct.enabled, fit_alpha_power,
%                  and transducer(1).freq_hz [Hz]
%   medium_masks - layer label map (values = medium index)
%   planimg      - planning image info from GRID_TISSUE_SETUP;
%                  may be empty for non-layered simulations
%   pseudoCT     - pseudo-CT Hounsfield values (same size as medium_masks);
%                  required when pct.enabled == 1 (optional, default: [])
%
% Output:
%   kwave_medium - struct with fields: sound_speed [m/s], density [kg/m^3],
%                  alpha_coeff [dB/(MHz cm)], alpha_power,
%                  thermal_conductivity [W/(m K)], specific_heat [J/(kg K)],
%                  perfusion_coeff [1/s], absorption_fraction, temp_0 [°C]
%
% See also: MEDIUM_PCT_DENSITY, MEDIUM_PCT_SOUNDSPEED, MEDIUM_PCT_ATTENUATION,
%   FITPOWERLAWPARAMSMULTI, MEDIUM_PROPERTIES_NIFTI

arguments
    parameters   (1,1) struct
    medium_masks {mustBeNumericOrLogical}
    planimg      (1,1) struct
    pseudoCT     {mustBeNumeric} = []
end

    % Set residual unmapped medium voxels to water
    medium_masks(medium_masks==0) = find(strcmp(fieldnames(parameters.medium_properties), 'water'));
    
    % Loads the medium settings from the config file
    medium = parameters.medium_properties;
    
    % Create empty matrices for medium properties
    empty_grid = NaN(parameters.grid.dims);
    sound_speed = empty_grid; 
    density = empty_grid;
    alpha_coeff = empty_grid;
    alpha_power = empty_grid;
    thermal_conductivity = empty_grid;
    specific_heat = empty_grid;
    perfusion = empty_grid;
    absorption_fraction = ones(parameters.grid.dims); % default: convert 100% of attenuation into absorption
    temp_0 = empty_grid;

    % Get layer and medium labels
    layer_labels = fieldnames(parameters.layers);
    medium_labels = fieldnames(parameters.medium_properties);

    % Iterate through each layer & assign medium ID to create medium mask
    for label_i = 1:length(layer_labels)
        label_name = layer_labels{label_i};
        medium_i = find(strcmp(medium_labels, label_name));
        if parameters.pct.enabled == 1 && strcmp(label_name, 'skull')
            skull_idx = find(ismember(medium_masks,medium_i));

            % set skull thermal conductivity
            thermal_conductivity(skull_idx) = medium.(label_name).thermal_conductivity;

            % set skull heat capacity     
            specific_heat(skull_idx) = medium.(label_name).specific_heat_capacity;

            % set skull perfusion
            if isfield(medium.(label_name), 'perfusion')
                perfusion(skull_idx) = medium.(label_name).perfusion;
            end

            % set skull absorption fraction
            if isfield(medium.(label_name), 'absorption_fraction')
                absorption_fraction(skull_idx) = medium.(label_name).absorption_fraction;
            end

            % set skull starting temperature
            if isfield(parameters.thermal.temp_0, label_name)
                temp_0(skull_idx) = parameters.thermal.temp_0.(label_name);
            end

            % map skull bone density with the desired algorithm
            if isfield(parameters, 'pct') && isfield(parameters.pct, 'mapping_density')
                pct_mapping_density = char(parameters.pct.mapping_density);
            else
                pct_mapping_density = 'none';
            end
            [density] = medium_pct_density(parameters, density, pseudoCT, skull_idx, pct_mapping_density);

            % map skull bone sound speed with the desired algorithm
            if isfield(parameters, 'pct') && isfield(parameters.pct, 'mapping_soundspeed')
                pct_mapping_soundspeed = char(parameters.pct.mapping_soundspeed);
            else
                pct_mapping_soundspeed = 'none';
            end
            [sound_speed] = medium_pct_soundspeed(parameters, sound_speed, density, pseudoCT, skull_idx, pct_mapping_soundspeed);

            % map skull bone attenuation with the desired algorithm
            if isfield(parameters, 'pct') && isfield(parameters.pct, 'mapping_attenuation')
                pct_mapping_attenuation = char(parameters.pct.mapping_attenuation);
            else
                pct_mapping_attenuation = 'none';
            end
            [alpha_coeff, alpha_power] = medium_pct_attenuation(parameters, alpha_coeff, alpha_power, pseudoCT, skull_idx, pct_mapping_attenuation);
            
            % [DEBUG] save pCT mapping overview
            if parameters.simulation.debug == 1
                h = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.25, 1]);
                subplot(4,1,1); 
                    hold on; histogram(pseudoCT(skull_idx)); xlabel("pseudo-HU")
                    title(['pseudoCT tissue property ranges']);
                subplot(4,1,2); 
                    hold on; histogram(density(skull_idx)); xlabel("Density [kg/m3]")
                    % add lines for the fixed parameters
                    xline(medium.skull_trabecular.density, 'r', 'LineWidth', 1);
                    xline(medium.skull_cortical.density, 'r', 'LineWidth', 1);
                    xline(medium.skull.density, 'k', 'LineWidth', 2);
                    title(sprintf('Density mapping: %s', pct_mapping_density))
                subplot(4,1,3); 
                    hold on; histogram(sound_speed(skull_idx)); xlabel("Sound speed [m/s]")
                    xline(medium.skull_trabecular.sound_speed, 'r', 'LineWidth', 1);
                    xline(medium.skull_cortical.sound_speed, 'r', 'LineWidth', 1);
                    xline(medium.skull.sound_speed, 'k', 'LineWidth', 2);
                    title(sprintf('Sound speed mapping: %s', pct_mapping_soundspeed))
                subplot(4,1,4); 
                    hold on; histogram(alpha_coeff(skull_idx)); xlabel("Attenuation [dB/(cm.MHzy)]")
                    xline(medium.skull_trabecular.alpha_coeff, 'r', 'LineWidth', 1);
                    xline(medium.skull_cortical.alpha_coeff, 'r', 'LineWidth', 1);
                    xline(medium.skull.alpha_coeff, 'k', 'LineWidth', 2);
                    title(sprintf('Attention mapping: %s', pct_mapping_attenuation))
                output_plot = fullfile(parameters.io.debug_dir_medium,...
                    sprintf('pCT_histograms%s.png',parameters.io.output_affix));
                exportgraphics(h, output_plot, 'Resolution', 150);
                close(h);
            end

            clear skull_idx

        else
            thermal_conductivity(medium_masks==medium_i) = medium.(label_name).thermal_conductivity;
            specific_heat(medium_masks==medium_i) = medium.(label_name).specific_heat_capacity;
            sound_speed(medium_masks==medium_i) = medium.(label_name).sound_speed; 
            density(medium_masks==medium_i) = medium.(label_name).density;
            alpha_coeff(medium_masks==medium_i) = medium.(label_name).alpha_coeff; 
            alpha_power(medium_masks==medium_i) = medium.(label_name).alpha_power;
            if isfield(medium.(label_name), 'perfusion')
                perfusion(medium_masks==medium_i) = medium.(label_name).perfusion;
            end
            if isfield(medium.(label_name), 'absorption_fraction')
                absorption_fraction(medium_masks==medium_i) = medium.(label_name).absorption_fraction;
            end
            if isfield(parameters.thermal.temp_0, label_name)
                temp_0(medium_masks==medium_i) = parameters.thermal.temp_0.(label_name);
            end
        end
    end

    % convert perfusion rate [mL/min/kg] into perfusion coefficient [1/s]
    perfusion_coeff = (perfusion ./ 60) .* density * 1e-6;

    %% account for k-Wave's actual attenuation behaviour
    % limits discrepancies for high attenuation estimates (see https://doi.org/10.1121/1.4894790).
    % 'alpha_coeff' is rescaled to match the specified alpha_coeff and alpa_power_true for the center frequency

    if ~isfield(parameters, 'fit_alpha_power')
        parameters.fit_alpha_power = 1;
    end

    if parameters.fit_alpha_power == 1
        alpha_power_fixed = 2;
    
        if parameters.simulation.debug == 1
            plot_fit = true;
        else
            plot_fit = false;
        end
    
        all_freqs = [parameters.transducer.freq_hz];
        source_freq = all_freqs(1);

        if numel(unique(all_freqs)) > 1
            warning('Transducers have different frequencies (%s Hz). Using %i Hz from first transducer.', ...
                    num2str(unique(all_freqs)), source_freq);
        end
        
        alpha_coeff_fixed = fitPowerLawParamsMulti(...
            alpha_coeff, ...
            alpha_power, ...
            sound_speed, ...
            source_freq, ...
            alpha_power_fixed, ...
            plot_fit);
    
        % DEBUG mode: save plot of fitted attenuation values
        if parameters.simulation.debug == 1
            fig_path = fullfile(parameters.io.debug_dir_medium,...
            ['attenuation_fit', char(parameters.io.output_affix), '.png']);
            saveas(gcf, fig_path);
            close(gcf);
        end
    else
        warning('Attenuation power-law is not remapped (as requested).')
        alpha_coeff_fixed = alpha_coeff;
        alpha_power_fixed = alpha_power;
    end

    %% smooth medium masks

    if isfield(parameters.headmodel, 'smooth_properties') && parameters.headmodel.smooth_properties == true
        disp("Smoothing acoustic property maps ...");

        tmp_density = density; % keep unsmoothed image for figure

        fwhm_mm = parameters.headmodel.smooth_fwhm_mm;
        grid_mm = parameters.grid.resolution_mm;
        smooth_method = parameters.headmodel.smooth_method;

        sound_speed = smooth_img(sound_speed, fwhm_mm, grid_mm, 0, smooth_method);
        density = smooth_img(density, fwhm_mm, grid_mm, 0, smooth_method);
        alpha_coeff_fixed = smooth_img(alpha_coeff_fixed, fwhm_mm, grid_mm, 0, smooth_method);
        thermal_conductivity = smooth_img(thermal_conductivity, fwhm_mm, grid_mm, 0, smooth_method);
        specific_heat = smooth_img(specific_heat, fwhm_mm, grid_mm, 0, smooth_method);
        perfusion_coeff = smooth_img(perfusion_coeff, fwhm_mm, grid_mm, 0, smooth_method);
        absorption_fraction = smooth_img(absorption_fraction, fwhm_mm, grid_mm, 0, smooth_method);
        temp_0 = smooth_img(temp_0, fwhm_mm, grid_mm, 0, smooth_method);
    
        % [DEBUG] Plot unsmoothed and smoothed density
        if parameters.simulation.debug == 1
            h = figure('Position', [100 100 800 400]);
            if numel(size(density))==3
                density_pre = squeeze(tmp_density(:,round(size(tmp_density,2)/2),:));
                density_post = squeeze(density(:,round(size(density,2)/2),:));
            else
                density_pre = squeeze(tmp_density);
                density_post = squeeze(density);
            end
            subplot(1,2,1); imagesc(density_pre); title('Original density')
            subplot(1,2,2); imagesc(density_post); title('Smoothed density')
            output_plot_filename = fullfile(parameters.io.debug_dir_medium,...
                sprintf('sub-%03d_%s_density_smoothing_changes%s.png', ...
                parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
            saveas(h, output_plot_filename, 'png')
            close(h);
            clear density_pre density_post;
        end; clear tmp_density;
    else
        disp('No smoothing applied to acoustic property maps ...')
    end


    %% specify the medium as a kWave-compatible structure 
    % Note: absorption_fraction and temp_0 need to be removed later)
    kwave_medium = struct('sound_speed', sound_speed, ...
                          'density', density, ...
                          'alpha_coeff', alpha_coeff_fixed,...
                          'alpha_power', alpha_power_fixed, ...
                          'thermal_conductivity', thermal_conductivity,...
                          'specific_heat', specific_heat,...
                          'perfusion_coeff', perfusion_coeff,...
                          'absorption_fraction', absorption_fraction,...
                          'temp_0', temp_0);
    
    %% Save acoustic property NIfTIs in T1 space (always for layered; raw matrices debug-only)
    %
    % T1-space maps (medium_properties_nifti): written unconditionally for
    % layered simulations so that they are available for QC and for the
    % uncertainty report without requiring debug mode.
    %
    % Raw simulation-grid matrices (matrix_*.nii.gz): debug-only because they
    % are large, anisotropic, and not directly interpretable without the grid
    % metadata.
    if contains(parameters.simulation.medium, {'layered'}) && exist('planimg','var') && ~isempty(planimg)
        try
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'sound_speed')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'density')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_coeff')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'alpha_power')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'thermal_conductivity')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'specific_heat')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'perfusion_coeff')
            medium_properties_nifti(parameters, kwave_medium, planimg.inv_transf, planimg.t1_header, 'absorption_fraction')
        catch ME
            warning('medium_setup:niftiWrite', ...
                'Could not write medium property NIfTIs: %s', ME.message);
        end
    end

    if parameters.simulation.debug == 1
        % [debug] save raw simulation-grid medium matrices as NIfTIs
        try
            filename_density = fullfile(parameters.io.debug_dir_medium,'matrix_density');
            niftiwrite(density, filename_density, 'Compressed',true); pause(0.1);
            filename_sound_speed = fullfile(parameters.io.debug_dir_medium,'matrix_sound_speed');
            niftiwrite(sound_speed, filename_sound_speed, 'Compressed',true); pause(0.1);
            filename_alpha_coeff = fullfile(parameters.io.debug_dir_medium,'matrix_alpha_coeff');
            niftiwrite(alpha_coeff, filename_alpha_coeff, 'Compressed',true); pause(0.1);
            filename_alpha_power = fullfile(parameters.io.debug_dir_medium,'matrix_alpha_power');
            niftiwrite(alpha_power, filename_alpha_power, 'Compressed',true); pause(0.1);
            filename_alpha_coeff_fixed = fullfile(parameters.io.debug_dir_medium,'matrix_alpha_coeff_fixed');
            niftiwrite(alpha_coeff_fixed, filename_alpha_coeff_fixed, 'Compressed',true); pause(0.1);
            filename_perfusion = fullfile(parameters.io.debug_dir_medium,'matrix_perfusion');
            niftiwrite(perfusion_coeff, filename_perfusion, 'Compressed',true); pause(0.1);
            filename_absorption = fullfile(parameters.io.debug_dir_medium,'matrix_absorption');
            niftiwrite(absorption_fraction, filename_absorption, 'Compressed',true); pause(0.1);
        catch
            warning('medium_setup:debugNifti', ...
                'Error saving debug grid matrices — may result from concurrent write attempts.');
        end

        % [debug] save pCT (if used)
        if parameters.pct.enabled == 1
            filename_pct = fullfile(parameters.io.debug_dir_medium,sprintf('pct%s', parameters.io.output_affix));
            niftiwrite(pseudoCT, filename_pct, 'Compressed',true);
        end
    end

end