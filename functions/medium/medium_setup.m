function [kwave_medium, medium_plus] = medium_setup(parameters, medium_masks, planimg, pseudoCT)
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
% saved only in debug mode. When pCT is enabled, empirical in-skull property
% distributions (pseudoCT HU, density, sound speed, attenuation) are cached
% unconditionally as JSON in io.dir_cache (sub-XXX_<medium>_medium_stats<affix>.json)
% for use by downstream tools (e.g. the Nexus GUI).
%
% Use as:
%   kwave_medium = medium_setup(parameters, medium_masks, planimg)
%   kwave_medium = medium_setup(parameters, medium_masks, planimg, pseudoCT)
%
% Input:
%   parameters   - PRESTUS config; must contain medium_properties, layers,
%                  grid.dims, grid.resolution_mm [mm], pct.enabled,
%                  transducer(1).freq_hz [Hz], and optionally fit_alpha_power (default 1)
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
%                  perfusion_coeff [1/s]
%                  (the absorption map and temp_0 are returned in medium_plus)
%   medium_plus  - struct with fields: temp_0 [°C], absorption_np [Np/m] at
%                  the source frequency, absorption_fraction [-] (absorption /
%                  attenuation), both effective for every tissue whatever its
%                  THERMAL_ABSORPTION_MODEL, and absorption_from_coeff (true
%                  when any tissue uses 'coefficient'; the thermal stage then
%                  uses absorption_np, otherwise absorption_fraction)
%                  (kept separate because k-Wave rejects these fields in
%                  kwave_medium during acoustic simulation)
%
% See also: MEDIUM_PCT_DENSITY, MEDIUM_PCT_SOUNDSPEED, MEDIUM_PCT_ATTENUATION,
%   FITPOWERLAWPARAMSMULTI, NIFTI_MEDIUM_PROPERTIES

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
    absorption_coeff = empty_grid; % [dB/(cm MHz^y)], tissues with absorption model 'coefficient'
    uses_absorption_coeff = false(parameters.grid.dims);
    temp_0 = empty_grid;

    % Get layer and medium labels
    layer_labels = fieldnames(parameters.layers);
    medium_labels = fieldnames(parameters.medium_properties);

    % Iterate through each layer & assign medium ID to create medium mask
    for label_i = 1:length(layer_labels)
        label_name = layer_labels{label_i};
        medium_i = find(strcmp(medium_labels, label_name));
        if pct_skull_mapped(parameters) && strcmp(label_name, 'skull')
            skull_idx = find(ismember(medium_masks,medium_i));

            % set skull thermal conductivity
            thermal_conductivity(skull_idx) = medium.(label_name).thermal_conductivity;

            % set skull heat capacity     
            specific_heat(skull_idx) = medium.(label_name).specific_heat_capacity;

            % set skull perfusion
            if isfield(medium.(label_name), 'perfusion')
                perfusion(skull_idx) = medium.(label_name).perfusion;
            end

            % set skull absorption fraction / coefficient
            [absorption_fraction, absorption_coeff, uses_absorption_coeff] = set_absorption( ...
                parameters, medium, label_name, skull_idx, absorption_fraction, absorption_coeff, uses_absorption_coeff);

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
            [alpha_coeff, alpha_power] = medium_pct_attenuation(parameters, alpha_coeff, alpha_power, pseudoCT, skull_idx, pct_mapping_attenuation, density);

            % [CACHE] save empirical in-skull property distributions (pseudoCT HU, density,
            % sound speed, attenuation) so downstream tools (e.g. the Nexus GUI) can show
            % them without a debug rerun. Written unconditionally, since it is cheap relative
            % to the per-voxel arrays already computed above; goes to dir_cache (regenerable
            % intermediates, always created) rather than the debug-only dir_debug_medium.
            medium_stats = struct();
            medium_stats.subject_id    = parameters.subject_id;
            medium_stats.medium        = parameters.simulation.medium;
            medium_stats.output_affix  = char(parameters.io.output_affix);
            medium_stats.layer         = label_name;
            medium_stats.n_voxels      = numel(skull_idx);
            medium_stats.mapping       = struct('density', pct_mapping_density, ...
                                                 'sound_speed', pct_mapping_soundspeed, ...
                                                 'attenuation', pct_mapping_attenuation);
            if strcmp(pct_mapping_density, 'k-plan')
                [~, ~, medium_stats.mapping.density_calibration] = medium_pct_density_calibration(parameters);
            end
            medium_stats.properties.pseudoCT_HU = distribution_stats(pseudoCT(skull_idx));
            medium_stats.properties.density     = distribution_stats(density(skull_idx));
            medium_stats.properties.sound_speed = distribution_stats(sound_speed(skull_idx));
            medium_stats.properties.alpha_coeff = distribution_stats(alpha_coeff(skull_idx));
            medium_stats.reference = struct();
            for ref_name = {'skull', 'skull_cortical', 'skull_trabecular'}
                rn = ref_name{1};
                if isfield(medium, rn)
                    medium_stats.reference.(rn) = struct( ...
                        'density', medium.(rn).density, ...
                        'sound_speed', medium.(rn).sound_speed, ...
                        'alpha_coeff', medium.(rn).alpha_coeff);
                end
            end
            try
                filename_medium_stats = fullfile(parameters.io.dir_cache, ...
                    sprintf('sub-%03d_%s_medium_stats%s.json', ...
                    parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix));
                fid = fopen(filename_medium_stats, 'w');
                fwrite(fid, jsonencode(medium_stats, 'PrettyPrint', true));
                fclose(fid);
            catch ME
                warn('medium_setup:statsCacheWrite', ...
                    'Could not write medium property stats cache: %s', ME.message);
            end
            clear medium_stats

            % [DEBUG] save pCT mapping overview
            if debug_on(parameters, 'medium', 'figures')
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
                    title(sprintf('Attenuation mapping: %s', pct_mapping_attenuation))
                debug_save_figure(parameters, 'medium', h, ...
                    sprintf('pCT_histograms%s', parameters.io.output_affix), ...
                    'CT skull voxels: HU and the mapped density, sound speed and attenuation', 'Resolution', 150);
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
            [absorption_fraction, absorption_coeff, uses_absorption_coeff] = set_absorption( ...
                parameters, medium, label_name, medium_masks==medium_i, absorption_fraction, absorption_coeff, uses_absorption_coeff);
            if isfield(parameters.thermal.temp_0, label_name)
                temp_0(medium_masks==medium_i) = parameters.thermal.temp_0.(label_name);
            end
        end
    end

    % convert perfusion rate [mL/min/kg] into perfusion coefficient [1/s]
    perfusion_coeff = (perfusion ./ 60) .* density * 1e-6;

    % effective absorption in [Np/m] at the source frequency and as a fraction
    % of attenuation, using each tissue's attenuation power law
    % (absorption_power = alpha_power): a coefficient tissue gets its fraction,
    % a fraction tissue its coefficient
    db2np_f = db_cm_mhz_to_np_m(1, parameters.transducer(1).freq_hz, alpha_power); % [dB/(cm MHz^y)] -> [Np/m] at f
    attenuation_np = alpha_coeff .* db2np_f;
    absorption_np = absorption_fraction .* attenuation_np;
    absorption_np(uses_absorption_coeff) = absorption_coeff(uses_absorption_coeff) .* db2np_f(uses_absorption_coeff);
    % fraction of a coefficient tissue; without attenuation it is 0 when nothing
    % is absorbed and 1 otherwise (the warning below flags that case), so the
    % map stays finite for smoothing and the NIfTI outputs
    from_coeff = uses_absorption_coeff & attenuation_np > 0;
    absorption_fraction(from_coeff) = absorption_np(from_coeff) ./ attenuation_np(from_coeff);
    no_attenuation = uses_absorption_coeff & ~(attenuation_np > 0);
    absorption_fraction(no_attenuation) = double(absorption_np(no_attenuation) > 0);
    clear from_coeff no_attenuation
    if any(absorption_np(:) > attenuation_np(:) * (1 + 1e-9))
        warn('medium_setup:absorptionExceedsAttenuation', ...
            'absorption exceeds attenuation in %d voxels; heating then exceeds the attenuated power.', ...
            nnz(absorption_np > attenuation_np * (1 + 1e-9)));
    end
    absorption_from_coeff = any(uses_absorption_coeff(:));
    absorption_log = absorption_summary(parameters, medium_masks, layer_labels, medium_labels, ...
        absorption_np, absorption_fraction, db2np_f);
    clear absorption_coeff attenuation_np db2np_f uses_absorption_coeff

    %% account for k-Wave's actual attenuation behaviour
    % limits discrepancies for high attenuation estimates (see https://doi.org/10.1121/1.4894790).
    % 'alpha_coeff' is rescaled to match the specified alpha_coeff and alpha_power at the center frequency

    if ~isfield(parameters, 'fit_alpha_power')
        parameters.fit_alpha_power = 1;
    end

    % jwave uses a fixed y=2 Stokes absorption law internally, so alpha_coeff
    % must always be remapped to alpha_power_fixed=2.  Override fit_alpha_power
    % when the jwave backend is active, regardless of the config value.
    if isfield(parameters, 'acoustic_backend') && strcmpi(parameters.acoustic_backend, 'jwave') ...
            && parameters.fit_alpha_power ~= 1
        warning('medium_setup:jwaveAlphaPower', ...
            'fit_alpha_power=%d overridden to 1 for jwave backend (requires alpha_power_fixed=2).', ...
            parameters.fit_alpha_power);
        parameters.fit_alpha_power = 1;
    end

    if parameters.fit_alpha_power == 1
        alpha_power_fixed = 2;
    
        plot_fit = debug_on(parameters, 'medium', 'figures');
    
        all_freqs = [parameters.transducer.freq_hz];
        source_freq = all_freqs(1);

        if numel(unique(all_freqs)) > 1
            warn('Transducers have different frequencies (%s Hz). Using %i Hz from first transducer.', ...
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
        if plot_fit
            debug_save_figure(parameters, 'medium', gcf, ...
                ['attenuation_fit', char(parameters.io.output_affix)], ...
                sprintf('Attenuation refitted to a single power law (alpha_power = %g)', alpha_power_fixed));
        end
    else
        warn('Attenuation power-law is not remapped (as requested).')
        alpha_coeff_fixed = alpha_coeff;
        % k-Wave accepts a single power law only (medium.alpha_power must be scalar)
        y = unique(alpha_power(~isnan(alpha_power)));
        if isscalar(y)
            alpha_power_fixed = y;
        else
            error('medium_setup:alphaPowerNotUniform', ...
                ['fit_alpha_power = 0 needs one alpha_power for all tissues (found %s); ' ...
                 'k-Wave supports a single power law. Use fit_alpha_power = 1 to remap.'], ...
                mat2str(y(:)'));
        end
    end

    %% smooth medium masks

    if isfield(parameters.headmodel, 'smooth_properties') && parameters.headmodel.smooth_properties == true
        disp("Smoothing acoustic property maps ...");

        tmp_density = density; % keep unsmoothed image for figure

        fwhm_mm = parameters.headmodel.smooth_fwhm_mm;
        grid_mm = parameters.grid.resolution_mm;
        smooth_method = parameters.headmodel.smooth_method;

        sigma_vox = fwhm_mm ./ (grid_mm * 2.35482);
        fprintf('Smoothing medium properties: FWHM=[%s] mm → σ=[%s] vox\n', num2str(fwhm_mm), num2str(sigma_vox));

        sound_speed = smooth_img(sound_speed, fwhm_mm, grid_mm, 0, smooth_method);
        density = smooth_img(density, fwhm_mm, grid_mm, 0, smooth_method);
        alpha_coeff_fixed = smooth_img(alpha_coeff_fixed, fwhm_mm, grid_mm, 0, smooth_method);
        thermal_conductivity = smooth_img(thermal_conductivity, fwhm_mm, grid_mm, 0, smooth_method);
        specific_heat = smooth_img(specific_heat, fwhm_mm, grid_mm, 0, smooth_method);
        perfusion_coeff = smooth_img(perfusion_coeff, fwhm_mm, grid_mm, 0, smooth_method);
        absorption_np = smooth_img(absorption_np, fwhm_mm, grid_mm, 0, smooth_method);
        absorption_fraction = smooth_img(absorption_fraction, fwhm_mm, grid_mm, 0, smooth_method);
        temp_0 = smooth_img(temp_0, fwhm_mm, grid_mm, 0, smooth_method);
    
        % [DEBUG] Plot unsmoothed and smoothed density
        if debug_on(parameters, 'medium', 'figures')
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
            debug_save_figure(parameters, 'medium', h, ...
                sprintf('sub-%03d_%s_density_smoothing_changes%s', ...
                parameters.subject_id, parameters.simulation.medium, parameters.io.output_affix), ...
                'Density before (left) and after (right) smoothing the property maps, centre y-slice');
            clear density_pre density_post;
        end; clear tmp_density;
    else
        disp('No smoothing applied to acoustic property maps ...')
    end


    %% specify the medium as a kWave-compatible structure
    kwave_medium = struct('sound_speed', sound_speed, ...
                          'density', density, ...
                          'alpha_coeff', alpha_coeff_fixed,...
                          'alpha_power', alpha_power_fixed, ...
                          'thermal_conductivity', thermal_conductivity,...
                          'specific_heat', specific_heat,...
                          'perfusion_coeff', perfusion_coeff,...
                          'temp_0', temp_0);

    % k-Wave's dispersion term is undefined at y = 1 (tan(pi*y/2)); a linear
    % law (fit_alpha_power = 0 with alpha_power 1, e.g. BabelBrain-matched
    % media) runs with absorption but without dispersion, as k-Wave advises.
    % The compiled k-Wave solvers (C++/CUDA) read the power law from their input
    % file and reject alpha_power = 1 whatever alpha_mode is set, so stop here
    % instead of inside k-Wave.
    if alpha_power_fixed == 1
        backend = lower(char(string(get_field(parameters, {'acoustic_backend'}, 'kwave'))));
        code_type = lower(char(string(get_field(parameters, {'simulation', 'code_type'}, ''))));
        runs_acoustics = logical(get_field(parameters, {'modules', 'run_acoustic_sims'}, 1));
        if strcmp(backend, 'kwave') && startsWith(code_type, 'cpp') && runs_acoustics
            error('medium_setup:alphaPowerOneCpp', ['alpha_power = 1 (fit_alpha_power = 0): the k-Wave C++/CUDA ' ...
                'solvers do not accept a linear attenuation law. Use fit_alpha_power = 1 (refit to y = 2), ' ...
                'a common alpha_power other than 1, or a MATLAB code_type (matlab_cpu / matlab_gpu).']);
        end
        kwave_medium.alpha_mode = 'no_dispersion';
        warn('alpha_power = 1: k-Wave runs with alpha_mode = ''no_dispersion'' (absorption kept, dispersion off).');
    end

    fprintf('[medium_setup] ---- medium summary ----\n');
    fprintf('[medium_setup] Medium type      : %s\n', parameters.simulation.medium);
    fprintf('[medium_setup] Tissue layers    : %d (%s)\n', numel(layer_labels), strjoin(layer_labels, ', '));
    fprintf('[medium_setup] sound_speed      : %.1f – %.1f m/s\n', min(sound_speed(:), [], 'omitnan'), max(sound_speed(:), [], 'omitnan'));
    fprintf('[medium_setup] density          : %.1f – %.1f kg/m3\n', min(density(:), [], 'omitnan'), max(density(:), [], 'omitnan'));
    fprintf('[medium_setup] alpha_coeff      : %.4f – %.4f\n', min(alpha_coeff_fixed(:), [], 'omitnan'), max(alpha_coeff_fixed(:), [], 'omitnan'));
    fprintf('%s', absorption_log);
    if pct_skull_mapped(parameters)
        pct_d = 'none'; pct_s = 'none'; pct_a = 'none';
        if isfield(parameters.pct, 'mapping_density'),     pct_d = char(parameters.pct.mapping_density);     end
        if isfield(parameters.pct, 'mapping_soundspeed'),  pct_s = char(parameters.pct.mapping_soundspeed);  end
        if isfield(parameters.pct, 'mapping_attenuation'), pct_a = char(parameters.pct.mapping_attenuation); end
        fprintf('[medium_setup] pCT mapping      : yes (density=%s, soundspeed=%s, attenuation=%s)\n', pct_d, pct_s, pct_a);
        if strcmp(pct_d, 'k-plan')
            [~, ~, cal_src] = medium_pct_density_calibration(parameters);
            fprintf('[medium_setup] CT calibration   : %s\n', cal_src);
        end
    else
        fprintf('[medium_setup] pCT mapping      : no\n');
    end
    if isfield(parameters.headmodel, 'smooth_properties') && parameters.headmodel.smooth_properties == true
        fprintf('[medium_setup] Smoothing        : yes\n');
    else
        fprintf('[medium_setup] Smoothing        : no\n');
    end
    fprintf('[medium_setup] ------------------------------------\n');

    %% Save raw simulation-grid medium matrices (debug only)
    %
    % NOTE: T1-space property maps (back-transformed, with output_affix, in
    % dir_output/ap/) are written by nifti_medium_properties when
    % io.save_property_maps is set. The matrices here are complementary:
    % they are in raw simulation-grid space, named matrix_*, and only
    % interpretable alongside the grid metadata. They are NOT duplicates.
    %
    % Raw simulation-grid matrices → debug_dir/medium/ (matrix_*.nii.gz):
    %   Written only with debug volumes on (DEBUG_ON). Large, anisotropic, and only
    %   interpretable with the grid metadata — not intended as pipeline outputs.
    if debug_on(parameters, 'medium', 'volumes')
        % [debug] save raw simulation-grid medium matrices as NIfTIs
        mats = {
            'matrix_density',           density,             'Density [kg/m^3], simulation grid'
            'matrix_sound_speed',       sound_speed,         'Sound speed [m/s], simulation grid'
            'matrix_alpha_coeff',       alpha_coeff,         'Attenuation coefficient [dB/(cm MHz^y)] before the power-law refit'
            'matrix_alpha_power',       alpha_power,         'Attenuation power-law exponent per tissue'
            'matrix_alpha_coeff_fixed', alpha_coeff_fixed,   'Attenuation coefficient after the power-law refit (used by k-Wave)'
            'matrix_perfusion',         perfusion_coeff,     'Perfusion coefficient, simulation grid'
            'matrix_absorption_np',     absorption_np,       'Absorption [Np/m], simulation grid'
            'matrix_absorption',        absorption_fraction, 'Absorption fraction of attenuation, simulation grid'
            };
        try
            for mi = 1:size(mats, 1)
                fname = fullfile(parameters.io.dir_debug_medium, mats{mi, 1});
                niftiwrite(mats{mi, 2}, fname, 'Compressed', true); pause(0.1);
                debug_register(parameters, 'medium', fname, mats{mi, 3}, 'volume');
            end
        catch
            warn('medium_setup:debugNifti', ...
                'Error saving debug grid matrices — may result from concurrent write attempts.');
        end
        clear mats

        % [debug] save pCT (if used)
        if parameters.pct.enabled == 1
            filename_pct = fullfile(parameters.io.dir_debug_medium,sprintf('pct%s', parameters.io.output_affix));
            niftiwrite(pseudoCT, filename_pct, 'Compressed',true);
            debug_register(parameters, 'medium', filename_pct, 'CT [HU], simulation grid', 'volume');
        end
    end

    %% Split fields incompatible with k-Wave acoustic validation
    medium_plus.temp_0 = kwave_medium.temp_0;
    medium_plus.absorption_np = absorption_np;
    medium_plus.absorption_fraction = absorption_fraction;
    medium_plus.absorption_from_coeff = absorption_from_coeff;
    kwave_medium = rmfield(kwave_medium, 'temp_0');

end

function [absorption_fraction, absorption_coeff, uses_absorption_coeff] = set_absorption( ...
        parameters, medium, label_name, idx, absorption_fraction, absorption_coeff, uses_absorption_coeff)
% Absorption of one tissue (voxels idx) by its THERMAL_ABSORPTION_MODEL
    a = thermal_absorption_properties(parameters, label_name);
    if a.forced
        warn('medium_setup:absorptionModelForced', ...
            ['the %s absorption_coeff has no pseudo-CT counterpart (pct.mapping_attenuation = ''%s''); ' ...
             'applying absorption_coeff / alpha_coeff = %.3g as its fraction of the mapped attenuation.'], ...
            label_name, char(parameters.pct.mapping_attenuation), a.absorption_fraction);
    end
    if strcmp(a.model, 'coefficient')
        absorption_coeff(idx) = get_absorption_coeff(medium, label_name);
        uses_absorption_coeff(idx) = true;
    else
        absorption_fraction(idx) = a.absorption_fraction;
    end
end

function s = absorption_summary(parameters, medium_masks, layer_labels, medium_labels, ...
        absorption_np, absorption_fraction, db2np_f)
% Log lines: per tissue the model and the effective coefficient (1 MHz
% prefactor and Np/m at the source frequency) and fraction; a range where
% the values vary (pseudo-CT skull)
    s = '';
    absorption_coeff = absorption_np ./ db2np_f; % [dB/(cm MHz^y)]
    for label_i = 1:numel(layer_labels)
        label_name = layer_labels{label_i};
        idx = medium_masks == find(strcmp(medium_labels, label_name));
        if ~any(idx(:))
            continue
        end
        [model, forced] = thermal_absorption_model(parameters, label_name);
        if forced
            model = 'fraction (from coeff)';
        end
        s = [s sprintf('[medium_setup] Absorption %-17s: %-17s coeff %s dB/(cm MHz^y) (%s Np/m), fraction %s\n', ...
            label_name, model, value_range(absorption_coeff(idx), '%.4g'), ...
            value_range(absorption_np(idx), '%.4g'), value_range(absorption_fraction(idx), '%.3g'))]; %#ok<AGROW>
    end
end

function s = value_range(x, fmt)
% 'v' for a constant, 'min – max' otherwise
    lo = min(x(:), [], 'omitnan');
    hi = max(x(:), [], 'omitnan');
    if isempty(lo) || isnan(lo)
        s = 'n/a';
    elseif abs(hi - lo) <= 1e-9 * max(abs(hi), eps)
        s = sprintf(fmt, lo);
    else
        s = sprintf([fmt ' – ' fmt], lo, hi);
    end
end

function a = get_absorption_coeff(medium, label_name)
% absorption_coeff of a tissue that sets it: a finite, non-negative scalar
    a = medium.(label_name).absorption_coeff;
    if ~isnumeric(a) || ~isscalar(a) || ~isfinite(a) || a < 0
        error('medium_setup:badAbsorptionCoeff', ...
            'medium_properties.%s.absorption_coeff must be a non-negative number [dB/(cm MHz^y)].', label_name);
    end
end
