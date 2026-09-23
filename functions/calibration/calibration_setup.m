function calibration_setup(parameters, equip_param)

% Set up calibration for multiple setup and depth requests.

% Create output folder (if not yet available)
if ~exist(parameters.calibration.path_output)
    mkdir(parameters.calibration.path_output);
end

% Iterate through equipment combinations
N_i = length(parameters.calibration.combinations);

% Coerce list-of-list fields to cell-of-cell if yaml ConvertToArray
% flattened them into numeric matrices.
%
% Ambiguous case: when every combo's list has exactly one value (e.g.
% N_i combos each requesting a single focal depth), ConvertToArray
% flattens "[[85],[85]]" into a plain 1xN_i (or N_ix1) vector that is
% structurally indistinguishable from "one combo's list of N_i values".
% Since N_i (from `combinations`) is now known, resolve the ambiguity by
% treating a flattened vector whose length equals N_i as one scalar per
% combo -- otherwise a genuinely single-combo multi-value list (N_i==1)
% would be mis-split, and a per-combo single value would otherwise get
% (via the old "replicate length-1 list across all combos" step below)
% incorrectly broadcast to every combo, silently duplicating work.
for coerce_field = {'focal_depths_wrt_exit_plane', 'desired_intensities'}
    f = coerce_field{1};
    val = parameters.calibration.(f);
    if ~iscell(val)
        if N_i > 1 && isvector(val) && numel(val) == N_i
            parameters.calibration.(f) = arrayfun(@(v) {v}, val(:)', 'UniformOutput', false);
        else
            parameters.calibration.(f) = ...
                cellfun(@num2cell, num2cell(val, 2), 'UniformOutput', false);
        end
    elseif ~iscell(val{1})
        parameters.calibration.(f) = ...
            cellfun(@num2cell, val, 'UniformOutput', false);
    end
end

% If a list-of-list field has only one entry, replicate it for all combinations
for coerce_field = {'focal_depths_wrt_exit_plane', 'desired_intensities'}
    f = coerce_field{1};
    if isfield(parameters.calibration, f) && iscell(parameters.calibration.(f)) && ...
            length(parameters.calibration.(f)) == 1 && N_i > 1
        parameters.calibration.(f) = repmat(parameters.calibration.(f), 1, N_i);
    end
end
for i = 1:N_i
    combo_name = parameters.calibration.combinations{i};
    combo = equip_param.combos.(combo_name);
    
    % Extract equipment details
    tran_serial = combo.tran_serial; % - Transducer serial number or identifier
    tran = equip_param.trans.(tran_serial);
    parameters.calibration.equipment_yaml_path = fullfile(equip_param.path, [char(tran_serial), '.yaml']);
    ds = combo.ds;                   % - Driving System metadata (inlined in combo)

    % Configure transducer parameters
    parameters.transducer = tran.transducer;
    parameters.transducer.annular.elem_phase_deg = zeros(1, tran.transducer.annular.elem_n); % initial phases
    parameters.transducer.annular.elem_p_amp = repmat(200000, 1, tran.transducer.annular.elem_n); % initial amplitude

    % Construct output directory for final profiles
    if ~isfield(parameters.calibration, 'path_output_profiles')
        parameters.calibration.path_output_profiles = ...
            fullfile(parameters.calibration.path_output, equip_param.gen.prestus_virt_folder);
    end

    % Extract name of current equipment
    [~, filename, ext] = fileparts(combo.char_data_path);
    equipment_name = erase(filename, equip_param.gen.axial_prof_name);

    % CSV file name of calibrated profile (prefix + equipment)
    parameters.calibration.filename_calibrated_CSV = ...
        strcat(equip_param.gen.prestus_virt_name, equipment_name, '.csv');

    % Load pre-calibrated per-element correction from the transducer equipment YAML.
    % Written there by save_elem_correction after a reference-depth calibration.
    % When present, activates geometric steering mode for all depths in this combo:
    %   phases(depth) = compute_phases(depth) + elem_phase_correction
    % Leave absent to use the default global-search path.
    parameters.calibration.elem_phase_correction_deg = [];
    if isfield(tran, 'elem_phase_correction') && ~isempty(tran.elem_phase_correction)
        parameters.calibration.elem_phase_correction_deg = tran.elem_phase_correction.deg(:)';
        fprintf('Loaded element correction (%d elements, ref depth %.2f mm) from transducer YAML.\n', ...
            numel(parameters.calibration.elem_phase_correction_deg), ...
            tran.elem_phase_correction.ref_depth_ep_mm);
    end

    % Load characterization data
    if ~isempty(parameters.calibration.path_input_axial)
        combo.char_data_path = fullfile(parameters.calibration.path_input_axial, combo.char_data_path);
    end
    charac_data = readmatrix(combo.char_data_path);
    available_foci_wrt_exit_plane = round(charac_data(2, 2:end), 2);
    dist_from_exit_plane = charac_data(3:end, 1);
    intens_data = charac_data(3:end, 2:end);

    % readmatrix pads ragged CSVs with NaN — drop those columns
    valid_cols                    = isfinite(available_foci_wrt_exit_plane);
    available_foci_wrt_exit_plane = available_foci_wrt_exit_plane(valid_cols);
    intens_data                   = intens_data(:, valid_cols);

    % Ensure focal depths are specified. 
    % If no focal depths, perform the simulations for all available focal depths
    if isempty(parameters.calibration.focal_depths_wrt_exit_plane{i})
        parameters.calibration.focal_depths_wrt_exit_plane{i} = available_foci_wrt_exit_plane;
    end

    fprintf('Equipment: transducer %s and driving system %s \n', [tran.name, ds.name])

    % Determine phase calibration mode.
    % 'single_depth' (default): per-depth global search + optional k-Wave validation via
    %   calibration_transducer.
    % 'multi_depth': joint analytical fit across all measured depths,
    %   followed by analytical amplitude calibration — no k-Wave needed.
    if isfield(parameters.calibration, 'phase_mode')
        phase_mode = parameters.calibration.phase_mode;
    else
        phase_mode = 'single_depth';
    end

    % =========================================================================
    % MULTI-DEPTH PATH
    %
    % Runs perform_joint_depth_fit once across all characterization depths,
    % saves the resulting per-element phase correction to the equipment YAML,
    % then calibrates amplitude analytically for each requested depth and
    % intensity. No k-Wave simulations are run (simulated_analytical_scaling=1).
    % Users who require simulation-validated amplitude should follow up with
    % single-depth calibration for specific depths of interest.
    % =========================================================================
    if strcmp(phase_mode, 'multi_depth')

        dist_bowl_offset = tran.transducer.annular.curv_radius_mm - ...
                           tran.transducer.annular.dist_geom_ep_mm;
        dist_bowl_all    = dist_from_exit_plane + dist_bowl_offset;

        desired_intensities_i = cell2mat(parameters.calibration.desired_intensities{i});
        ref_intensity = desired_intensities_i(1);

        fprintf('=== multi_depth: joint fit across %d characterization depths ===\n', ...
            numel(available_foci_wrt_exit_plane));

        [opt_delta_rad, opt_velocity_ref, ~] = perform_joint_depth_fit(...
            parameters, intens_data, available_foci_wrt_exit_plane, ...
            dist_bowl_all, dist_bowl_offset, ref_intensity, tran);

        % Persist correction to the library YAML (not the equipment YAML, which
        % is a read-only hardware spec). Stored as global_model.elem_phase_correction_deg
        % so load_transducer_from_library can apply it on top of geometric phases
        % at any requested depth without the correction being baked into the
        % stored phase table.
        lib_path = fullfile(get_prestus_path(), 'config', 'calibration');
        if isfield(parameters.calibration, 'library_path') && ...
                ~isempty(parameters.calibration.library_path)
            lib_path = parameters.calibration.library_path;
        end
        lib_yaml = fullfile(lib_path, [combo_name, '.yaml']);
        if isfile(lib_yaml)
            lib_data = yaml.loadFile(lib_yaml, 'ConvertToArray', true);
        else
            lib_data = struct('meta', struct('tran_serial', tran_serial, ...
                'created_at', datestr(now, 'yyyy-mm-dd'), ...  %#ok<TNOW1,DATST>
                'calibration_type', 'multi_depth'));
        end
        if ~isfield(lib_data, 'global_model')
            lib_data.global_model = struct();
        end
        lib_data.global_model.elem_phase_correction_deg = opt_delta_rad * 180/pi;
        lib_data.meta.calibration_type = 'multi_depth';
        yaml.dumpFile(lib_yaml, lib_data);
        fprintf('Phase correction (%d elements) written to library: %s\n', ...
            numel(opt_delta_rad), lib_yaml);

        % Analytical amplitude calibration for each requested depth × intensity
        parameters.calibration.equipment_name = equipment_name;
        N_j_md = length(parameters.calibration.focal_depths_wrt_exit_plane{i});
        for j = 1:N_j_md
            focal_distance_ep = round(parameters.calibration.focal_depths_wrt_exit_plane{i}{j}, 2);

            if focal_distance_ep < tran.min_foc || focal_distance_ep > tran.max_foc
                warning('multi_depth: focus %.2f mm outside range [%.2f, %.2f] mm. Skipping.', ...
                    focal_distance_ep, tran.min_foc, tran.max_foc);
                continue;
            end

            % Phases at this depth: geometric + fitted delta
            params_j = parameters;
            params_j.transducer.type = 'annular';
            params_j.transducer.focal_distance_ep   = focal_distance_ep;
            params_j.transducer.focal_distance_bowl = focal_distance_ep + dist_bowl_offset;
            params_j.calibration.desired_focal_distance_ep = focal_distance_ep;
            params_j = load_parameters(params_j);

            geo_phases_j_deg  = set_real_phases(tran, focal_distance_ep, params_j);
            opt_phases_j_deg  = mod(geo_phases_j_deg + opt_delta_rad * 180/pi, 360);
            opt_phases_j_rad  = opt_phases_j_deg * pi / 180;

            params_j.transducer.annular.elem_phase_deg = opt_phases_j_deg;
            params_j.transducer.annular.elem_phase_rad = opt_phases_j_rad;

            % Analytical bootstrap profile at this depth (provides velocity baseline)
            profile_bootstrap_j = extract_analytical_profile(params_j);

            for k = 1:numel(desired_intensities_i)
                desired_intensity_k = desired_intensities_i(k);
                params_j.calibration.desired_intensity = desired_intensity_k;

                % Fit velocity to target intensity (I ∝ v²)
                [profile_focus_ep, ~] = extract_real_intensity_profile(...
                    params_j, available_foci_wrt_exit_plane, focal_distance_ep, ...
                    intens_data, equipment_name, dist_from_exit_plane);
                profile_focus_ep(profile_focus_ep < 0) = 0;
                profile_empirical_jk.axial_intensity    = profile_focus_ep(:);
                profile_empirical_jk.axial_distance_bowl = dist_bowl_all;

                [opt_velocity_k, ~, ~] = fit_velocity_to_intensity(...
                    params_j, profile_bootstrap_j, opt_phases_j_rad, ...
                    opt_velocity_ref, desired_intensity_k);

                % Source amplitude: velocity scaling from bootstrap (no k-Wave correction)
                velocity_scaling = opt_velocity_k / profile_bootstrap_j.velocity;
                opt_amp_k = round(velocity_scaling * params_j.transducer.annular.elem_p_amp(1));

                fprintf('multi_depth: depth=%.1f mm  I=%.0f W/cm²  amp=%d\n', ...
                    focal_distance_ep, desired_intensity_k, opt_amp_k);

                save_params = params_j;
                save_params.transducer.annular.elem_phase_deg = opt_phases_j_deg;
                save_params.transducer.annular.elem_phase_rad = opt_phases_j_rad;
                save_params.transducer.annular.elem_p_amp = repmat(opt_amp_k, 1, ...
                    tran.transducer.annular.elem_n);
                save_params.calibration.desired_focal_distance_ep = focal_distance_ep;
                save_params.calibration.desired_intensity         = desired_intensity_k;
                save_optimized_values(save_params);
            end
        end

        continue;  % skip single-depth loop for this combo
    end

    % Iterate across focal depths (single_depth path)
    N_j = length(parameters.calibration.focal_depths_wrt_exit_plane{i});
    for j = 1:N_j

        focal_distance_ep = round(parameters.calibration.focal_depths_wrt_exit_plane{i}{j}, 2);
        fprintf('Focus from exit plane: %.2f (%s) \n', focal_distance_ep, ...
            focal_reference_label(parameters))

        parameters.transducer.type = 'annular';
        parameters.transducer.focal_distance_ep = focal_distance_ep;

        % Load default parameters incl. additional chosen parameters like transducer
        parameters = load_parameters(parameters);

        % Verify focal range
        if focal_distance_ep < tran.min_foc || focal_distance_ep > tran.max_foc
            warning('Focus %.2f mm is outside the range [%.2f, %.2f] mm. Skipping.\n', ...
                focal_distance_ep, tran.min_foc, tran.max_foc);
            continue;
        end

        % [SIM] Set the manufacturer-specified phases of the transducer for the focal distance
        elem_phase_deg = set_real_phases(tran, focal_distance_ep, parameters);
        parameters.transducer.annular.elem_phase_deg = elem_phase_deg; clear elem_phase_deg;

        % Interpolate or select axial profile
        [profile_focus_ep, max_intens] = extract_real_intensity_profile(...
            parameters,...
            available_foci_wrt_exit_plane, ...
            focal_distance_ep, ...
            intens_data, ...
            equipment_name, ...
            dist_from_exit_plane);

        % figure; plot(dist_from_exit_plane, profile_focus_ep); hold on; xline(focal_distance_ep)

        % Collect all desired intensities for this combo as a vector
        desired_intensities = cell2mat(parameters.calibration.desired_intensities{i});
        desired_focal_distance_ep = focal_distance_ep;

        % assign a loop-specific simulation id (one per focal depth)
        sim_id = (i-1)*N_j + (j-1) + 1;

        % By default, values are reported from the exit plane of the transducer
        % For simulations, we need to implement the distance from the bowl

        % If the profile is taken from the bowl: add the distance from the exit plane (if requested)
        if isfield(parameters.calibration, 'add_FDO') && parameters.calibration.add_FDO == 1
            dist_bowl_exit_plane = parameters.transducer.annular.curv_radius_mm - ...
                parameters.transducer.annular.dist_geom_ep_mm;
        else
            dist_bowl_exit_plane = 0;
        end
        dist_bowl_focus = dist_from_exit_plane + dist_bowl_exit_plane;

        % Set the expected focal distance to exit plane
        % Account for the potential distance between bowl (actual focal distance) and exit plane (axial profile definition)
        parameters.transducer.focal_distance_ep = focal_distance_ep;
        parameters.transducer.focal_distance_bowl = focal_distance_ep + dist_bowl_exit_plane;

        % if the original profiles are measured from the exit plane,
        % we do not model the space until the exit plane; this can lead
        % to suboptimal solutions; assume that the amplitude is zero in
        % that space (i.e., no strong near-field interference)

        % add distance values prior to the exit plane
        if isfield(parameters.calibration, 'add_FDO') && parameters.calibration.add_FDO == 1
            Nvals = round(dist_bowl_focus(1)/(dist_bowl_focus(2)-dist_bowl_focus(1)));
            dist_bowl_focus = cat(1, linspace(0, dist_bowl_focus(1), Nvals)',dist_bowl_focus);
            % set those to zero in amplitude profile
            profile_focus = cat(1, repmat(0, Nvals,1), profile_focus_ep');
        else
            profile_focus = profile_focus_ep;
        end

        % regularize lower values to 0 [possible after interpolation]
        profile_focus(profile_focus<0) = 0;

        % collect data on empirical profile
        profile_empirical.axial_intensity = profile_focus;
        profile_empirical.axial_distance_bowl = dist_bowl_focus;

        % figure; plot(profile_empirical.axial_distance_bowl, profile_empirical.profile_focus);
        % hold on; xline(parameters.focal_distance_bowl)

        % Show current iteration
        disp(['Profiling: ', num2str(sim_id), ' - ', equipment_name{1}, ...
            ' F ', num2str(desired_focal_distance_ep), ' I ', mat2str(desired_intensities)])

        % perform the calibration for all intensities in one call
        % calibration_pipeline_start dispatches to HPC or runs locally
        % depending on parameters.platform ('matlab'/'slurm'/'qsub'/'auto')
        calibration_pipeline_start(...
            profile_empirical, ...
            equipment_name, ...
            desired_intensities, ...
            desired_focal_distance_ep,...
            parameters,...
            sim_id)
    end
end
end

% -------------------------------------------------------------------------
function label = focal_reference_label(parameters)
% Human-readable tag for the depth-setting convention, for console output.
if isfield(parameters, 'calibration') && isfield(parameters.calibration, 'focal_reference')
    focal_reference = parameters.calibration.focal_reference;
else
    focal_reference = 'flhm_3db';
end
switch focal_reference
    case 'flhm_3db'
        label = 'FLHM -3dB';
    case 'flhm_6db'
        label = 'FLHM -6dB';
    case 'peak'
        label = 'peak';
    otherwise
        label = focal_reference;
end
end