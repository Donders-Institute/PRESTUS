function calibration_setup(parameters, equip_param)

% Set up calibration for multiple setup and depth requests.

% TO DO: Make independent from Donders equipment info.

% Create output folder (if not yet available)
if ~exist(parameters.calibration.path_output)
    mkdir(parameters.calibration.path_output);
end

% Iterate through equipment combinations
N_i = length(parameters.calibration.combinations);
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
    parameters.transducer.annular.elem_amp = repmat(200000, 1, tran.transducer.annular.elem_n); % initial amplitude

    % Construct output directory for final profiles
    if ~isfield(parameters.calibration, 'path_output_profiles') || isempty(parameters.calibration)
        parameters.calibration.path_output_profiles = ...
            fullfile(parameters.calibration.path_output, equip_param.gen.prestus_virt_folder);
    end

    % Extract name of current equipment
    [~, filename, ext] = fileparts(combo.char_data_path);
    equipment_name = erase(filename, equip_param.gen.axial_prof_name);

    % CSV file name of calibrated profile (prefix + equipment)
    parameters.calibration.filename_calibrated_CSV = ...
        strcat(equip_param.gen.prestus_virt_name, equipment_name, '.csv');

    % Load manufacturer-provided phase data.
    % Sonic Concepts and Imasonic require an external phase-table file.
    % For any other manufacturer ('generic' or custom), phases are computed
    % geometrically from the ring dimensions in the YAML — no file needed.
    if isequal(tran.manufact, 'Sonic Concepts')
        combo.phase_table = fullfile(parameters.calibration.path_input_phase, combo.phase_table);
        phase_table = readtable(combo.phase_table);
    elseif isequal(tran.manufact, 'Imasonic')
        combo.phase_table = fullfile(parameters.calibration.path_input_phase, combo.phase_table);
        phase_table = read_ini_file(combo.phase_table);
    else
        phase_table = struct();   % unused; set_real_phases calls generate_tran_ini_from_geometry
    end

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
    combo.char_data_path = fullfile(parameters.calibration.path_input_axial, combo.char_data_path);
    charac_data = readmatrix(combo.char_data_path);
    available_foci_wrt_exit_plane = round(charac_data(2, 2:end), 2);
    dist_from_exit_plane = charac_data(3:end, 1);
    intens_data = charac_data(3:end, 2:end);

    % Ensure focal depths are specified. 
    % If no focal depths, perform the simulations for all available focal depths
    if isempty(parameters.calibration.focal_depths_wrt_exit_plane{i})
        parameters.calibration.focal_depths_wrt_exit_plane{i} = available_foci_wrt_exit_plane;
    end

    fprintf('Equipment: transducer %s and driving system %s \n', [tran.name, ds.name])

    % Iterate across focal depths
    N_j = length(parameters.calibration.focal_depths_wrt_exit_plane{i});
    for j = 1:N_j

        focal_distance_ep = round(parameters.calibration.focal_depths_wrt_exit_plane{i}{j}, 2);
        fprintf('Focus from exit plane: %.2f \n', focal_distance_ep)

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
        elem_phase_deg = set_real_phases(phase_table, tran, focal_distance_ep, parameters);
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