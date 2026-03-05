%% Main script to deploy transducer calibrations
%
% This script configures and executes calibration simulations 
% using equipment configurations and user-defined input parameters. 
% By performing focal adjustments, intensity scaling, and interpolating
% acoustic profiles, it determines the real axial profile. This profile
% is then used to optimize the virtual transducer's amplitude and phases
% to replicate the behavior of the actual transducer measured in a water tank.
%
% Due to model constraints, additional virtual elements might be required 
% to mimic the actual transducer behavior accurately. The script
% calc_virtual_elem.m in functions/calibration can be used for this purpose.
% 
% This preprocessing step generates optimized amplitude and phase values, 
% which serve as inputs for subsequent acoustic and/or heating simulations.

% Clear workspace and close all figures
close all; clear;
format long; % Ensures accurate display of intensity values

%% Initialization

% Set up paths
% Determine the current and main folder paths
func_path = fullfile(fileparts(mfilename('fullpath')), '..', '..');
main_folder = fileparts(func_path);
cd(main_folder); % Change directory to the main folder

% Add necessary paths for functions and toolboxes
addpath(genpath('functions'));
addpath(genpath('configs'));
addpath(genpath('toolboxes'));

%% Load configuration settings

% Load equipment parameters from the YAML configuration file
equip_param = yaml.loadFile('equipment_config.yaml', 'ConvertToArray', true);
% Equipment information (by default Donders-specific)
parameters = yaml.loadFile('default_config.yaml', 'ConvertToArray', true);
% User-defined calibration parameters
parameters.calibration = yaml.loadFile('calibration_config.yaml');

% Display available equipment combinations (Donders equipment)
available_combos = fieldnames(equip_param.combos);
disp('Available Equipment Combinations:');
disp(available_combos);

%% Iterate through equipment combinations
N_i = length(parameters.calibration.combinations);
for i = 1:N_i
    combo_name = parameters.calibration.combinations{i};
    combo = equip_param.combos.(combo_name);
    
    % Extract equipment details
    ds_serial = combo.ds_serial; % - Driving System serial number or identifier
    ds = equip_param.ds.(ds_serial);
    tran_serial = combo.tran_serial; % - Transducer serial number or identifier
    tran = equip_param.trans.(tran_serial);

    % Configure transducer parameters
    parameters.transducer = tran.prestus.transducer;
    parameters.transducer.source_phase_deg = zeros(1, tran.prestus.transducer.n_elements); % initial phases
    parameters.transducer.source_amp = repmat(200000, 1, tran.prestus.transducer.n_elements); % initial amplitude

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

    % Load manufacturer-provided phase data (only possible if phase tables are available)
    combo.phase_table = fullfile(parameters.calibration.path_input_phase, combo.phase_table);
    if isequal(tran.manufact, "Sonic Concepts")
        phase_table = readtable(combo.phase_table);
    elseif isequal(tran.manufact, "Imasonic")
        phase_table = read_ini_file(combo.phase_table);
    else
        error('Unsupported transducer manufacturer: %s | Please reach out to the developers.', tran.manufact);
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

        parameters.expected_focal_distance_ep = focal_distance_ep;

        % Load default parameters incl. additional chosen parameters like transducer
        parameters = load_parameters(parameters);

        % Verify focal range
        if focal_distance_ep < tran.min_foc || focal_distance_ep > tran.max_foc
            warning('Focus %.2f mm is outside the range [%.2f, %.2f] mm. Skipping.\n', ...
                focal_distance_ep, tran.min_foc, tran.max_foc);
            continue;
        end

        % [SIM] Set the manufacturer-specified phases of the transducer for the focal distance
        source_phase_deg = set_real_phases(phase_table, tran, focal_distance_ep, parameters);
        parameters.transducer.source_phase_deg = source_phase_deg; clear source_phase_deg;

        % Interpolate or select axial profile
        [profile_focus_ep, max_intens] = extract_real_intensity_profile(...
            parameters,...
            available_foci_wrt_exit_plane, ...
            focal_distance_ep, ...
            intens_data, ...
            equipment_name, ...
            dist_from_exit_plane);

        % figure; plot(dist_from_exit_plane, profile_focus_ep); hold on; xline(focal_distance_ep)

        % Iterate across intensities
        N_k = length(parameters.calibration.desired_intensities{i});
        for k = 1:N_k
            desired_intensity = parameters.calibration.desired_intensities{i}{k};
            desired_focal_distance_ep = focal_distance_ep;
            
            % assign a loop-specific simulation id
            sim_id = (i-1)*N_j*N_k + (j-1)*N_k + (k-1) + 1;

            % By default, values are reported from the exit plane of the transducer
            % For simulations, we need to implement the distance from the bowl

            % If the profile is taken from the bowl: add the distance from the exit plane (if requested)
            if isfield(parameters.calibration, 'addEPdistance') && parameters.calibration.addEPdistance == 1
                dist_bowl_exit_plane = parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm;
            else
                dist_bowl_exit_plane = 0;
            end
            dist_bowl_focus = dist_from_exit_plane + dist_bowl_exit_plane;

            % Set the expected focal distance to exit plane 
            % Account for the potential distance between bowl (actual focal distance) and exit plane (axial profile definition)
            parameters.expected_focal_distance_ep = focal_distance_ep;
            parameters.expected_focal_distance_bowl = focal_distance_ep + dist_bowl_exit_plane;

            % if the original profiles are measured from the exit plane,
            % we do not model the space until the exit plane; this can lead
            % to suboptimal solutions; assume that the amplitude is zero in
            % that space (i.e., no strong near-field interference)

            % add distance values prior to the exit plane
            if isfield(parameters.calibration, 'addEPdistance') && parameters.calibration.addEPdistance == 1
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
            % hold on; xline(parameters.expected_focal_distance_bowl)

            % convert from default 3D to 2D axisymmetric simulation (if requested)
            if isfield(parameters.calibration, 'axisymmetric2D') && ...
                    parameters.calibration.axisymmetric2D == 1
                parameters.n_sim_dims = 2;
                parameters.axisymmetric = 1;
                if numel(parameters.default_grid_dims)==3
                    parameters.default_grid_dims(2) = [];
                end
            end

            % Show current iteration
            disp(['Profiling: ', num2str(sim_id), ' - ', equipment_name{1}, ...
                ' F ', num2str(desired_focal_distance_ep), ' I ', num2str(desired_intensity)])

            % perform the calibration
            calibration_transducer(...
                profile_empirical, ...
                equipment_name, ...
                desired_intensity, ...
                desired_focal_distance_ep,...
                parameters,...
                sim_id)

        end
    end
end