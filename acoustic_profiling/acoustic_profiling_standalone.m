%% Main script to deploy acoustic profiling simulations
% This script configures and executes acoustic profiling simulations 
% using equipment configurations and user-defined input parameters. 
% By performing focal adjustments, intensity scaling, and interpolating
% acoustic profiles, it determines the real axial profile. This profile
% is then used to optimize the virtual transducer's amplitude and phases
% to replicate the behavior of the actual transducer measured in a water tank.
% Due to model constraints, additional virtual elements might be required 
% to mimic the actual transducer behavior accurately. The script
% calc_virtual_elem.m in acoustic_profiling/functions can be used for this purpose.
% 
% This preprocessing step generates optimized amplitude and phase values, 
% which serve as inputs for subsequent acoustic and/or heating simulations.

% Clear workspace and close all figures
close all; clear;
format long; % Ensures accurate display of intensity values

%% Initialization

% Set up paths
% Determine the current and main folder paths
func_path = fileparts(mfilename('fullpath'));
main_folder = fileparts(func_path);
cd(main_folder); % Change directory to the main folder

% Add necessary paths for functions and toolboxes
addpath('functions');
addpath(genpath('acoustic_profiling'));
addpath(genpath('configs'));
addpath(genpath('toolboxes'));

% Load configuration
% Load equipment parameters from the YAML configuration file
equip_param = yaml.loadFile('equipment_config.yaml', 'ConvertToArray', true);
parameters = yaml.loadFile('default_config.yaml', 'ConvertToArray', true);
parameters.calibration = yaml.loadFile('calibration_config.yaml', 'ConvertToArray', true);

% Display available equipment combinations (Donders equipment)
available_combos = fieldnames(equip_param.combos);
disp('Available Equipment Combinations:');
disp(available_combos);

%% User-defined input parameters
% Define data path and simulation medium
% Path containing 'Axial profiles' and 'PRESTUS virtual parameter' folders 
% of /home/common/matlab/PRESTUS.
% data_path = '\\ru.nl\WrkGrp\FUS_Researchers\';

% Specify equipment combinations and simulation parameters
parameters.calibration.combinations = {...
    'IS_PCD15287_01001_IGT_128_ch_comb_10_ch', ...
    'IS_PCD15287_01002_IGT_128_ch_comb_10_ch'};

% Every [] will be performed for same index equipment combination. Example:
% all foci within the first [] will be executed for the first equipment 
% combination. If array is empty ([]), focal depths of available 
% characterization data will be used.
parameters.calibration.focal_depths_wrt_exit_plane = {[85], [30, 50]}; % Focal depths in mm

% Every [] will be performed for same index equipment combination and for 
% each foci in same index focal depth set []. Example: for the first
% equipment combination, all foci within the first [] will be executed per
% intensity within the first desired intensity set [].
parameters.calibration.desired_intensities = {[30], [30, 60]}; % Desired intensities in W/cm^2

%% Iterate through equipment combinations
for i = 1:length(parameters.calibration.combinations)
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
    combo.phase_table = fullfile(parameters.calibration.path_input_profiles, combo.phase_table);
    if isequal(tran.manufact, "Sonic Concepts")
        phase_table = readtable(combo.phase_table);
    elseif isequal(tran.manufact, "Imasonic")
        phase_table = read_ini_file(combo.phase_table);
    else
        error('Unsupported transducer manufacturer: %s | Please reach out to the developers.', tran.manufact);
    end

    % Load characterization data
    combo.char_data_path = fullfile(parameters.calibration.path_input_profiles, combo.char_data_path);
    charac_data = readmatrix(combo.char_data_path);
    available_foci_wrt_exit_plane = round(charac_data(2, 2:end), 2);
    dist_from_exit_plane = charac_data(3:end, 1);
    intens_data = charac_data(3:end, 2:end);

    % Convert exit plane reference to mid-bowl reference
    % Simulation results are w.r.t. mid-bowl of transducer
    % Profiles are referenced to 0 = exit plane.
    % If the profile is taken from the bowl: add the distance from the exit plane (if requested)
    if isfield(parameters.calibration, 'addEPdistance') && parameters.calibration.addEPdistance == 1
        dist_tran_exit_plane = parameters.transducer.curv_radius_mm - parameters.transducer.dist_to_plane_mm;
    else
        dist_tran_exit_plane = 0;
    end
    dist_from_tran = dist_from_exit_plane + dist_tran_exit_plane;

    % Ensure focal depths are specified. 
    % If no focal depths, perform the simulations for all available focal depths
    if isempty(parameters.calibration.focal_depths_wrt_exit_plane{i})
        parameters.calibration.focal_depths_wrt_exit_plane{i} = available_foci_wrt_exit_plane;
    end

    % Load default parameters incl. additional chosen parameters like transducer
    parameters = load_parameters(parameters);

    fprintf('Equipment: transducer %s and driving system %s \n', [tran.name, ds.name])

    % Iterate across focal depths
    for j = 1:length(parameters.calibration.focal_depths_wrt_exit_plane{i})
        focus_wrt_exit_plane = round(parameters.calibration.focal_depths_wrt_exit_plane{i}(j), 2);
        fprintf('Focus: %.2f \n', focus_wrt_exit_plane)

        % Verify focal range
        if focus_wrt_exit_plane < tran.min_foc || focus_wrt_exit_plane > tran.max_foc
            warning('Focus %.2f mm is outside the range [%.2f, %.2f] mm. Skipping.\n', ...
                focus_wrt_exit_plane, tran.min_foc, tran.max_foc);
            continue;
        end

        % [SIM] Set the expected focal distance to exit plane distance 
        % Account for the potential distance between bowl (actual focal
        % distance) and exit plane (axial profile definition)
        parameters.expected_focal_distance_mm = focus_wrt_exit_plane + dist_tran_exit_plane;
        parameters.expected_focal_distance_EP_mm = focus_wrt_exit_plane;

        % [SIM] Set the manufacturer-specified phases of the transducer for the focal distance
        source_phase_deg = set_real_phases(phase_table, tran, focus_wrt_exit_plane, parameters);
        parameters.transducer.source_phase_deg = source_phase_deg; clear source_phase_deg;

        % Interpolate or select axial profile
        [profile_focus, max_intens] = extract_real_intensity_profile(...
            parameters,...
            available_foci_wrt_exit_plane, ...
            focus_wrt_exit_plane, ...
            intens_data, ...
            equipment_name, ...
            dist_from_tran);

        % Iterate across intensities
        for k = 1:length(parameters.calibration.desired_intensities{i})
            desired_intensity = parameters.calibration.desired_intensities{i}(k);
            
            % assign a loop-specific simulation id
            sim_id = i*j*k;

            % if the original profiles are measured from the exit plane,
            % we do not model the space until the exit plane; this can lead
            % to suboptimal solutions; assume that the amplitude is zero in
            % that space (i.e., no strong near-field interference)

            % add distance values prior to the exit plane
            Nvals = round(dist_from_tran(1)/(dist_from_tran(2)-dist_from_tran(1)));
            dist_from_tran = cat(1, ...
                linspace(0, dist_from_tran(1), Nvals)',...
                dist_from_tran);
            % set those to initial value in amplitude profile
            profile_focus = cat(1, repmat(profile_focus(1),Nvals,1), profile_focus');

            % collect data on empirical profile
            profile_empirical.profile_focus = profile_focus;
            profile_empirical.dist_from_tran = dist_from_tran;
            profile_empirical.focus_wrt_exit_plane = focus_wrt_exit_plane;

            % perform the acoustic profiling
            acoustic_profiling(...
                profile_empirical, ...
                equipment_name, ...
                desired_intensity, ...
                parameters,...
                sim_id)

        end
    end
end