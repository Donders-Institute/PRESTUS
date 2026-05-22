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

% Add necessary paths for functions and external
addpath(genpath('functions'));
addpath(genpath('config'));
addpath(genpath('external'));

%% Telemetry consent
if ~exist('~/.prestus/telemetry.json', 'file')
    telemetry_func();
end

%% Load configuration settings

config_folder = ""; % [optional] specify a application-specific config folder
if strcmp(config_folder, "")
    warning("Using configurations in PRESTUS config folder. It is recommended that you work with a local copy instead...");
else
    cd(config_folder) % move to local folder containing configs
end

% Load equipment parameters from individual YAML files under config/equipment/
equip_param = load_equipment_config();
% Equipment information (by default Donders-specific)
parameters = yaml.loadFile('config_default.yaml', 'ConvertToArray', true);
% User-defined calibration parameters
parameters.calibration = yaml.loadFile('config_calibration.yaml');

% Display available equipment combinations (Donders equipment)
available_combos = fieldnames(equip_param.combos);
disp('Available Equipment Combinations:');
disp(available_combos);

% Additional user-specific parameters
% parameters.xxx = '...';

% Set up calibration deployment
calibration_setup(parameters, equip_param);