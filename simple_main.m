%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MAIN FUNCTION                                %
% This script is used to perform and to debug the simulation.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set medium, submit medium, subject number and configuration file
clear; close all;

% add paths
addpath(genpath('functions'))
addpath(genpath('config'))
addpath(genpath('external'))

parameters = load_parameters('config_tutorial.yaml'); % load the configuration file

parameters.subject_id = 1;                      % subject number
parameters.simulation.medium = 'layered';       % water or layered
parameters.simulation.code_type = 'matlab_cpu';
parameters.simulation.interactive = 1;          % for interactive debugging
parameters.platform = 'matlab';                 % or 'auto'
parameters.io.overwrite_files = 'always';       % overwrite?

% -------------------------------------------------------------------------
% Optional: load transducer phases and amplitude from the calibration
% library instead of using the values in the YAML config.
% -------------------------------------------------------------------------
% combo_name    = 'CTX_250_001_SC_105_010';
% focal_dist_ep = 50;   % mm from exit plane
% intensity     = 30;   % W/cm²
% equip_param   = load_equipment_config();
% lib_params    = load_transducer_from_library(combo_name, focal_dist_ep, intensity, equip_param);
% parameters.transducer = lib_params.transducer;

prestus_pipeline_start(parameters)