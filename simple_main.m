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

prestus_pipeline_start(parameters)