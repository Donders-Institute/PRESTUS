%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MAIN FUNCTION                                %
% This script is used to perform and to debug the simulation.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set medium, submit medium, subject number and configuration file 
clear; close all;

subject_id = 1; % subject number, if none, choose 1

% add paths
addpath(genpath('functions'))
addpath(genpath('toolboxes')) 

parameters = load_parameters('tutorial_config.yaml'); % load the configuration file

parameters.simulation_medium = 'layered';   % water or layered
parameters.code_type = 'matlab_cpu';
parameters.hpc_submit_medium = 'matlab';    % or 'auto'
parameters.interactive = 1;                 % for interactive debugging
parameters.overwrite_files = 'always';      % overwrite?

prestus_pipeline_start(subject_id, parameters)