%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MAIN FUNCTION                                %
% This script is used to perform and to debug the simulation.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set medium, submit medium, subject number and configuration file 
clear; close all;
medium = 'layered'; % water or layered
submit = 'slurm'; % run scripts via 'matlab' (debugging) or via a job using 'slurm' (recommended) or 'qsub'
subject_id = 1; % subject number, if none, choose 1

% add paths
addpath(genpath('functions'))
addpath(genpath('toolboxes')) 
addpath(genpath('configs')) 

parameters = load_parameters('tutorial_config.yaml'); % load the configuration file

parameters.simulation_medium = medium;
parameters.hpc_submit_medium = submit;

if strcmp(parameters.hpc_submit_medium, 'matlab') == true
    parameters.code_type = 'matlab_cpu';
    single_subject_pipeline(subject_id, parameters);
elseif strcmp(parameters.hpc_submit_medium, 'qsub') == true
    parameters.interactive = 0;
    parameters.overwrite_files = 'always';
    single_subject_pipeline_with_qsub(subject_id, parameters);
elseif strcmp(parameters.hpc_submit_medium, 'slurm') == true
    parameters.interactive = 0;
    parameters.overwrite_files = 'always';
    single_subject_pipeline_with_slurm(subject_id, parameters);
else
    error('Submit medium does not correspond to available options.')
end