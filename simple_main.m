%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MAIN FUNCTION                                %
% This script is used to perform and to debug the simulation.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set medium, submit medium, subject number and configuration file 
clear; close all;
medium = 'layered'; % water or layered
submit = 'slurm'; % run scripts via 'matlab' (debugging) or via a job using 'slurm' (recommended) or 'qsub'
subject_id = 9; % subject number, if none, choose 1

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 

parameters = load_parameters('tmp_test_3-transducers.yaml'); % load the configuration file

parameters.simulation_medium = medium;
if strcmp(submit, 'matlab') == true
    parameters.code_type = 'matlab_cpu';
    parameters.using_donders_hpc = 0;
    single_subject_pipeline(subject_id, parameters);
elseif strcmp(submit, 'qsub') == true
    parameters.interactive = 0;
    parameters.overwrite_files = 'always';
    single_subject_pipeline_with_qsub(subject_id, parameters);
elseif strcmp(submit, 'slurm') == true
    parameters.interactive = 0;
    parameters.overwrite_files = 'always';
    single_subject_pipeline_with_slurm(subject_id, parameters);
else
    error('Submit medium does not correspond to available options.')
end