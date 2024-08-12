%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            MAIN FUNCTION                                %
% This script is used to perform and to debug the simulation.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set medium, qsub and path parameter 
clear; close all;
medium = 'layered'; % water or layered
qsub = 0; %run scripts via matlab (0) or via a job (1)
cd '/home/neuromod/marcorn/Documents/simnibs4_examples/m2m_ernie/' % change path to demo data here

% add paths
addpath('functions')
addpath(genpath('toolboxes')) 
addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

if gpuDeviceCount==0 && ~exist('/home/common/matlab/fieldtrip/qsub','dir')
    error('Many of the examples in this tutorial assume that you have a GPU available for computations or that you''re using the Donders HPC cluster. It looks like this is not the case. You can still run the tutorial but you''ll need to switch to using CPU (see matlab_code parameter in the config) and it would be slow.')
end

parameters = load_parameters('clus_ernie_dataset_config.yaml'); % load the configuration file

parameters.simulation_medium = medium;
subject_id = 1;
if strcmp(medium, 'layered')
    parameters.transducer.source_amp = [116947	116947	116947	116947];
    parameters.transducer.source_phase_rad = [0 6.28309051155650	0.748390467392312	0.000141417558751905];
    parameters.transducer.source_phase_deg = [0 6.28309051155650	0.748390467392312	0.000141417558751905]/pi*180;
    parameters.results_filename_affix = '_optimized';
end

if qsub == 0
    single_subject_pipeline(subject_id, parameters);
else
    parameters.interactive = 0;
    parameters.overwrite_files = 'always';
    single_subject_pipeline_with_qsub(subject_id, parameters);
end