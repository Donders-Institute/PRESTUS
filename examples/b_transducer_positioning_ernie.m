%% Pipeline to obtain sensible transducer positions.
% It assumes that you have already did the segmentation for your subjects in SimNIBS (can be done through the main pipeline)
% The script uses gpuArray so run it on a CUDA-enabled machine. 

restoredefaultpath; clear all; clc

% get root path (script must be run)
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile); 
cd(fullfile(pathstr,'..'))
rootpath = pwd;

pn.tuSIM = fullfile(rootpath, 'tools', 'PRESTUS');
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions'); addpath(pn.tuSIM_fun);
pn.tuSIM_tools = fullfile(pn.tuSIM, 'toolboxes'); addpath(genpath(pn.tuSIM_tools));
pn.kwave = fullfile(rootpath, 'tools', 'k-wave-toolbox-version-1.4', 'k-Wave'); addpath(pn.kwave);
pn.simnibs = fullfile(rootpath, '..', '.conda', 'envs', 'simnibs_env', 'lib', 'python3.9', 'site-packages', 'simnibs');
pn.simnibs_bin = fullfile(rootpath, '..', 'SimNIBS', 'bin');
    addpath(fullfile(pn.simnibs, 'matlab_tools'))
pn.seg_path = fullfile(rootpath, 'data', 'ernie', 'simnibs');

addpath('/home/common/matlab/fieldtrip/qsub') % uncomment if you are using Donders HPC

cd(fullfile(rootpath, 'data')) % needs to contain a folder called "configs" in this setup

transducer_list = {'CTX250-001_64.5mm';...
    'CTX250-011_64.5mm';...
    'CTX500-026_73.5mm';...
    'CTX500-024_72.6mm'};
% the specification should actually not vary by depth

mni_targets = struct('right_mediodorsal_thalamus', [6.5 -18 9], ...
    'left_mediodorsal_thalamus', [-6.5 -18 9],...
    'right_amygdala',[26 -4 -20], ...
    'right_pulvinar', [15 -27 6.5]);

all_targets = fieldnames(mni_targets);

all_subjs = 1;

for i_transducer = 1:length(transducer_list)
    transducer_name = transducer_list{i_transducer};
    pn.sim_path = fullfile(rootpath, 'data', 'ernie', 'tussim', transducer_name);
        if ~isfolder(pn.sim_path); mkdir(pn.sim_path); end
    parameters = load_parameters(['config_thal_',transducer_name,'.yaml']);
    % Note: load_parameters always requires a default config that will then be overwritten
    if ismac
        % reset paths to local mac install; NOTE: horrible not to have this be dynamic
        parameters.simnibs_bin_path = '/Users/julian.kosciessa/SimNIBS/bin/';
    else
        parameters.simnibs_bin_path = fullfile('/home', 'neuromod', 'julkos', 'SimNIBS', 'bin');
    end
    parameters.ld_library_path ="/opt/gcc/7.2.0/lib64";
    parameters.paths_to_add = {pn.kwave};
    parameters.dist_close = 80; % set distance in mm that is close enough
    
    for subject_id = all_subjs

        fprintf('Current subject: %03i\n', subject_id)

        % Make subfolder (if enabled) and check if directory exists
        if isfield(parameters,'subject_subfolder') && parameters.subject_subfolder == 1
            parameters.output_dir = fullfile(pn.sim_path, sprintf('sub-%03d', subject_id));
        else 
            parameters.output_dir = pn.sim_path;
        end
        
        if ~isfolder(parameters.output_dir)
            mkdir(parameters.output_dir);
        end  
        for i_target = 1:length(all_targets)
            target_name = all_targets{i_target};
            timelimit = 60*60*0.5; % time limit for a job in seconds
            memorylimit = 8; % memory limit for a job in Gb
            parameters.overwrite_files = "always";
            tpos_output_file = fullfile(parameters.output_dir, sprintf('tpars_sub-%03i_%s.csv', subject_id, target_name));
            if ~isfile(tpos_output_file) || strcmp(parameters.overwrite_files, 'always')
                transducer_positioning_with_qsub(subject_id, parameters, pn, target_name, mni_targets, timelimit, memorylimit)
                %transducer_positioning(parameters, pn, subject_id, target_name, mni_targets)
            end
        end
    end
end