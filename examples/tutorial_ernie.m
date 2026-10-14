%% TUTORIAL_ERNIE  Head simulation on SimNIBS' Ernie example (see doc/doc_tutorial_head.md)
%
% Run the three sections one after the other. Edit the paths in
% tutorial_ernie.yaml first. Each section can also be submitted to a cluster
% by setting platform: slurm (or qsub) in the YAML.

prestus_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(prestus_root, 'functions')));
addpath(genpath(fullfile(prestus_root, 'config')));
addpath(genpath(fullfile(prestus_root, 'external')));

cfg_dir = fullfile(prestus_root, 'examples');

%% 1. Segmentation only (SimNIBS charm)
parameters = load_parameters('tutorial_ernie.yaml', cfg_dir, ...
    {'subject_id', 1; 'modules.segmentation_only', 1});
prestus_pipeline_start(parameters);

%% 2. Placement and acoustic simulation
parameters = load_parameters('tutorial_ernie.yaml', cfg_dir, {'subject_id', 1});
prestus_pipeline_start(parameters);

%% 3. Add the thermal simulation
% Repeats the acoustic simulation (io.overwrite_files: always, the default)
% and then simulates heating over the timing protocol in the YAML.
parameters = load_parameters('tutorial_ernie.yaml', cfg_dir, ...
    {'subject_id', 1; 'modules.run_heating_sims', 1});
prestus_pipeline_start(parameters);
