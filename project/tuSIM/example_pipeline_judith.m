%% Initialization

% add paths
addpath('.');
addpath('functions')
addpath('toolboxes/kwave') % set your kwave path here
addpath('toolboxes/Colormaps') % set your path to Colormaps files here
addpath('toolboxes/export_fig') % set your path to export_fig files here
addpath('toolboxes/yaml') % set your path to yaml files here

parameters = load_parameters('judith_config.yaml');

parameters.simulation_medium = 'water_and_skull';
parameters.overwrite_files = 'always';
parameters.overwrite_simnibs = 0;

parameters.run_posthoc_water_sims = 1;

subject_list = [1]; % each subject should have T1 & T2 files named sub-%03d_orig_T1w.nii.gz and sub-%03d_orig_T2w.nii.gz in the data folder

%% Run the pipeline for each subject
for subject_id = subject_list
    single_subject_pipeline(subject_id, parameters)
end