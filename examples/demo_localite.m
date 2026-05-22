% DEMO_LOCALITE - Demo extraction of Localite coordinates

% The goal is to obtain transducer and target positions as they have been
% captured by localite trackers. They will be provided in the format
% expected by PRESTUS (native image space), and can optionally be output in
% MNI space (experimental).

% Dependency:
% The demo will require a matching input of localite files and the subject T1w scan. 
% For optional MNI analysis, a matching simnibs segmentation and simbins binary path (parameters.simnibs_bin_path) must be specified.

% Clear environment to ensure a clean state before running the demo
restoredefaultpath;  % Restore MATLAB's default search path to avoid conflicts
clear all;           % Clear all variables from workspace
close all;           % Close all figure windows
clc;                 % Clear the command window

% --- Hard-coded demo inputs ---
demo_sub_id = 'sub-003';                         % Subject identifier for demo processing
demo_session = 'ses-01';                         % Session identifier to process (localite files are expected to follow a sub-xxx/ses-xx organization)
demo_setup = 'CTX500-026-010_79.6mm';            % Transducer/setup configuration name (suffix of 'config_XXX')
demo_positions = {'lVS', 'rVS'};                 % Labels for multiple positions (left and right hemisphere)
demo_markertype = 'GUMMarkers';                  % TriggerMarkers/GUMMarkers

rootpath = '/project/2425122.01/kplan_example/';  % Root directory of data and tools, adjust to your environment

% Indicate whether to compute and include MNI-space coordinates (optional)
mni_coords_requested = 0;

% Define key path structure for data and external
pn.tuSIM = fullfile('/project/2425122.01/v05/', 'tools', 'PRESTUS');           % PRESTUS toolbox path
addpath(pn.tuSIM);
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions');              % PRESTUS functions path
addpath(genpath(pn.tuSIM_fun));
pn.tuSIM_tools = fullfile(pn.tuSIM, 'external');            % PRESTUS external path
addpath(genpath(pn.tuSIM_tools));

% [Optional, MNI transform] Set user/environment-specific paths and parameters
% parameters.io.simnibs_bin_path = fullfile('/home', 'neuromod', 'julkos', '.conda', 'envs', 'simnibs_env', 'bin');
% parameters.io.data_path = fullfile(pn.data_seg, sprintf('m2m_%s', demo_sub_id));
% parameters.io.seg_path = pn.data_seg;
% parameters.hpc.ld_library_path = "/opt/gcc/7.2.0/lib64";

% The expected folder structure for key data and intermediate files is:
% (1a) [TriggerMarkers] Localite trigger position files in: 
%     <pn.data_postlocalite>/sub-XXX/ses-XX/localite/Session_XXXX/TMSTrigger
% (1b) [GUMMarkers] Localite General positioning markers in:
%     <pn.data_postlocalite>/sub-XXX/ses-XX/localite/
% (2) A planning image aligned with localite coordinate expectations:
%     <pn.data_prelocalite>/sub-XXX_T1(_forneuronav).nii(.gz)
% (3) SIMNIBS segmentation folder for the subject:
%     <pn.data_seg>/m2m_sub-XXX/
% (4) PRESTUS config files, for example in:
%     <pn.configs>/

pn.data_path = fullfile(rootpath, 'data');                      % Base data path
pn.data_postlocalite = fullfile(pn.data_path, 'localite');      % Folder [1] where localite XML are stored
pn.data_prelocalite = fullfile(pn.data_path, 'localite');       % Folder [2] containing T1 planning images
pn.data_seg = fullfile(pn.data_path, 'simnibs');                % Folder [3] Segmentation files for SIMNIBS
pn.configs = fullfile(pn.data_path, 'config');                 % Folder [4] containing config YAMLs

% Load stimulation parameters from the transducer YAML config
cd(pn.configs);
parameters = load_parameters(['config_', demo_setup, '.yaml']); 

% Select the latest valid Localite XML file for the subject and session
localite = neuronav_select_localite(pn, demo_sub_id, demo_session, demo_markertype);

% Calculate statistics over stimulus series (e.g., average position, variability)
% Parameters: voxel size of planning image  (1 mm), expected localite trigger train length (2)
[results] = neuronav_compute_series_statistics(localite, 1, 2, demo_markertype);

% Fix demo position indices as 1 (left) and 2 (right), matching demo_positions
positions = 1:numel(demo_positions);

% Convert averaged triggers to voxel indices and RAS coordinates in native space
[trans_ras, trans_pos, target_ras, target_pos, ~] = ...
    neuronav_convert_trigger_to_voxels(demo_sub_id, positions, results, parameters, pn);

clear results;

% Conditional logic to convert native space coordinates to MNI space
if mni_coords_requested == 1
    [trans_ras_seg, target_ras_seg, trans_mni_pos, target_mni_pos, trans_mni_ras, targ_mni_ras] = ...
        neuronav_convert_native_to_MNI(demo_sub_id, parameters, pn, trans_ras, target_ras, demo_positions);
else
    % If MNI coordinates not requested, fill outputs with NaNs to maintain output shape consistency
    trans_ras_seg(positions, 1:3) = NaN;
    target_ras_seg(positions, 1:3) = NaN;
    trans_mni_pos(positions, 1:3) = NaN;
    target_mni_pos(positions, 1:3) = NaN;
    trans_mni_ras(positions, 1:3) = NaN;
    targ_mni_ras(positions, 1:3) = NaN;
end

% Flag indicating if interpolation was used (0 = no)
interpolated = 0;

% Export all computed data and statistics to a CSV file within the output folder
neuronav_export_session_csv(...
    demo_sub_id, ...
    demo_session, ...
    demo_positions, ...
    pn.data_postlocalite, ...
    trans_pos, ...
    target_pos, ...
    trans_ras_seg, ...
    target_ras_seg, ...
    trans_mni_pos, ...
    target_mni_pos, ...
    trans_mni_ras, ...
    targ_mni_ras, ...
    interpolated);

% Inform user of processing completion and output location
fprintf('Demo processing completed. CSV saved to:\n%s\n', pn.data_postlocalite);

% If interpolation across a group of subjects would be desirable, the
% following additional steps could be used:
% [INT1] Get GROUP MNI coordinates
% [trans_mni_ras, targ_mni_ras, trans_mni_pos, target_mni_pos] = neuronav_get_group_mean_mni(session_target, parameters);
% [INT2] TRANSFORM MNI averages to subject space
% [target_ras_seg, trans_ras_seg, target_pos, trans_pos] = ...
%             neuronav_convert_MNI_to_native(sub_id, parameters, pn, trans_mni_ras, targ_mni_ras);
% [INT3] Write output CSV
% interpolated = 1;
% neuronav_export_session_csv(...)