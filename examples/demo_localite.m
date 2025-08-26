% DEMO_THALSTIM_LOCALITE - Demo extraction of Localite coordinates

% The goal is to obtain transducer and target positions as they have been
% captured by localite trackers. They will be provided in the format
% expected by PRESTUS (native image space), and can optionally be output in
% MNI space (experimental).

% Dependency:
% The demo will require a matching input of localite files and the subject T1w scan. 
% For optional MNI analysis, a matching simnibs segmentation and simbins binary path (parameters.simnibs_bin_path) must be specified.

% Clear environment to ensure a clean state before running the demo
restoredefaultpath;   % Restore MATLAB's default search path to avoid conflicts
clear all;           % Clear all variables from workspace
close all;           % Close all figure windows
clc;                 % Clear the command window

% --- Hard-coded demo inputs ---
demo_sub_id = 'sub-001';                         % Subject identifier for demo processing
demo_session = 'ses-02';                         % Session identifier to process (localite files are expected to follow a sub-xxx/ses-xx organization)
demo_setup = 'CTX500-026-010_79.6mm';            % Transducer/setup configuration name (suffix of 'config_XXX')
demo_sim_dir = 'CTX500-026-010_79.6mm_60W_post'; % Simulation output folder name based on setup and intensity
demo_positions = {'left X', 'right X'};          % Labels for multiple positions (left and right hemisphere)

rootpath = '';  % Root directory of data and tools, adjust to your environment

% Indicate whether to compute and include MNI-space coordinates (optional)
mni_coords_requested = 0;

% Define key path structure for data and toolboxes
pn.tuSIM = fullfile(rootpath, 'tools', 'PRESTUS');           % PRESTUS toolbox path
addpath(pn.tuSIM);
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions');              % PRESTUS functions path
addpath(genpath(pn.tuSIM_fun));
pn.tuSIM_tools = fullfile(pn.tuSIM, 'toolboxes');            % PRESTUS toolboxes path
addpath(genpath(pn.tuSIM_tools));

% The expected folder structure for key data and intermediate files is:
% (1) Localite trigger position files in: 
%     <pn.data_postlocalite>/sub-XXX/ses-XX/localite/Session_XXXX/TMSTrigger
% (2) A planning image aligned with localite coordinate expectations:
%     <pn.data_prelocalite>/sub-XXX_T1_forneuronav.nii.gz
% (3) SIMNIBS segmentation folder for the subject:
%     <pn.data_seg>/m2m_sub-XXX/
% (4) PRESTUS config files, for example in:
%     <pn.configs>/

pn.data_path = fullfile(rootpath, 'data');                      % Base data path
pn.data_postlocalite = fullfile(pn.data_path, 'demo_neuronav'); % Folder [1] where localite XML are stored
pn.data_prelocalite = fullfile(pn.data_path, 'localite');       % Folder [2] containing T1 planning images
pn.data_seg = fullfile(pn.data_path, 'simnibs');                % Folder [3] Segmentation files for SIMNIBS
pn.configs = fullfile(pn.data_path, 'configs');                 % Folder [4] containing config YAMLs

% Move to data root for relative file handling
cd(fullfile(pn.configs, '..'));

% Load stimulation parameters from the transducer YAML config
parameters = load_parameters(['config_', demo_setup, '.yaml']); 

% Set user/environment-specific paths and parameters
parameters.simnibs_bin_path = fullfile('/home', 'neuromod', 'julkos', '.conda', 'envs', 'simnibs_env', 'bin'); 
parameters.data_path = fullfile(pn.data_seg, sprintf('m2m_%s', demo_sub_id));
parameters.seg_path = pn.data_seg;
parameters.ld_library_path = "/opt/gcc/7.2.0/lib64";
parameters.sim_path = fullfile(pn.data_path, 'tussim', demo_sim_dir);

% Select the latest valid Localite trigger XML file for the subject and session
localite = neuronav_select_and_average_localite(demo_sub_id, demo_session, pn);

% Calculate statistics over stimulus series (e.g., average position, variability)
% Parameters: voxel size of planning image  (0.9 mm), expected localite trigger train length (80)
[results, n_series] = neuronav_compute_series_statistics(localite, parameters, 0.9, 80);

% Generate averaged trigger marker structure for subsequent conversions and visualization
outputStruct = neuronav_create_marker_average(localite, results);

% Fix demo position indices as 1 (left) and 2 (right), matching demo_positions
positions = 1:numel(demo_positions);

% Convert averaged triggers to voxel indices and RAS coordinates in native space
[trans_ras, trans_pos, target_ras, target_pos, ~] = ...
    neuronav_convert_trigger_to_voxels(demo_sub_id, positions, outputStruct, parameters, pn);

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

% Set output folder for CSV as the localite data folder for the subject/session
csv_outfile = pn.data_postlocalite;

% Flag indicating if interpolation was used (0 = no)
interpolated = 0;

% Export all computed data and statistics to a CSV file within the output folder
neuronav_export_session_csv(...
    demo_sub_id, ...
    demo_session, ...
    demo_positions, ...
    n_series, ...
    csv_outfile, ...
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
fprintf('Demo processing completed. CSV saved to:\n%s\n', csv_outfile);

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
