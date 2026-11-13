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
demo_session = 'ses-01';                         % Session identifier (inputs follow sub-XXX/ses-XX, see below)
demo_setup = 'CTX500-026-010_79.6mm';            % Transducer/setup configuration name (suffix of 'config_XXX')
demo_positions = {'lVS', 'rVS'};                 % Target names of the recorded series 1 and 2

rootpath = '/path/to/localite_example/';  % Root directory of data and tools, adjust to your environment

% Indicate whether to compute and include MNI-space coordinates (optional)
mni_coords_requested = 0;

% Define key path structure for data and external
pn.tuSIM = fullfile('/path/to', 'tools', 'PRESTUS');           % PRESTUS toolbox path
addpath(pn.tuSIM);
pn.tuSIM_fun = fullfile(pn.tuSIM, 'functions');              % PRESTUS functions path
addpath(genpath(pn.tuSIM_fun));
pn.tuSIM_tools = fullfile(pn.tuSIM, 'external');            % PRESTUS external path
addpath(genpath(pn.tuSIM_tools));

% The expected folder structure (see doc_placement_neuronav.md, "Localite I/O"):
% (1) Raw Localite exports, the Localite folder structure below sub/ses:
%     <path.localite_raw>/sub-XXX/ses-XX/<loc_folder>/Session_XXXX/TMSTrigger/TriggerMarkers_Coil0_<ts>.xml
%     <path.localite_raw>/sub-XXX/ses-XX/<loc_folder>/GUMMarkers<ts>.xml
% (2) The planning image loaded into Localite:
%     <path.localite_t1>/[sub-XXX/]sub-XXX_*T1*.nii(.gz)   (path.localite_t1 defaults to path.localite_raw)
% (3) SimNIBS segmentation folder for the subject (MNI step only):
%     <path.seg>/m2m_sub-XXX/
% (4) PRESTUS config files, for example in:
%     <pn.configs>/
% Ingestion writes <path.localite_post>/sub-XXX/ses-XX/InstrumentMarker_sub-XXX_ses-XX.xml/.json

pn.data_path = fullfile(rootpath, 'data');                      % Base data path
pn.configs = fullfile(pn.data_path, 'config');                 % Folder [4] containing config YAMLs

% Load stimulation parameters from the transducer YAML config
cd(pn.configs);
parameters = load_parameters(['config_', demo_setup, '.yaml']);
parameters.path.localite_raw  = fullfile(pn.data_path, 'localite');            % [1]
parameters.path.localite_t1   = fullfile(pn.data_path, 'localite');            % [2]
parameters.path.localite_post = fullfile(pn.data_path, 'localite_ingested');   % outputs
parameters.path.seg           = fullfile(pn.data_path, 'simnibs');             % [3]
parameters.placement.localite.loc_folder = 'localite';

% Stage 1: select the newest valid export, average each series and write the
% named targets (series 1 -> lVS, series 2 -> rVS). neuronav_gui does the same
% interactively, with QC plots and suggested names.
positions = 1:numel(demo_positions);
target_map = struct('series_index', num2cell(positions), 'transducer_id', 0, ...
    'target_name', demo_positions);
[~, ~, results] = neuronav_ingest_markers(parameters, demo_sub_id, demo_session, target_map);

% Convert averaged series to voxel indices and RAS coordinates in native space
[trans_ras, trans_pos, target_ras, target_pos, ~] = ...
    neuronav_convert_trigger_to_voxels(demo_sub_id, positions, results, parameters);

clear results;

% Conditional logic to convert native space coordinates to MNI space
if mni_coords_requested == 1
    [trans_ras_seg, target_ras_seg, trans_mni_pos, target_mni_pos, trans_mni_ras, targ_mni_ras] = ...
        neuronav_convert_native_to_MNI(demo_sub_id, parameters, trans_ras, target_ras, demo_positions);
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
    parameters.path.localite_post, ...
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
fprintf('Demo processing completed. CSV saved to:\n%s\n', parameters.path.localite_post);

% If interpolation across a group of subjects would be desirable, the
% following additional steps could be used:
% [INT1] Get GROUP MNI coordinates
% [trans_mni_ras, targ_mni_ras, trans_mni_pos, target_mni_pos] = neuronav_get_group_mean_mni(session_target, parameters);
% [INT2] TRANSFORM MNI averages to subject space
% [target_ras_seg, trans_ras_seg, target_pos, trans_pos] = ...
%             neuronav_convert_MNI_to_native(sub_id, parameters, trans_mni_ras, targ_mni_ras);
% [INT3] Write output CSV
% interpolated = 1;
% neuronav_export_session_csv(...)