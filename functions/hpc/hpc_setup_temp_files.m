function [log_dir, prestus_path, temp_data_path, temp_m_path, temp_m_file] = ...
    hpc_setup_temp_files(parameters)
% HPC_SETUP_TEMP_FILES  Create output/log directories and timestamped temp files
%
% Creates the simulation output directory and a hpc_log/ subdirectory,
% then generates uniquely named temporary files for passing parameters to the
% HPC worker (a .mat data file) and as the submitted MATLAB script (a .m file).
% File names are timestamped and randomised to avoid collisions between
% concurrent job submissions.
%
% Use as:
%   [log_dir, prestus_path, temp_data_path, temp_m_path, temp_m_file] = ...
%       hpc_setup_temp_files(parameters)
%
% Input:
%   parameters - PRESTUS config struct; must contain subject_id and io fields
%                used by get_output_dir
%
% Output:
%   log_dir        - path to <output_dir>/hpc_log/
%   prestus_path   - absolute path to the PRESTUS root (from get_prestus_path)
%   temp_data_path - full path to the temporary .mat file (written by caller)
%   temp_m_path    - full path to the temporary .m script (written by caller)
%   temp_m_file    - basename of the .m script without extension (used in
%                    scheduler script as the MATLAB entry point)
%
% See also: HPC_SUBMIT_JOB, GET_OUTPUT_DIR, GET_PRESTUS_PATH

% Setup output directory
output_dir = get_output_dir(parameters);
if ~isfolder(output_dir), mkdir(output_dir); end

% Setup HPC log directory
log_dir = fullfile(output_dir, 'hpc_log');
if ~isfolder(log_dir), mkdir(log_dir); end

prestus_path = get_prestus_path;

% Generate temp files
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
temp_base = tempname(log_dir);
[~, temp_base_name] = fileparts(temp_base);
temp_base_name = temp_base_name(end-7:end);

temp_data_path = fullfile(log_dir, sprintf('temp_data_%s_%s.mat', timestamp, temp_base_name));
temp_m_file = sprintf('temp_matlab_%s_%s', timestamp, temp_base_name);
temp_m_path = fullfile(log_dir, [temp_m_file, '.m']);

end
