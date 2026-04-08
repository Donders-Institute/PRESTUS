function [log_dir, prestus_path, temp_data_path, temp_m_path, temp_m_file] = ...
    hpc_setup_temp_files(parameters)
%% HPC_SETUP_TEMP_FILES  Setup directories and generate temporary files
%
%   Creates output directory, log directory, and timestamped temporary files
%   for MATLAB data and script.
%   Subject ID is read from parameters.subject_id.
%
%   Outputs:
%     log_dir         - Path to batch_job_logs directory
%     prestus_path    - PRESTUS path
%     temp_data_path  - Path for temporary .mat data file
%     temp_m_path     - Path for temporary .m script file
%     temp_m_file     - Basename of MATLAB script (no path)
%
%   See also HPC_SUBMIT_JOB.

% Setup output directory
if isfield(parameters.path, 'subject_subfolder') && parameters.path.subject_subfolder
    output_dir = fullfile(parameters.path.sim, sprintf('sub-%03d', parameters.subject_id));
else
    output_dir = parameters.path.sim;
end
if ~isfolder(output_dir), mkdir(output_dir); end

% Setup batch/log directory
log_dir = fullfile(output_dir, 'batch_job_logs');
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
