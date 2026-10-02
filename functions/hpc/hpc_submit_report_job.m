function report_job_id = hpc_submit_report_job(main_job_id, platform, parameters, log_dir, prestus_path)
% HPC_SUBMIT_REPORT_JOB  Submit a lightweight dependent job to generate the post-run report
%
% Submits a minimal SLURM job (1 CPU, 2 GB, 10 min) that depends on the
% completion of main_job_id and calls hpc_job_report to write the resource
% usage, carbon, and cost report to log_dir. For qsub clusters this is a
% no-op (sacct is SLURM-only) and returns [].
%
% Use as:
%   report_job_id = hpc_submit_report_job(main_job_id, platform, parameters, log_dir, prestus_path)
%
% Input:
%   main_job_id  - numeric SLURM job ID of the main simulation job
%   platform     - 'slurm' (only SLURM is supported; qsub returns [])
%   parameters   - PRESTUS config struct (saved to file for the report job)
%   log_dir      - directory for the report file and job logs
%   prestus_path - root path of the PRESTUS installation
%
% Output:
%   report_job_id - SLURM job ID of the submitted report job, or [] if skipped
%
% See also: HPC_JOB_REPORT, HPC_SUBMIT_JOB

report_job_id = [];

if ~strcmpi(platform, 'slurm')
    return;
end

main_job_id_str = sprintf('%.0f', main_job_id);

% Save a minimal copy of parameters for the report job
params_file = fullfile(log_dir, sprintf('report_params_job%s.mat', main_job_id_str));
save(params_file, 'parameters');

% Write the MATLAB report script
m_file = fullfile(log_dir, sprintf('report_job%s.m', main_job_id_str));
fid = fopen(m_file, 'w+');
fprintf(fid, 'load(''%s'', ''parameters'');\n', matlab_escape_quotes(params_file));
fprintf(fid, 'addpath(''%s'');\n', matlab_escape_quotes(fullfile(prestus_path, 'functions', 'helper')));
fprintf(fid, 'safe_addpath(''%s'');\n', matlab_escape_quotes(fullfile(prestus_path, 'functions')));
fprintf(fid, 'hpc_job_report(%s, ''slurm'', parameters, ''%s'');\n', main_job_id_str, log_dir);
fprintf(fid, 'delete(''%s'');\n', matlab_escape_quotes(params_file));
fprintf(fid, 'delete(''%s'');\n', matlab_escape_quotes(m_file));
fclose(fid);

% Determine which MATLAB module to load (match main job)
if isfield(parameters, 'hpc') && isfield(parameters.hpc, 'name') && ...
        strcmpi(parameters.hpc.name, 'snellius')
    matlab_module = 'MATLAB/2024b';
    module_year   = 'module load 2024';
else
    matlab_module = 'matlab/R2023b';
    module_year   = '';
end

% Write the SLURM batch script
sh_file = fullfile(log_dir, sprintf('report_slurm_job%s.sh', main_job_id_str));
subj_str = sprintf('sub-%03d', parameters.subject_id);
fid = fopen(sh_file, 'w+');
fprintf(fid, '#!/bin/bash\n');
fprintf(fid, '#SBATCH --job-name=PS_report_%s\n', main_job_id_str);
fprintf(fid, '#SBATCH --dependency=afterany:%s\n', main_job_id_str);
fprintf(fid, '#SBATCH --kill-on-invalid-dep=yes\n');
fprintf(fid, '#SBATCH --cpus-per-task=1\n');
fprintf(fid, '#SBATCH --mem=2G\n');
fprintf(fid, '#SBATCH --time=00:10:00\n');
fprintf(fid, '#SBATCH --output=%s_report_output_%%j.log\n', subj_str);
fprintf(fid, '#SBATCH --error=%s_report_error_%%j.log\n',  subj_str);
fprintf(fid, '#SBATCH --chdir=%s\n', log_dir);
if ~isempty(module_year)
    fprintf(fid, '%s\n', module_year);
end
fprintf(fid, 'module load %s\n', matlab_module);
[~, m_name] = fileparts(m_file);
fprintf(fid, 'matlab -batch "%s"\n', m_name);
fclose(fid);

% Submit. run_command invokes sbatch directly via argv, no shell at all.
[status, out] = run_command({'sbatch', '--parsable', sh_file}, 'WorkingDir', log_dir);
if status ~= 0
    warning('PRESTUS:hpc_submit_report_job:submitFailed', ...
        'Report job submission failed: %s', out);
    return;
end
report_job_id = parse_sbatch_job_id(out);
if isnan(report_job_id)
    report_job_id = [];
else
    fprintf('  Report job %d submitted (depends on job %s)\n', report_job_id, main_job_id_str);
end
end
