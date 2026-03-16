function job_id = hpc_submit_job(hpc_type, temp_m_file, parameters, subject_id, log_dir)
%% HPC_SUBMIT_JOB  Submit HPC batch job (SLURM or qsub)
%
%   Generates scheduler script and submits job. Supports SLURM (sbatch) and
%   PBS (qsub) with automatic GPU/partition detection.
%
%   Inputs:
%     hpc_type     - 'slurm' or 'qsub'
%     temp_m_file  - MATLAB script basename
%     parameters   - Job parameters (hpc_partition, hpc_gpu, etc.)
%     subject_id   - Subject number for naming
%     parameters.hpc_timelimit    - Walltime in format 'HH:MM:SS' or minutes
%     parameters.hpc_memorylimit  - Memory in GB
%     log_dir      - Log directory path
%
%   Output:
%     job_id     - Job ID string/number
%
%   See also HPC_DETECT_SYSTEM, HPC_WAIT_FOR_COMPLETION.

subj_id_string = sprintf('sub-%03d', subject_id);

switch hpc_type
    case 'slurm'
        temp_slurm_path = fullfile(log_dir, sprintf('temp_slurm_%s.sh', datestr(now, 'yyyymmdd_HHMMSS')));
        write_slurm_script(temp_slurm_path, parameters, subject_id, temp_m_file, log_dir);
        job_id = submit_slurm_job(temp_slurm_path, log_dir);
        
    case 'qsub'
        temp_qsub_path = fullfile(log_dir, sprintf('temp_qsub_%s.sh', datestr(now, 'yyyymmdd_HHMMSS')));
        write_qsub_script(temp_qsub_path, parameters, subject_id, temp_m_file, log_dir);
        job_id = submit_qsub_job(temp_qsub_path, log_dir);
        
    otherwise
        error('Unsupported HPC type: %s', hpc_type);
end

fprintf('Job "%s" (ID: %s) submitted successfully\n', hpc_job_name(parameters, subject_id), sprintf('%d', job_id));

% ========== LOCAL FUNCTIONS ==========
function write_slurm_script(temp_slurm_path, parameters, subject_id, temp_m_file, log_dir)
    subj_id_string = sprintf('sub-%03d', subject_id);
    job_name = hpc_job_name(parameters, subject_id);
    
    fid = fopen(temp_slurm_path, 'w+');
    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '#SBATCH --job-name=%s\n', job_name);
    
    % Partition & GPU detection
    needs_gpu = (isfield(parameters, 'code_type') && any(strcmp(parameters.code_type, {'matlab_gpu', 'cpp_gpu'})));
    
    if isfield(parameters, 'hpc_partition') && ~isempty(strtrim(char(parameters.hpc_partition)))
        fprintf(fid, '#SBATCH --partition=%s\n', strtrim(char(parameters.hpc_partition)));
    elseif needs_gpu
        fprintf(fid, '#SBATCH --partition=gpu\n');
    end
    
    if isfield(parameters, 'hpc_gpu') && ~isempty(strtrim(char(parameters.hpc_gpu)))
        fprintf(fid, '#SBATCH --gres=%s\n', strtrim(char(parameters.hpc_gpu)));
    elseif needs_gpu
        fprintf(fid, '#SBATCH --gres=gpu:1\n');
    end

    if isfield(parameters, 'hpc_reservation') && ~isempty(strtrim(char(parameters.hpc_reservation)))
        fprintf(fid, '#SBATCH --reservation=%s\n', strtrim(char(parameters.hpc_reservation)));
    end
    
    fprintf(fid, '#SBATCH --mem=%iG\n', parameters.hpc_memorylimit);
    fprintf(fid, '#SBATCH --time=%s\n', parameters.hpc_timelimit);
    fprintf(fid, '#SBATCH --output=%s_slurm_output_%%j.log\n', subj_id_string);
    fprintf(fid, '#SBATCH --error=%s_slurm_error_%%j.log\n', subj_id_string);
    fprintf(fid, '#SBATCH --chdir=%s\n', log_dir);
    
    if needs_gpu, fprintf(fid, 'nvidia-smi\n'); end
    fprintf(fid, 'module load matlab/R2023b\n');
    fprintf(fid, 'matlab -batch "%s"\n', temp_m_file);
    fclose(fid);
end

function job_id = submit_slurm_job(temp_slurm_path, log_dir)
    sbatch_call = sprintf('sbatch %s', temp_slurm_path);
    full_cmd = sprintf('cd %s; %s', log_dir, sbatch_call);
    
    fprintf('SLURM command: %s\n', full_cmd);
    [status, out] = system(full_cmd);
    if status ~= 0, error('SLURM submission failed: %s', out); end
    
    job_ids = regexp(out, '\d+', 'match');
    if isempty(job_ids), error('No SLURM job ID returned'); end
    job_id = str2double(job_ids{1});
end

function write_qsub_script(temp_qsub_path, parameters, subject_id, temp_m_file, log_dir)
    subj_id_string = sprintf('sub-%03d', subject_id);
    job_name = hpc_job_name(parameters, subject_id);
    
    fid = fopen(temp_qsub_path, 'w+');
    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '#PBS -N %s\n', job_name);
    fprintf(fid, '#PBS -l nodes=1:gpus=1,feature=cuda,reqattr=cudacap>=8.0,mem=%iGb,walltime=%i\n', ...
        parameters.hpc_memorylimit, parameters.hpc_timelimit);
    fprintf(fid, '#PBS -o %s_qsub_output_%%j.log\n', subj_id_string);
    fprintf(fid, '#PBS -e %s_qsub_error_%%j.log\n', subj_id_string);
    fprintf(fid, '#PBS -d %s\n', log_dir);
    fprintf(fid, 'module load matlab/R2023b\n');
    fprintf(fid, 'matlab -batch "%s"\n', temp_m_file);
    fclose(fid);
end

function job_id = submit_qsub_job(temp_qsub_path, log_dir)
    qsub_call = sprintf('qsub %s', temp_qsub_path);
    full_cmd = sprintf('cd %s; %s', log_dir, qsub_call);
    
    fprintf('qsub command: %s\n', full_cmd);
    [status, out] = system(full_cmd);
    if status ~= 0, error('qsub submission failed: %s', out); end
    job_id = strtrim(out);
end
end
