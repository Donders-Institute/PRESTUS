function hpc_matlab_pipeline(temp_m_path, temp_data_path, prestus_path, options)
%% HPC_MATLAB_PIPELINE  Generate temporary MATLAB batch script
%
%   Creates self-deleting MATLAB batch script that loads parameters and runs
%   prestus_pipeline. Handles optional sequential_configs.
%
%   Inputs:
%     temp_m_path     - Path for output .m script
%     temp_data_path  - Path to .mat parameter file
%     prestus_path    - Directory containing prestus_pipeline.m
%     options         - Optional struct with sequential_configs field
%
%   See also HPC_SETUP_TEMP_FILES, HPC_SUBMIT_JOB.

fid = fopen(temp_m_path, 'w+');
fprintf(fid, 'load(''%s'');\n', temp_data_path);
fprintf(fid, 'addpath(genpath(''%s''));\n', prestus_path);
if ismember(fieldnames(options), 'sequential_configs')
    sequential_configs = options.sequential_configs;
    save(temp_data_path, 'sequential_configs', '-append');
    fprintf(fid, 'prestus_pipeline(subject_id, parameters, options);\n');
else
    fprintf(fid, 'prestus_pipeline(subject_id, parameters);\n');
end
fprintf(fid, 'delete(''%s'');\n', temp_data_path);
fprintf(fid, 'delete(''%s'');\n', temp_m_path);
fclose(fid);
end