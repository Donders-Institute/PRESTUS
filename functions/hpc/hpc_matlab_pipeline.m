function hpc_matlab_pipeline(temp_m_path, temp_data_path, path_to_pipeline, options)
%% HPC_MATLAB_PIPELINE  Generate temporary MATLAB batch script
%
%   Creates self-deleting MATLAB batch script that loads parameters and runs
%   prestus_pipeline. Handles optional sequential_configs.
%
%   Inputs:
%     temp_m_path     - Path for output .m script
%     temp_data_path  - Path to .mat parameter file
%     path_to_pipeline - Directory containing prestus_pipeline.m
%     options         - Optional struct with sequential_configs field
%
%   See also HPC_SETUP_TEMP_FILES, HPC_SUBMIT_JOB.

fid = fopen(temp_m_path, 'w+');
if ismember(fieldnames(options), 'sequential_configs')
    sequential_configs = options.sequential_configs;
    save(temp_data_path, 'sequential_configs', '-append');
    fprintf(fid, "load '%s'; cd '%s'; prestus_pipeline(subject_id, parameters, options); delete '%s'; delete '%s';", ...
        temp_data_path, path_to_pipeline, temp_data_path, temp_m_path);
else
    fprintf(fid, "load '%s'; cd '%s'; prestus_pipeline(subject_id, parameters); delete '%s'; delete '%s';", ...
        temp_data_path, path_to_pipeline, temp_data_path, temp_m_path);
end
fclose(fid);
end