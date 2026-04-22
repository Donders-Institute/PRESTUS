function output_dir = get_output_dir(parameters)
% GET_OUTPUT_DIR  Derive the simulation output directory from parameters
%
%   Single source of truth for output_dir derivation, used by
%   path_log_setup, hpc_setup_temp_files, and uncertainty_pipeline.
%   Mirrors the logic in path_log_setup without any side effects.
%
% Use as:
%   output_dir = get_output_dir(parameters)
%
% Input:
%   parameters - PRESTUS config with path.sim and subject_id
%
% Output:
%   output_dir - resolved output directory path
%
% See also: PATH_LOG_SETUP, HPC_SETUP_TEMP_FILES

arguments
    parameters (1,1) struct
end

    if ~isfield(parameters, 'path') || ~isfield(parameters.path, 'sim') || ...
            isempty(parameters.path.sim) || ...
            (ischar(parameters.path.sim)   && strcmp(parameters.path.sim, '')) || ...
            (isstring(parameters.path.sim) && parameters.path.sim == "")
        error('get_output_dir:missingPathSim', ...
            ['parameters.path.sim is empty or missing. ' ...
             'Set parameters.path.sim to the simulation output directory before calling the pipeline.']);
    end

    output_dir = fullfile(parameters.path.sim, sprintf('sub-%03d', parameters.subject_id));
end
