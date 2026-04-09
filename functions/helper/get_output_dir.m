function output_dir = get_output_dir(parameters)
%GET_OUTPUT_DIR  Derive the simulation output directory from parameters.
%
%   output_dir = get_output_dir(parameters)
%
%   Single source of truth for output_dir derivation, used by
%   path_log_setup, hpc_setup_temp_files, and uncertainty_pipeline.
%   Mirrors the logic in path_log_setup without any side effects.
%
%   Requires parameters.path.sim (non-empty) and parameters.subject_id.

    if ~isfield(parameters, 'path') || ~isfield(parameters.path, 'sim') || ...
            isempty(parameters.path.sim) || ...
            (ischar(parameters.path.sim)   && strcmp(parameters.path.sim, '')) || ...
            (isstring(parameters.path.sim) && parameters.path.sim == "")
        error('get_output_dir:missingPathSim', ...
            ['parameters.path.sim is empty or missing. ' ...
             'Set parameters.path.sim to the simulation output directory before calling the pipeline.']);
    end

    if parameters.path.subject_subfolder == 1
        output_dir = fullfile(parameters.path.sim, ...
            sprintf('sub-%03d', parameters.subject_id));
    else
        output_dir = parameters.path.sim;
    end
end
