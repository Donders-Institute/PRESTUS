function gui_run_worker(yaml_path, prestus_path)
%GUI_RUN_WORKER  Helper to run the pipeline from a given YAML config.
%   Can be called directly or via parfeval with a process-based pool.
%   Note: backgroundPool (thread workers) do not support addpath or diary,
%   so this function should not be used with backgroundPool.

if nargin > 1 && ~isempty(prestus_path)
    addpath(genpath(prestus_path));
end

parameters = load_parameters([], yaml_path);
prestus_pipeline_start(parameters);
end
