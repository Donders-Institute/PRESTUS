function gui_run_worker(yaml_path, prestus_path)
% GUI_RUN_WORKER  Run the PRESTUS pipeline from a YAML config file
%
% Loads parameters from the given YAML and calls prestus_pipeline_start.
% Intended to be submitted via parfeval with a process-based parallel pool
% so the GUI remains responsive during a simulation run.
%
% Use as:
%   gui_run_worker(yaml_path)
%   gui_run_worker(yaml_path, prestus_path)
%
% Input:
%   yaml_path    - path to the YAML config file to run
%   prestus_path - (optional) PRESTUS root path to add to the MATLAB search
%                  path in the worker process
%
% Note:
%   backgroundPool (thread-based workers) does not support addpath or diary.
%   Always use a process-based pool (parpool('Processes')) with this function.
%
% See also: PRESTUS_GUI, PRESTUS_PIPELINE_START, LOAD_PARAMETERS

if nargin > 1 && ~isempty(prestus_path)
    addpath(genpath(prestus_path));
end

parameters = load_parameters([], yaml_path);
prestus_pipeline_start(parameters);
end
