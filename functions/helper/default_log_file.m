function filename_log = default_log_file(parameters, dir_logs)
% DEFAULT_LOG_FILE  Timestamped pipeline log (diary) path for a run
%
% Single definition of the default log name, shared by path_log_setup (local
% runs) and prestus_pipeline_start (HPC submission). Assigning io.log_file at
% submission lets the post-job report find the log without guessing the
% timestamp the job would otherwise pick.
%
% Use as:
%   parameters.io.log_file = default_log_file(parameters, parameters.io.dir_logs)
%
% Input:
%   parameters - PRESTUS config; uses subject_id, simulation.medium, io.output_affix
%   dir_logs   - (1,:) char, log directory (<output_dir>/log)
%
% Output:
%   filename_log - (1,:) char, <dir_logs>/sub-XXX_<medium>_<yyMMdd_HHmm><affix>.txt
%
% See also: PATH_LOG_SETUP, PRESTUS_PIPELINE_START

filename_log = fullfile(dir_logs, ...
    sprintf('sub-%03d_%s_%s%s.txt', ...
    parameters.subject_id, parameters.simulation.medium, ...
    char(string(datetime('now'), 'yyMMdd_HHmm')), parameters.io.output_affix));
end
