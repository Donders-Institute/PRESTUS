function tf = telemetry_is_headless_cluster()
% TELEMETRY_IS_HEADLESS_CLUSTER  True for compiled apps and SLURM/PBS jobs.
%
% These sessions can never show a prompt of any kind, console or dialog:
% a compiled/deployed application has no attached terminal, and a SLURM or
% PBS/Torque allocation typically runs on a compute node with no display.
%
% One of the conditions of TELEMETRY_IS_NONINTERACTIVE, which also covers
% `matlab -batch` sessions and sessions with neither desktop nor display.
%
% Use as:
%   tf = telemetry_is_headless_cluster()
%
% Output:
%   tf - logical
%
% See also: TELEMETRY_IS_NONINTERACTIVE, TELEMETRY_SETUP

    tf = true;
    try
        tf = isdeployed || ~isempty(getenv('SLURM_JOB_ID')) || ~isempty(getenv('PBS_JOBID'));
    catch
        tf = true;
    end
end
