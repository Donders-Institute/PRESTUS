function single_subject_pipeline_with_slurm(subject_id, parameters, wait_for_job, timelimit, memorylimit, options)
% SINGLE_SUBJECT_PIPELINE_WITH_SLURM  SLURM entry point of PRESTUS v0.5.0 and earlier; submits via PRESTUS_PIPELINE_START
%
% Kept so that scripts written for v0.5.0 keep running. Warns
% (prestus:legacyFunction) with the current call. Old parameter names are
% converted (LEGACY_BRIDGE_STRUCT); the arguments set hpc.wait_for_job,
% hpc.timelimit and hpc.memorylimit, with the v0.5.0 defaults (no wait,
% 04:00:00, 40 GB) when omitted, as v0.5.0 did.
%
% Use as:
%   single_subject_pipeline_with_slurm(subject_id, parameters)
%   single_subject_pipeline_with_slurm(subject_id, parameters, wait_for_job, timelimit, memorylimit)
%   single_subject_pipeline_with_slurm(..., sequential_configs=s)
%
% See also: PRESTUS_PIPELINE_START, SINGLE_SUBJECT_PIPELINE, LEGACY_BRIDGE_ENTRY

    arguments
        subject_id (1,1) double
        parameters struct
        wait_for_job logical = false
        timelimit string = "04:00:00"
        memorylimit (1,1) double = 40
        options.sequential_configs struct = struct()
    end
    legacy_hpc_submit('slurm', 'single_subject_pipeline_with_slurm', subject_id, parameters, ...
        wait_for_job, char(timelimit), memorylimit, options.sequential_configs);
end
