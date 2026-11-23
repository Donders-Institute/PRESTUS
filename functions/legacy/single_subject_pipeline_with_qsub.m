function single_subject_pipeline_with_qsub(subject_id, parameters, wait_for_job, timelimit, memorylimit)
% SINGLE_SUBJECT_PIPELINE_WITH_QSUB  qsub entry point of PRESTUS v0.5.0 and earlier; submits via PRESTUS_PIPELINE_START
%
% Kept so that scripts written for v0.5.0 keep running. Warns
% (prestus:legacyFunction) with the current call. Old parameter names are
% converted (LEGACY_BRIDGE_STRUCT); the arguments set hpc.wait_for_job,
% hpc.timelimit (v0.5.0 took seconds; converted to HH:MM:SS) and
% hpc.memorylimit, with the v0.5.0 defaults (no wait, 4 h, 40 GB) when omitted.
%
% Use as:
%   single_subject_pipeline_with_qsub(subject_id, parameters)
%   single_subject_pipeline_with_qsub(subject_id, parameters, wait_for_job, timelimit_s, memorylimit)
%
% See also: PRESTUS_PIPELINE_START, SINGLE_SUBJECT_PIPELINE, LEGACY_BRIDGE_ENTRY

    arguments
        subject_id (1,1) double
        parameters struct
        wait_for_job logical = false
        timelimit (1,1) double = 60*60*4
        memorylimit (1,1) double = 40
    end
    hms = sprintf('%02d:%02d:%02d', floor(timelimit/3600), floor(mod(timelimit, 3600)/60), round(mod(timelimit, 60)));
    legacy_hpc_submit('qsub', 'single_subject_pipeline_with_qsub', subject_id, parameters, ...
        wait_for_job, hms, memorylimit, struct());
end
