function legacy_hpc_submit(platform, caller, subject_id, parameters, wait_for_job, timelimit, memorylimit, sequential_configs)
% LEGACY_HPC_SUBMIT  Shared body of the v0.5.0 HPC entry points (SINGLE_SUBJECT_PIPELINE_WITH_SLURM / _QSUB)
%
% Converts old parameter names (LEGACY_BRIDGE_STRUCT), applies the platform and
% the job settings the v0.5.0 function took as arguments, warns
% (prestus:legacyFunction) with the equivalent current call, and submits
% through PRESTUS_PIPELINE_START.
%
% Use as:
%   legacy_hpc_submit(platform, caller, subject_id, parameters, wait_for_job, timelimit, memorylimit, sequential_configs)
%
% Input:
%   platform           - 'slurm' or 'qsub'
%   caller             - name of the v0.5.0 function (for messages)
%   subject_id         - subject number
%   parameters         - struct from LOAD_PARAMETERS (old key names allowed)
%   wait_for_job       - logical
%   timelimit          - 'HH:MM:SS'
%   memorylimit        - GB
%   sequential_configs - struct (empty: no chain)
%
% See also: PRESTUS_PIPELINE_START, LEGACY_BRIDGE_STRUCT

    seq = '';
    if ~isempty(fieldnames(sequential_configs))
        seq = ', struct(''sequential_configs'', sequential_configs)';
    end
    warning('prestus:legacyFunction', ...
        ['%s is an entry point of PRESTUS v0.5.0; it still works. Update the script to:\n' ...
         '  parameters.subject_id = %d;\n  parameters.platform = ''%s'';\n' ...
         '  parameters.hpc.wait_for_job = %d;\n  parameters.hpc.timelimit = ''%s'';\n' ...
         '  parameters.hpc.memorylimit = %g;\n  prestus_pipeline_start(parameters%s)'], ...
        caller, subject_id, platform, wait_for_job, timelimit, memorylimit, seq);

    parameters.subject_id = subject_id;
    parameters = legacy_bridge_struct(parameters, caller);
    parameters.platform = platform;
    parameters.hpc.wait_for_job = logical(wait_for_job);
    parameters.hpc.timelimit = timelimit;
    parameters.hpc.memorylimit = memorylimit;
    options = struct();
    if ~isempty(fieldnames(sequential_configs))
        options.sequential_configs = sequential_configs;
    end
    prestus_pipeline_start(parameters, options);
end
