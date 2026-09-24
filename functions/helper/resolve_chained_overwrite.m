function ow = resolve_chained_overwrite(parameters)
% RESOLVE_CHAINED_OVERWRITE  Map io.overwrite_files onto the stages of a
%                            chained pipeline (async, multi-ISPPA, uncertainty).
%
%   Chained pipelines split a run into intermediate stages (preprocessing,
%   per-transducer / stage-1 acoustics, field combination) whose caches feed
%   output stages (thermal, CSV, NIfTI, HTML report).
%
%     overwrite_files   intermediate stages          output stages
%     'never'           resume from existing cache    skip if output exists
%     'always'          resume from existing cache    always re-run
%     'force'           always re-run                 always re-run
%
%   Any other value (e.g. 'ask', or unset) behaves like 'never'; sub-jobs are
%   non-interactive and cannot prompt.
%
% Use as:
%   ow = resolve_chained_overwrite(parameters)
%
% Input:
%   parameters - PRESTUS config; reads io.overwrite_files
%
% Output:
%   ow - struct with fields
%        .rerun_outputs        logical, re-run output stages even if outputs exist
%        .rerun_intermediates  logical, re-run intermediate stages even if cached
%        .outputs              overwrite_files for output-stage sub-jobs
%        .intermediates        overwrite_files for intermediate-stage sub-jobs
%
% See also: CONFIRM_OVERWRITING, REUSE_UPSTREAM_CACHE

    mode = 'never';
    if isfield(parameters, 'io') && isfield(parameters.io, 'overwrite_files') && ...
            ~isempty(parameters.io.overwrite_files)
        mode = char(parameters.io.overwrite_files);
    end

    ow.rerun_outputs       = any(strcmp(mode, {'always', 'force'}));
    ow.rerun_intermediates = strcmp(mode, 'force');
    ow.outputs             = 'never';
    ow.intermediates       = 'never';
    if ow.rerun_outputs;       ow.outputs       = 'always'; end
    if ow.rerun_intermediates; ow.intermediates = 'always'; end
end
