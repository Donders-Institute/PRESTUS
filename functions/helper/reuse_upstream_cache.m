function reuse = reuse_upstream_cache(filename, parameters, kind)
% REUSE_UPSTREAM_CACHE  True when an existing upstream cache must be loaded
%                       regardless of io.overwrite_files.
%
%   Chained pipelines (async, multi-ISPPA, uncertainty) set
%   io.reuse_upstream_caches on output-stage sub-jobs. Those jobs honour
%   io.overwrite_files = 'always' for their own products (acoustic/heating
%   results, NIfTIs, CSV, report) but must load the caches that earlier
%   stages built, since upstream results were computed on them.
%
% Use as:
%   if reuse_upstream_cache(f, parameters, 'grid') || ~confirm_overwriting(f, parameters)
%       % load cache
%   end
%
% Input:
%   filename   - path of the upstream cache file
%   parameters - PRESTUS config; reads io.reuse_upstream_caches, a cellstr
%                of cache kinds to reuse (e.g. {'head','grid','medium'}; default {})
%   kind       - cache kind of filename: 'head' | 'grid' | 'medium'
%
% Output:
%   reuse      - [1x1] logical, true if the file exists and its kind is listed
%
% See also: CONFIRM_OVERWRITING

    reuse = isfield(parameters, 'io') && ...
        isfield(parameters.io, 'reuse_upstream_caches') && ...
        any(strcmp(cellstr(parameters.io.reuse_upstream_caches), kind)) && ...
        isfile(filename);
end
