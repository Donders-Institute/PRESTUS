function ok = check_cache_signature(filename, stored, current, parameters)
% CHECK_CACHE_SIGNATURE  Refuse to reuse a cache built from different settings
%
% Compares the signature stored with a grid/medium cache (CACHE_SIGNATURE) to
% the one of the current run. io.cache_signature_check selects the response to
% a mismatch: 'error' (default) stops with the differing keys, 'warn' reports
% them and reuses the cache anyway, 'off' skips the check. A cache written
% before signatures existed cannot be verified and is reused with a warning.
%
% Use as:
%   ok = check_cache_signature(filename, stored, current, parameters)
%
% Input:
%   filename   - cache file (for messages)
%   stored     - signature saved with the cache ([] for legacy caches)
%   current    - signature of the current run
%   parameters - PRESTUS parameters (io.cache_signature_check)
%
% Output:
%   ok - true when the signatures match (or the check is off)
%
% See also: CACHE_SIGNATURE, GRID_STAGE, MEDIUM_STAGE

    mode = 'error';
    if isfield(parameters, 'io') && isfield(parameters.io, 'cache_signature_check') && ...
            ~isempty(parameters.io.cache_signature_check)
        mode = parameters.io.cache_signature_check;
        if islogical(mode) && isscalar(mode)
            mode = 'off';                % unquoted YAML 'off' loads as false
            if parameters.io.cache_signature_check; mode = 'error'; end
        end
        mode = lower(char(mode));
    end
    if ~any(strcmp(mode, {'error', 'warn', 'off'}))
        error('prestus:cache:badCheckMode', ...
            'io.cache_signature_check must be ''error'', ''warn'' or ''off'' (got ''%s'').', mode);
    end

    ok = true;
    if strcmp(mode, 'off')
        return
    end

    if isempty(stored) || ~isstruct(stored)
        warning('prestus:cache:unverified', ...
            ['Reusing %s, which was written before cache signatures existed, so it cannot be ' ...
             'checked against the current settings. Set io.overwrite_files: always to rebuild it.'], ...
            filename);
        return
    end

    diffs = struct_diff_keys(stored, current, '');
    if isempty(diffs)
        return
    end
    ok = false;
    shown = strjoin(diffs(1:min(end, 10)), ', ');
    if numel(diffs) > 10
        shown = sprintf('%s, ... (%d more)', shown, numel(diffs) - 10);
    end
    msg = sprintf(['Cache %s was built with different settings than this run (%s). ' ...
        'Reusing it would ignore those changes. Set io.overwrite_files: always to rebuild it ' ...
        '(chained pipelines - async, multi-ISPPA, uncertainty - reuse intermediate caches ' ...
        'under always: use force), or io.cache_signature_check: warn to reuse it anyway.'], ...
        filename, shown);
    if strcmp(mode, 'error')
        error('prestus:cache:signatureMismatch', '%s', msg);
    else
        warning('prestus:cache:signatureMismatch', '%s', msg);
    end
end

function keys = struct_diff_keys(a, b, prefix)
% Dotted keys whose values differ between two structs (fields of either).
    keys = {};
    if isstruct(a) && isstruct(b) && numel(a) == numel(b)
        fa = fieldnames(a);
        fb = fieldnames(b);
        all_fields = union(fa, fb, 'stable');
        for i = 1:numel(all_fields)
            f = all_fields{i};
            path = f;
            if ~isempty(prefix); path = [prefix '.' f]; end
            if ~isfield(a, f) || ~isfield(b, f)
                keys{end+1} = path; %#ok<AGROW>
                continue
            end
            for k = 1:numel(a)
                sub = struct_diff_keys(a(k).(f), b(k).(f), path);
                keys = [keys, sub(~ismember(sub, keys))]; %#ok<AGROW>
            end
        end
    elseif ~isequaln(a, b)
        keys = {prefix};
    end
end
