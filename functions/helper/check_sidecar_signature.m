function [ok, diffs] = check_sidecar_signature(sidecar_file, stored, current, parameters, remedy)
% CHECK_SIDECAR_SIGNATURE  Compare a sidecar's recorded settings with the current run's
%
% Like CHECK_CACHE_SIGNATURE (and with the same io.cache_signature_check
% response: 'error' | 'warn' | 'off'), for signatures read back from a JSON
% sidecar. Both signatures are first passed through JSON, so that the
% comparison does not depend on how JSON represents arrays (rows come back
% as columns), empty values or numeric precision.
%
% Use as:
%   ok = check_sidecar_signature(sidecar_file, stored, current, parameters, remedy)
%
% Input:
%   sidecar_file - sidecar the stored signature came from (for messages)
%   stored       - signature read from the sidecar ([] when it has none)
%   current      - signature of the current run (struct)
%   parameters   - PRESTUS parameters (io.cache_signature_check)
%   remedy       - (optional) how to rebuild the output, for the mismatch message
%
% Output:
%   ok    - true when the signatures match (or the check is off)
%   diffs - dotted keys whose values differ
%
% See also: CHECK_CACHE_SIGNATURE, WRITE_SIDECAR, READ_SIDECAR

arguments
    sidecar_file {mustBeTextScalar}
    stored
    current      (1,1) struct
    parameters   (1,1) struct
    remedy       {mustBeTextScalar} = ''
end

    if ~isempty(stored)
        stored = normalise(stored);
    end
    [ok, diffs] = check_cache_signature(sidecar_file, stored, normalise(current), parameters, remedy);
end

function s = normalise(s)
    s = jsondecode(jsonencode(s));
end
