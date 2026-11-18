function ok = check_pct_sidecar(pct_dir, parameters)
% CHECK_PCT_SIDECAR  Refuse an existing pseudo-CT made with other settings
%
% An existing pseudoCT.nii.gz is reused instead of generated again
% (PREPROC_SEGMENTATION); it must have been made with this run's skull
% mapping (PCT_SIGNATURE). Reads pseudoCT.json (PCT_CREATE_PCT_RUN); a sidecar
% without a signature but with SkullMapping has that compared. A mismatch is handled per
% io.cache_signature_check (error | warn | off); a pseudo-CT without a
% sidecar is used with a warning.
%
% Use as:
%   ok = check_pct_sidecar(pct_dir, parameters)
%
% Input:
%   pct_dir    - folder holding pseudoCT.nii.gz (FIND_PSEUDOCT_DIR)
%   parameters - PRESTUS parameters (pct.*, io.cache_signature_check)
%
% Output:
%   ok - true when the recorded settings match (or cannot be checked)
%
% See also: PCT_SIGNATURE, CHECK_SIDECAR_SIGNATURE, FIND_PSEUDOCT_DIR

arguments
    pct_dir    {mustBeTextScalar}
    parameters (1,1) struct
end

    pct_dir = char(pct_dir);
    sidecar = find_sidecar(fullfile(pct_dir, 'pseudoCT.nii.gz'));
    stored = [];
    if ~isempty(sidecar)
        sc = read_sidecar(sidecar);
        if isfield(sc, 'signature')
            stored = sc.signature;
        elseif isfield(sc, 'SkullMapping')          % sidecar without a signature
            stored = struct('signature_version', 1, 'skull_mapping', sc.SkullMapping);
        end
    else
        sidecar = fullfile(pct_dir, 'pseudoCT.nii.gz');
    end
    ok = check_sidecar_signature(sidecar, stored, pct_signature(parameters), parameters, ...
        sprintf(['Move %s aside (or set path.pct to another folder) to generate a pseudo-CT ' ...
                 'with the current settings'], pct_dir));
end
