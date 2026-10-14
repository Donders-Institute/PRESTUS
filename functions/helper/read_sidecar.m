function sc = read_sidecar(sidecar_file)
% READ_SIDECAR  Read a stage-output JSON sidecar
%
% Reads a sidecar written by WRITE_SIDECAR (or an older PRESTUS JSON sidecar
% such as pseudoCT.json with only SkullMapping). JSON arrays come back as
% column vectors or struct arrays; compare signatures with
% CHECK_SIDECAR_SIGNATURE, which normalises both sides, rather than isequal.
%
% Use as:
%   sc = read_sidecar(sidecar_file)
%
% Input:
%   sidecar_file - path to a .json sidecar
%
% Output:
%   sc - struct with the sidecar's fields
%
% See also: WRITE_SIDECAR, FIND_SIDECAR, CHECK_SIDECAR_SIGNATURE

arguments
    sidecar_file {mustBeTextScalar}
end

    sidecar_file = char(sidecar_file);
    if ~isfile(sidecar_file)
        error('prestus:sidecar:notFound', 'Sidecar not found: %s', sidecar_file);
    end
    sc = jsondecode(fileread(sidecar_file));
    if ~isstruct(sc) || ~isscalar(sc)
        error('prestus:sidecar:badContent', 'Sidecar does not hold a JSON object: %s', sidecar_file);
    end
end
