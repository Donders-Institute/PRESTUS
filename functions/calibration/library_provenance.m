function [version_str, commit] = library_provenance()
% LIBRARY_PROVENANCE  PRESTUS version string and commit hash for library entry meta
%
% Every calibration writer stamps meta.prestus_version (version string) and
% meta.prestus_commit (short git hash) from this function, so all entries
% carry the same fields. Both are 'unknown' when they cannot be resolved.
%
% Use as:
%   [lib.meta.prestus_version, lib.meta.prestus_commit] = library_provenance()
%
% Output:
%   version_str - version tag (git describe or VERSION file), or 'unknown'
%   commit      - short git commit hash, or 'unknown'
%
% See also: PRESTUS_VERSION, SAVE_LIBRARY_ENTRY, CALIBRATION_SAVE_ENTRY

    try
        [commit, version_str] = prestus_version();
    catch
        commit = 'unknown';
        version_str = 'unknown';
    end
    commit = char(commit);
    version_str = char(version_str);
end
