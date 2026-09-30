function paths = neuronav_resolve_paths(parameters)
% NEURONAV_RESOLVE_PATHS  Resolve the Localite folders of a PRESTUS config.
%
% Each Localite path key has one meaning:
%
%   path.localite_raw   raw Localite exports, read-only
%                       <raw>/sub-XXX/ses-XX/[<placement.localite.loc_folder>/]…
%   path.localite_post  ingested outputs written by NEURONAV_INGEST_MARKERS
%                       <post>/sub-XXX/ses-XX/InstrumentMarker_sub-XXX_ses-XX.xml/.json
%   path.localite_t1    folder with the neuronavigation planning T1s
%                       (<t1>/[sub-XXX/]<path.localite_t1_pattern>); defaults to
%                       the raw folder
%
% path.localite is a deprecated alias: when path.localite_raw is empty it is
% used as the raw folder (with a warning). Its other legacy use, as the output
% folder for heuristic-placement exports, is now path.localite_export
% (resolved per subject by PATH_LOG_SETUP into io.dir_localite_export).
%
% Use as:
%   paths = neuronav_resolve_paths(parameters)
%
% Output:
%   paths - struct with char fields .raw, .post, .t1 ('' when not set)
%
% See also: NEURONAV_INGEST_MARKERS, NEURONAV_SELECT_LOCALITE, NEURONAV_LOCALITE_POSE

    p = struct();
    if isfield(parameters, 'path') && isstruct(parameters.path)
        p = parameters.path;
    end

    paths.raw  = get_char(p, 'localite_raw');
    paths.post = get_char(p, 'localite_post');
    paths.t1   = get_char(p, 'localite_t1');

    legacy = get_char(p, 'localite');
    if isempty(paths.raw) && ~isempty(legacy)
        warning('PRESTUS:localite:deprecatedPath', ...
            ['path.localite is deprecated as the Localite data folder; set ' ...
             'path.localite_raw (raw exports) and path.localite_post (ingested ' ...
             'outputs) instead. Using path.localite = %s as path.localite_raw.'], legacy);
        paths.raw = legacy;
    end
    if isempty(paths.t1)
        paths.t1 = paths.raw;
    end
end

function v = get_char(s, f)
    v = '';
    if isfield(s, f) && ~isempty(s.(f))
        v = char(s.(f));
    end
end
