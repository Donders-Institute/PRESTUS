function t1_path = neuronav_find_planning_t1(sub_id, parameters)
% NEURONAV_FIND_PLANNING_T1  Locate the neuronavigation planning T1 of a subject.
%
% Localite coordinates refer to the T1 loaded into Localite for planning, so
% voxel conversions of Localite positions need that image's header. It is
% searched in <path.localite_t1>/<sub_id>/ and then <path.localite_t1>/
% (path.localite_t1 defaults to the raw Localite folder, see
% NEURONAV_RESOLVE_PATHS) with the pattern path.localite_t1_pattern
% (sprintf template taking sub_id; default '%s_*T1*.nii*'; the older key
% path.neuronav_t1_pattern is still read).
%
% Use as:
%   t1_path = neuronav_find_planning_t1(sub_id, parameters)
%
% See also: NEURONAV_RESOLVE_PATHS, NEURONAV_CONVERT_TRIGGER_TO_VOXELS

    paths = neuronav_resolve_paths(parameters);
    if isempty(paths.t1)
        error('neuronav_find_planning_t1:noFolder', ...
            'Set path.localite_t1 (or path.localite_raw) to find the planning T1 of %s.', sub_id);
    end

    pat = '%s_*T1*.nii*';
    for key = {'localite_t1_pattern', 'neuronav_t1_pattern'}
        if isfield(parameters.path, key{1}) && ~isempty(parameters.path.(key{1}))
            pat = char(parameters.path.(key{1}));
            break;
        end
    end
    t1_pat = sprintf(pat, sub_id);

    hits = dir(fullfile(paths.t1, sub_id, t1_pat));
    if isempty(hits)
        hits = dir(fullfile(paths.t1, t1_pat));
    end
    if isempty(hits)
        error('neuronav_find_planning_t1:notFound', ...
            'No planning T1 for %s.\nSearched:\n  %s\n  %s\nPattern: %s', ...
            sub_id, fullfile(paths.t1, sub_id), paths.t1, t1_pat);
    end
    t1_path = fullfile(hits(1).folder, hits(1).name);
end
