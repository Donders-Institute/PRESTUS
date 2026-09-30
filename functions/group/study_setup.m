function study = study_setup(study_config)
% STUDY_SETUP  Load a study config and resolve its manifest into runs
%
% Shared by PRESTUS_STUDY_START (launching) and PRESTUS_GROUP_REPORT_START
% (reporting) so both resolve the same runs and output names.
%
% Study outputs (group report, group maps, run logs) go to
%   <path.sim>/group/<study.name>/
% where study.name defaults to the manifest file name without extension.
%
% Use as:
%   study = study_setup(study_config)
%
% Input:
%   study_config - path to the study YAML, or a parameters struct already
%                  loaded with LOAD_PARAMETERS (a relative study.manifest
%                  then resolves against the current folder)
%
% Output:
%   study - struct with fields
%             config        absolute path of the study YAML ('' for a struct)
%             parameters    loaded study config (no row applied)
%             manifest_path absolute manifest path
%             manifest      table from STUDY_READ_MANIFEST
%             runs          table from STUDY_RUN_TABLE (all rows, see .include)
%             name          study name
%             out_dir       <path.sim>/group/<name>
%
% See also: STUDY_READ_MANIFEST, STUDY_RUN_TABLE, PRESTUS_STUDY_START,
%           PRESTUS_GROUP_REPORT_START

    if ischar(study_config) || isstring(study_config)
        cfg = char(study_config);
        if ~isfile(cfg)
            error('study_setup:configNotFound', 'Study config not found: %s', cfg);
        end
        cfg = absolute_path(cfg);
        [cdir, cname, cext] = fileparts(cfg);
        parameters = load_parameters([cname cext], cdir);
        base_dir = cdir;
    elseif isstruct(study_config)
        cfg = '';
        parameters = study_config;
        base_dir = pwd;
    else
        error('study_setup:badConfig', 'study_config must be a YAML path or a parameters struct.');
    end

    if ~isfield(parameters, 'study') || ~isfield(parameters.study, 'manifest') || ...
            isempty(parameters.study.manifest)
        error('study_setup:noManifest', ...
            'The study config does not set study.manifest (see doc_group.md, "Study setup").');
    end
    mpath = char(parameters.study.manifest);
    if ~startsWith(mpath, filesep) && isempty(regexp(mpath, '^[A-Za-z]:[\\/]', 'once'))
        mpath = fullfile(base_dir, mpath);
    end

    [manifest, ~] = study_read_manifest(mpath);
    runs = study_run_table(parameters, manifest);

    name = '';
    if isfield(parameters.study, 'name') && ~isempty(parameters.study.name)
        name = char(parameters.study.name);
    end
    if isempty(name)
        [~, name] = fileparts(mpath);
    end
    name = regexprep(name, '[^a-zA-Z0-9_-]', '_');

    if ~isfield(parameters, 'path') || ~isfield(parameters.path, 'sim') || isempty(parameters.path.sim)
        error('study_setup:noSimPath', 'The study config must set path.sim.');
    end

    study = struct('config', cfg, 'parameters', parameters, 'manifest_path', mpath, ...
        'manifest', manifest, 'runs', runs, 'name', name, ...
        'out_dir', fullfile(char(parameters.path.sim), 'group', name));
end

function p = absolute_path(p)
    if ~startsWith(p, filesep) && isempty(regexp(p, '^[A-Za-z]:[\\/]', 'once'))
        p = fullfile(pwd, p);
    end
end
