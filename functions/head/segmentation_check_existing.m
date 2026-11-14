function diffs = segmentation_check_existing(segmentation_folder, parameters)
% SEGMENTATION_CHECK_EXISTING  How an existing segmentation differs from the selected SimNIBS
%
% Compares an existing m2m folder with what the selected SimNIBS install and
% charm settings would produce: the SimNIBS version (first line of
% charm_log.html vs. simnibs_install_info) and the charm settings
% (settings.ini, which charm copies into the m2m folder, vs.
% segmentation.settings_file or the install's charm.ini). Comments, blank
% lines and spacing are ignored. Mirrors nexus/prestus_nexus/simnibs.py.
%
% Use as:
%   diffs = segmentation_check_existing(segmentation_folder, parameters)
%
% Input:
%   segmentation_folder - m2m_sub-XXX folder
%   parameters          - PRESTUS config struct
%
% Output:
%   diffs - cellstr, one entry per difference; empty when they match or
%           cannot be compared
%
% See also: PREPROC_SEGMENTATION, SIMNIBS_INSTALL_INFO, SIMNIBS_VERSION

arguments
    segmentation_folder (1,:) char
    parameters          (1,1) struct
end

    diffs = {};
    if ~isfile(fullfile(segmentation_folder, 'final_tissues.nii.gz'))
        return
    end

    bin_path = '';
    if isfield(parameters, 'startup') && isfield(parameters.startup, 'simnibs_bin_path')
        bin_path = char(parameters.startup.simnibs_bin_path);
    end
    info = struct('version', '', 'default_settings', '');
    if ~isempty(bin_path)
        info = simnibs_install_info(bin_path);
    end

    % version
    made_with = simnibs_version(segmentation_folder, struct());
    made_with = made_with.simnibs_version;
    if ~isempty(info.version) && ~any(strcmp(made_with, {'unknown', 'error'})) ...
            && ~strcmp(made_with, info.version)
        diffs{end+1} = sprintf('made with SimNIBS %s; the selected SimNIBS is %s', made_with, info.version);
    end

    % charm settings
    settings_file = '';
    if isfield(parameters, 'segmentation') && isfield(parameters.segmentation, 'settings_file')
        settings_file = char(parameters.segmentation.settings_file);
    end
    requested = settings_file;
    if isempty(requested)
        requested = info.default_settings;
    end
    used = fullfile(segmentation_folder, 'settings.ini');
    if ~isempty(requested) && isfile(requested) && isfile(used) ...
            && ~isequal(normalised_settings(used), normalised_settings(requested))
        if isempty(settings_file)
            what = 'the selected install''s defaults';
        else
            what = 'segmentation.settings_file';
        end
        diffs{end+1} = sprintf('made with other charm settings than %s', what);
    end
end


function lines = normalised_settings(file)
% Setting lines without comments, blank lines and spacing (order kept).
    raw = splitlines(fileread(file));
    lines = {};
    for i = 1:numel(raw)
        s = regexprep(raw{i}, '#.*$', '');
        s = strjoin(strsplit(strtrim(s)), ' ');
        s = regexprep(s, '\s*=\s*', '=');
        if ~isempty(s)
            lines{end+1} = s; %#ok<AGROW>
        end
    end
end
