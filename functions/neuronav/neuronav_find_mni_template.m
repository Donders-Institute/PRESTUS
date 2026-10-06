function f = neuronav_find_mni_template(parameters)
% NEURONAV_FIND_MNI_TEMPLATE  SimNIBS' MNI152_T1_1mm.nii.gz, for head outlines in MNI views
%
% Looks in the SimNIBS Python package resources of startup.simnibs_bin_path,
% the PlanTUS env_path and the default install folders.
%
% Use as:
%   f = neuronav_find_mni_template(parameters)
%
% Input:
%   parameters - config struct (startup.simnibs_bin_path, placement.plantus.env_path)
%
% Output:
%   f - template path, or '' when no SimNIBS installation is found
%
% See also: NEURONAV_HEAD_SURFACE, NEURONAV_SUMMARY_TO_MNI

    f = '';
    roots = {};
    bin = char(get_nested(parameters, {'startup', 'simnibs_bin_path'}, ''));
    if ~isempty(bin), roots = {fileparts(regexprep(bin, '[\\/]+$', '')), bin}; end
    env = char(get_nested(parameters, {'placement', 'plantus', 'env_path'}, ''));
    if ~isempty(env), roots{end+1} = env; end
    home = getenv('HOME');
    if isempty(home), home = getenv('USERPROFILE'); end
    roots = [roots, {fullfile(home, 'Applications', 'SimNIBS*'), fullfile(home, 'SimNIBS*'), '/opt/SimNIBS*'}];
    pkgs = {fullfile('simnibs_env', 'lib', 'python*', 'site-packages'), ...
            fullfile('simnibs_env', 'Lib', 'site-packages')};
    for r = roots
        for k = 1:numel(pkgs)
            hit = dir(fullfile(r{1}, pkgs{k}, 'simnibs', 'resources', 'templates', 'MNI152_T1_1mm.nii.gz'));
            if ~isempty(hit), f = fullfile(hit(1).folder, hit(1).name); return; end
        end
    end
end

function val = get_nested(s, fields, default)
    val = s;
    for i = 1:numel(fields)
        if isstruct(val) && isscalar(val) && isfield(val, fields{i})
            val = val.(fields{i});
        else
            val = default; return;
        end
    end
end
