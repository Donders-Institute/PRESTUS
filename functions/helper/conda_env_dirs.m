function dirs = conda_env_dirs()
% CONDA_ENV_DIRS  Folders of the conda environments on this computer
%
% Lists the environments under the usual conda installs in the home folder
% (.conda, miniforge3, miniconda3, anaconda3, mambaforge) and the ones conda
% records in ~/.conda/environments.txt, plus $CONDA_PREFIX. Mirrors
% _conda_env_bins in nexus/prestus_nexus/simnibs.py.
%
% Use as:
%   dirs = conda_env_dirs()
%
% Output:
%   dirs - (1,:) cell of char environment folders (no bin/ or Scripts\),
%          existing folders only, without duplicates
%
% See also: FIND_NEXUS_PYTHON, PRESTUS_CHECK_INSTALL

    home = char(java.lang.System.getProperty('user.home'));
    dirs = {};
    roots = {'.conda', 'miniforge3', 'miniconda3', 'anaconda3', 'mambaforge'};
    for k = 1:numel(roots)
        listing = dir(fullfile(home, roots{k}, 'envs'));
        for j = 1:numel(listing)
            if listing(j).isdir && ~startsWith(listing(j).name, '.')
                dirs{end + 1} = fullfile(listing(j).folder, listing(j).name); %#ok<AGROW>
            end
        end
    end

    txt = fullfile(home, '.conda', 'environments.txt');
    if isfile(txt)
        lines = strtrim(splitlines(fileread(txt)));
        dirs = [dirs, lines(~cellfun(@isempty, lines))'];
    end

    prefix = getenv('CONDA_PREFIX');
    if ~isempty(prefix)
        dirs{end + 1} = prefix;
    end

    dirs = dirs(cellfun(@isfolder, dirs));
    [~, keep] = unique(dirs, 'stable');
    dirs = dirs(sort(keep));
end
