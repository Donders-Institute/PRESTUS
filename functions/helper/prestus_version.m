function [hash, version_str] = prestus_version(prestus_path)
% PRESTUS_VERSION  Print PRESTUS version info and optionally return the git hash.
%
% Resolves the PRESTUS root (via get_prestus_path if not supplied), then uses
% git to retrieve the short commit hash, branch name, and commit date.
%
% Use as:
%   prestus_version()                    % print only
%   hash = prestus_version()             % print and return short commit hash
%   [hash, ver] = prestus_version()      % also return version string (tag/VERSION)
%   hash = prestus_version(prestus_path)
%
% Input:
%   prestus_path - path to the PRESTUS root directory (optional)
%
% Output:
%   hash        - short git commit hash string, or 'unknown' if git is unavailable
%   version_str - version tag string (from git describe or VERSION file), or 'unknown'
%
% See also: KWAVE_VERSION, SIMNIBS_VERSION, PATH_LOG_SETUP

    if nargin == 0
        prestus_path = get_prestus_path;
    end

    hash        = 'unknown';
    version_str = 'unknown';

    try
        [status, git_info] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null', prestus_path));
        [~,      branch]   = system(sprintf('git -C "%s" branch --show-current 2>/dev/null', prestus_path));
        [~,      date]     = system(sprintf('git -C "%s" --no-pager log -1 --format="%%cd" --date=short 2>/dev/null', prestus_path));

        if status == 0 && ~isempty(strtrim(git_info))
            hash = strtrim(git_info);
        end

        [tag_status, tag_str] = system(sprintf('git -C "%s" describe --tags --always 2>/dev/null', prestus_path));
        if tag_status == 0 && ~isempty(strtrim(tag_str))
            version_str = strtrim(tag_str);
        else
            version_file = fullfile(prestus_path, 'VERSION');
            if isfile(version_file)
                fid = fopen(version_file, 'r');
                version_str = strtrim(fgetl(fid));
                fclose(fid);
            end
        end

        if nargout == 0
            fprintf('============================================================ \n');
            fprintf('PRESTUS %s \n', version_str);
            fprintf('Commit: %s (branch: %s)  \n', hash, strtrim(branch));
            fprintf('Commit date: %s \n', strtrim(date));
            fprintf('============================================================ \n');
        end
    catch
        if nargout == 0
            fprintf('Could not retrieve git repo status\n');
        end
    end

end