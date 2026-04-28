function hash = prestus_version(prestus_path)
% PRESTUS_VERSION  Print PRESTUS version info and optionally return the git hash.
%
% Resolves the PRESTUS root (via GET_PRESTUSPATH if not supplied), then uses
% git to retrieve the short commit hash, branch name, and commit date.
%
% Use as:
%   prestus_version()           % print only
%   hash = prestus_version()    % print and return short commit hash
%   hash = prestus_version(prestus_path)
%
% Input:
%   prestus_path - path to the PRESTUS root directory (optional)
%
% Output:
%   hash - short git commit hash string, or 'unknown' if git is unavailable
%
% See also: KWAVE_VERSION, SIMNIBS_VERSION, PATH_LOG_SETUP

    if nargin == 0
        prestus_path = get_prestus_path;
    end

    hash = 'unknown';

    try
        [status, git_info] = system(sprintf('git -C "%s" rev-parse --short HEAD 2>/dev/null', prestus_path));
        [~,      branch]   = system(sprintf('git -C "%s" branch --show-current 2>/dev/null', prestus_path));
        [~,      date]     = system(sprintf('git -C "%s" --no-pager log -1 --format="%%cd" --date=short 2>/dev/null', prestus_path));

        if status == 0 && ~isempty(strtrim(git_info))
            hash = strtrim(git_info);
        end

        if nargout == 0
            fprintf('============================================================ \n');
            fprintf('PRESTUS 0.4.0 + \n');
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