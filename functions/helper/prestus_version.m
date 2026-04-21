function prestus_version(prestus_path)
% PRESTUS_VERSION  Print PRESTUS version, git commit hash, and branch to console
%
% Resolves the PRESTUS root (via GET_PRESTUSPATH if not supplied), then uses
% git to retrieve the short commit hash, branch name, and commit date.
%
% Use as:
%   prestus_version()
%   prestus_version(prestus_path)
%
% Input:
%   prestus_path - path to the PRESTUS root directory (optional, default: get_prestus_path)
%
% See also: KWAVE_VERSION, SIMNIBS_VERSION, PATH_LOG_SETUP

    if nargin == 0
        prestus_path = get_prestus_path;
    end
    
    try
        cd(prestus_path);
        
        % Shell-based git queries
        [status, git_info] = system('git rev-parse --short HEAD');
        [~, branch] = system('git branch --show-current');
        branch = strtrim(branch);
        [~, date] = system('git --no-pager log -1 --format="%cd" --date=short');
        
        fprintf('============================================================ \n');
        fprintf('PRESTUS 0.4.0 + \n');
        fprintf('Commit: %s (branch: %s)  \n', ...
            strtrim(git_info), branch);
        fprintf('Commit date: %s \n', strtrim(date));
        fprintf('============================================================ \n');
    catch
        fprintf('Could not retrieve git repo status\n');
        return;
    end

end