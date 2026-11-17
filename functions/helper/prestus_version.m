function [hash, version_str, info] = prestus_version(prestus_path)
% PRESTUS_VERSION  Print PRESTUS version info and optionally return the git hash.
%
% Resolves the PRESTUS root (via get_prestus_path if not supplied), then uses
% git to retrieve the short commit hash, branch name, and commit date.
%
% Use as:
%   prestus_version()                    % print only
%   hash = prestus_version()             % print and return short commit hash
%   [hash, ver] = prestus_version()      % also return version string (tag/VERSION)
%   [hash, ver, info] = prestus_version()   % also return a provenance struct
%   hash = prestus_version(prestus_path)
%
% Input:
%   prestus_path - path to the PRESTUS root directory (optional)
%
% Output:
%   hash        - short git commit hash string, or 'unknown' if git is unavailable
%   version_str - version tag string (from git describe or VERSION file), or 'unknown'
%   info        - struct with fields version, commit, commit_date (yyyy-mm-dd),
%                 branch and dirty ('yes' when tracked files differ from the
%                 commit, 'no', or 'unknown'). The dirty check runs only when
%                 info is requested.
%
% See also: KWAVE_VERSION, SIMNIBS_VERSION, PATH_LOG_SETUP

    if nargin == 0
        prestus_path = get_prestus_path;
    end

    hash        = 'unknown';
    version_str = 'unknown';
    branch      = '';
    date        = '';
    dirty       = 'unknown';

    try
        % run_command invokes git directly via argv, no shell at all --
        % prestus_path reaches the process as a literal argument
        % regardless of its contents.
        [status, git_info] = run_command({'git', '-C', prestus_path, 'rev-parse', '--short', 'HEAD'});
        [b_status, b_out]  = run_command({'git', '-C', prestus_path, 'branch', '--show-current'});
        [d_status, d_out]  = run_command({'git', '-C', prestus_path, '--no-pager', 'log', '-1', '--format=%cd', '--date=short'});
        if b_status == 0, branch = strtrim(b_out); end   % else git's error text, not a branch
        if d_status == 0, date   = strtrim(d_out); end

        if status == 0 && ~isempty(strtrim(git_info))
            hash = strtrim(git_info);
        end

        [tag_status, tag_str] = run_command({'git', '-C', prestus_path, 'describe', '--tags', '--always'});
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

        if nargout >= 3 && ~strcmp(hash, 'unknown')
            % Tracked files only: untracked files (local scratch) do not
            % change what the toolbox computes.
            [st_status, st_out] = run_command({'git', '-C', prestus_path, 'status', ...
                '--porcelain', '--untracked-files=no'});
            if st_status == 0
                if isempty(strtrim(st_out)), dirty = 'no'; else, dirty = 'yes'; end
            end
        end

        if nargout == 0
            fprintf('============================================================ \n');
            fprintf('PRESTUS %s \n', version_str);
            fprintf('Commit: %s (branch: %s)  \n', hash, branch);
            fprintf('Commit date: %s \n', date);
            fprintf('============================================================ \n');
        end
    catch
        if nargout == 0
            fprintf('Could not retrieve git repo status\n');
        end
    end

    info = struct('version', version_str, 'commit', hash, 'commit_date', date, ...
        'branch', branch, 'dirty', dirty);

end
