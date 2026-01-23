function prestus_version(prestus_path)
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