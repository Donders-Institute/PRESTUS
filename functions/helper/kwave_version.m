function kwave_version(kpath)
    if nargin == 0
        kpath = getkWavePath;
    end
    % establish general version of kWave
    setupInfo = getComputerInfo;
    fprintf('k-Wave toolbox: %s\n', setupInfo.kwave_version);
    
    try
        cd(kpath);
        
        % Shell-based git queries
        [~, hash] = system('git rev-parse --short HEAD');
        hash = hash(1:end-1);
        [~, date] = system('git --no-pager log -1 --format="%cd" --date=short');
        date = date(1:end-1);
        
        fprintf('Git commit: %s\n', strtrim(hash));
        fprintf('Commit date: %s\n', strtrim(date));
    catch
        fprintf('Could not retrieve git repo status\n');
        return;
    end

end
