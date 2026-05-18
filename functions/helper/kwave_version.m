function kwave_version(kpath)
% KWAVE_VERSION  Print k-Wave toolbox version, git commit, and CUDA binary version
%
% Queries getComputerInfo for the k-Wave version string, the kpath git
% repository for the short commit hash and date, and extracts the version
% string embedded in the CUDA binary (kspaceFirstOrder-CUDA) via `strings`.
%
% Use as:
%   kwave_version()
%   kwave_version(kpath)
%
% Input:
%   kpath - path to the k-Wave root directory (optional, default: getkWavePath)
%
% See also: PRESTUS_VERSION, SIMNIBS_VERSION

    if nargin == 0
        kpath = getkWavePath;
    end

    try
        % Use git tag if kwave is a git repo; otherwise fall back to getkWaveVersion
        [is_git] = system(sprintf('git -C "%s" rev-parse --git-dir > /dev/null 2>&1', kpath));
        if is_git == 0
            [status, tag] = system(sprintf('git -C "%s" describe --tags --exact-match 2>/dev/null', kpath));
            if status == 0 && ~isempty(strtrim(tag))
                version_str = strtrim(tag);
                if startsWith(version_str, 'v')
                    version_str = version_str(2:end);
                end
            else
                setupInfo = getComputerInfo;
                version_str = setupInfo.kwave_version;
            end
        else
            setupInfo = getComputerInfo;
            version_str = setupInfo.kwave_version;
        end
        fprintf('k-Wave toolbox: %s\n', version_str);

        [~, hash] = system(sprintf('git -C "%s" rev-parse --short HEAD', kpath));
        hash = strtrim(hash);
        [~, date] = system(sprintf('git -C "%s" --no-pager log -1 --format="%%cd" --date=short', kpath));
        date = strtrim(date);

        fprintf('Git commit: %s\n', hash);
        fprintf('Commit date: %s\n', date);
    catch
        fprintf('Could not retrieve git repo status\n');
        return;
    end

    % Report CUDA binary version (embedded string, may differ from toolbox version)
    if isunix
        cuda_binary = fullfile(getkWavePath('binaries'), 'kspaceFirstOrder-CUDA');
    else
        cuda_binary = fullfile(getkWavePath('binaries'), 'kspaceFirstOrder-CUDA.exe');
    end
    if exist(cuda_binary, 'file')
        [status, out] = system(sprintf('strings "%s" 2>/dev/null | grep -m1 "kspaceFirstOrder-CUDA v"', cuda_binary));
        if status == 0 && ~isempty(strtrim(out))
            fprintf('CUDA binary:    %s\n', strtrim(out));
        else
            fprintf('CUDA binary:    version string not found in %s\n', cuda_binary);
        end
    end

end
