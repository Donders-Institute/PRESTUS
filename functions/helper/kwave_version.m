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
        % Use git tag if kwave is a git repo; otherwise fall back to
        % getkWaveVersion. run_command invokes git directly via argv, no
        % shell at all -- kpath reaches the process as a literal argument
        % regardless of its contents.
        [is_git] = run_command({'git', '-C', kpath, 'rev-parse', '--git-dir'});
        if is_git == 0
            [status, tag] = run_command({'git', '-C', kpath, 'describe', '--tags', '--exact-match'});
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

        [~, hash] = run_command({'git', '-C', kpath, 'rev-parse', '--short', 'HEAD'});
        hash = strtrim(hash);
        [~, date] = run_command({'git', '-C', kpath, '--no-pager', 'log', '-1', '--format=%cd', '--date=short'});
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
        % run_command invokes `strings` directly via argv; the `| grep -m1
        % "..."` first-match selection is replicated in MATLAB instead of
        % a shell pipe.
        [status, strings_out] = run_command({'strings', cuda_binary});
        match_line = '';
        if status == 0
            str_lines = strsplit(strings_out, newline);
            for li = 1:numel(str_lines)
                if contains(str_lines{li}, 'kspaceFirstOrder-CUDA v')
                    match_line = str_lines{li};
                    break;
                end
            end
        end
        out = match_line;
        if status == 0 && ~isempty(strtrim(out))
            fprintf('CUDA binary:    %s\n', strtrim(out));
        else
            fprintf('CUDA binary:    version string not found in %s\n', cuda_binary);
        end
    end

end
