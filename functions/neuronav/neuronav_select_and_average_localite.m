function localite = neuronav_select_and_average_localite(sub_id, session, pn)
% NEURONAV_SELECT_AND_AVERAGE_LOCALITE - Select most recent Localite XML for a session.
% Compares against Localite files from the immediately preceding session (if available) to skip duplicates.
%
% INPUT:
%   sub_id  - e.g., 'sub-010'
%   session - e.g., 'ses-03' (format must be 'ses-XX' with XX numeric)
%   pn      - struct with .data_postlocalite field
%
% OUTPUT:
%   localite - struct parsed from selected Localite XML file ([] if none match)

    localite = [];
    excludeFile = struct('name', {});

    % --- Parse session number from input ---
    ses_num = regexp(session, 'ses-(\d{2})', 'tokens', 'once');
    if isempty(ses_num)
        error('Session input must be in the format "ses-XX" where XX is two digits');
    end
    ses_num = str2double(ses_num{1});

    % Define previous session string (with zero padding)
    prev_ses_num = ses_num - 1;
    if prev_ses_num < 1
        session_earlier = ''; % No earlier session available
    else
        session_earlier = sprintf('ses-%02d', prev_ses_num);
    end

    % --- Helper function to list valid Localite files ---
    function [files, datetimes] = list_valid_localite_files(basepath)
        files = dir(fullfile(basepath, 'TriggerMarkers_Coil0*.xml'));
        files = files([files.bytes] > 20480);
        datetimes = datetime.empty;
        for i = 1:length(files)
            tokens = regexp(files(i).name, 'TriggerMarkers_Coil0_(\d{17})', 'tokens', 'once');
            if ~isempty(tokens)
                try
                    t = datetime(tokens{1}, 'InputFormat', 'yyyyMMddHHmmssSSS');
                    datetimes(end+1) = t;
                catch
                    datetimes(end+1) = NaT;
                end
            else
                datetimes(end+1) = NaT;
            end
        end
        valid_idx = ~isnat(datetimes);
        files = files(valid_idx);
        datetimes = datetimes(valid_idx);
    end

    % --- List files for current session ---
    localite_path_later = fullfile(pn.data_postlocalite, sub_id, session, 'localite', '*', 'TMSTrigger');
    [localite_files, dt_later] = list_valid_localite_files(localite_path_later);

    if isempty(localite_files)
        warning("⚠ No Localite files >20kB with valid timestamps found for %s %s", sub_id, session);
        return;
    end

    % --- List files for previous session (if it exists) ---
    compare_files = [];
    if ~isempty(session_earlier)
        localite_path_earlier = fullfile(pn.data_postlocalite, sub_id, session_earlier, 'localite', '*', 'TMSTrigger');
        [compare_files, dt_earlier] = list_valid_localite_files(localite_path_earlier);
        % Sort previous session files by date desc
        [~, sidx] = sort(dt_earlier, 'descend');
        compare_files = compare_files(sidx);
    end

    % Sort later session files by date desc (most recent first)
    [~, sidx] = sort(dt_later, 'descend');
    localite_files = localite_files(sidx);

    % --- Main selection loop ---
    status = 0;
    while status == 0 && ~isempty(localite_files)
        current_file = [];

        % Select first unused valid file (most recent)
        for k = 1:length(localite_files)
            if ~ismember(localite_files(k).name, {excludeFile.name})
                current_file = localite_files(k);
                break;
            end
        end

        if isempty(current_file)
            warning("⚠ No Localite file could be selected for %s %s", sub_id, session);
            return;
        end

        sel_file_path = fullfile(current_file.folder, current_file.name);

        try
            current_localite = readstruct(sel_file_path);
        catch
            warning("⚠ Failed to read file: %s", sel_file_path);
            excludeFile(end+1).name = current_file.name;
            continue;
        end

        % If previous session files exist, compare files to skip duplicates
        if ~isempty(compare_files)
            is_duplicate = false;
            for cfi = 1:length(compare_files)
                earlier_file_path = fullfile(compare_files(cfi).folder, compare_files(cfi).name);
                cmd = sprintf('diff "%s" "%s"', sel_file_path, earlier_file_path);
                [cmp_status, ~] = system(cmd);
                if cmp_status == 0
                    disp("⚠ File is identical to previous session — skipping.");
                    excludeFile(end+1).name = current_file.name;
                    is_duplicate = true;
                    break;
                end
            end
            if is_duplicate
                continue; % skip to next later session file
            else
                disp("File differs from all previous session files — selecting.");
            end
        end

        % Select this Localite file
        localite = current_localite;
        break;
    end

    if isempty(localite)
        warning("⚠ No valid Localite XML could be selected for %s %s", sub_id, session);
    end
end
