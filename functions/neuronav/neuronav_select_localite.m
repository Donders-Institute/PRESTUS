function localite = neuronav_select_localite(pn, sub_id, ses_id, markertype)
% NEURONAV_SELECT_LOCALITE - Select most recent Localite XML for a session.
% Compares against Localite files from the immediately preceding session (if available) to skip duplicates.
%
% INPUT:
%   pn         - struct with .data_postlocalite field (FIRST input)
%   sub_id     - e.g., 'sub-010'
%   ses_id     - e.g., 3 (numeric) or 'ses-03'
%   markertype - 'TriggerMarkers' (default) or 'GUMMarkers'
%
% OUTPUT:
%   localite - struct parsed from selected Localite XML file ([] if none match)

    if nargin < 4 || isempty(markertype)
        markertype = 'TriggerMarkers';  % Default
    end
    
    localite = [];
    excludeFile = struct('name', {});

    % --- Harmonize ses_id to session string ---
    if isnumeric(ses_id)
        ses_num = double(ses_id);
        session = sprintf('ses-%02d', ses_num);
    else
        ses_num = regexp(ses_id, 'ses-(\d{2})', 'tokens', 'once');
        if isempty(ses_num)
            error('ses_id must be numeric or "ses-XX" format');
        end
        ses_num = str2double(ses_num{1});
        session = char(ses_id);
    end

    % Define previous session string (with zero padding)
    prev_ses_num = ses_num - 1;
    if prev_ses_num < 1
        session_earlier = ''; % No earlier session available
    else
        session_earlier = sprintf('ses-%02d', prev_ses_num);
    end

    % --- Helper function to list valid Localite files ---
    function [files, datetimes] = list_valid_localite_files(basepath, markertype)
        datetimes = datetime.empty;
        
        if strcmp(markertype, 'TriggerMarkers')
            pattern = 'TriggerMarkers_Coil0*.xml';
            search_path = fullfile(basepath, 'localite', '*', 'TMSTrigger');
            files = dir(fullfile(search_path, pattern));
            files = files([files.bytes] > 10000);  % 10kB ONLY for TriggerMarkers
        else  % GUMMarkers
            pattern = 'GUMMarkers*.xml';  % Matches GUMMarkers_sub-003_VS.xml
            search_path = fullfile(basepath, 'localite');
            files = dir(fullfile(search_path, pattern));
            % No size criterion for GUMMarkers
        end
        
        for i = 1:length(files)
            if strcmp(markertype, 'TriggerMarkers')
                % Extract filename timestamp for TriggerMarkers
                tokens = regexp(files(i).name, '_(\d{17})', 'tokens', 'once');
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
            else
                % Use file modification timestamp for GUMMarkers
                datetimes(end+1) = datetime(files(i).datenum, 'ConvertFrom', 'datenum');
            end
        end
        
        valid_idx = ~isnat(datetimes);
        files = files(valid_idx);
        datetimes = datetimes(valid_idx);
    end

    % --- List files for current session ---
    if strcmp(markertype, 'TriggerMarkers')
        localite_path_later = fullfile(pn.data_postlocalite, sub_id, session, 'localite', '*', 'TMSTrigger');
    else
        localite_path_later = fullfile(pn.data_postlocalite, sub_id, session);  % Session level for GUMMarkers
    end
    [localite_files, dt_later] = list_valid_localite_files(localite_path_later, markertype);

    if isempty(localite_files)
        if strcmp(markertype, 'TriggerMarkers')
            warning("⚠ No %s files >10kB with valid timestamps found for %s %s", markertype, sub_id, session);
        else
            warning("⚠ No %s files found for %s %s", markertype, sub_id, session);
        end
        return;
    end

    % --- List files for previous session (if it exists) ---
    compare_files = [];
    if ~isempty(session_earlier)
        if strcmp(markertype, 'TriggerMarkers')
            localite_path_earlier = fullfile(pn.data_postlocalite, sub_id, session_earlier, 'localite', '*', 'TMSTrigger');
        else
            localite_path_earlier = fullfile(pn.data_postlocalite, sub_id, session_earlier);
        end
        [compare_files, dt_earlier] = list_valid_localite_files(localite_path_earlier, markertype);
        if ~isempty(compare_files)
            [~, sidx] = sort(dt_earlier, 'descend');
            compare_files = compare_files(sidx);
        end
    end

    % --- Main selection loop ---
    [~, sidx] = sort(dt_later, 'descend');
    localite_files = localite_files(sidx);

    status = 0;
    while status == 0 && ~isempty(localite_files)
        current_file = [];
        for k = 1:length(localite_files)
            if ~ismember(localite_files(k).name, {excludeFile.name})
                current_file = localite_files(k);
                break;
            end
        end

        if isempty(current_file)
            warning("⚠ No %s file could be selected for %s %s", markertype, sub_id, session);
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
                continue;
            else
                disp("File differs from all previous session files — selecting.");
            end
        end

        localite = current_localite;
        break;
    end

    if isempty(localite)
        warning("⚠ No valid %s XML could be selected for %s %s", markertype, sub_id, session);
    end
end
