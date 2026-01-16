function check_availability(files_to_check)
    % Loop through each filename passed
    for i = 1:length(files_to_check)
        filename = files_to_check{i};
        if contains(filename, '*')
            matching_files = dir(filename);
            if length(matching_files) > 1
                error('More than 1 file matches the template %s', filename);
            elseif isempty(matching_files)
                error('No files match the template %s', filename);
            else
                filename = fullfile(matching_files.folder, matching_files.name);
            end
        end

        if ~isfile(filename)
            warning('File does not exist: \r\n%s', filename);
        end
    end
end
