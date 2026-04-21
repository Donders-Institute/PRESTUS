function check_availability(files_to_check)
% CHECK_AVAILABILITY  Verify that a list of files or glob patterns exist on disk
%
% Checks each entry in files_to_check. If an entry contains a wildcard (*),
% exactly one matching file must exist; otherwise the literal path must be
% a regular file. Errors if a wildcard matches zero or multiple files;
% warns if a literal path does not exist.
%
% Use as:
%   check_availability(files_to_check)
%
% Input:
%   files_to_check - [1xN] cell array of file paths or glob patterns
%
% See also: DIR, ISFILE

arguments
    files_to_check (1,:) cell
end

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
