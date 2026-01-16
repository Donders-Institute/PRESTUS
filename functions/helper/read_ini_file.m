function iniData = read_ini_file(filename)
    % Reads an INI file into a MATLAB struct.
    % 
    % Arguments:
    % - filename: Path to the INI file.

    % Returns:
    % - iniData: Struct containing the parsed data.

    % Open the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open the file: %s', filename);
    end

    iniData = struct();
    currentSection = '';
    
    % Read the file line by line
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        
        % Skip empty lines or comments
        if isempty(line) || startsWith(line, ';') || startsWith(line, '#')
            continue;
        end
        
        % Check if the line is a section header
        if startsWith(line, '[') && endsWith(line, ']')
            currentSection = matlab.lang.makeValidName(line(2:end-1)); % Remove [ and ]
            iniData.(currentSection) = struct(); % Create a new section
            continue;
        end
        
        % Parse key-value pairs
        if contains(line, '=')
            tokens = strsplit(line, '=');
            key = strtrim(tokens{1});
            value = strtrim(tokens{2});
            
            % Convert numeric values if possible
            numericValue = str2double(value);
            if ~isnan(numericValue)
                value = numericValue;
            end
            
            % Assign the key-value pair to the current section
            % Convert key to a valid field name
            validKey = matlab.lang.makeValidName(key);
            if isempty(currentSection)
                iniData.(validKey) = value;
            else
                iniData.(currentSection).(validKey) = value;
            end
        end
    end
    
    % Close the file
    fclose(fid);
    
end