function parameters = simnibs_version(segmentation_folder, parameters)
% SIMNIBS_VERSION  Extract SimNIBS version from charm_log.html and store in parameters
%
% Reads the first line of charm_log.html, parses the SimNIBS version string,
% and stores it in parameters.simnibs_version. Falls back to 'unknown' if
% the file is missing or the version cannot be parsed.
%
% Use as:
%   parameters = simnibs_version(segmentation_folder, parameters)
%
% Input:
%   segmentation_folder - path to folder containing charm_log.html
%   parameters          - PRESTUS config struct to update
%
% Output:
%   parameters - updated with parameters.simnibs_version field
%
% See also: KWAVE_VERSION, PRESTUS_VERSION

arguments
    segmentation_folder (1,:) char
    parameters          (1,1) struct
end

    log_file = fullfile(segmentation_folder, 'charm_log.html');
    
    if ~exist(log_file, 'file')
        warning('charm_log.html not found in %s', segmentation_folder);
        parameters.simnibs_version = 'unknown';
        fprintf('LOG: SimNIBS version: UNKNOWN (no charm_log.html)\n');
        return;
    end
    
    try
        % Read first few lines (fast)
        fid = fopen(log_file, 'r');
        first_line = fgetl(fid);
        fclose(fid);
        
        % Parse "INFO: simnibs version X.Y.Z"
        pattern = 'INFO:\s*simnibs\s*version\s*([\d.]+)';
        tokens = regexp(first_line, pattern, 'tokens');
        
        if ~isempty(tokens)
            version = tokens{1}{1};
            parameters.simnibs_version = version;
            fprintf('SimNIBS version %s detected from %s\n', version, log_file);
        else
            parameters.simnibs_version = 'unknown';
            fprintf('SimNIBS version UNKNOWN (could not parse "%s")\n', first_line);
        end
        
    catch ME
        warning('Error parsing %s: %s', log_file, ME.message);
        parameters.simnibs_version = 'error';
        fprintf('SimNIBS version: ERROR\n');
    end
end
