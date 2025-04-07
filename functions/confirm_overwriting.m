function overwrite = confirm_overwriting(filename, parameters)

% CONFIRM_OVERWRITING Handles file overwriting based on user preferences and interactivity.
%
% This function checks whether a file should be overwritten based on its existence, 
% user-defined parameters (`parameters.overwrite_files`), and interactive mode settings. 
% If interactive mode is enabled, the user is prompted to confirm overwriting.
%
% Input:
%   filename    - String specifying the name of the file to check.
%   parameters  - Struct containing overwrite preferences and interactivity settings:
%                 * overwrite_files: 'always', 'never', or interactive confirmation.
%                 * interactive: 1 for interactive mode, 0 for non-interactive mode.
%
% Output:
%   overwrite   - Boolean flag indicating whether the file should be overwritten (1 = yes, 0 = no).

    % Check if the file does not exist or if overwriting is set to 'always'
    if ~exist(filename, 'file') || strcmp(parameters.overwrite_files, 'always')
        overwrite = 1; % Overwrite allowed
        return;
    % If overwriting is set to 'never', do not overwrite
    elseif strcmp(parameters.overwrite_files, 'never')
        overwrite = 0; % Do not overwrite
        return;
    end

    % Ensure valid behavior in non-interactive mode
    if parameters.interactive == 0
        error("In non-interactive mode, parameters.overwrite_files should be 'always' or 'never'");
    end

    % Interactive mode: Prompt user for confirmation to overwrite the file
    promptMessage = sprintf('%s already exists. \nDo you want to overwrite it?', filename);
    overwrite = confirmation_dlg(promptMessage, 'Overwrite', 'Cancel'); % User confirmation dialog
end
