function overwrite = confirm_overwriting(filename, parameters)
% CONFIRM_OVERWRITING  Decide whether an existing file should be overwritten
%
%   Returns 1 (overwrite) if the file does not exist or if
%   parameters.io.overwrite_files is 'always'. Returns 0 if it is 'never'.
%   In interactive mode, prompts the user via a dialog. Raises an error in
%   non-interactive mode if overwrite_files is neither 'always' nor 'never'.
%
% Use as:
%   overwrite = confirm_overwriting(filename, parameters)
%
% Input:
%   filename   - path of the file to check
%   parameters - PRESTUS config with io.overwrite_files ('always'|'never'|'')
%                and simulation.interactive
%
% Output:
%   overwrite  - [1x1] 1 = proceed with writing, 0 = skip
%
% See also: CONFIRMATION_DLG

    % Check if the file does not exist or if overwriting is set to 'always'
    if ~exist(filename, 'file') || strcmp(parameters.io.overwrite_files, 'always')
        overwrite = 1; % Overwrite allowed
        return;
    % If overwriting is set to 'never', do not overwrite
    elseif strcmp(parameters.io.overwrite_files, 'never')
        overwrite = 0; % Do not overwrite
        return;
    end

    % Ensure valid behavior in non-interactive mode
    if parameters.simulation.interactive == 0
        error("In non-interactive mode, parameters.io.overwrite_files should be 'always' or 'never'");
    end

    % Interactive mode: Prompt user for confirmation to overwrite the file
    promptMessage = sprintf('%s already exists. \nDo you want to overwrite it?', filename);
    overwrite = confirmation_dlg(promptMessage, 'Overwrite', 'Cancel'); % User confirmation dialog
end
