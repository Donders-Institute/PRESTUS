function overwrite = confirm_overwriting(filename, parameters)
% CONFIRM_OVERWRITING  Decide whether an existing file should be overwritten
%
%   Returns 1 (overwrite) if the file does not exist or if
%   parameters.io.overwrite_files is 'always' or 'force'. Returns 0 if it is
%   'never'. In interactive mode, prompts the user via a dialog. Raises an
%   error in non-interactive mode for any other value.
%
%   'always' and 'force' differ only in chained pipelines (async, multi-ISPPA,
%   uncertainty), where 'force' also re-runs cached intermediate stages; see
%   RESOLVE_CHAINED_OVERWRITE.
%
% Use as:
%   overwrite = confirm_overwriting(filename, parameters)
%
% Input:
%   filename   - path of the file to check
%   parameters - PRESTUS config with io.overwrite_files
%                ('always'|'force'|'never'|'ask')
%                and simulation.interactive
%
% Output:
%   overwrite  - [1x1] 1 = proceed with writing, 0 = skip
%
% See also: CONFIRMATION_DLG, RESOLVE_CHAINED_OVERWRITE

    if ~exist(filename, 'file') || any(strcmp(parameters.io.overwrite_files, {'always', 'force'}))
        overwrite = 1; % Overwrite allowed
        return;
    % If overwriting is set to 'never', do not overwrite
    elseif strcmp(parameters.io.overwrite_files, 'never')
        overwrite = 0; % Do not overwrite
        return;
    end

    % Ensure valid behavior in non-interactive mode
    if parameters.simulation.interactive == 0
        error("In non-interactive mode, parameters.io.overwrite_files should be 'always', 'force' or 'never'");
    end

    % Interactive mode: Prompt user for confirmation to overwrite the file
    promptMessage = sprintf('%s already exists. \nDo you want to overwrite it?', filename);
    overwrite = confirmation_dlg(promptMessage, 'Overwrite', 'Cancel'); % User confirmation dialog
end
