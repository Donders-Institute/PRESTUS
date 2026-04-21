function response = confirmation_dlg(promptMessage, yes_text, no_text)
% CONFIRMATION_DLG  Display a two-button modal dialog and return the user's choice
%
%   Creates a questdlg with two customizable button labels and returns 1 if
%   the user selects the affirmative option, 0 otherwise.
%
% Use as:
%   response = confirmation_dlg(promptMessage, yes_text, no_text)
%
% Input:
%   promptMessage - message to display in the dialog
%   yes_text      - label for the affirmative button
%   no_text       - label for the negative button (also the default)
%
% Output:
%   response      - [1x1] 1 = affirmative, 0 = negative or dialog closed
%
% See also: QUESTDLG, CONFIRM_OVERWRITING

arguments
    promptMessage (1,:) char
    yes_text      (1,:) char
    no_text       (1,:) char
end

    % Display a confirmation dialog box with custom button labels
    selection = questdlg(promptMessage, 'Confirmation needed', ...
               yes_text, no_text, no_text);

    % Check user selection and return corresponding response
    if strcmp(selection, yes_text)
        response = 1; % User selected "Yes"
    else 
        response = 0; % User selected "No"
    end
end
