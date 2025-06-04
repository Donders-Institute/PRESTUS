function response = confirmation_dlg(promptMessage, yes_text, no_text)

% CONFIRMATION_DLG Displays a dialog box to confirm user action.
%
% This function creates a modal dialog box with customizable options for 
% confirmation. The user can choose between two options, and the function 
% returns a binary response based on the selection.
%
% Input:
%   promptMessage - String specifying the message displayed in the dialog box.
%   yes_text      - String specifying the text for the "Yes" button.
%   no_text       - String specifying the text for the "No" button.
%
% Output:
%   response      - Boolean flag indicating the user's choice (1 = Yes, 0 = No).

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
