function response = confirmation_dlg(promptMessage, yes_text, no_text)
        selection = questdlg(promptMessage , 'Confirmation needed', ...
               yes_text, no_text, no_text);

        if strcmp(selection, yes_text)
          response = 1; 
        else 
          response = 0;
        end
end
