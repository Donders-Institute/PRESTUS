function overwrite = confirm_overwriting(filename, parameters)
    if ~exist(filename,'file') || strcmp(parameters.overwrite_files, 'always')
        overwrite = 1;
        return;
    elseif strcmp(parameters.overwrite_files, 'never')
        overwrite = 0;
        return
    end
    if parameters.interactive == 0 
        error("In non-interactive mode, parameters.overwrite_files should be 'always' or 'never'")
    end
    
    promptMessage = sprintf('%s already exists. \nDo you want to overwrite it?',filename);
    overwrite = confirmation_dlg(promptMessage , 'Overwrite','Cancel');
end