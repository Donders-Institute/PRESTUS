function telemetry_setup_reset()
% TELEMETRY_SETUP_RESET  Remove the telemetry consent file so the opt-in
% prompt is shown again on the next pipeline run.
%
% Use as:
%   telemetry_setup_reset()

    cfg_file = fullfile(prefdir_prestus(), 'telemetry.json');
    if isfile(cfg_file)
        delete(cfg_file);
        fprintf('Telemetry consent reset. You will be prompted again on the next run.\n');
    else
        fprintf('No telemetry consent file found (already unset).\n');
    end
end
