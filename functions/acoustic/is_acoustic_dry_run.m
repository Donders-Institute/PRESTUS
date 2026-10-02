function dry = is_acoustic_dry_run(parameters)
% IS_ACOUSTIC_DRY_RUN  True when the configured acoustic backend only writes inputs
%
% jwave.dry_run / babelbrain.dry_run make the backend write its solver inputs
% and return no field (same test as JWAVE_WRAPPER / BABELBRAIN_WRAPPER). The
% pipeline then stops after the acoustic stage, caches nothing and skips the
% free-water baseline. k-Wave has no dry run.
%
% Use as:
%   dry = is_acoustic_dry_run(parameters)
%
% Input:
%   parameters - pipeline parameters struct
%
% Output:
%   dry - logical scalar
%
% See also: ACOUSTIC_STAGE, ACOUSTIC_WATER_BASELINE

arguments
    parameters (1,1) struct
end

    dry = false;
    backend = 'kwave';
    if isfield(parameters, 'acoustic_backend') && ~isempty(parameters.acoustic_backend)
        backend = lower(char(parameters.acoustic_backend));
    end
    if ~any(strcmp(backend, {'jwave', 'babelbrain'}))
        return;
    end
    if isfield(parameters, backend) && isfield(parameters.(backend), 'dry_run')
        flag = parameters.(backend).dry_run;
        dry  = ~isempty(flag) && logical(flag(1));
    end
end
