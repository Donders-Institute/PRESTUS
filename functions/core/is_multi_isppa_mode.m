function tf = is_multi_isppa_mode(parameters)
% IS_MULTI_ISPPA_MODE  True when any transducer carries more than one ISPPA target.
%
% Multi-ISPPA mode is active when at least one entry in parameters.transducer
% has a target_isppa_wcm2 field containing more than one finite value.  Each
% unique target spawns a separate thermal job (see MULTI_ISPPA_PIPELINE).
%
% See also: MULTI_ISPPA_PIPELINE, ASYNC_TRANSDUCER_PIPELINE

    tf = false;
    if ~isfield(parameters, 'transducer')
        return;
    end
    for ti = 1:numel(parameters.transducer)
        t = parameters.transducer(ti);
        if isfield(t, 'target_isppa_wcm2') && ...
                numel(t.target_isppa_wcm2) > 1 && ...
                any(isfinite(t.target_isppa_wcm2))
            tf = true;
            return;
        end
    end
end
