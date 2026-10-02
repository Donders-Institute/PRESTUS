function sig = cache_signature(parameters, kind)
% CACHE_SIGNATURE  Settings that determine a cached pipeline stage's output
%
% Grid and medium caches are named only by subject, medium and affix, so a
% cache written before a settings change would otherwise be reloaded as if
% nothing had changed. The signature collects the settings each stage's output
% depends on; it is stored with the cache and compared on reuse
% (CHECK_CACHE_SIGNATURE).
%
% Paths, I/O, HPC and module flags are excluded (moving data or changing what
% is run downstream does not change the cached result), as are the drive level
% (target_isppa_wcm2, tx_calibration: they scale the acoustic field, not the
% grid or medium).
%
% Both the writing and the reusing run compute the signature at the same point
% of the pipeline, so settings that earlier stages derive (e.g. transducer
% positions from placement) are compared like the configured ones.
%
% Use as:
%   sig = cache_signature(parameters, 'grid')
%   sig = cache_signature(parameters, 'medium')
%
% Input:
%   parameters - PRESTUS parameter struct
%   kind       - 'grid' or 'medium'
%
% Output:
%   sig - struct of the relevant settings
%
% See also: CHECK_CACHE_SIGNATURE, GRID_STAGE, MEDIUM_STAGE

    sig = struct();
    sig.signature_version = 1;
    sig.medium = get_field(parameters, {'simulation', 'medium'});
    for f = {'grid', 'layers', 'headmodel', 'pct'}
        sig.(f{1}) = get_field(parameters, f);
    end
    sig.transducer = strip_drive(get_field(parameters, {'transducer'}));

    switch kind
        case 'grid'
        case 'medium'
            sig.medium_properties = get_field(parameters, {'medium_properties'});
        otherwise
            error('cache_signature:badKind', 'Unknown cache kind ''%s'' (grid | medium).', kind);
    end
end

function v = get_field(s, path)
    v = [];
    for i = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{i})
            v = [];
            return
        end
        s = s(1).(path{i});
    end
    v = s;
end

function t = strip_drive(t)
    if ~isstruct(t)
        return
    end
    for f = {'target_isppa_wcm2', 'tx_calibration'}
        if isfield(t, f{1})
            t = rmfield(t, f{1});
        end
    end
end
