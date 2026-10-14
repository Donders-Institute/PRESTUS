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
% grid or medium) and, for grid and medium caches, the partial acoustic FOV
% (grid.acoustic_fov_*: cropped after these stages).
%
% Both the writing and the reusing run compute the signature at the same point
% of the pipeline, so settings that earlier stages derive (e.g. transducer
% positions from placement) are compared like the configured ones.
%
% Use as:
%   sig = cache_signature(parameters, 'grid')
%   sig = cache_signature(parameters, 'medium')
%   sig = cache_signature(parameters, 'acoustic')
%
% Input:
%   parameters - PRESTUS parameter struct
%   kind       - 'grid', 'medium' or 'acoustic' (sidecar of a reusable pressure map)
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
            sig.grid = drop_acoustic_fov(sig.grid);
        case 'medium'
            sig.grid = drop_acoustic_fov(sig.grid);
            sig.medium_properties = get_field(parameters, {'medium_properties'});
        case 'acoustic'
            % A pressure map reused by a later run (path.acoustic_input): its
            % grid and medium, but not the thermal-only grid settings, so the
            % thermal stage can be re-run on another grid.
            sig.subject_id = get_field(parameters, {'subject_id'});
            sig.medium_properties = get_field(parameters, {'medium_properties'});
            if isstruct(sig.grid)
                g = fieldnames(sig.grid);
                sig.grid = rmfield(sig.grid, g(startsWith(g, 'thermal_')));
            end
        otherwise
            error('cache_signature:badKind', 'Unknown cache kind ''%s'' (grid | medium | acoustic).', kind);
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

function g = drop_acoustic_fov(g)
% The partial acoustic FOV is cropped after the grid and medium stages
% (ACOUSTIC_GRID_FOV), so their caches serve runs with or without one.
    if ~isstruct(g)
        return
    end
    names = fieldnames(g);
    g = rmfield(g, names(startsWith(names, 'acoustic_fov_')));
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
