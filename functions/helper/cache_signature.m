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
% is run downstream does not change the cached result); a phantom's file is
% identified by subject, size and modification date, a head segmentation by its
% content (SEGMENTATION_FINGERPRINT), so a re-run of charm is detected. Excluded are the drive level
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
%   sig = cache_signature(parameters, 'headmodel')
%
% Input:
%   parameters - PRESTUS parameter struct
%   kind       - 'grid', 'medium', 'acoustic' (sidecar of a reusable pressure map)
%                or 'headmodel' (HEAD_MODEL_LOAD_OR_BUILD)
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
    % A measured CT (path.ct_pattern) replaces the pseudo-CT: identify the file
    % (field only when set, so caches of pseudo-CT and label runs keep their signature)
    ct = measured_ct_identity(parameters);
    if ~isempty(ct)
        sig.measured_ct = ct;
    end
    % A phantom is read from path.seg/m2m_sub-*/final_tissues.nii.gz: a
    % regenerated phantom changes the grid, medium and field built from it
    if strcmp(sig.medium, 'phantom') && ~strcmp(kind, 'headmodel')
        sig.phantom_file = segmentation_identity(parameters);
    end
    % A head is built from the SimNIBS segmentation: a re-run of charm (another
    % SimNIBS version or other settings) changes everything built from it
    if strcmp(sig.medium, 'layered') || strcmp(kind, 'headmodel')
        fp = segmentation_fingerprint(parameters);
        if ~isempty(fp)
            sig.segmentation = fp;
        end
    end

    switch kind
        case 'grid'
            sig.grid = drop_acoustic_fov(sig.grid);
        case 'medium'
            sig.grid = drop_acoustic_fov(sig.grid);
            % an uncertainty bound variant of a pseudo-CT run: homogeneous skull
            % (field only when set, so other caches keep their signature)
            if ~pct_skull_mapped(parameters) && isstruct(sig.pct) && isfield(sig.pct, 'enabled') && ...
                    isequal(sig.pct.enabled, 1)
                sig.homogeneous_skull = true;
            end
            % a tissue's medium uses absorption_coeff or absorption_fraction, whichever is set
            sig.medium_properties = drop_absorption_coeff(get_field(parameters, {'medium_properties'}), parameters);
        case 'acoustic'
            % A pressure map reused by a later run (path.acoustic_input): its
            % grid and medium, but not the thermal-only grid settings, so the
            % thermal stage can be re-run on another grid.
            sig.subject_id = get_field(parameters, {'subject_id'});
            sig.medium_properties = drop_absorption_coeff(get_field(parameters, {'medium_properties'}), []);
            if isstruct(sig.grid)
                g = fieldnames(sig.grid);
                sig.grid = rmfield(sig.grid, g(startsWith(g, 'thermal_')));
            end
        case 'headmodel'
            % The pose-independent head model (HEAD_MODEL_LOAD_OR_BUILD): its
            % resolution, smoothing / hole filling, labels and pseudo-CT, and
            % the segmentation it was made from; not the transducer or the crop
            sig = rmfield(sig, {'transducer', 'grid', 'headmodel'});
            sig.subject_id    = get_field(parameters, {'subject_id'});
            sig.resolution_mm = get_field(parameters, {'grid', 'resolution_mm'});
            hmod = get_field(parameters, {'headmodel'});
            sig.headmodel = struct();
            for f = {'smooth_method', 'smooth_threshold_skull', 'smooth_threshold_other', 'smooth_fwhm_mm', ...
                     'skull_fill_method', 'skull_wrap_radius'}
                if isstruct(hmod) && isfield(hmod, f{1}), sig.headmodel.(f{1}) = hmod.(f{1}); end
            end
            mp = get_field(parameters, {'medium_properties'});
            if isstruct(mp), sig.tissue_labels = fieldnames(mp); end
        otherwise
            error('cache_signature:badKind', 'Unknown cache kind ''%s'' (grid | medium | acoustic | headmodel).', kind);
    end
end

function id = segmentation_identity(parameters)
% Subject, size and modification date of the segmentation (final_tissues)
    id = struct('subject_id', get_field(parameters, {'subject_id'}), 'bytes', [], 'date', []);
    try
        f = dir(fullfile(parameters.path.seg, sprintf('m2m_sub-%03d', parameters.subject_id), 'final_tissues.nii.gz'));
        if ~isempty(f), id.bytes = f(1).bytes; id.date = round(f(1).datenum * 86400); end   % [s]
    catch
    end
end

function id = measured_ct_identity(parameters)
% Name, size and modification date of the measured CT; [] without one or with pct off
    id = [];
    pct = get_field(parameters, {'pct'});
    if ~(isstruct(pct) && isfield(pct, 'enabled') && isequal(pct.enabled, 1))
        return
    end
    pattern = get_field(parameters, {'path', 'ct_pattern'});
    if ~((ischar(pattern) || isstring(pattern)) && strlength(string(pattern)) > 0)
        return
    end
    id = struct('name', char(pattern), 'bytes', [], 'date', []);
    try
        f = dir(pct_ct_file(parameters));
        if ~isempty(f), id.name = f(1).name; id.bytes = f(1).bytes; id.date = round(f(1).datenum * 86400); end   % [s]
    catch
    end
end

function mp = drop_absorption_coeff(mp, parameters)
% absorption_coeff is thermal-only and unused by tissues that set
% absorption_fraction, as is an emptied absorption_fraction by tissues that
% set absorption_coeff; parameters = [] drops absorption_coeff from every tissue
    if ~isstruct(mp)
        return
    end
    for t = fieldnames(mp)'
        if ~isstruct(mp.(t{1}))
            continue
        end
        by_fraction = isempty(parameters) || strcmp(thermal_absorption_model(parameters, t{1}), 'fraction') && ...
            ~thermal_absorption_forced(parameters, t{1});
        if by_fraction && isfield(mp.(t{1}), 'absorption_coeff')
            mp.(t{1}) = rmfield(mp.(t{1}), 'absorption_coeff');
        elseif ~by_fraction && isfield(mp.(t{1}), 'absorption_fraction')
            mp.(t{1}) = rmfield(mp.(t{1}), 'absorption_fraction');
        end
    end
end

function tf = thermal_absorption_forced(parameters, tissue)
    [~, tf] = thermal_absorption_model(parameters, tissue);
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
