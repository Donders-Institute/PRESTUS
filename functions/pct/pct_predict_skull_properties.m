function stats = pct_predict_skull_properties(parameters, options)
% PCT_PREDICT_SKULL_PROPERTIES  Skull properties a pseudo-CT run will model, without simulating
%
% Applies the configured pct.mapping_density / _soundspeed / _attenuation
% (MEDIUM_PCT_DENSITY, MEDIUM_PCT_SOUNDSPEED, MEDIUM_PCT_ATTENUATION, as in
% MEDIUM_SETUP) to the skull voxels of the pseudo-CT and summarises the
% mapped distributions. The skull voxels and HU values come from
%   'grid'    the head-preprocessing cache (cache/sub-XXX_<medium>_cropped_smoothed.mat),
%             i.e. exactly what MEDIUM_SETUP maps, when it exists; else from
%   'native'  pseudoCT.nii.gz and tissues_mask.nii.gz in the pCT folder
%             (FIND_PSEUDOCT_DIR), skull = layers.skull labels, before
%             cropping, resampling and smoothing.
% Alongside, the homogeneous skull values of the uncertainty bound variants
% (UNCERTAINTY_VARIANTS) are listed, so the mapped distribution can be
% compared with the range the bound variants model.
%
% Use as:
%   stats = pct_predict_skull_properties(parameters)
%   stats = pct_predict_skull_properties(parameters, write=false)
%
% Input:
%   parameters - PRESTUS config (LOAD_PARAMETERS); pct.enabled is not required
%
% Name-value options:
%   write  - write sub-XXX_<medium>_pct_prediction<affix>.json to the
%            output's cache folder (default true)
%   source - 'auto' (default) | 'grid' | 'native'
%
% Output:
%   stats - struct in the layout of the medium_stats cache of MEDIUM_SETUP
%           (.mapping, .n_voxels, .properties.<pseudoCT_HU|density|sound_speed|
%           alpha_coeff>, .reference.skull), plus .predicted = true, .source,
%           .file (written JSON or '') and .bounds.<variant> (skull values of
%           each bound variant's medium file)
%
% See also: MEDIUM_SETUP, PCT_SKULL_MAPPED, UNCERTAINTY_PIPELINE

arguments
    parameters (1,1) struct
    options.write (1,1) logical = true
    options.source (1,:) char {mustBeMember(options.source, {'auto', 'grid', 'native'})} = 'auto'
end

    subject_id = parameters.subject_id;
    medium     = parameters.simulation.medium;
    cache_dir  = fullfile(get_output_dir(parameters), 'cache');
    affix      = get_field(parameters, {'io', 'output_affix'}, '');

    [pseudoCT, skull_idx, source] = load_skull_hu(parameters, cache_dir, options.source);

    % Same calls as MEDIUM_SETUP; arrays outside the skull are not read
    mapping = struct('density',     mapping_name(parameters, 'mapping_density'), ...
                     'sound_speed', mapping_name(parameters, 'mapping_soundspeed'), ...
                     'attenuation', mapping_name(parameters, 'mapping_attenuation'));
    p = parameters;
    if ~isfield(p, 'simulation') || ~isfield(p.simulation, 'debug'); p.simulation.debug = 0; end
    skull = p.medium_properties.skull;
    density     = repmat(skull.density,     size(pseudoCT));
    sound_speed = repmat(skull.sound_speed, size(pseudoCT));
    alpha_coeff = repmat(skull.alpha_coeff, size(pseudoCT));
    alpha_power = repmat(skull.alpha_power, size(pseudoCT));
    density = medium_pct_density(p, density, pseudoCT, skull_idx, mapping.density);
    sound_speed = medium_pct_soundspeed(p, sound_speed, density, pseudoCT, skull_idx, mapping.sound_speed);
    alpha_coeff = medium_pct_attenuation(p, alpha_coeff, alpha_power, pseudoCT, skull_idx, ...
                                         mapping.attenuation, density);

    stats = struct();
    stats.subject_id   = subject_id;
    stats.medium       = medium;
    stats.output_affix = char(affix);
    stats.layer        = 'skull';
    stats.predicted    = true;
    stats.source       = source;
    stats.n_voxels     = numel(skull_idx);
    stats.mapping      = mapping;
    if strcmp(mapping.density, 'k-plan')
        [~, ~, stats.mapping.density_calibration] = medium_pct_density_calibration(p);
    end
    stats.properties.pseudoCT_HU = distribution_stats(pseudoCT(skull_idx));
    stats.properties.density     = distribution_stats(density(skull_idx));
    stats.properties.sound_speed = distribution_stats(sound_speed(skull_idx));
    stats.properties.alpha_coeff = distribution_stats(alpha_coeff(skull_idx));
    stats.reference.skull = struct('density', skull.density, 'sound_speed', skull.sound_speed, ...
                                   'alpha_coeff', skull.alpha_coeff);
    stats.bounds = bound_values(parameters);

    stats.file = '';
    if options.write
        if ~isfolder(cache_dir); mkdir(cache_dir); end
        stats.file = fullfile(cache_dir, sprintf('sub-%03d_%s_pct_prediction%s.json', subject_id, medium, affix));
        fid = fopen(stats.file, 'w');
        if fid == -1
            error('pct_predict_skull_properties:write', 'Cannot write %s.', stats.file);
        end
        fwrite(fid, jsonencode(stats, 'PrettyPrint', true));
        fclose(fid);
    end
end

function [pseudoCT, skull_idx, source] = load_skull_hu(parameters, cache_dir, source)
    medium_labels = fieldnames(parameters.medium_properties);
    grid_file = fullfile(cache_dir, sprintf('sub-%03d_%s_cropped_smoothed.mat', ...
        parameters.subject_id, parameters.simulation.medium));
    if strcmp(source, 'auto')
        source = 'native';
        if isfile(grid_file) && ~isempty(whos('-file', grid_file, 'pseudoCT_crop'))
            source = 'grid';
        end
    end
    switch source
        case 'grid'
            c = load(grid_file, 'medium_masks', 'pseudoCT_crop');
            if ~isfield(c, 'pseudoCT_crop') || isempty(c.pseudoCT_crop)
                error('pct_predict_skull_properties:noGridPct', ...
                    '%s holds no pseudo-CT (preprocessed with pct.enabled = 0?).', grid_file);
            end
            pseudoCT  = double(c.pseudoCT_crop);
            skull_idx = find(c.medium_masks == find(strcmp(medium_labels, 'skull')));
        case 'native'
            pct_dir = find_pseudoct_dir(parameters);
            if isempty(pct_dir)
                error('pct_predict_skull_properties:notFound', ...
                    'No pseudoCT.nii.gz found in:\n  %s', get_pct_dir(parameters));
            end
            pseudoCT = double(niftiread(fullfile(pct_dir, 'pseudoCT.nii.gz')));
            mask = niftiread(fullfile(pct_dir, 'tissues_mask.nii.gz'));
            skull_idx = find(ismember(mask, parameters.layers.skull));
    end
    if isempty(skull_idx)
        error('pct_predict_skull_properties:noSkull', 'No skull voxels found (%s source).', source);
    end
end

function name = mapping_name(parameters, key)
    name = 'none';
    if isfield(parameters, 'pct') && isfield(parameters.pct, key) && ~isempty(parameters.pct.(key))
        name = char(parameters.pct.(key));
    end
end

function bounds = bound_values(parameters)
% Skull values of the uncertainty bound variants' medium files (unreadable files are skipped).
    bounds = struct();
    variants = uncertainty_variants();
    for v = variants(2:end)
        try
            file = yaml.loadFile(uncertainty_medium_config(parameters, v.name), 'ConvertToArray', true);
            s = file.medium_properties.skull;
            bounds.(v.name) = struct('label', v.label, 'density', s.density, ...
                                     'sound_speed', s.sound_speed, 'alpha_coeff', s.alpha_coeff);
        catch
        end
    end
end

function v = get_field(s, path, default)
    v = default;
    for k = 1:numel(path)
        if ~isstruct(s) || ~isfield(s, path{k}); return; end
        s = s.(path{k});
    end
    v = s;
end
