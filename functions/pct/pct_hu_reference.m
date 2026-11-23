function [hu_ref, info] = pct_hu_reference(parameters)
% PCT_HU_REFERENCE  HU that max-normalised skull mappings take as the densest bone
%
% Several HU mappings scale by the densest bone (density 'marsac-babelbrain'
% and 'aubry', sound speed 'aubry', attenuation 'mueller' and 'aubry'). The
% reference is a percentile of the CT within layers of the segmentation over
% the whole head (pct.hu_reference_layers, default the skull layers;
% pct.hu_reference_percentile), on the segmentation grid, after the HU limits
% (pct.limits.hu_*). It does not depend on the simulation grid, its crop or
% the interpolation onto it. A percentile below 100 keeps single bright voxels
% (metal, artefacts) from setting the scale of the whole skull; HU above the
% reference map as the densest bone (MEDIUM_PCT_POROSITY), in the BabelBrain
% bridge too (run_babelbrain_case.py --ct_hu_max_reference). 100 takes the
% largest HU, so nothing is limited.
%
% With pct.hu_reference_resolution_mm set, the CT is first resampled over
% the whole head onto an isotropic grid of that spacing with
% pct.interpolation (labels nearest neighbour), then limited, and the
% percentile is taken there. This reproduces tools that take the densest bone
% after resampling (BabelBrain: cubic spline onto its simulation grid). Cubic
% and spline interpolation overshoot at the skull surface, so the reference
% then rises with the overshoot; it is meant for comparisons with such tools.
%
% The CT and its tissue labels are read with PCT_LOAD_CT; the result is kept
% for the MATLAB session, keyed by the CT file, its date and the settings.
%
% Use as:
%   [hu_ref, info] = pct_hu_reference(parameters)
%
% Input:
%   parameters - PRESTUS config with pct.enabled and the CT inputs (PCT_LOAD_CT)
%
% Output:
%   hu_ref - reference HU
%   info   - struct: layers, labels, percentile, n_voxels, hu_max (largest HU
%            in the layers), source (CT file)
%
% See also: MEDIUM_PCT_HU_MAX, PCT_MAP_PROPERTIES, PCT_LOAD_CT

    arguments
        parameters (1,1) struct
    end
    persistent memo
    if isempty(memo), memo = containers.Map(); end

    layers = get_field(parameters, {'pct', 'hu_reference_layers'}, {'skull'});
    layers = cellstr(layers);
    pctl = double(get_field(parameters, {'pct', 'hu_reference_percentile'}, 99.9));   % as config_default.yaml
    if ~isscalar(pctl) || pctl <= 0 || pctl > 100
        error('pct_hu_reference:badPercentile', 'pct.hu_reference_percentile must be in (0, 100].');
    end
    labels = [];
    for k = 1:numel(layers)
        if ~isfield(parameters.layers, layers{k})
            error('pct_hu_reference:badLayer', 'pct.hu_reference_layers: ''%s'' is not a layer of this run.', layers{k});
        end
        labels = [labels, getidx(parameters.layers, layers{k})]; %#ok<AGROW>
    end
    labels = unique(labels);

    ct_file = pct_ct_file(parameters);
    if isempty(ct_file), ct_file = fullfile(find_pseudoct_dir(parameters), 'pseudoCT.nii.gz'); end
    d = dir(ct_file);
    lim = get_field(parameters, {'pct', 'limits'}, struct());
    res = get_field(parameters, {'pct', 'hu_reference_resolution_mm'}, []);
    if ~isempty(res) && ~(isscalar(res) && isnumeric(res) && res > 0)
        error('pct_hu_reference:badResolution', 'pct.hu_reference_resolution_mm must be [] or a positive spacing [mm].');
    end
    method = char(get_field(parameters, {'pct', 'interpolation'}, 'linear'));
    key = jsonencode(struct('f', ct_file, 'd', datestr_safe(d), 'l', labels, 'p', pctl, ...
        'lim', lim, 'i', method, 'r', res));
    if isKey(memo, key)
        s = memo(key); hu_ref = s.hu_ref; info = s.info;
        return
    end

    [ct, ct_hdr, mask] = pct_load_ct(parameters);
    in_layers = ismember(mask, labels);
    if ~any(in_layers(:))
        error('pct_hu_reference:empty', 'No voxel of the segmentation is in pct.hu_reference_layers (labels %s).', ...
            mat2str(labels));
    end
    if isempty(res)
        hu = ct(in_layers);
    else
        hu = resampled_in_layers(ct, in_layers, double(ct_hdr.PixelDimensions(1:3)), res, method);
    end
    hu = sort(double(pct_apply_limits(parameters, 'hu', hu)));
    hu_ref = hu(max(1, ceil(pctl / 100 * numel(hu))));        % nearest rank; no Statistics Toolbox
    info = struct('layers', {layers}, 'labels', labels, 'percentile', pctl, 'n_voxels', numel(hu), ...
        'hu_max', hu(end), 'hu_ref', hu_ref, 'resolution_mm', res, 'source', ct_file);
    grid_txt = 'segmentation grid';
    if ~isempty(res), grid_txt = sprintf('%g mm grid (%s)', res, method); end
    fprintf('[pct_hu_reference] %g th percentile of %d voxels in layers %s, %s: %.1f HU (largest %.1f HU)\n', ...
        pctl, numel(hu), strjoin(layers, ', '), grid_txt, hu_ref, hu(end));
    memo(key) = struct('hu_ref', hu_ref, 'info', info);
end

function hu = resampled_in_layers(ct, in_layers, vox_mm, res, method)
% CT values at the points of an isotropic res-mm grid over the layer voxels
% (nearest label in the layers), interpolated with method; slab by slab, so
% fine grids over the whole head fit in memory
    F = griddedInterpolant(double(ct), method, 'nearest');
    L = griddedInterpolant(single(in_layers), 'nearest', 'nearest');
    [i, j, k] = ind2sub(size(in_layers), find(in_layers));
    ax = cell(1, 3);
    lo = [min(i) min(j) min(k)]; hi = [max(i) max(j) max(k)];
    for a = 1:3
        ax{a} = lo(a) - 1 : res / vox_mm(a) : hi(a) + 1;    % voxel coordinates
    end
    hu = cell(numel(ax{3}), 1);
    for z = 1:numel(ax{3})
        sel = L({ax{1}, ax{2}, ax{3}(z)}) > 0.5;
        if ~any(sel(:)), continue; end
        v = F({ax{1}, ax{2}, ax{3}(z)});
        hu{z} = single(v(sel));
    end
    hu = vertcat(hu{:});
end

function s = datestr_safe(d)
    s = '';
    if ~isempty(d), s = sprintf('%.6f_%d', d(1).datenum, d(1).bytes); end
end
